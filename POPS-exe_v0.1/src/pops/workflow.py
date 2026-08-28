from __future__ import annotations

import concurrent.futures
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from . import __version__
from .command import CommandRunner, shell_join
from .config import ToolConfig, load_config
from .fasta import fasta_lengths, filter_fasta_by_length, write_ranked_fasta
from .manifests import Library, Sample, read_manifest, validate_unique_across
from .metrics import merge_metric_matrix, write_contig_bed, write_sample_metrics
from .ranking import calculate_rankings, read_numeric_matrix, write_ranking_tables


class WorkflowError(RuntimeError):
    pass


TOOLS = ["megahit", "mmseqs", "bowtie2", "bowtie2-build", "samtools", "bedtools"]


def _sha256(path: Path, block: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(block)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _file_signature(path: Path) -> dict:
    st = path.stat()
    return {"path": str(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def _input_signature(group_manifest: Path, control_manifest: Path, cfg_path: Path,
                     group: list[Sample], control: list[Sample], x: int, y: int) -> dict:
    return {
        "group_manifest": {"path": str(group_manifest), "sha256": _sha256(group_manifest)},
        "control_manifest": {"path": str(control_manifest), "sha256": _sha256(control_manifest)},
        "config": {"path": str(cfg_path), "sha256": _sha256(cfg_path)},
        "x": x,
        "y": y,
        "fastqs": [
            _file_signature(p)
            for s in group + control
            for lib in s.libraries
            for p in ([lib.read1] + ([lib.read2] if lib.read2 else []))
        ],
    }


def _require_tools() -> None:
    missing = [tool for tool in TOOLS if shutil.which(tool) is None]
    if missing:
        raise WorkflowError("Required executable(s) not found on PATH: " + ", ".join(missing))


def _tool_version(tool: str) -> str:
    variants = {
        "megahit": ["megahit", "--version"],
        "mmseqs": ["mmseqs", "version"],
        "bowtie2": ["bowtie2", "--version"],
        "bowtie2-build": ["bowtie2-build", "--version"],
        "samtools": ["samtools", "--version"],
        "bedtools": ["bedtools", "--version"],
    }
    try:
        p = subprocess.run(variants[tool], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=20)
        return (p.stdout or "").splitlines()[0].strip()
    except Exception as exc:
        return f"version query failed: {exc}"


def _write_versions(path: Path) -> None:
    with path.open("w") as out:
        out.write("software\tversion\n")
        out.write(f"POPS\t{__version__}\n")
        out.write(f"Python\t{sys.version.split()[0]}\n")
        for t in TOOLS:
            out.write(f"{t}\t{_tool_version(t)}\n")


def _stage_done(stage_dir: Path, name: str) -> bool:
    return (stage_dir / f"{name}.done").is_file()


def _mark_done(stage_dir: Path, name: str) -> None:
    (stage_dir / f"{name}.done").write_text("ok\n")


def _megahit_inputs(samples: list[Sample]) -> list[str]:
    libraries = [lib for sample in samples for lib in sample.libraries]
    paired_r1 = [str(lib.read1) for lib in libraries if lib.read2]
    paired_r2 = [str(lib.read2) for lib in libraries if lib.read2]
    singles = [str(lib.read1) for lib in libraries if not lib.read2]
    args: list[str] = []
    if paired_r1:
        args += ["-1", ",".join(paired_r1), "-2", ",".join(paired_r2)]
    if singles:
        args += ["-r", ",".join(singles)]
    return args


def _safe_component(value: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in value)
    return safe or "library"


def _map_library(library: Library, sample_id: str, index_prefix: Path, out_bam: Path,
                 cfg: ToolConfig, runner: CommandRunner) -> None:
    """Map one sequencing library and tag alignments with its library id.

    The RG tag is used downstream to distinguish reads with identical QNAMEs that
    originated from different libraries of the same biological sample.
    """
    bt = [
        "bowtie2", f"--{cfg.bowtie2_preset}", "-a", "-x", str(index_prefix),
        "-p", str(cfg.threads_bowtie2),
        "--rg-id", library.library_id,
        "--rg", f"SM:{sample_id}",
        "--rg", f"LB:{library.library_id}",
    ]
    if library.read2:
        bt += ["-1", str(library.read1), "-2", str(library.read2)]
    else:
        bt += ["-U", str(library.read1)]
    pipeline = (
        f"{shell_join(bt)} | samtools view -b - | "
        f"samtools sort -@ {cfg.threads_bowtie2} -o {shell_join([out_bam])} -"
    )
    runner.run_pipeline(pipeline)


def _run_mapping(sample: Sample, rep_fasta: Path, index_prefix: Path, contig_bed: Path,
                 sample_dir: Path, cfg: ToolConfig, runner: CommandRunner, keep_bam: bool) -> Path:
    sample_work = sample_dir / _safe_component(sample.sample_id)
    sample_work.mkdir(parents=True, exist_ok=True)
    metric_tsv = sample_work / f"{_safe_component(sample.sample_id)}.metrics.tsv"
    if metric_tsv.is_file():
        return metric_tsv

    library_bams: list[Path] = []
    for i, lib in enumerate(sample.libraries, start=1):
        lib_bam = sample_work / f"lib{i:03d}.{_safe_component(lib.library_id)}.bam"
        _map_library(lib, sample.sample_id, index_prefix, lib_bam, cfg, runner)
        library_bams.append(lib_bam)

    bam = sample_work / f"{_safe_component(sample.sample_id)}.bam"
    if len(library_bams) == 1:
        shutil.copy2(library_bams[0], bam)
    else:
        runner.run(["samtools", "merge", "-f", "-@", str(cfg.threads_bowtie2), bam, *library_bams])
    runner.run(["samtools", "index", bam])

    breadth = sample_work / f"{_safe_component(sample.sample_id)}.breadth.tsv"
    depth = sample_work / f"{_safe_component(sample.sample_id)}.depth.tsv"
    with breadth.open("w") as fh:
        runner.run(["bedtools", "coverage", "-a", contig_bed, "-b", bam, "-sorted"], stdout=fh)
    with depth.open("w") as fh:
        runner.run(["bedtools", "coverage", "-a", contig_bed, "-b", bam, "-mean", "-sorted"], stdout=fh)
    write_sample_metrics(sample.sample_id, rep_fasta, bam, breadth, depth, metric_tsv,
                         cfg.detection_max_overlap_fraction)
    breadth.unlink(missing_ok=True)
    depth.unlink(missing_ok=True)
    for lib_bam in library_bams:
        lib_bam.unlink(missing_ok=True)
    if not keep_bam:
        bam.unlink(missing_ok=True)
        Path(str(bam) + ".bai").unlink(missing_ok=True)
    return metric_tsv

def run_workflow(group_manifest: str, control_manifest: str, x: int, y: int,
                 outprefix: str, cfg_path: str, resume: bool = False, keep_bam: bool = False) -> None:
    if x < 1:
        raise WorkflowError("-x must be >=1")
    if y < 0:
        raise WorkflowError("-y must be >=0")

    group_manifest_p = Path(group_manifest).resolve()
    control_manifest_p = Path(control_manifest).resolve()
    cfg_path_p = Path(cfg_path).resolve()
    group, group_fields = read_manifest(group_manifest_p)
    control, control_fields = read_manifest(control_manifest_p)
    validate_unique_across(group, control)
    if x > len(group):
        raise WorkflowError(f"-x ({x}) exceeds number of group samples ({len(group)})")
    cfg, resolved_cfg = load_config(cfg_path_p)
    _require_tools()

    prefix = Path(outprefix).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    workdir = Path(str(prefix) + ".work")
    manifest_path = Path(str(prefix) + ".run_manifest.json")
    commands_path = Path(str(prefix) + ".commands.jsonl")
    output_sentinels = [manifest_path, Path(str(prefix) + ".all_contigs.tsv"), workdir]
    sig = _input_signature(group_manifest_p, control_manifest_p, cfg_path_p, group, control, x, y)

    if any(p.exists() for p in output_sentinels):
        if not resume:
            raise WorkflowError(f"Output prefix already exists: {prefix}. Use --resume only for an interrupted identical run.")
        if not manifest_path.is_file():
            raise WorkflowError("Cannot resume: run manifest is missing")
        previous = json.loads(manifest_path.read_text())
        if previous.get("status") == "complete":
            raise WorkflowError("Cannot resume a completed run; use a new outprefix")
        if previous.get("input_signature") != sig:
            raise WorkflowError("Cannot resume: manifests, FASTQ file metadata, configuration, X, or Y changed")
    else:
        workdir.mkdir(parents=True)

    stages = workdir / "stages"
    stages.mkdir(exist_ok=True)
    runner = CommandRunner(commands_path)
    run_manifest = {
        "pops_version": __version__,
        "status": "running",
        "started_unix": time.time(),
        "group_manifest": str(group_manifest_p),
        "control_manifest": str(control_manifest_p),
        "group_samples": [
            {
                "sample_id": s.sample_id,
                "libraries": [
                    {
                        "library_id": lib.library_id,
                        "read1": str(lib.read1),
                        "read2": str(lib.read2) if lib.read2 else None,
                        "metadata": lib.metadata,
                    }
                    for lib in s.libraries
                ],
            }
            for s in group
        ],
        "control_samples": [
            {
                "sample_id": s.sample_id,
                "libraries": [
                    {
                        "library_id": lib.library_id,
                        "read1": str(lib.read1),
                        "read2": str(lib.read2) if lib.read2 else None,
                        "metadata": lib.metadata,
                    }
                    for lib in s.libraries
                ],
            }
            for s in control
        ],
        "manifest_columns": {"group": group_fields, "control": control_fields},
        "parameters": {"x": x, "y": y, "keep_bam": keep_bam},
        "resolved_configuration": resolved_cfg,
        "input_signature": sig,
        "outputs": {},
    }
    manifest_path.write_text(json.dumps(run_manifest, indent=2, default=str) + "\n")

    try:
        versions = Path(str(prefix) + ".software_versions.tsv")
        _write_versions(versions)

        assembly_dir = workdir / "megahit"
        raw_contigs = assembly_dir / "final.contigs.fa"
        if not _stage_done(stages, "assembly"):
            if assembly_dir.exists():
                shutil.rmtree(assembly_dir)
            cmd = ["megahit", *_megahit_inputs(group), "-o", assembly_dir, "-t", str(cfg.threads_megahit)]
            runner.run(cmd)
            if not raw_contigs.is_file():
                raise WorkflowError("MEGAHIT completed but final.contigs.fa was not produced")
            _mark_done(stages, "assembly")

        filtered_fa = workdir / "length_filtered.fasta"
        length_report = Path(str(prefix) + ".contig_length_filter.tsv")
        if not _stage_done(stages, "length_filter"):
            n = filter_fasta_by_length(raw_contigs, filtered_fa, length_report, cfg.min_contig_length)
            if n == 0:
                raise WorkflowError("No contigs remain after the minimum-length filter")
            _mark_done(stages, "length_filter")

        mmseqs_dir = workdir / "mmseqs"
        mmseqs_dir.mkdir(exist_ok=True)
        mm_prefix = mmseqs_dir / "cluster"
        reps = workdir / "representatives.fasta"
        cluster_out = Path(str(prefix) + ".mmseqs_clusters.tsv")
        if not _stage_done(stages, "mmseqs"):
            tmp = mmseqs_dir / "tmp"
            runner.run([
                "mmseqs", "easy-cluster", filtered_fa, mm_prefix, tmp,
                "--min-seq-id", str(cfg.mmseqs_min_seq_id),
                "-c", str(cfg.mmseqs_min_coverage),
                "--cov-mode", str(cfg.mmseqs_cov_mode),
                "--threads", str(cfg.threads_mmseqs),
            ])
            generated_reps = Path(str(mm_prefix) + "_rep_seq.fasta")
            generated_cluster = Path(str(mm_prefix) + "_cluster.tsv")
            if not generated_reps.is_file() or not generated_cluster.is_file():
                raise WorkflowError("MMseqs2 output files were not produced as expected")
            shutil.copy2(generated_reps, reps)
            shutil.copy2(generated_cluster, cluster_out)
            _mark_done(stages, "mmseqs")

        contig_bed = workdir / "representatives.bed"
        write_contig_bed(reps, contig_bed)
        index_dir = workdir / "bowtie2_index"
        index_dir.mkdir(exist_ok=True)
        index_prefix = index_dir / "representatives"
        if not _stage_done(stages, "bowtie_index"):
            runner.run(["bowtie2-build", reps, index_prefix])
            _mark_done(stages, "bowtie_index")

        sample_dir = workdir / "samples"
        sample_dir.mkdir(exist_ok=True)
        all_samples = group + control
        metric_paths: dict[str, Path] = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=cfg.mapping_jobs) as ex:
            futs = {
                ex.submit(_run_mapping, s, reps, index_prefix, contig_bed, sample_dir, cfg, runner, keep_bam): s
                for s in all_samples
            }
            for fut in concurrent.futures.as_completed(futs):
                s = futs[fut]
                metric_paths[s.sample_id] = fut.result()
        _mark_done(stages, "mapping_metrics")

        contig_ids = list(fasta_lengths(reps))
        group_ids = [s.sample_id for s in group]
        control_ids = [s.sample_id for s in control]
        matrices = {
            "group.coverage": (group_ids, "breadth"),
            "group.depth": (group_ids, "mean_depth"),
            "group.detected": (group_ids, "detected"),
            "control.coverage": (control_ids, "breadth"),
            "control.depth": (control_ids, "mean_depth"),
            "control.detected": (control_ids, "detected"),
        }
        matrix_paths: dict[str, Path] = {}
        for stem, (sids, field) in matrices.items():
            out = Path(str(prefix) + f".{stem}.tsv")
            merge_metric_matrix(contig_ids, sids, metric_paths, field, out, cfg.merge_batch_size)
            matrix_paths[stem] = out

        _, gdepth = read_numeric_matrix(matrix_paths["group.depth"])
        _, cdepth = read_numeric_matrix(matrix_paths["control.depth"])
        _, gdet = read_numeric_matrix(matrix_paths["group.detected"])
        _, cdet = read_numeric_matrix(matrix_paths["control.detected"])
        rows = calculate_rankings(fasta_lengths(reps), gdepth, cdepth, gdet, cdet, x, y)
        all_tsv = Path(str(prefix) + ".all_contigs.tsv")
        ranked_tsv = Path(str(prefix) + ".ranked_contigs.tsv")
        ordered_ids = write_ranking_tables(rows, all_tsv, ranked_tsv)
        ranked_fa = Path(str(prefix) + ".ranked_contigs.fasta")
        write_ranked_fasta(reps, ranked_fa, ordered_ids)

        run_manifest["status"] = "complete"
        run_manifest["finished_unix"] = time.time()
        run_manifest["outputs"] = {
            "all_contigs": str(all_tsv),
            "ranked_contigs": str(ranked_tsv),
            "ranked_fasta": str(ranked_fa),
            "mmseqs_clusters": str(cluster_out),
            "length_filter": str(length_report),
            "software_versions": str(versions),
            "commands": str(commands_path),
            **{k: str(v) for k, v in matrix_paths.items()},
        }
        manifest_path.write_text(json.dumps(run_manifest, indent=2, default=str) + "\n")
    except Exception as exc:
        run_manifest["status"] = "failed"
        run_manifest["finished_unix"] = time.time()
        run_manifest["error"] = f"{type(exc).__name__}: {exc}"
        manifest_path.write_text(json.dumps(run_manifest, indent=2, default=str) + "\n")
        raise
