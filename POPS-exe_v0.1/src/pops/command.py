from __future__ import annotations

import json
import shlex
import subprocess
import time
from pathlib import Path
from typing import Sequence


class CommandError(RuntimeError):
    pass


def shell_join(parts: Sequence[str | Path]) -> str:
    return " ".join(shlex.quote(str(x)) for x in parts)


class CommandRunner:
    def __init__(self, log_path: str | Path):
        self.log_path = Path(log_path)
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

    def run(self, command: Sequence[str | Path], *, cwd: str | Path | None = None, stdout=None) -> None:
        cmd = [str(x) for x in command]
        start = time.time()
        proc = subprocess.run(cmd, cwd=cwd, stdout=stdout, stderr=subprocess.PIPE, text=True)
        self._record(cmd, proc.returncode, start, time.time(), proc.stderr)
        if proc.returncode != 0:
            tail = (proc.stderr or "")[-4000:]
            raise CommandError(f"Command failed ({proc.returncode}): {shell_join(cmd)}\n{tail}")

    def run_pipeline(self, command: str, *, cwd: str | Path | None = None) -> None:
        start = time.time()
        proc = subprocess.run(["bash", "-o", "pipefail", "-c", command], cwd=cwd, stderr=subprocess.PIPE, text=True)
        self._record(command, proc.returncode, start, time.time(), proc.stderr)
        if proc.returncode != 0:
            tail = (proc.stderr or "")[-4000:]
            raise CommandError(f"Pipeline failed ({proc.returncode}): {command}\n{tail}")

    def _record(self, command, returncode: int, start: float, end: float, stderr: str | None) -> None:
        record = {
            "command": command,
            "started_unix": start,
            "finished_unix": end,
            "elapsed_seconds": round(end - start, 6),
            "exit_code": returncode,
            "stderr_tail": (stderr or "")[-2000:],
        }
        with self.log_path.open("a") as fh:
            fh.write(json.dumps(record, ensure_ascii=False) + "\n")
