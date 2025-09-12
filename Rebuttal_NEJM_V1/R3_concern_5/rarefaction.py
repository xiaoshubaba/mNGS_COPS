#!/usr/bin/env python3
# usage:
#   python rarefy_from_matrix.py -in k8.matrix.tsv.gz -out rarefaction.csv \
#     --points 60 --min-frac 0.01 --max-frac 0.95 --ci 0.95
#   # optional parametric bootstrap under independence:
#   # ... --ci-method mc --mc-sims 500 --seed 42
import argparse, os, sys, gzip, math, csv
import numpy as np

def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def parse_header_and_M(path):
    with open_maybe_gzip(path, "rt") as f:
        header = f.readline().rstrip("\n\r")
    if not header:
        raise ValueError("Empty header line.")
    parts = header.split("\t")
    if len(parts) < 2 or parts[0] != "":
        raise ValueError("Header must start with a leading tab then sample names.")
    return parts[1:], len(parts) - 1

def build_row_sum_histogram(path, M):
    n_f = np.zeros(M + 1, dtype=np.int64)
    rows = 0
    with open_maybe_gzip(path, "rt") as f:
        _ = f.readline()
        for line in f:
            line = line.rstrip("\n\r")
            if not line: continue
            parts = line.split("\t")
            vals = parts[1:]
            if len(vals) < M: vals += ["0"] * (M - len(vals))
            elif len(vals) > M: vals = vals[:M]
            s = 0
            for v in vals:
                s += (v != "0")
            n_f[s] += 1
            rows += 1
    if rows == 0: raise ValueError("No data rows found.")
    return n_f, rows

def precompute_lgamma(M):
    L = np.empty(M + 1, dtype=np.float64)
    for n in range(M + 1):
        L[n] = math.lgamma(n + 1.0)
    return L

def logC(L, n, k):
    if k < 0 or k > n: return -np.inf
    return L[n] - L[k] - L[n - k]

def expected_and_var(n_f, M, m_grid):
    f_vals = np.nonzero(n_f[1:])[0] + 1
    counts = n_f[f_vals].astype(np.float64)
    L = precompute_lgamma(M)
    means, vars_ = [], []
    for m in m_grid:
        log_den = logC(L, M, m)
        log_num = np.array([logC(L, M - int(f), m) for f in f_vals], dtype=np.float64)
        absent = np.exp(log_num - log_den)
        p = 1.0 - absent
        mu = float(np.sum(counts * p))
        # independence approximation for variance:
        var = float(np.sum(counts * p * (1.0 - p)))
        means.append(mu)
        vars_.append(var)
    return np.array(means), np.array(vars_)

# Normal quantile (Acklam’s approximation); stable without SciPy
def ndtri(p):
    if not (0.0 < p < 1.0):
        if p == 0.0: return -np.inf
        if p == 1.0: return  np.inf
        raise ValueError("p must be in (0,1)")
    a = (-3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,
          1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00)
    b = (-5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,
          6.680131188771972e+01, -1.328068155288572e+01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00)
    d = ( 7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,
          3.754408661907416e+00)
    plow  = 0.02425
    phigh = 1 - plow
    if p < plow:
        q = math.sqrt(-2*math.log(p))
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    if phigh < p:
        q = math.sqrt(-2*math.log(1-p))
        return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                 ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    q = p - 0.5
    r = q*q
    return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q / \
           (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-in",  dest="inp",  required=True)
    ap.add_argument("-out", dest="out",  required=True)
    ap.add_argument("--points", type=int, default=60)
    ap.add_argument("--min-frac", type=float, default=0.01)
    ap.add_argument("--max-frac", type=float, default=0.95)
    ap.add_argument("--ci", type=float, default=0.95, help="confidence level, e.g. 0.95")
    ap.add_argument("--ci-method", choices=["normal","mc"], default="normal")
    ap.add_argument("--mc-sims", type=int, default=0, help="MC draws if --ci-method mc (suggest 200–1000)")
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()

    label = os.path.splitext(os.path.basename(args.inp))[0]
    sample_names, M = parse_header_and_M(args.inp)

    # grid
    lo = max(0.0, float(args.min_frac)); hi = min(1.0, float(args.max_frac))
    if not (0 < lo < hi <= 1.0): raise ValueError("Require 0 < --min-frac < --max-frac ≤ 1")
    m_grid = np.linspace(lo, hi, args.points) * M
    m_grid = np.clip(np.rint(m_grid), 1, M).astype(int)
    m_grid = np.unique(m_grid)

    # data
    n_f, total_fragments = build_row_sum_histogram(args.inp, M)
    observed_fragments = int(n_f[1:].sum())

    means, vars_ = expected_and_var(n_f, M, m_grid)

    # CIs
    ci_lo = np.empty_like(means); ci_hi = np.empty_like(means)
    if args.ci_method == "normal":
        z = ndtri(0.5 + 0.5*args.ci)
        sd = np.sqrt(np.maximum(vars_, 0.0))
        ci_lo = np.maximum(0.0, means - z*sd)
        ci_hi = np.minimum(observed_fragments, means + z*sd)
        ci_note = "normal"
    else:
        rng = np.random.default_rng(args.seed)
        ci_alpha = (1 - args.ci) / 2.0
        ci_lo = np.empty_like(means); ci_hi = np.empty_like(means)
        # precompute unique f, counts:
        f_vals = np.nonzero(n_f[1:])[0] + 1
        counts = n_f[f_vals]
        L = precompute_lgamma(M)
        for i, m in enumerate(m_grid):
            log_den = logC(L, M, m)
            log_num = np.array([logC(L, M - int(f), m) for f in f_vals], dtype=np.float64)
            p = 1.0 - np.exp(log_num - log_den)
            # parametric bootstrap under independence: S = sum Binomial(n_f, p_f)
            draws = rng.binomial(counts, p, size=(max(1, args.mc_sims), counts.size)).sum(axis=1)
            ci_lo[i] = np.quantile(draws, ci_alpha)
            ci_hi[i] = np.quantile(draws, 1 - ci_alpha)
        ci_note = f"mc{args.mc_sims}"

    write_header = not os.path.exists(args.out) or os.path.getsize(args.out) == 0
    with open(args.out, "a", newline="") as fh:
        w = csv.writer(fh)
        if write_header:
            w.writerow(["label","n_samples","frac_of_samples","expected_uniques",
                        "ci_low","ci_high","ci_method",
                        "total_samples","total_fragments","observed_fragments"])
        for m, mu, lo, hi in zip(m_grid, means, ci_lo, ci_hi):
            w.writerow([label, int(m), float(m)/M, float(mu),
                        float(lo), float(hi), ci_note,
                        M, total_fragments, observed_fragments])

    print(f"[{label}] samples={M}, fragments={total_fragments} (nonzero={observed_fragments}); "
          f"wrote {len(m_grid)} points with {args.ci}-CI ({ci_note}) to {args.out}")

if __name__ == "__main__":
    main()

