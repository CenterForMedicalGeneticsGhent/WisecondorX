import logging
from pathlib import Path

import numpy as np

MINREFBINS = 150
OUTLIER_N_SIGMA = 3


def _get_gender_suffixes(ref):
    keys = list(ref.keys())
    out = []
    if "bins_per_chr.F" in keys:
        out.append(".F")
    if "bins_per_chr.M" in keys:
        out.append(".M")
    if "bins_per_chr" in keys and not out:
        out.append("")
    return out


def _compute_per_bin_stats(indexes, distances):
    n = len(indexes)
    mean_d = np.zeros(n, dtype=float)
    max_d = np.zeros(n, dtype=float)
    n_refs = np.zeros(n, dtype=int)
    for i in range(n):
        d = np.atleast_1d(distances[i]).ravel()
        idx = np.atleast_1d(indexes[i]).ravel()
        if len(d) == 0:
            mean_d[i] = np.nan
            max_d[i] = np.nan
            n_refs[i] = 0
        else:
            mean_d[i] = np.mean(d)
            max_d[i] = np.max(d)
            n_refs[i] = len(idx)
    return mean_d, max_d, n_refs


def _chrY_metrics(ref, suf, mean_d, n_refs, cutoff_outlier):
    if suf != ".M":
        return None
    key = "masked_bins_per_chr_cum" + suf
    if key not in ref:
        return None
    mbpcc = np.atleast_1d(ref[key][...])
    if len(mbpcc) < 24:
        return None
    start, end = int(mbpcc[22]), int(mbpcc[23])
    if start >= end:
        return {"n_bins": 0}
    m = mean_d[start:end]
    r = n_refs[start:end]
    valid = np.isfinite(m)
    n_valid = int(np.sum(valid))
    if n_valid == 0:
        return {"n_bins": end - start, "n_valid": 0, "mean_of_means": np.nan}
    return {
        "n_bins": end - start,
        "n_valid": n_valid,
        "mean_of_means": float(np.mean(m[valid])),
        "std_of_means": float(np.std(m[valid])),
        "n_mean_outlier": int(np.sum(m[valid] >= cutoff_outlier)),
        "n_low_refs": int(np.sum(r < MINREFBINS)),
    }


def _compute_metrics(ref, suf):
    idx_key = "indexes" + suf
    dist_key = "distances" + suf
    if idx_key not in ref or dist_key not in ref:
        return None
    indexes, distances = ref[idx_key], ref[dist_key]
    n_bins = len(indexes)
    if n_bins == 0:
        return {"n_bins": 0}

    mean_d, max_d, n_refs = _compute_per_bin_stats(indexes, distances)
    valid = np.isfinite(mean_d)
    n_valid = int(np.sum(valid))
    if n_valid == 0:
        return {"n_bins": n_bins, "n_valid": 0}

    mean_of_means = float(np.mean(mean_d[valid]))
    std_of_means = float(np.std(mean_d[valid]))
    cutoff_outlier = mean_of_means + OUTLIER_N_SIGMA * float(np.std(mean_d[valid]))
    n_mean_outlier = int(np.sum(mean_d[valid] >= cutoff_outlier))
    n_low_refs = int(np.sum(n_refs < MINREFBINS))
    outlier_pct = 100.0 * n_mean_outlier / n_valid

    chrY = _chrY_metrics(ref, suf, mean_d, n_refs, cutoff_outlier)
    return {
        "n_bins": n_bins,
        "n_valid": n_valid,
        "mean_of_means": mean_of_means,
        "std_of_means": std_of_means,
        "n_mean_outlier": n_mean_outlier,
        "outlier_pct": outlier_pct,
        "n_low_refs": n_low_refs,
        "chrY": chrY,
    }


def _verdict_f(m):
    if m is None or m.get("n_valid", 0) == 0:
        return "FAIL", "no data"
    if m["n_low_refs"] > 0:
        return "WARN", f"n_refs<{MINREFBINS} in {m['n_low_refs']} bins"
    if m["std_of_means"] > 10:
        return "FAIL", f"std(per-bin mean dist) = {m['std_of_means']:.2f} (high)"
    if m["std_of_means"] > 2:
        return "WARN", f"std(per-bin mean dist) = {m['std_of_means']:.2f}"
    if m["outlier_pct"] > 1:
        return "WARN", f"outlier bins = {m['outlier_pct']:.2f}%"
    return "PASS", ""


def _verdict_m(m):
    if m is None or m.get("n_valid", 0) == 0:
        return "FAIL", "no data"
    if m["n_low_refs"] > 0:
        return "WARN", f"n_refs<{MINREFBINS} in {m['n_low_refs']} bins"
    if m["mean_of_means"] > 10:
        return "FAIL", f"mean(per-bin mean dist) = {m['mean_of_means']:.2f} (heavy tail)"
    if m["mean_of_means"] > 2:
        return "WARN", f"mean(per-bin mean dist) = {m['mean_of_means']:.2f}"
    cy = m.get("chrY")
    if cy and cy.get("n_valid", 0) > 0 and np.isfinite(cy.get("mean_of_means", np.nan)):
        ym = cy["mean_of_means"]
        if ym > 100:
            return "FAIL", f"chrY mean distance = {ym:.1f} (very poor chrY)"
        if ym > 5:
            return "WARN", f"chrY mean distance = {ym:.1f}"
    if m["outlier_pct"] > 1:
        return "WARN", f"outlier bins = {m['outlier_pct']:.2f}%"
    return "PASS", ""


def qc_reference(npz_path: Path):
    """
    Reads the given numpy array representing a WisecondorX reference
    and checks it for common quality issues using predefined heuristics.

    Returns the worst severity code found: 0 (PASS), 1 (WARN), 2 (FAIL).
    """
    npz = Path(npz_path).resolve()
    if not npz.exists():
        logging.error(f"QC check skipped: file not found: {npz}")
        return 2

    ref = np.load(npz, encoding="latin1", allow_pickle=True)
    try:
        binsize = int(np.atleast_1d(ref["binsize"])[0])
    except Exception:
        binsize = None

    suffixes = _get_gender_suffixes(ref)
    if not suffixes:
        logging.error("QC failed: no bins_per_chr / bins_per_chr.F / bins_per_chr.M in npz")
        ref.close()
        return 2

    logging.info("Starting ref-QC for file: {}".format(npz))
    if binsize:
        logging.info("Reference binsize: {} bp".format(binsize))
    else:
        logging.info("Reference binsize: (unknown)")

    worst = 0  # 0 pass, 1 warn, 2 fail
    for suf in suffixes:
        label = "F" if suf == ".F" else "M" if suf == ".M" else "A"
        m = _compute_metrics(ref, suf)
        if m is None:
            logging.warning(f"[{label}] no indexes/distances — skip")
            continue
        if m.get("n_valid", 0) == 0:
            logging.error(f"[{label}] n_bins={m['n_bins']}, n_valid=0 — FAIL")
            worst = max(worst, 2)
            continue

        verdict_fn = _verdict_m if label == "M" else _verdict_f
        verdict, msg = verdict_fn(m)

        if verdict == "FAIL":
            worst = max(worst, 2)
            logger_func = logging.error
        elif verdict == "WARN":
            worst = max(worst, 1)
            logger_func = logging.warning
        else:
            logger_func = logging.info

        logger_func(
            f"[{label}] n_bins={m['n_bins']}, mean(dist)={m['mean_of_means']:.4f}, std(dist)={m['std_of_means']:.4f}, "
            f"outliers={m['n_mean_outlier']} ({m['outlier_pct']:.2f}%), n_refs<{MINREFBINS}={m['n_low_refs']}"
        )
        if m.get("chrY") and m["chrY"].get("n_valid", 0) > 0:
            cy = m["chrY"]
            logger_func(
                f"       chrY: n_bins={cy['n_bins']}, mean={cy['mean_of_means']:.4f}, std={cy['std_of_means']:.4f}, "
                f"outliers={cy['n_mean_outlier']}, n_refs<{MINREFBINS}={cy['n_low_refs']}"
            )

        logger_func(f"         -> {verdict}" + (f": {msg}" if msg else ""))

    ref.close()

    if worst == 0:
        logging.info("QC Overall Verdict: PASS")
    elif worst == 1:
        logging.warning("QC Overall Verdict: WARN (review metrics above)")
    else:
        logging.error(
            "QC Overall Verdict: FAIL (ref may cause poor predictions; consider rebuilding or more samples)"
        )

    return worst
