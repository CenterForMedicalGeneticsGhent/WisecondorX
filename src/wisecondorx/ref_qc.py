import json
import logging
from pathlib import Path

import numpy as np

MINREFBINS = 150
OUTLIER_N_SIGMA = 3
CHRY_MIN_USABLE_PCT_WARN = 50.0
CHRY_MIN_USABLE_PCT_FAIL = 20.0


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
    usable = valid & (r >= MINREFBINS)
    n_usable = int(np.sum(usable))
    n_bins = end - start
    usable_pct = 100.0 * n_usable / max(n_bins, 1)
    return {
        "n_bins": n_bins,
        "n_valid": n_valid,
        "mean_of_means": float(np.mean(m[valid])),
        "std_of_means": float(np.std(m[valid])),
        "n_mean_outlier": int(np.sum(m[valid] >= cutoff_outlier)),
        "n_low_refs": int(np.sum(r < MINREFBINS)),
        "n_usable": n_usable,
        "usable_pct": usable_pct,
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
    cutoff_outlier = mean_of_means + OUTLIER_N_SIGMA * float(
        np.std(mean_d[valid])
    )
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


def _legacy_autosomal_prefix_compat_issues(ref):
    issues = []
    if (
        "mask" not in ref
        or "bins_per_chr" not in ref
        or "masked_bins_per_chr_cum" not in ref
    ):
        return issues

    a_bins = np.atleast_1d(ref["bins_per_chr"][...])
    a_cum = np.atleast_1d(ref["masked_bins_per_chr_cum"][...])
    if len(a_bins) < 22 or len(a_cum) < 22:
        return issues

    a_autosomal_span = int(np.sum(a_bins[:22]))
    a_autosomal_masked = int(a_cum[21])
    a_autosomal_prefix = np.atleast_1d(ref["mask"][...])[:a_autosomal_span]

    for suf in [".F", ".M"]:
        bins_key = "bins_per_chr{}".format(suf)
        cum_key = "masked_bins_per_chr_cum{}".format(suf)
        mask_key = "mask{}".format(suf)
        if bins_key not in ref or cum_key not in ref or mask_key not in ref:
            continue

        s_bins = np.atleast_1d(ref[bins_key][...])
        s_cum = np.atleast_1d(ref[cum_key][...])
        if len(s_bins) < 22 or len(s_cum) < 22:
            continue

        s_autosomal_span = int(np.sum(s_bins[:22]))
        s_autosomal_masked = int(s_cum[21])
        if s_autosomal_masked != a_autosomal_masked:
            issues.append(
                "{} autosomal masked bins ({}) differ from A ({})".format(
                    suf, s_autosomal_masked, a_autosomal_masked
                )
            )

        s_autosomal_prefix = np.atleast_1d(ref[mask_key][...])[
            :s_autosomal_span
        ]
        if s_autosomal_span != a_autosomal_span or not np.array_equal(
            s_autosomal_prefix, a_autosomal_prefix
        ):
            issues.append(
                "{} autosomal mask prefix differs from A (legacy predict incompatible)".format(
                    suf
                )
            )
    return issues


def _verdict_f(m):
    if m is None or m.get("n_valid", 0) == 0:
        return "FAIL", "no data"
    if m["n_low_refs"] > 0:
        return "WARN", f"n_refs<{MINREFBINS} in {m['n_low_refs']} bins"
    if m["std_of_means"] > 10:
        return (
            "FAIL",
            f"std(per-bin mean dist) = {m['std_of_means']:.2f} (high)",
        )
    if m["std_of_means"] > 2:
        return "WARN", f"std(per-bin mean dist) = {m['std_of_means']:.2f}"
    if m["outlier_pct"] > 1:
        return "WARN", f"outlier bins = {m['outlier_pct']:.2f}%"
    return "PASS", ""


def _verdict_m(m):
    if m is None or m.get("n_valid", 0) == 0:
        return "FAIL", "no data"

    level = "PASS"
    msg = ""

    def update(new_level, new_msg):
        nonlocal level, msg
        if new_level == "FAIL" or (new_level == "WARN" and level == "PASS"):
            level = new_level
            msg = new_msg

    if m["n_low_refs"] > 0:
        update("WARN", f"n_refs<{MINREFBINS} in {m['n_low_refs']} bins")
    if m["mean_of_means"] > 10:
        update(
            "FAIL",
            f"mean(per-bin mean dist) = {m['mean_of_means']:.2f} (heavy tail)",
        )
    elif m["mean_of_means"] > 2:
        update("WARN", f"mean(per-bin mean dist) = {m['mean_of_means']:.2f}")
    cy = m.get("chrY")
    if cy and cy.get("n_valid", 0) > 0:
        usable_pct = cy.get("usable_pct", np.nan)
        if np.isfinite(usable_pct):
            if usable_pct < CHRY_MIN_USABLE_PCT_FAIL:
                update(
                    "FAIL",
                    f"usable chrY bins = {usable_pct:.1f}% (<{CHRY_MIN_USABLE_PCT_FAIL:.0f}%)",
                )
            elif usable_pct < CHRY_MIN_USABLE_PCT_WARN:
                update(
                    "WARN",
                    f"usable chrY bins = {usable_pct:.1f}% (<{CHRY_MIN_USABLE_PCT_WARN:.0f}%)",
                )
    if (
        cy
        and cy.get("n_valid", 0) > 0
        and np.isfinite(cy.get("mean_of_means", np.nan))
    ):
        ym = cy["mean_of_means"]
        if ym > 100:
            update("FAIL", f"chrY mean distance = {ym:.1f} (very poor chrY)")
        elif ym > 5:
            update("WARN", f"chrY mean distance = {ym:.1f}")
    if m["outlier_pct"] > 1:
        update("WARN", f"outlier bins = {m['outlier_pct']:.2f}%")
    return level, msg


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

    if npz.suffix == ".npz":
        out_json = npz.with_name(npz.stem + "_qc.json")
    else:
        out_json = npz.with_name(npz.name + "_qc.json")

    ref = np.load(npz, encoding="latin1", allow_pickle=True)
    try:
        binsize = int(np.atleast_1d(ref["binsize"])[0])
    except Exception:
        binsize = None

    suffixes = _get_gender_suffixes(ref)
    if not suffixes:
        logging.error(
            "QC failed: no bins_per_chr / bins_per_chr.F / bins_per_chr.M in npz"
        )
        ref.close()
        return 2

    logging.info("Starting ref-QC for file: {}".format(npz))
    if binsize:
        logging.info("Reference binsize: {} bp".format(binsize))
    else:
        logging.info("Reference binsize: (unknown)")

    worst = 0  # 0 pass, 1 warn, 2 fail
    compat_issues = _legacy_autosomal_prefix_compat_issues(ref)

    qc_data = {
        "overall_verdict": "PASS",
        "worst_severity": 0,
        "compat_issues": compat_issues,
        "metrics": {},
    }

    for issue in compat_issues:
        logging.error("[COMPAT] {}".format(issue))
    if compat_issues:
        worst = max(worst, 2)

    for suf in suffixes:
        label = "F" if suf == ".F" else "M" if suf == ".M" else "A"
        m = _compute_metrics(ref, suf)
        if m is None:
            logging.warning(f"[{label}] no indexes/distances — skip")
            continue

        qc_data["metrics"][label] = {"n_bins": m.get("n_bins")}

        if m.get("n_valid", 0) == 0:
            logging.error(f"[{label}] n_bins={m['n_bins']}, n_valid=0 — FAIL")
            worst = max(worst, 2)
            qc_data["metrics"][label]["n_valid"] = 0
            qc_data["metrics"][label]["verdict"] = "FAIL"
            continue

        verdict_fn = _verdict_m if label == "M" else _verdict_f
        verdict, msg = verdict_fn(m)

        qc_data["metrics"][label].update(
            {
                "n_valid": m.get("n_valid", 0),
                "verdict": verdict,
                "message": msg,
                "mean_of_means": m["mean_of_means"],
                "std_of_means": m["std_of_means"],
                "n_mean_outlier": m["n_mean_outlier"],
                "outlier_pct": m["outlier_pct"],
                "n_low_refs": m["n_low_refs"],
            }
        )

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
            qc_data["metrics"][label]["chrY"] = cy
            logger_func(
                f"       chrY: n_bins={cy['n_bins']}, mean={cy['mean_of_means']:.4f}, std={cy['std_of_means']:.4f}, "
                f"outliers={cy['n_mean_outlier']}, n_refs<{MINREFBINS}={cy['n_low_refs']}, "
                f"usable={cy.get('n_usable', 0)} ({cy.get('usable_pct', np.nan):.2f}%)"
            )

        logger_func(f"         -> {verdict}" + (f": {msg}" if msg else ""))

    ref.close()

    if worst == 0:
        logging.info("QC Overall Verdict: PASS")
        qc_data["overall_verdict"] = "PASS"
    elif worst == 1:
        logging.warning("QC Overall Verdict: WARN (review metrics above)")
        qc_data["overall_verdict"] = "WARN"
    else:
        logging.error(
            "QC Overall Verdict: FAIL (ref may cause poor predictions; consider rebuilding or more samples)"
        )
        qc_data["overall_verdict"] = "FAIL"

    qc_data["worst_severity"] = worst

    try:
        with open(out_json, "w") as f:
            json.dump(qc_data, f, indent=4)
        logging.info(f"Wrote structured QC data to {out_json}")
    except Exception as e:
        logging.error(f"Failed to write structured QC data to {out_json}: {e}")

    return worst
