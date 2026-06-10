import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
import re
import pandas as pd

# Color scheme
BLACK = "#3f3f3f"
LIGHTER_GREY = "#e0e0e0"

COLOR_A = (84 / 255, 84 / 255, 84 / 255)
COLOR_B = (227 / 255, 200 / 255, 138 / 255)
COLOR_C = (141 / 255, 209 / 255, 198 / 255)
COLOR_D = (150 / 255, 80 / 255, 33 / 255)

COLOR_AA = (84 / 255, 84 / 255, 84 / 255, 80 / 255)
COLOR_BB = (227 / 255, 200 / 255, 138 / 255, 80 / 255)
COLOR_CC = (141 / 255, 209 / 255, 198 / 255, 80 / 255)


def get_aberration_cutoff(beta, ploidy):
    loss_cutoff = np.log2((ploidy - (beta / 2)) / ploidy)
    gain_cutoff = np.log2((ploidy + (beta / 2)) / ploidy)
    return loss_cutoff, gain_cutoff


def format_n_reads(n_reads):
    s = str(int(n_reads))
    parts = []
    while s:
        parts.append(s[-3:])
        s = s[:-3]
    return ".".join(reversed(parts))


def plot_constitutionals(ax, ploidy, start, end):
    ax.hlines(
        np.log2(1 / ploidy),
        start,
        end,
        colors=[COLOR_B],
        linestyles="dotted",
        linewidth=2,
    )
    ax.hlines(
        np.log2(2 / ploidy),
        start,
        end,
        colors=[COLOR_A],
        linestyles="dotted",
        linewidth=2,
    )
    ax.hlines(
        np.log2(3 / ploidy),
        start,
        end,
        colors=[COLOR_C],
        linestyles="dotted",
        linewidth=2,
    )


def get_boxplot_stats(data):
    if len(data) == 0:
        return [np.nan, np.nan]
    q1, q3 = np.nanpercentile(data, [25, 75])
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    # Clip to actual data range as R's boxplot does
    lower = max(lower, np.nanmin(data))
    upper = min(upper, np.nanmax(data))
    return [lower, upper]


def _normalize_chrom(chrom):
    chrom_str = str(chrom).strip()
    chrom_str = (
        chrom_str[3:] if chrom_str.lower().startswith("chr") else chrom_str
    )
    chrom_upper = chrom_str.upper()

    if chrom_upper == "X":
        return 23
    if chrom_upper == "Y":
        return 24

    try:
        return int(chrom_str)
    except ValueError:
        return chrom_str


def _chrom_label(chrom):
    normalized = _normalize_chrom(chrom)
    if normalized == 23:
        return "chrX"
    if normalized == 24:
        return "chrY"
    return f"chr{normalized}"


def _chrom_sort_key(chrom):
    normalized = _normalize_chrom(chrom)
    if isinstance(normalized, int):
        return (0, normalized)
    return (1, str(normalized))


def _safe_numeric(series):
    return pd.to_numeric(series, errors="coerce").to_numpy(dtype=float).copy()


def write_plots(data):
    plt.switch_backend("Agg")
    out_dir = data["out_dir"]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    gender = data["ref_gender"]
    beta = float(data["beta"]) if data["beta"] != "None" else None
    zcutoff = float(data["zscore"])
    binsize = int(data["binsize"])
    n_reads_raw = int(data["n_reads"])
    n_reads_formatted = format_n_reads(n_reads_raw)
    ylim_arg = data["ylim"]
    regions_path = data["regions"]
    plot_title = data.get("plot_title")

    cna_df = data["cna_df"].copy()
    segments_df = data["segments_df"].copy()

    cna_df = cna_df.sort_values(by=["chrom", "maploc"])
    chrom_keys = sorted(
        cna_df["chrom"].drop_duplicates().tolist(), key=_chrom_sort_key
    )

    chrom_offsets = {}
    chrom_lengths = {}
    chrom_labels = []
    ratio_parts = []
    weight_parts = []
    box_list = []
    l_whis_per_chr = []
    h_whis_per_chr = []

    running_offset = 0
    for chrom in chrom_keys:
        chrom_df = cna_df[cna_df["chrom"] == chrom].sort_values(by="maploc")
        chrom_offsets[chrom] = running_offset
        chrom_lengths[chrom] = len(chrom_df)
        chrom_labels.append(_chrom_label(chrom))

        chrom_ratio = _safe_numeric(chrom_df["ratios"])
        chrom_ratio[chrom_ratio == 0] = np.nan
        chrom_weights = _safe_numeric(chrom_df["weights"])
        chrom_weights[chrom_weights == 0] = np.nan

        ratio_parts.append(chrom_ratio)
        weight_parts.append(chrom_weights)
        box_list.append(chrom_ratio[~np.isnan(chrom_ratio)])

        stats = get_boxplot_stats(chrom_ratio)
        l_whis_per_chr.append(stats[0])
        h_whis_per_chr.append(stats[1])

        running_offset += len(chrom_df)

    ratio = np.concatenate(ratio_parts) if ratio_parts else np.array([])
    weights = np.concatenate(weight_parts) if weight_parts else np.array([])
    num_chrs = len(chrom_keys)
    chr_ends = [0] + np.cumsum(
        [chrom_lengths[chrom] for chrom in chrom_keys]
    ).tolist()
    chr_mids = [
        chr_ends[i] + chrom_lengths[chrom_keys[i]] / 2 for i in range(num_chrs)
    ]
    labels = chrom_labels

    # Boxplot data per chromosome
    chr_wide_upper_limit = (
        max(0.65, np.nanmax(h_whis_per_chr) if h_whis_per_chr else 0.65) * 1.25
    )
    chr_wide_lower_limit = (
        min(-0.95, np.nanmin(l_whis_per_chr) if l_whis_per_chr else -0.95)
        * 1.25
    )

    if ylim_arg != "def":
        ylim_parts = re.findall(r"[-+]?\d*\.\d+|\d+", ylim_arg)
        if len(ylim_parts) == 2:
            chr_wide_lower_limit = float(ylim_parts[0])
            chr_wide_upper_limit = float(ylim_parts[1])

    segment_starts = _safe_numeric(segments_df["start"])
    segment_ends = _safe_numeric(segments_df["end"])
    segment_means = _safe_numeric(segments_df["ratio"])
    segment_zscores = _safe_numeric(segments_df["zscore"])

    # Determine dot colors
    dot_cols: list[tuple[float, float, float] | str] = [COLOR_A] * len(ratio)
    for idx, segment in segments_df.reset_index(drop=True).iterrows():
        chrom_key = _normalize_chrom(segment["chrom"])
        if chrom_key not in chrom_offsets:
            continue

        chr_offset = chrom_offsets[chrom_key]
        chr_len = chrom_lengths[chrom_key]
        start_idx = chr_offset + int(np.floor(segment_starts[idx] / binsize))
        end_idx = chr_offset + max(
            int(np.ceil(segment_ends[idx] / binsize)) - 1,
            int(np.floor(segment_starts[idx] / binsize)),
        )
        end_idx = min(end_idx, chr_offset + chr_len - 1)
        z_val = segment_zscores[idx]
        height = segment_means[idx]
        ploidy = 1 if chrom_key in [23, 24] and gender == "M" else 2

        color = COLOR_A
        if beta is not None:
            l_cut, g_cut = get_aberration_cutoff(beta, ploidy)
            if height < l_cut:
                color = COLOR_B
            elif height > g_cut:
                color = COLOR_C
        else:
            if not np.isnan(z_val):
                if z_val < -zcutoff:
                    color = COLOR_B
                elif z_val > zcutoff:
                    color = COLOR_C
            else:
                color = "grey"

        for dot_idx in range(start_idx, min(end_idx + 1, len(dot_cols))):
            dot_cols[dot_idx] = color

    # Parse regions
    gene_labels = []
    if regions_path and os.path.exists(regions_path):
        with open(regions_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                chrom_key = _normalize_chrom(parts[0])
                if chrom_key not in chrom_offsets:
                    continue

                try:
                    r_start = int(parts[1])
                    r_end = int(parts[2])
                except ValueError:
                    continue

                start_bin = r_start // binsize + chrom_offsets[chrom_key]
                end_bin = r_end // binsize + chrom_offsets[chrom_key]

                if start_bin >= end_bin:
                    continue

                reg_ratio = np.nanmean(ratio[start_bin : end_bin + 1])
                if reg_ratio > 0:
                    label_pos = np.nanmax(ratio[start_bin : end_bin + 1]) + 0.2
                    label_adj = 0  # top
                else:
                    label_pos = np.nanmin(ratio[start_bin : end_bin + 1]) - 0.2
                    label_adj = 1  # bottom

                gene_labels.append(
                    {
                        "start_bin": start_bin,
                        "end_bin": end_bin,
                        "label": parts[3],
                        "pos": label_pos,
                        "adj": label_adj,
                    }
                )

    # Genome-wide plot
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(17, 25)
    ax_main = fig.add_subplot(gs[:10, :])
    ax_box_auto = fig.add_subplot(gs[10:, :22])
    ax_box_sex = fig.add_subplot(gs[10:, 22:])

    # Main plot
    ax_main.set_facecolor("white")
    ax_main.set_xlim(0, chr_ends[num_chrs])
    ax_main.set_ylim(chr_wide_lower_limit, chr_wide_upper_limit)
    ax_main.set_ylabel(r"$\log_2(\text{ratio})$")

    # Constitutional lines
    if gender == "F":
        plot_constitutionals(ax_main, 2, 0, chr_ends[num_chrs])
    else:
        plot_constitutionals(ax_main, 2, 0, chr_ends[22])
        plot_constitutionals(ax_main, 1, chr_ends[22], chr_ends[num_chrs])

    # Masked regions
    masked_indices = np.where(np.isnan(ratio))[0]
    for idx in masked_indices:
        ax_main.axvline(idx, color=LIGHTER_GREY, linewidth=0.1, alpha=0.5)

    # Scatter
    dot_sizes = (weights / np.pi) ** 0.5 * 20  # Adjusted multiplier
    ax_main.scatter(
        np.arange(len(ratio)),
        ratio,
        c=dot_cols,
        s=dot_sizes,
        edgecolors="none",
        alpha=0.8,
    )

    # Gene labels
    for gl in gene_labels:
        ax_main.scatter(
            range(gl["start_bin"], gl["end_bin"] + 1),
            ratio[gl["start_bin"] : gl["end_bin"] + 1],
            facecolors="none",
            edgecolors=COLOR_D,
            s=dot_sizes[gl["start_bin"] : gl["end_bin"] + 1] * 1.5,
            linewidths=1.5,
        )
        ax_main.text(
            (gl["start_bin"] + gl["end_bin"]) / 2,
            gl["pos"],
            gl["label"],
            color=COLOR_D,
            rotation=90,
            va="bottom" if gl["adj"] == 0 else "top",
            ha="center",
            fontsize=8,
        )

    # Segments
    for idx, segment in segments_df.reset_index(drop=True).iterrows():
        chrom_key = _normalize_chrom(segment["chrom"])
        if chrom_key not in chrom_offsets:
            continue

        start = chrom_offsets[chrom_key] + segment_starts[idx] / binsize
        end = chrom_offsets[chrom_key] + segment_ends[idx] / binsize
        height = segment_means[idx]
        color = dot_cols[int(np.floor(start))]

        # In R it used AA, BB, CC with alpha. Here we use the base color with alpha
        face_color = list(color) + [0.3] if isinstance(color, tuple) else color
        ax_main.add_patch(
            Rectangle(
                (start, 0),
                end - start,
                height,
                facecolor=face_color,
                edgecolor=face_color,
                linewidth=0.1,
            )
        )
        # Draw the segment line
        ax_main.hlines(height, start, end, colors=[LIGHTER_GREY], linewidth=2)

    # X-axis ticks
    ax_main.set_xticks(chr_mids)
    ax_main.set_xticklabels(labels, rotation=45)
    for end in chr_ends:
        ax_main.axvline(end, color=BLACK, linestyle="dotted", linewidth=1)

    # Legend
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D(
            [0],
            [0],
            color=COLOR_C,
            linestyle="dotted",
            lw=2,
            label="Constitutional 3n",
        ),
        Line2D(
            [0],
            [0],
            color=COLOR_A,
            linestyle="dotted",
            lw=2,
            label="Constitutional 2n",
        ),
        Line2D(
            [0],
            [0],
            color=COLOR_B,
            linestyle="dotted",
            lw=2,
            label="Constitutional 1n",
        ),
    ]
    ax_main.legend(
        handles=legend_elements,
        loc="upper right",
        bbox_to_anchor=(1, 1.25),
        frameon=False,
        ncol=3,
    )

    # Info legend
    info_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Gain",
            markerfacecolor=COLOR_C,
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Loss",
            markerfacecolor=COLOR_B,
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="",
            color="w",
            label=f"Number of reads: {n_reads_formatted}",
        ),
    ]
    ax_main.add_artist(
        ax_main.legend(
            handles=info_elements,
            loc="upper left",
            bbox_to_anchor=(0, 1.25),
            frameon=False,
            ncol=3,
        )
    )

    if plot_title:
        fig.suptitle(plot_title, color=COLOR_A, fontsize=16)

    # Boxplots
    ax_box_auto.boxplot(
        box_list[:22],
        tick_labels=labels[:22],
        flierprops=dict(marker="o", markersize=2, markerfacecolor=BLACK),
    )
    ax_box_auto.set_ylim(
        min(l_whis_per_chr[:22]) if l_whis_per_chr[:22] else -1,
        max(h_whis_per_chr[:22]) if h_whis_per_chr[:22] else 1,
    )
    ax_box_auto.set_ylabel(r"$\log_2(\text{ratio})$")
    ax_box_auto.tick_params(axis="x", rotation=45)
    plot_constitutionals(ax_box_auto, 2, 0.5, 22.5)

    if num_chrs > 22:
        ax_box_sex.boxplot(
            box_list[22:],
            tick_labels=labels[22:],
            flierprops=dict(marker="o", markersize=2, markerfacecolor=BLACK),
        )
        y_sex_down = min(l_whis_per_chr[22:]) if l_whis_per_chr[22:] else -1
        y_sex_up = max(h_whis_per_chr[22:]) if h_whis_per_chr[22:] else 1
        ax_box_sex.set_ylim(y_sex_down, y_sex_up)
        ax_box_sex.tick_params(axis="x", rotation=45)
        if gender == "F":
            plot_constitutionals(ax_box_sex, 2, 0.5, len(labels[22:]) + 0.5)
        else:
            plot_constitutionals(ax_box_sex, 1, 0.5, len(labels[22:]) + 0.5)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "genome_wide.png"), dpi=300)
    plt.close()

    # Chromosome-specific plots
    for c, chrom in enumerate(chrom_keys):
        fig_chr = plt.figure(figsize=(10, 6))
        ax_chr = fig_chr.add_subplot(111)

        chr_start = chr_ends[c]
        chr_end = chr_ends[c + 1]

        chr_ratio = ratio[chr_start:chr_end]
        chr_weights = weights[chr_start:chr_end]
        chr_dot_cols = dot_cols[chr_start:chr_end]
        chr_dot_sizes = (chr_weights / np.pi) ** 0.5 * 30

        stats = get_boxplot_stats(chr_ratio)
        if np.isnan(stats[0]):
            plt.close()
            continue

        upper_limit = 0.6 + stats[1]
        lower_limit = -1.05 + stats[0]
        upper_limit = max(upper_limit, np.nanmax(chr_ratio))
        lower_limit = min(lower_limit, np.nanmin(chr_ratio))

        if ylim_arg != "def":
            ylim_parts = re.findall(r"[-+]?\d*\.\d+|\d+", ylim_arg)
            if len(ylim_parts) == 2:
                lower_limit = float(ylim_parts[0])
                upper_limit = float(ylim_parts[1])

        ax_chr.set_xlim(chr_start, chr_end)
        ax_chr.set_ylim(lower_limit, upper_limit)
        ax_chr.set_title(labels[c])
        ax_chr.set_ylabel(r"$\log_2(\text{ratio})$")

        if gender == "F":
            plot_constitutionals(ax_chr, 2, chr_start, chr_end)
        else:
            if c >= 22:
                plot_constitutionals(ax_chr, 1, chr_start, chr_end)
            else:
                plot_constitutionals(ax_chr, 2, chr_start, chr_end)

        # Masked regions
        chr_masked = np.where(np.isnan(chr_ratio))[0]
        for idx in chr_masked:
            ax_chr.axvline(
                chr_start + idx, color=COLOR_A, linewidth=1, alpha=0.3
            )

        # Scatter
        ax_chr.scatter(
            range(chr_start, chr_end),
            chr_ratio,
            c=chr_dot_cols,
            s=chr_dot_sizes,
            edgecolors="none",
            alpha=0.8,
        )

        # Segments
        for idx, segment in segments_df.reset_index(drop=True).iterrows():
            if _normalize_chrom(segment["chrom"]) != _normalize_chrom(chrom):
                continue
            seg_start = chr_start + segment_starts[idx] / binsize
            seg_end = chr_start + segment_ends[idx] / binsize
            height = segment_means[idx]
            color = dot_cols[int(np.floor(seg_start))]
            face_color = (
                list(color) + [0.3] if isinstance(color, tuple) else color
            )
            ax_chr.add_patch(
                Rectangle(
                    (seg_start, 0),
                    seg_end - seg_start,
                    height,
                    facecolor=face_color,
                    edgecolor=face_color,
                    linewidth=0.1,
                )
            )
            ax_chr.hlines(
                height, seg_start, seg_end, colors=[LIGHTER_GREY], linewidth=3
            )

        # Gene labels
        for gl in gene_labels:
            if gl["start_bin"] >= chr_start and gl["end_bin"] <= chr_end:
                ax_chr.scatter(
                    range(gl["start_bin"], gl["end_bin"] + 1),
                    ratio[gl["start_bin"] : gl["end_bin"] + 1],
                    facecolors="none",
                    edgecolors=COLOR_D,
                    s=chr_dot_sizes[
                        gl["start_bin"] - chr_start : gl["end_bin"]
                        - chr_start
                        + 1
                    ]
                    * 1.5,
                    linewidths=1.5,
                )
                ax_chr.text(
                    (gl["start_bin"] + gl["end_bin"]) / 2,
                    gl["pos"],
                    gl["label"],
                    color=COLOR_D,
                    rotation=90,
                    va="bottom" if gl["adj"] == 0 else "top",
                    ha="center",
                    fontsize=9,
                )

        # X-axis (MB)
        x_ticks = np.linspace(chr_start, chr_end, 11)
        x_labels = [
            f"{((t - chr_start) * binsize) / 1e6:.1f}" for t in x_ticks
        ]
        ax_chr.set_xticks(x_ticks)
        ax_chr.set_xticklabels(x_labels)
        ax_chr.set_xlabel("Position (MB)")

        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{labels[c]}.png"), dpi=200)
        plt.close()
