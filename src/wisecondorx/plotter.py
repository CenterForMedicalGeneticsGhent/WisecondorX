import os
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

# Define colors
COLOR_BLACK = "#3f3f3f"
COLOR_LIGHTER_GREY = "#e0e0e0"
COLOR_A = "#545454"
COLOR_B = "#E3C88A"
COLOR_C = "#8DD1C6"
COLOR_D = "#965021"
COLOR_AA = "#545454CC"  # Alpha ~80/255 -> 0.31 -> maybe just use native matplotlib alpha instead of hex
COLOR_BB = "#E3C88ACC"
COLOR_CC = "#8DD1C6CC"


def get_aberration_cutoff(beta: float, ploidy: int) -> tuple[float, float]:
    loss_cutoff = math.log2((ploidy - (beta / 2)) / ploidy)
    gain_cutoff = math.log2((ploidy + (beta / 2)) / ploidy)
    return loss_cutoff, gain_cutoff


def format_n_reads(n_reads: int) -> str:
    # formatted with dots thousands separator
    return f"{n_reads:,}".replace(",", ".")


def create_plots(
    out_dir: str,
    ref_gender: str,
    beta: float,
    zscore: float,
    binsize: int,
    n_reads: int,
    cairo_flag: bool,
    ylim_str: str,
    regions_file: str,
    plot_title: str,
    results_r: list[list[float]],
    results_w: list[list[float]],
    results_c: list[list[float]],
):
    if cairo_flag:
        matplotlib.use("cairo")
    else:
        matplotlib.use("Agg")

    os.makedirs(out_dir, exist_ok=True)

    n_reads_str = format_n_reads(n_reads)

    # Flatten ratios and weights (R code treats them as flat arrays + NAs for 0s)
    # R script: ratio <- unlist(input$results_r); ratio[which(ratio == 0)] <- NA
    ratios = []
    weights = []

    chrs = list(range(1, 25)) if ref_gender == "M" else list(range(1, 24))

    bins_per_chr = []
    for c in chrs:
        r_c_arr = np.array(results_r[c - 1])
        r_c_arr[r_c_arr == 0] = np.nan
        r_c_arr[np.isinf(r_c_arr)] = np.nan

        w_c_arr = np.array(results_w[c - 1])
        w_c_arr[w_c_arr == 0] = np.nan
        w_c_arr[np.isinf(w_c_arr)] = np.nan

        ratios.extend(r_c_arr.tolist())
        weights.extend(w_c_arr.tolist())
        bins_per_chr.append(len(r_c_arr))

    ratios = np.array(ratios)
    weights = np.array(weights)

    labels = [f"chr{c}" for c in chrs]
    labels = ["chrX" if label == "chr23" else label for label in labels]
    labels = ["chrY" if label == "chr24" else label for label in labels]

    chr_ends = np.concatenate(([0], np.cumsum(bins_per_chr)))
    chr_mids = chr_ends[1:] - np.array(bins_per_chr) / 2

    # Get margins and whiskers
    box_list = []
    l_whis_per_chr = []
    h_whis_per_chr = []

    for i, c in enumerate(chrs):
        start = chr_ends[i]
        end = chr_ends[i + 1]
        data = ratios[start:end]
        data_clean = data[~np.isnan(data)]
        box_list.append(data_clean)

        if len(data_clean) > 0:
            q1, q3 = np.percentile(data_clean, [25, 75])
            iqr = q3 - q1
            l_whis = (
                data_clean[data_clean >= q1 - 1.5 * iqr].min()
                if len(data_clean[data_clean >= q1 - 1.5 * iqr]) > 0
                else q1
            )
            h_whis = (
                data_clean[data_clean <= q3 + 1.5 * iqr].max()
                if len(data_clean[data_clean <= q3 + 1.5 * iqr]) > 0
                else q3
            )
            l_whis_per_chr.append(l_whis)
            h_whis_per_chr.append(h_whis)
        else:
            l_whis_per_chr.append(np.nan)
            h_whis_per_chr.append(np.nan)

    # Genome-wide limits
    if ylim_str == "def":
        max_h_whis = (
            np.nanmax(h_whis_per_chr)
            if not np.all(np.isnan(h_whis_per_chr))
            else 0
        )
        min_l_whis = (
            np.nanmin(l_whis_per_chr)
            if not np.all(np.isnan(l_whis_per_chr))
            else 0
        )
        chr_wide_upper_limit = max(0.65, max_h_whis) * 1.25
        chr_wide_lower_limit = min(-0.95, min_l_whis) * 1.25
    else:
        # e.g., '[-2,2]'
        cleaned_ylim = ylim_str.replace("[", "").replace("]", "")
        parts = cleaned_ylim.split(",")
        chr_wide_lower_limit = float(parts[0])
        chr_wide_upper_limit = float(parts[1])

    dot_sizes = (
        (weights / np.pi) ** 0.5 * 0.8 * 10
    )  # adjust multiplier for matplotlib scatter size

    # Default colors
    dot_cols = np.array([COLOR_A] * len(ratios), dtype=object)

    for ab in results_c:
        c_idx = int(ab[0]) + 1
        start = int(ab[1]) + chr_ends[c_idx - 1]
        end = int(ab[2]) + chr_ends[c_idx - 1]
        z = float(ab[3]) if ab[3] is not None else np.nan
        height = float(ab[4])

        ploidy = 2
        if (c_idx == 23 or c_idx == 24) and ref_gender == "M":
            ploidy = 1

        if beta is not None:
            l_cut, g_cut = get_aberration_cutoff(beta, ploidy)
            if height < l_cut:
                dot_cols[start:end] = COLOR_B
            if height > g_cut:
                dot_cols[start:end] = COLOR_C
        else:
            if np.isnan(z):
                dot_cols[start:end] = "grey"
                continue
            if z < -zscore:
                dot_cols[start:end] = COLOR_B
            if z > zscore:
                dot_cols[start:end] = COLOR_C

    # Read regions
    gene_labels = []
    if regions_file and os.path.exists(regions_file):
        with open(regions_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    chr_name = parts[0].replace("chr", "")
                    if chr_name == "X":
                        chr_name = "23"
                    if chr_name == "Y":
                        chr_name = "24"
                    try:
                        c_idx = int(chr_name)
                        if c_idx < 1 or c_idx > 24:
                            continue
                        r_start = int(parts[1])
                        r_end = int(parts[2])
                        start_bin = (
                            math.ceil(r_start / binsize) + chr_ends[c_idx - 1]
                        )
                        end_bin = (
                            math.ceil(r_end / binsize) + chr_ends[c_idx - 1]
                        )

                        region_ratio = np.nanmean(ratios[start_bin:end_bin])
                        if region_ratio > 0:
                            label_pos = (
                                np.nanmax(ratios[start_bin:end_bin]) + 0.2
                            )
                            label_adj = 0  # bottom aligned
                        else:
                            label_pos = (
                                np.nanmin(ratios[start_bin:end_bin]) - 0.2
                            )
                            label_adj = 1  # top aligned

                        gene_labels.append(
                            (
                                start_bin,
                                end_bin,
                                parts[3],
                                label_pos,
                                label_adj,
                            )
                        )
                    except ValueError:
                        continue

    # ------------------
    # Genome-wide plot
    # ------------------
    plt.ioff()
    fig = plt.figure(figsize=(14, 10), dpi=256)
    gs = gridspec.GridSpec(10, 1, height_ratios=[7] + [1] * 9)
    # The R script had layout of 10x25 where top plot takes 7 rows,
    # bottom takes 3 rows.
    gs = gridspec.GridSpec(10, 1)
    ax_main = fig.add_subplot(gs[0:7, 0])
    ax_box = fig.add_subplot(gs[7:10, 0], sharex=ax_main)

    # Plot components for main ax...

    # Plot constitutionals
    def plot_constitutionals(ax, ploidy, start_x, end_x):
        ax.plot(
            [start_x, end_x],
            [math.log2(1 / ploidy), math.log2(1 / ploidy)],
            color=COLOR_B,
            lw=2,
            ls=":",
        )
        ax.plot(
            [start_x, end_x],
            [math.log2(2 / ploidy), math.log2(2 / ploidy)],
            color=COLOR_A,
            lw=2,
            ls=":",
        )
        ax.plot(
            [start_x, end_x],
            [math.log2(3 / ploidy), math.log2(3 / ploidy)],
            color=COLOR_C,
            lw=2,
            ls=":",
        )

    genome_len = chr_ends[-1]
    autosome_len = chr_ends[22] if len(chr_ends) > 22 else genome_len

    if ref_gender == "F":
        plot_constitutionals(
            ax_main, 2, -genome_len * 0.025, genome_len * 1.025
        )
    else:
        plot_constitutionals(ax_main, 2, -genome_len * 0.025, autosome_len)
        plot_constitutionals(ax_main, 1, autosome_len, genome_len * 1.025)

    # Missing regions lines
    nas = np.where(np.isnan(ratios))[0]
    for i in nas:
        ax_main.axvline(i, color=COLOR_LIGHTER_GREY, lw=0.1)

    # Scatter
    x_coords = np.arange(len(ratios))
    mask = ~np.isnan(ratios)
    ax_main.scatter(
        x_coords[mask],
        ratios[mask],
        c=dot_cols[mask],
        s=dot_sizes[mask],
        zorder=3,
        edgecolors="none",
        alpha=0.9,
    )

    # Gene labels highlight
    for start_bin, end_bin, label, label_pos, label_adj in gene_labels:
        ax_main.scatter(
            x_coords[start_bin:end_bin],
            ratios[start_bin:end_bin],
            facecolors="none",
            edgecolors=COLOR_D,
            s=dot_sizes[start_bin:end_bin] * 1.5,
            lw=1.5,
            zorder=4,
        )
        va = "bottom" if label_adj == 0 else "top"
        ax_main.text(
            start_bin + (end_bin - start_bin) / 2,
            label_pos,
            label,
            color=COLOR_D,
            fontsize=8,
            rotation=90,
            va=va,
            ha="center",
        )

    # CBS segments
    for ab in results_c:
        c_idx = int(ab[0]) + 1
        start = int(ab[1]) + chr_ends[c_idx - 1]
        end = int(ab[2]) + chr_ends[c_idx - 1]
        height = float(ab[4])

        c = dot_cols[start]
        fill_alpha = 0.3
        ax_main.add_patch(
            patches.Rectangle(
                (start, 0),
                end - start,
                height,
                linewidth=0.1,
                edgecolor=c,
                facecolor=c,
                alpha=fill_alpha,
                zorder=1,
            )
        )

        # horizontal segment
        lw = (
            5 * np.nanmean(dot_sizes[start:end]) / 10
            if not np.isnan(np.nanmean(dot_sizes[start:end]))
            else 2
        )
        ax_main.plot(
            [start, end],
            [height, height],
            color=COLOR_LIGHTER_GREY,
            lw=lw,
            zorder=2,
        )

    ax_main.set_xlim(0, genome_len)
    ax_main.set_ylim(chr_wide_lower_limit, chr_wide_upper_limit)
    ax_main.set_ylabel(r"$\log_2$(ratio)", fontsize=14)
    if plot_title:
        ax_main.set_title(plot_title, loc="right", color=COLOR_A, fontsize=14)

    ax_main.set_xticks(chr_mids)
    ax_main.set_xticklabels(labels, rotation=45, ha="right", fontsize=12)
    ax_main.spines["top"].set_visible(False)
    ax_main.spines["right"].set_visible(False)
    ax_main.spines["bottom"].set_visible(False)
    ax_main.tick_params(axis="x", length=0)

    # Draw chromosome dividers
    for x in chr_ends:
        ax_main.axvline(x, color=COLOR_BLACK, lw=1.2, ls=":")

    # Boxplots
    ax_box.boxplot(
        box_list[:22],
        positions=chr_mids[:22],
        widths=np.diff(chr_ends[:23]) * 0.6,
        patch_artist=True,
        boxprops=dict(facecolor=COLOR_BLACK),
        medianprops=dict(color="white"),
        flierprops=dict(
            marker="o",
            markerfacecolor=COLOR_BLACK,
            markersize=3,
            markeredgecolor="none",
        ),
        showfliers=False,
    )  # R's outpch=16 means it shows outliers but we can just let matplotlib do it. Let's do showfliers=True

    if len(chrs) > 22:
        ax_box.boxplot(
            box_list[22:],
            positions=chr_mids[22:],
            widths=np.diff(chr_ends[22:]) * 0.6,
            patch_artist=True,
            boxprops=dict(facecolor=COLOR_BLACK),
            medianprops=dict(color="white"),
            flierprops=dict(
                marker="o",
                markerfacecolor=COLOR_BLACK,
                markersize=3,
                markeredgecolor="none",
            ),
        )

    ax_box.set_xlim(0, genome_len)
    ax_box.get_ylim()
    # constitutionals for boxplots
    plot_constitutionals(ax_box, 2, 0, autosome_len)
    if ref_gender == "F":
        plot_constitutionals(ax_box, 2, autosome_len, genome_len)
    else:
        plot_constitutionals(ax_box, 1, autosome_len, genome_len)

    ax_box.set_xticks(chr_mids)
    ax_box.set_xticklabels(labels, rotation=45, ha="right", fontsize=12)
    ax_box.spines["top"].set_visible(False)
    ax_box.spines["right"].set_visible(False)
    ax_box.spines["bottom"].set_visible(False)
    ax_box.tick_params(axis="x", length=0)

    # Legend
    legend_elements = [
        plt.Line2D(
            [0], [0], color=COLOR_B, lw=2, ls=":", label="Constitutional 1n"
        ),
        plt.Line2D(
            [0], [0], color=COLOR_A, lw=2, ls=":", label="Constitutional 2n"
        ),
        plt.Line2D(
            [0], [0], color=COLOR_C, lw=2, ls=":", label="Constitutional 3n"
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=COLOR_C,
            markersize=8,
            label="Gain",
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=COLOR_B,
            markersize=8,
            label="Loss",
        ),
        plt.Line2D(
            [0],
            [0],
            marker="None",
            color="w",
            label=f"Number of reads: {n_reads_str}",
        ),
    ]
    ax_main.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.15),
        ncol=3,
        frameon=False,
        fontsize=12,
    )

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "genome_wide.png"))
    plt.close()

    # ------------------
    # Chromosome wide plots
    # ------------------
    for i, c in enumerate(chrs):
        start = chr_ends[i]
        end = chr_ends[i + 1]
        if start == end:
            continue

        fig = plt.figure(figsize=(14, 10), dpi=256)
        ax = fig.add_subplot(111)

        whis_lower = l_whis_per_chr[i]
        whis_upper = h_whis_per_chr[i]
        if np.isnan(whis_lower):
            continue

        upper_limit = 0.6 + whis_upper
        lower_limit = -1.05 + whis_lower

        c_ratios = ratios[start:end]
        valid_ratios = c_ratios[~np.isnan(c_ratios)]
        if len(valid_ratios) > 0:
            upper_limit = max(upper_limit, np.max(valid_ratios))
            lower_limit = min(lower_limit, np.min(valid_ratios))

        if ylim_str != "def":
            lower_limit = chr_wide_lower_limit
            upper_limit = chr_wide_upper_limit

        x_c_coords = np.arange(start, end)

        ax.set_xlim(start, end)
        ax.set_ylim(lower_limit, upper_limit)
        ax.set_title(labels[i], fontsize=16)
        ax.set_ylabel(r"$\log_2$(ratio)", fontsize=14)

        # Constitutionals
        ploidy_c = 2
        if (c == 23 or c == 24) and ref_gender == "M":
            ploidy_c = 1
        plot_constitutionals(
            ax,
            ploidy_c,
            start - bins_per_chr[i] * 0.02,
            end + bins_per_chr[i] * 0.02,
        )

        # Missing data lines
        c_nas = np.where(np.isnan(c_ratios))[0] + start
        for na_idx in c_nas:
            ax.axvline(
                na_idx, color=COLOR_A, lw=(1 / (end - start)) * 200, alpha=0.5
            )

        # Scatter
        c_mask = ~np.isnan(c_ratios)
        ax.scatter(
            x_c_coords[c_mask],
            c_ratios[c_mask],
            c=dot_cols[start:end][c_mask],
            s=dot_sizes[start:end][c_mask] * 1.5,
            zorder=3,
            edgecolors="none",
        )

        # Segments
        for ab in results_c:
            c_idx = int(ab[0]) + 1
            if c_idx != c:
                continue

            s = int(ab[1]) + start
            e = int(ab[2]) + start
            height = float(ab[4])
            col = dot_cols[s]

            ax.add_patch(
                patches.Rectangle(
                    (s, 0),
                    e - s,
                    height,
                    linewidth=0.1,
                    edgecolor=col,
                    facecolor=col,
                    alpha=0.3,
                    zorder=1,
                )
            )

            lw = (
                6 * np.nanmean(dot_sizes[s:e]) / 10
                if not np.isnan(np.nanmean(dot_sizes[s:e]))
                else 2
            )
            ax.plot(
                [s, e],
                [height, height],
                color=COLOR_LIGHTER_GREY,
                lw=lw,
                zorder=2,
            )

        # Highlight regions
        for start_bin, end_bin, label, label_pos, label_adj in gene_labels:
            if start_bin >= start and start_bin < end:
                s_b = max(start_bin, start)
                e_b = min(end_bin, end)
                ax.scatter(
                    x_coords[s_b:e_b],
                    ratios[s_b:e_b],
                    facecolors="none",
                    edgecolors=COLOR_D,
                    s=dot_sizes[s_b:e_b] * 2,
                    lw=1.5,
                    zorder=4,
                )
                va = "bottom" if label_adj == 0 else "top"
                ax.text(
                    s_b + (e_b - s_b) / 2,
                    label_pos,
                    label,
                    color=COLOR_D,
                    fontsize=10,
                    rotation=90,
                    va=va,
                    ha="center",
                )

        # X ticks in binsize steps
        num_ticks = 10
        step = bins_per_chr[i] // num_ticks
        if step == 0:
            step = 1
        x_ticks_at = np.arange(0, bins_per_chr[i], step)
        x_tick_labels = (x_ticks_at * binsize).astype(int)

        ax.set_xticks(x_ticks_at + start)
        ax.set_xticklabels(x_tick_labels, rotation=45, ha="right")

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{labels[i]}.png"))
        plt.close()
