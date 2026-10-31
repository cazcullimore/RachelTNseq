#!/usr/bin/env python3
"""Run PyDESeq2 on a TnSeq GENE count table (raw counts only)."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

MIN_LIBRARY_COUNTS = 100_000
MIN_GENE_COUNTS = 10
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 1.0
SAMPLE_RE = re.compile(r"^(M108|M1)_(HIS|NHS)_R(\d)_")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run PyDESeq2 on raw (non-normalized) sample columns from a "
            "tab-delimited GENE count table."
        )
    )
    parser.add_argument("input", type=Path, help="GENE.txt count table")
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        required=True,
        help="Directory for results tables and plots",
    )
    parser.add_argument(
        "--min-library-counts",
        type=int,
        default=MIN_LIBRARY_COUNTS,
        help=(
            "Drop samples with fewer total raw counts than this "
            f"(default: {MIN_LIBRARY_COUNTS})"
        ),
    )
    parser.add_argument(
        "--min-gene-counts",
        type=int,
        default=MIN_GENE_COUNTS,
        help=f"Keep genes with at least this many total counts (default: {MIN_GENE_COUNTS})",
    )
    parser.add_argument("--alpha", type=float, default=PADJ_CUTOFF)
    parser.add_argument("--n-cpus", type=int, default=4)
    parser.add_argument(
        "--only-r1-vs-r2",
        action="store_true",
        help=(
            "Only run R1 vs R2 contrasts (mislabeling check). "
            "Does not overwrite PCA or HIS-vs-NHS result tables."
        ),
    )
    parser.add_argument(
        "--diagnostic-ratios",
        action="store_true",
        help="Write matched replicate HIS/NHS enrichment ratios without statistical testing",
    )
    parser.add_argument(
        "--group-a",
        help="Comma-separated sample prefixes for a custom A-vs-B comparison",
    )
    parser.add_argument(
        "--group-b",
        help="Comma-separated sample prefixes for a custom A-vs-B comparison",
    )
    parser.add_argument(
        "--custom-stem",
        default="custom_A_vs_B",
        help="Output filename stem for a custom comparison",
    )
    return parser.parse_args()


def raw_count_columns(columns: list[str]) -> list[str]:
    return [
        name
        for name in columns
        if name.endswith("_filtered") and not name.endswith("(Normalized)")
    ]


def parse_sample(column: str) -> dict[str, str]:
    match = SAMPLE_RE.match(column)
    if not match:
        raise ValueError(f"Cannot parse strain/condition/replicate from {column!r}")
    strain, condition, replicate = match.groups()
    return {
        "sample": column,
        "strain": strain,
        "condition": condition,
        "replicate": f"R{replicate}",
        "group": f"{strain}_{condition}",
    }


def load_table(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    table = pd.read_csv(path, sep="\t")
    if "GeneID" not in table.columns:
        raise SystemExit("Input is missing a GeneID column")
    if table["GeneID"].duplicated().any():
        raise SystemExit("GeneID values are not unique")
    sample_cols = raw_count_columns(table.columns.tolist())
    if not sample_cols:
        raise SystemExit("No raw count columns ending in '_filtered' were found")
    annot_cols = [
        name
        for name in (
            "GeneID",
            "Gene",
            "Function",
            "Start",
            "Stop",
            "Order",
            "Strand",
            "Length",
            "PID",
        )
        if name in table.columns
    ]
    annot = table[annot_cols].set_index("GeneID")
    counts = table.set_index("GeneID")[sample_cols].apply(pd.to_numeric, errors="coerce")
    if counts.isna().any().any():
        bad = counts.columns[counts.isna().any()].tolist()
        raise SystemExit(f"Non-numeric raw counts in columns: {bad}")
    counts = counts.fillna(0).clip(lower=0).round().astype(int)
    return counts, annot, sample_cols


def build_metadata(sample_cols: list[str]) -> pd.DataFrame:
    meta = pd.DataFrame([parse_sample(col) for col in sample_cols]).set_index("sample")
    meta["strain"] = pd.Categorical(meta["strain"], categories=["M1", "M108"])
    meta["condition"] = pd.Categorical(meta["condition"], categories=["NHS", "HIS"])
    meta["replicate"] = pd.Categorical(meta["replicate"], categories=["R1", "R2"])
    return meta


def qc_drop_low_depth(
    counts: pd.DataFrame, metadata: pd.DataFrame, min_library_counts: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, list[str]]:
    library_sizes = counts.sum(axis=0)
    keep = library_sizes >= min_library_counts
    dropped = library_sizes.index[~keep].tolist()
    return counts.loc[:, keep], metadata.loc[keep], library_sizes, dropped


def subset_samples(counts: pd.DataFrame, metadata: pd.DataFrame, mask: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = metadata.loc[mask].copy()
    meta["strain"] = meta["strain"].cat.remove_unused_categories()
    meta["condition"] = meta["condition"].cat.remove_unused_categories()
    meta["replicate"] = meta["replicate"].cat.remove_unused_categories()
    return counts.loc[:, meta.index], meta


def filter_genes(counts: pd.DataFrame, min_gene_counts: int) -> pd.DataFrame:
    keep = counts.sum(axis=1) >= min_gene_counts
    return counts.loc[keep]


def write_diagnostic_ratios(
    counts: pd.DataFrame, metadata: pd.DataFrame, annot: pd.DataFrame, outdir: Path
) -> None:
    library_sizes = counts.sum(axis=0)
    pseudocount = 0.5
    for strain in ("M1", "M108"):
        for replicate in ("R1", "R2"):
            mask = (
                metadata["strain"].astype(str).eq(strain)
                & metadata["replicate"].astype(str).eq(replicate)
            )
            samples = metadata.loc[mask]
            if set(samples["condition"].astype(str)) != {"HIS", "NHS"}:
                raise SystemExit(f"Missing matched HIS/NHS samples for {strain} {replicate}")
            nhs = samples.index[samples["condition"].astype(str).eq("NHS")][0]
            his = samples.index[samples["condition"].astype(str).eq("HIS")][0]
            nhs_cpm = (counts[nhs] / library_sizes[nhs]) * 1_000_000
            his_cpm = (counts[his] / library_sizes[his]) * 1_000_000
            result = pd.DataFrame(
                {
                    "NHS_sample": nhs,
                    "HIS_sample": his,
                    "NHS_raw_count": counts[nhs],
                    "HIS_raw_count": counts[his],
                    "NHS_CPM": nhs_cpm,
                    "HIS_CPM": his_cpm,
                    "NHS_over_HIS_ratio": (counts[nhs] + pseudocount)
                    / library_sizes[nhs]
                    / ((counts[his] + pseudocount) / library_sizes[his]),
                    "HIS_over_NHS_ratio": (counts[his] + pseudocount)
                    / library_sizes[his]
                    / ((counts[nhs] + pseudocount) / library_sizes[nhs]),
                },
                index=counts.index,
            )
            result["log2_HIS_over_NHS"] = np.log2(result["HIS_over_NHS_ratio"])
            result = result.join(annot, how="left")
            result.to_csv(
                outdir / f"{strain}_{replicate}_NHS_vs_HIS_diagnostic.tsv",
                sep="\t",
            )
            print(
                f"{strain} {replicate}: NHS={nhs} HIS={his} "
                f"(library sizes {library_sizes[nhs]:,}, {library_sizes[his]:,})"
            )


def run_wald(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    contrast: list[str],
    inference: DefaultInference,
    alpha: float,
    min_gene_counts: int,
) -> tuple[DeseqDataSet, pd.DataFrame]:
    counts = filter_genes(counts, min_gene_counts)
    # PyDESeq2 wants samples as rows, genes as columns.
    dds = DeseqDataSet(
        counts=counts.T,
        metadata=metadata,
        design=design,
        size_factors_fit_type="poscounts",
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=contrast,
        alpha=alpha,
        cooks_filter=True,
        independent_filter=True,
        inference=inference,
        quiet=True,
    )
    stats.summary()
    result = stats.results_df.copy()
    result.index.name = "GeneID"
    return dds, result


def annotate_results(result: pd.DataFrame, annot: pd.DataFrame, alpha: float) -> pd.DataFrame:
    out = result.join(annot, how="left")
    out["significant"] = (out["padj"] < alpha) & out["padj"].notna()
    out["direction"] = "NS"
    up = out["significant"] & (out["log2FoldChange"] > 0)
    down = out["significant"] & (out["log2FoldChange"] < 0)
    out.loc[up, "direction"] = "up"
    out.loc[down, "direction"] = "down"
    out = out.sort_values(["padj", "pvalue"], na_position="last")
    return out


def label_for_gene(row: pd.Series) -> str:
    gene = str(row.get("Gene", "")).strip()
    if gene and gene not in {"-", "nan", "None"}:
        return gene
    return str(row.name)


def volcano(result: pd.DataFrame, title: str, outfile: Path, alpha: float) -> None:
    df = result.copy()
    df["neglog10p"] = -np.log10(df["pvalue"].clip(lower=np.nextafter(0, 1)))
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ns = df["direction"].eq("NS") | df["direction"].isna()
    ax.scatter(
        df.loc[ns, "log2FoldChange"],
        df.loc[ns, "neglog10p"],
        s=8,
        c="#b0b0b0",
        alpha=0.6,
        linewidths=0,
        label="NS",
    )
    up = df["direction"].eq("up")
    down = df["direction"].eq("down")
    ax.scatter(
        df.loc[up, "log2FoldChange"],
        df.loc[up, "neglog10p"],
        s=14,
        c="#d62728",
        alpha=0.85,
        linewidths=0,
        label=f"padj < {alpha}, LFC > 0",
    )
    ax.scatter(
        df.loc[down, "log2FoldChange"],
        df.loc[down, "neglog10p"],
        s=14,
        c="#1f77b4",
        alpha=0.85,
        linewidths=0,
        label=f"padj < {alpha}, LFC < 0",
    )
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.axhline(-np.log10(alpha), color="black", linewidth=0.8, linestyle=":")
    labeled = df.loc[df["significant"]].head(18)
    for _, row in labeled.iterrows():
        ax.annotate(
            label_for_gene(row),
            (row["log2FoldChange"], row["neglog10p"]),
            fontsize=7,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 p-value")
    ax.set_title(title)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)


def ma_plot(result: pd.DataFrame, title: str, outfile: Path, alpha: float) -> None:
    df = result.copy()
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ns = ~df["significant"].fillna(False)
    ax.scatter(
        df.loc[ns, "baseMean"].clip(lower=0.1),
        df.loc[ns, "log2FoldChange"],
        s=8,
        c="#b0b0b0",
        alpha=0.5,
        linewidths=0,
        label="NS",
    )
    sig = df["significant"].fillna(False)
    ax.scatter(
        df.loc[sig, "baseMean"].clip(lower=0.1),
        df.loc[sig, "log2FoldChange"],
        s=14,
        c="#d62728",
        alpha=0.85,
        linewidths=0,
        label=f"padj < {alpha}",
    )
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlabel("mean of normalized counts")
    ax.set_ylabel("log2 fold change")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)


def pca_plot(dds: DeseqDataSet, metadata: pd.DataFrame, outfile: Path, title: str) -> pd.DataFrame:
    dds.vst(use_design=True)
    vst = pd.DataFrame(
        dds.layers["vst_counts"],
        index=dds.obs_names,
        columns=dds.var_names,
    )
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(vst.values)
    pcs = pd.DataFrame(coords, index=vst.index, columns=["PC1", "PC2"])
    pcs = pcs.join(metadata)
    fig, ax = plt.subplots(figsize=(7.5, 6))
    markers = {"HIS": "o", "NHS": "s"}
    colors = {"M1": "#1f77b4", "M108": "#ff7f0e"}
    for strain in pcs["strain"].unique():
        for condition in pcs["condition"].unique():
            sub = pcs[(pcs["strain"] == strain) & (pcs["condition"] == condition)]
            if sub.empty:
                continue
            ax.scatter(
                sub["PC1"],
                sub["PC2"],
                c=colors.get(str(strain), "gray"),
                marker=markers.get(str(condition), "o"),
                s=90,
                edgecolors="black",
                linewidths=0.6,
                label=f"{strain} {condition}",
            )
            for sample, row in sub.iterrows():
                ax.annotate(str(sample).split("_")[2], (row["PC1"], row["PC2"]), fontsize=7, xytext=(5, 5), textcoords="offset points")
    var = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)
    return pcs


def write_result(
    result: pd.DataFrame,
    annot: pd.DataFrame,
    outdir: Path,
    stem: str,
    title: str,
    contrast_note: str,
    alpha: float,
) -> dict[str, int]:
    annotated = annotate_results(result, annot, alpha)
    annotated.to_csv(outdir / f"{stem}.tsv", sep="\t")
    sig = annotated.loc[annotated["significant"]].copy()
    sig.to_csv(outdir / f"{stem}_significant.tsv", sep="\t")
    volcano(annotated, title, outdir / f"{stem}_volcano.png", alpha)
    ma_plot(annotated, title, outdir / f"{stem}_MA.png", alpha)
    n = len(annotated)
    n_sig = int(annotated["significant"].sum())
    n_up = int((annotated["direction"] == "up").sum())
    n_down = int((annotated["direction"] == "down").sum())
    print(f"{stem}: {n_sig}/{n} significant (up={n_up}, down={n_down})  [{contrast_note}]")
    return {"genes": n, "significant": n_sig, "up": n_up, "down": n_down}


def main() -> int:
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    counts, annot, sample_cols = load_table(args.input)
    metadata = build_metadata(sample_cols)
    if args.diagnostic_ratios:
        write_diagnostic_ratios(counts, metadata, annot, outdir)
        print(f"\nWrote diagnostic ratios to {outdir}")
        return 0
    counts, metadata, library_sizes, dropped = qc_drop_low_depth(
        counts, metadata, args.min_library_counts
    )
    if counts.shape[1] < 3:
        raise SystemExit(f"Only {counts.shape[1]} samples passed QC; need at least 3")

    if bool(args.group_a) != bool(args.group_b):
        raise SystemExit("--group-a and --group-b must be provided together")

    if args.group_a and args.group_b:
        selected: dict[str, str] = {}
        for label, prefixes in (("A", args.group_a), ("B", args.group_b)):
            for prefix in prefixes.split(","):
                matches = [sample for sample in metadata.index if sample.startswith(prefix + "_")]
                if len(matches) != 1:
                    raise SystemExit(f"{prefix!r} matched {len(matches)} samples: {matches}")
                if matches[0] in selected:
                    raise SystemExit(f"Sample selected more than once: {matches[0]}")
                selected[matches[0]] = label
        if (
            len(selected) < 4
            or list(selected.values()).count("A") < 2
            or list(selected.values()).count("B") < 2
        ):
            raise SystemExit("Custom comparison needs at least two samples in each group")
        counts = counts.loc[:, list(selected)]
        metadata = metadata.loc[counts.columns].copy()
        metadata["comparison_group"] = pd.Categorical(
            [selected[sample] for sample in metadata.index], categories=["A", "B"]
        )
        library_sizes = counts.sum(axis=0)
        counts.to_csv(outdir / "raw_counts_used.tsv", sep="\t")
        metadata.to_csv(outdir / "sample_metadata.tsv", sep="\t")
        library_sizes.to_csv(outdir / "library_sizes.tsv", sep="\t", header=["total_raw_counts"])
        inference = DefaultInference(n_cpus=args.n_cpus)
        design = "~comparison_group"
        contrast = ["comparison_group", "A", "B"]
        print(f"\nRunning {args.custom_stem}: design={design} contrast={contrast} n={metadata.shape[0]}")
        print(metadata[["strain", "condition", "replicate", "comparison_group"]].to_string())
        _dds, result = run_wald(
            counts,
            metadata,
            design,
            contrast,
            inference,
            args.alpha,
            args.min_gene_counts,
        )
        stats = write_result(
            result,
            annot,
            outdir,
            args.custom_stem,
            "Custom group A vs B",
            "LFC>0: higher insertion counts in group A than group B",
            args.alpha,
        )
        summary = [
            f"input\t{args.input}",
            f"samples_in_file\t{len(sample_cols)}",
            f"samples_used\t{counts.shape[1]}",
            f"dropped_low_depth\t{','.join(dropped) if dropped else 'none'}",
            f"group_A\t{args.group_a}",
            f"group_B\t{args.group_b}",
            f"design\t{design}",
            f"contrast\tA vs B (comparison_group)",
            f"genes_tested\t{stats['genes']}",
            f"significant\t{stats['significant']}",
            f"up\t{stats['up']}",
            f"down\t{stats['down']}",
        ]
        (outdir / "summary.tsv").write_text("\n".join(summary) + "\n")
        print(f"\nWrote results to {outdir}")
        return 0

    library_sizes.to_csv(outdir / "library_sizes.tsv", sep="\t", header=["total_raw_counts"])
    metadata.to_csv(outdir / "sample_metadata.tsv", sep="\t")
    counts.to_csv(outdir / "raw_counts_used.tsv", sep="\t")

    inference = DefaultInference(n_cpus=args.n_cpus)
    summaries: list[str] = []
    summaries.append(f"input\t{args.input}")
    summaries.append(f"samples_in_file\t{len(sample_cols)}")
    summaries.append(f"samples_used\t{counts.shape[1]}")
    summaries.append(f"dropped_low_depth\t{','.join(dropped) if dropped else 'none'}")
    summaries.append(f"genes_in_file\t{counts.shape[0]}")
    summaries.append(
        "note\tDESeq2 uses RAW integer counts only; normalized columns were ignored"
    )
    summaries.append(
        "note\tposcounts size factors were used because TnSeq tables are zero-heavy"
    )
    summaries.append(
        f"note\tpositive log2FoldChange means higher insertion counts in the first contrast level"
    )
    if dropped:
        print(
            "Dropped low-depth samples "
            f"(<{args.min_library_counts} total counts): {dropped}"
        )
        print("Library sizes:")
        print(library_sizes.to_string())

    comparisons: list[tuple[str, pd.Series, str, list[str], str, str]] = []

    m1 = metadata["strain"].astype(str).eq("M1")
    m108 = metadata["strain"].astype(str).eq("M108")
    nhs = metadata["condition"].astype(str).eq("NHS")
    his = metadata["condition"].astype(str).eq("HIS")

    if not args.only_r1_vs_r2:
        # HIS vs NHS within M1 (2 vs 2)
        if m1.sum() >= 3 and metadata.loc[m1, "condition"].nunique() == 2:
            comparisons.append(
                (
                    "M1_HIS_vs_NHS",
                    m1,
                    "~condition",
                    ["condition", "HIS", "NHS"],
                    "M1: HIS vs NHS",
                    "LFC>0: more insertions in HIS than NHS (M1 only)",
                )
            )
            # PCA shows R1 vs R2 is the main axis, so also estimate HIS vs NHS
            # after blocking on replicate (paired HIS_R1/NHS_R1 and HIS_R2/NHS_R2).
            if metadata.loc[m1, "replicate"].nunique() >= 2:
                comparisons.append(
                    (
                        "M1_HIS_vs_NHS_paired",
                        m1,
                        "~replicate + condition",
                        ["condition", "HIS", "NHS"],
                        "M1: HIS vs NHS (paired by replicate)",
                        "LFC>0: more insertions in HIS than NHS (M1, controlling for R1/R2)",
                    )
                )

        # HIS vs NHS within M108 (may be 1 vs 2 if HIS_R1 was dropped)
        if m108.sum() >= 3 and metadata.loc[m108, "condition"].nunique() == 2:
            n_his = int(metadata.loc[m108, "condition"].astype(str).eq("HIS").sum())
            note = "LFC>0: more insertions in HIS than NHS (M108 only)"
            if n_his < 2:
                note += f"; HIS has only {n_his} replicate after QC"
            comparisons.append(
                (
                    "M108_HIS_vs_NHS",
                    m108,
                    "~condition",
                    ["condition", "HIS", "NHS"],
                    "M108: HIS vs NHS",
                    note,
                )
            )

        # All samples, HIS vs NHS controlling for strain
        if metadata["condition"].nunique() == 2 and metadata["strain"].nunique() >= 1:
            design = "~strain + condition" if metadata["strain"].nunique() == 2 else "~condition"
            comparisons.append(
                (
                    "all_samples_HIS_vs_NHS",
                    pd.Series(True, index=metadata.index),
                    design,
                    ["condition", "HIS", "NHS"],
                    "All samples: HIS vs NHS (strain as covariate)"
                    if metadata["strain"].nunique() == 2
                    else "All samples: HIS vs NHS",
                    "LFC>0: more insertions in HIS than NHS, controlling for strain",
                )
            )

        # Strain comparison within NHS (2 vs 2)
        if nhs.sum() >= 3 and metadata.loc[nhs, "strain"].nunique() == 2:
            comparisons.append(
                (
                    "NHS_M108_vs_M1",
                    nhs,
                    "~strain",
                    ["strain", "M108", "M1"],
                    "NHS: M108 vs M1",
                    "LFC>0: more insertions in M108 than M1 (NHS only)",
                )
            )

        # Strain comparison within HIS
        if his.sum() >= 3 and metadata.loc[his, "strain"].nunique() == 2:
            n_m108_his = int(metadata.loc[his, "strain"].astype(str).eq("M108").sum())
            note = "LFC>0: more insertions in M108 than M1 (HIS only)"
            if n_m108_his < 2:
                note += f"; M108 HIS has only {n_m108_his} replicate after QC"
            comparisons.append(
                (
                    "HIS_M108_vs_M1",
                    his,
                    "~strain",
                    ["strain", "M108", "M1"],
                    "HIS: M108 vs M1",
                    note,
                )
            )

    # R1 vs R2 mislabeling checks. True biological replicates should yield far
    # fewer DE genes than HIS vs NHS; a large R1/R2 effect suggests mix-up or batch.
    if m1.sum() >= 3 and metadata.loc[m1, "replicate"].nunique() == 2:
        comparisons.append(
            (
                "M1_R1_vs_R2",
                m1,
                "~condition + replicate",
                ["replicate", "R1", "R2"],
                "M1: R1 vs R2 (controlling for HIS/NHS)",
                "LFC>0: more insertions in R1 than R2 (M1, controlling for condition)",
            )
        )
        comparisons.append(
            (
                "M1_R1_vs_R2_unpaired",
                m1,
                "~replicate",
                ["replicate", "R1", "R2"],
                "M1: R1 vs R2 (ignoring HIS/NHS)",
                "LFC>0: more insertions in R1 than R2 (M1, not controlling for condition)",
            )
        )

    if nhs.sum() >= 3 and metadata.loc[nhs, "replicate"].nunique() == 2:
        nhs_design = (
            "~strain + replicate" if metadata.loc[nhs, "strain"].nunique() == 2 else "~replicate"
        )
        comparisons.append(
            (
                "NHS_R1_vs_R2",
                nhs,
                nhs_design,
                ["replicate", "R1", "R2"],
                "NHS: R1 vs R2 (controlling for strain)"
                if metadata.loc[nhs, "strain"].nunique() == 2
                else "NHS: R1 vs R2",
                "LFC>0: more insertions in R1 than R2 (NHS only, controlling for strain)",
            )
        )

    if his.sum() >= 3 and metadata.loc[his, "replicate"].nunique() == 2:
        his_design = (
            "~strain + replicate"
            if metadata.loc[his, "strain"].nunique() == 2 and his.sum() > 3
            else "~replicate"
        )
        his_note = "LFC>0: more insertions in R1 than R2 (HIS only)"
        if his_design == "~replicate" and metadata.loc[his, "strain"].nunique() == 2:
            his_note += "; strain is confounded with replicate after QC"
        comparisons.append(
            (
                "HIS_R1_vs_R2",
                his,
                his_design,
                ["replicate", "R1", "R2"],
                "HIS: R1 vs R2 (controlling for strain)"
                if his_design == "~strain + replicate"
                else "HIS: R1 vs R2",
                his_note,
            )
        )

    # Drop strain×condition groups that only have one replicate (M108 HIS after QC)
    present_groups = set(
        zip(metadata["strain"].astype(str), metadata["condition"].astype(str), metadata["replicate"].astype(str))
    )
    both_reps = {
        (strain, condition)
        for strain, condition, _rep in present_groups
        if (strain, condition, "R1") in present_groups and (strain, condition, "R2") in present_groups
    }
    balanced = metadata.apply(
        lambda row: (str(row["strain"]), str(row["condition"])) in both_reps, axis=1
    )
    if balanced.sum() >= 3 and metadata.loc[balanced, "replicate"].nunique() == 2:
        covars = []
        if metadata.loc[balanced, "strain"].nunique() == 2:
            covars.append("strain")
        if metadata.loc[balanced, "condition"].nunique() == 2:
            covars.append("condition")
        design = "~" + " + ".join(covars + ["replicate"]) if covars else "~replicate"
        comparisons.append(
            (
                "all_samples_R1_vs_R2",
                balanced,
                design,
                ["replicate", "R1", "R2"],
                "Balanced samples: R1 vs R2",
                "LFC>0: more insertions in R1 than R2; groups missing a replicate were dropped",
            )
        )

    if not comparisons:
        raise SystemExit("No valid DESeq2 comparisons could be constructed from the samples")

    if not args.only_r1_vs_r2:
        # PCA from the broadest model (all remaining samples)
        all_counts = filter_genes(counts, args.min_gene_counts)
        pca_design = "~strain + condition" if metadata["strain"].nunique() == 2 else "~condition"
        pca_dds = DeseqDataSet(
            counts=all_counts.T,
            metadata=metadata,
            design=pca_design,
            size_factors_fit_type="poscounts",
            inference=inference,
            quiet=True,
        )
        pca_dds.deseq2()
        size_factors = pca_dds.obs["size_factors"].rename("size_factor")
        size_factors.to_csv(outdir / "size_factors.tsv", sep="\t")
        norm_counts = pd.DataFrame(
            pca_dds.layers["normed_counts"],
            index=pca_dds.obs_names,
            columns=pca_dds.var_names,
        ).T
        norm_counts.to_csv(outdir / "deseq2_normalized_counts.tsv", sep="\t")
        pcs = pca_plot(pca_dds, metadata, outdir / "pca.png", "VST PCA of samples used in DESeq2")
        pcs.to_csv(outdir / "pca_coordinates.tsv", sep="\t")

    for stem, mask, design, contrast, title, note in comparisons:
        sub_counts, sub_meta = subset_samples(counts, metadata, mask)
        print(f"\nRunning {stem}: design={design} contrast={contrast} n={sub_meta.shape[0]}")
        print(sub_meta[["strain", "condition", "replicate"]].to_string())
        dds, result = run_wald(
            sub_counts,
            sub_meta,
            design,
            contrast,
            inference,
            args.alpha,
            args.min_gene_counts,
        )
        stats = write_result(result, annot, outdir, stem, title, note, args.alpha)
        summaries.append(f"{stem}_design\t{design}")
        summaries.append(f"{stem}_contrast\t{contrast[1]} vs {contrast[2]} ({contrast[0]})")
        summaries.append(f"{stem}_n_samples\t{sub_meta.shape[0]}")
        summaries.append(f"{stem}_genes_tested\t{stats['genes']}")
        summaries.append(f"{stem}_significant\t{stats['significant']}")
        summaries.append(f"{stem}_up\t{stats['up']}")
        summaries.append(f"{stem}_down\t{stats['down']}")
        summaries.append(f"{stem}_note\t{note}")

    summary_name = "r1_vs_r2_summary.tsv" if args.only_r1_vs_r2 else "summary.tsv"
    (outdir / summary_name).write_text("\n".join(summaries) + "\n")
    print(f"\nWrote results to {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
