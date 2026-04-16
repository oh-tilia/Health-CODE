"""
PCA pipeline: log2(counts+1) -> transpose (samples x features) -> drop zero-variance
columns -> PCA -> explained variance ratios.

Optional outputs with --out-prefix: CSV tables and PCA figures (matplotlib):
scree, sample maps (PC1–PC2, PC1–PC3), top variable loadings,
biplots PC1–PC2 / PC1–PC3 / PC1–PC4 (labels lisibles).
Coloration: détail (lignée|condition), KD vs non-KD (NT), ou les deux (--color-by both, défaut).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def parse_sample_group(sample_id: str) -> str:
    """Derive a short group label from names like counts.CellLine.condition...."""
    s = str(sample_id)
    if s.startswith("counts."):
        parts = s.split(".")
        if len(parts) >= 3:
            return f"{parts[1]} | {parts[2]}"
    return "other"


def parse_cell_line(sample_id: str) -> str:
    """Cell line token after 'counts.' (e.g. HS578T, MDAMB231, SUM159)."""
    s = str(sample_id)
    if s.startswith("counts."):
        parts = s.split(".")
        if len(parts) >= 2:
            return parts[1]
    return "other"


def parse_kd_group(sample_id: str) -> str:
    """Knockdown vs negative control from sample name (sh/si NRP1 vs sh/si NT)."""
    s = str(sample_id)
    if not s.startswith("counts."):
        return "autre"
    parts = s.split(".")
    if len(parts) < 3:
        return "autre"
    cond = parts[2]
    if cond in ("shNT", "siNT"):
        return "non-KD (NT)"
    if cond.startswith("shNRP1") or cond.startswith("siNRP1"):
        return "KD (NRP1)"
    return f"autre ({cond})"


def sample_hue(sample_id: str, hue_mode: str) -> str:
    if hue_mode == "kd":
        return parse_kd_group(sample_id)
    if hue_mode == "cellline":
        return parse_cell_line(sample_id)
    return parse_sample_group(sample_id)


CELL_LINE_MARKERS: dict[str, str] = {
    "HS578T": "o",
    "MDAMB231": "s",
    "SUM159": "^",
    "other": "D",
}


def _kd_group_color(label: str) -> tuple[float, float, float]:
    if label.startswith("KD"):
        return (0.75, 0.2, 0.18)
    if label.startswith("non-KD"):
        return (0.15, 0.45, 0.72)
    return (0.45, 0.45, 0.48)


def _scatter_scores_by_group(
    ax: plt.Axes,
    scores: np.ndarray,
    pc_i: int,
    pc_j: int,
    ratios: np.ndarray,
    sample_index: pd.Index,
    title: str,
    annotate_samples: bool,
    hue_mode: str,
    markers_by_cellline: bool,
) -> None:
    hues = [sample_hue(s, hue_mode) for s in sample_index]
    unique_hues = list(dict.fromkeys(hues))
    use_markers = markers_by_cellline and hue_mode == "kd"
    cells = [parse_cell_line(s) for s in sample_index] if use_markers else None

    if not use_markers:
        if hue_mode == "kd":
            color_of = {h: _kd_group_color(h) for h in unique_hues}
        else:
            cmap = colormaps["tab20" if hue_mode == "detail" else "tab10"]
            nmax = 20 if hue_mode == "detail" else 10
            color_of = {h: cmap(i % nmax) for i, h in enumerate(unique_hues)}
        for h in unique_hues:
            idx = [k for k, g in enumerate(hues) if g == h]
            ax.scatter(
                scores[idx, pc_i],
                scores[idx, pc_j],
                label=h,
                alpha=0.9,
                s=50,
                color=color_of[h],
                edgecolors="white",
                linewidths=0.45,
            )
    else:
        assert cells is not None
        unique_cells = list(dict.fromkeys(cells))
        for h in unique_hues:
            col = _kd_group_color(h)
            for cl in unique_cells:
                idx = [k for k in range(len(sample_index)) if hues[k] == h and cells[k] == cl]
                if not idx:
                    continue
                mk = CELL_LINE_MARKERS.get(cl, "o")
                ax.scatter(
                    scores[idx, pc_i],
                    scores[idx, pc_j],
                    label=f"{h} — {cl}",
                    alpha=0.9,
                    s=52,
                    color=col,
                    marker=mk,
                    edgecolors="white",
                    linewidths=0.45,
                )

    ax.axhline(0, color="grey", lw=0.6, zorder=0)
    ax.axvline(0, color="grey", lw=0.6, zorder=0)
    ri, rj = float(ratios[pc_i]), float(ratios[pc_j])
    ax.set_xlabel(f"PC{pc_i + 1} ({ri:.1%} variance expliquée)")
    ax.set_ylabel(f"PC{pc_j + 1} ({rj:.1%} variance expliquée)")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7.5, framealpha=0.92)

    if annotate_samples:
        for k, sid in enumerate(sample_index):
            short = str(sid).removeprefix("counts.")[:40]
            ax.annotate(
                short,
                (scores[k, pc_i], scores[k, pc_j]),
                fontsize=5,
                alpha=0.85,
                xytext=(3, 2),
                textcoords="offset points",
            )


def _biplot_label_offsets(n: int) -> list[tuple[int, int]]:
    """Cyclic pixel offsets so variable captions overlap less."""
    base = [
        (14, 10),
        (-14, 10),
        (14, -10),
        (-14, -10),
        (22, 0),
        (-22, 0),
        (0, 16),
        (0, -16),
        (18, 14),
        (-18, -14),
        (18, -14),
        (-18, 14),
    ]
    return [base[i % len(base)] for i in range(n)]


def save_biplot(
    prefix: Path,
    name_suffix: str,
    plot_suffix: str,
    scores: np.ndarray,
    comp: np.ndarray,
    ratios: np.ndarray,
    sample_index: pd.Index,
    names: np.ndarray,
    pc_i: int,
    pc_j: int,
    top_n: int,
    n_feat: int,
    hue_mode: str,
    markers_by_cellline: bool,
) -> None:
    """Biplot: samples on (PC_i+1, PC_j+1) + top variable loadings; readable transcript labels."""
    top_n = min(max(top_n, 1), n_feat)
    contrib = comp[pc_i] ** 2 + comp[pc_j] ** 2
    top_idx = np.argsort(-contrib)[:top_n]
    L = np.column_stack([comp[pc_i, top_idx], comp[pc_j, top_idx]])

    si, sj = scores[:, pc_i], scores[:, pc_j]
    max_score = max(float(np.max(np.abs(si))), float(np.max(np.abs(sj))), 1e-12)
    max_load = max(float(np.max(np.abs(L[:, 0]))), float(np.max(np.abs(L[:, 1]))), 1e-12)
    scale = 0.42 * max_score / max_load
    vx = L[:, 0] * scale
    vy = L[:, 1] * scale

    fig, ax = plt.subplots(figsize=(11.5, 9.5))
    hues = [sample_hue(s, hue_mode) for s in sample_index]
    unique_hues = list(dict.fromkeys(hues))
    use_markers = markers_by_cellline and hue_mode == "kd"
    cells = [parse_cell_line(s) for s in sample_index] if use_markers else None

    if not use_markers:
        if hue_mode == "kd":
            color_of = {h: _kd_group_color(h) for h in unique_hues}
        else:
            cmap = colormaps["tab20" if hue_mode == "detail" else "tab10"]
            nmax = 20 if hue_mode == "detail" else 10
            color_of = {h: cmap(i % nmax) for i, h in enumerate(unique_hues)}
        for h in unique_hues:
            idx = [k for k, g in enumerate(hues) if g == h]
            ax.scatter(
                si[idx],
                sj[idx],
                label=h,
                alpha=0.9,
                s=54,
                color=color_of[h],
                edgecolors="white",
                linewidths=0.45,
                zorder=3,
            )
    else:
        assert cells is not None
        unique_cells = list(dict.fromkeys(cells))
        for h in unique_hues:
            col = _kd_group_color(h)
            for cl in unique_cells:
                idx = [k for k in range(len(sample_index)) if hues[k] == h and cells[k] == cl]
                if not idx:
                    continue
                mk = CELL_LINE_MARKERS.get(cl, "o")
                ax.scatter(
                    si[idx],
                    sj[idx],
                    label=f"{h} — {cl}",
                    alpha=0.9,
                    s=56,
                    color=col,
                    marker=mk,
                    edgecolors="white",
                    linewidths=0.45,
                    zorder=3,
                )

    for j, fi in enumerate(top_idx):
        ax.annotate(
            "",
            xy=(vx[j], vy[j]),
            xytext=(0.0, 0.0),
            arrowprops=dict(arrowstyle="->", color="firebrick", lw=0.9, alpha=0.78, shrinkA=0, shrinkB=0),
            zorder=2,
        )

    offsets = _biplot_label_offsets(len(top_idx))
    for j, fi in enumerate(top_idx):
        label = str(names[fi])
        fs = max(6.0, min(8.5, 11.0 - len(label) * 0.12))
        ox, oy = offsets[j]
        ax.annotate(
            label,
            xy=(vx[j], vy[j]),
            xytext=(ox, oy),
            textcoords="offset points",
            fontsize=fs,
            ha="center",
            va="center",
            color="darkred",
            zorder=5,
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="white",
                edgecolor="firebrick",
                linewidth=0.55,
                alpha=0.96,
            ),
            arrowprops=dict(
                arrowstyle="-",
                color="firebrick",
                lw=0.55,
                alpha=0.75,
                shrinkA=2,
                shrinkB=2,
                connectionstyle="arc3,rad=0.12",
            ),
        )

    ax.axhline(0, color="grey", lw=0.55, zorder=1)
    ax.axvline(0, color="grey", lw=0.55, zorder=1)
    ri, rj = float(ratios[pc_i]), float(ratios[pc_j])
    ax.set_xlabel(f"PC{pc_i + 1} ({ri:.1%} variance expliquée)")
    ax.set_ylabel(f"PC{pc_j + 1} ({rj:.1%} variance expliquée)")
    hue_note = {"kd": "couleur KD / non-KD", "detail": "lignée | condition", "cellline": "lignée"}.get(
        hue_mode, hue_mode
    )
    ax.set_title(
        f"Biplot PC{pc_i + 1} vs PC{pc_j + 1} — {hue_note} + top {top_n} variables "
        f"(flèches = loadings redimensionnés)"
    )
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7.5, framealpha=0.92)
    fig.tight_layout()
    out = prefix.parent / f"{prefix.name}_pca_biplot_{name_suffix}{plot_suffix}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def save_pca_plots(
    prefix: Path,
    pca: PCA,
    ratios: np.ndarray,
    scores: np.ndarray,
    sample_index: pd.Index,
    feature_names: list[str],
    plot_n_pcs: int,
    plot_n_vars: int,
    biplot_n_vars: int,
    annotate_samples: bool,
    color_by: str,
    markers_by_cellline: bool,
) -> None:
    """Write scree, sample maps, top loading bars, and biplots (PNG)."""

    def hue_runs() -> list[tuple[str, str]]:
        if color_by == "both":
            return [("detail", ""), ("kd", "_kd")]
        return [(color_by, "")]
    k = min(plot_n_pcs, len(ratios))
    ratios_k = ratios[:k]
    cumsum_k = np.cumsum(ratios_k)
    x = np.arange(1, k + 1)
    pc_labels = [f"PC{i}" for i in x]

    # 1) Scree + cumulative
    fig1, ax1 = plt.subplots(figsize=(9, 4.5))
    ax1.bar(x, ratios_k, color="steelblue", alpha=0.85, label="Variance expliquée")
    ax1.set_xlabel("Composante principale")
    ax1.set_ylabel("Variance expliquée", color="steelblue")
    ax1.set_xticks(x)
    ax1.set_xticklabels(pc_labels, rotation=45, ha="right")
    ax1.tick_params(axis="y", labelcolor="steelblue")

    ax1b = ax1.twinx()
    ax1b.plot(x, cumsum_k, color="darkorange", marker="o", lw=2, label="Variance cumulée")
    ax1b.set_ylabel("Variance cumulée", color="darkorange")
    ax1b.set_ylim(0, 1.05)
    ax1b.tick_params(axis="y", labelcolor="darkorange")
    fig1.suptitle("Éboulis (scree) et variance cumulée", y=1.02)
    fig1.tight_layout()
    scree_path = prefix.parent / f"{prefix.name}_pca_scree.png"
    fig1.savefig(scree_path, dpi=150, bbox_inches="tight")
    plt.close(fig1)
    print(f"Wrote {scree_path}")

    # 2–3) Score maps (one file per hue_mode when color_by=both)
    for hue_mode, file_tag in hue_runs():
        sub = ""
        if hue_mode == "kd":
            sub = " — KD vs non-KD (NT)"
        elif hue_mode == "cellline":
            sub = " — par lignée"

        fig2, ax = plt.subplots(figsize=(8, 6))
        _scatter_scores_by_group(
            ax,
            scores,
            0,
            1,
            ratios,
            sample_index,
            f"Scores PCA — PC1 vs PC2{sub}",
            annotate_samples,
            hue_mode=hue_mode,
            markers_by_cellline=markers_by_cellline,
        )
        fig2.tight_layout()
        pc12_path = prefix.parent / f"{prefix.name}_pca_pc1_pc2{file_tag}.png"
        fig2.savefig(pc12_path, dpi=150, bbox_inches="tight")
        plt.close(fig2)
        print(f"Wrote {pc12_path}")

        if scores.shape[1] >= 3:
            fig3, ax = plt.subplots(figsize=(8, 6))
            _scatter_scores_by_group(
                ax,
                scores,
                0,
                2,
                ratios,
                sample_index,
                f"Scores PCA — PC1 vs PC3{sub}",
                annotate_samples,
                hue_mode=hue_mode,
                markers_by_cellline=markers_by_cellline,
            )
            fig3.tight_layout()
            pc13_path = prefix.parent / f"{prefix.name}_pca_pc1_pc3{file_tag}.png"
            fig3.savefig(pc13_path, dpi=150, bbox_inches="tight")
            plt.close(fig3)
            print(f"Wrote {pc13_path}")

    # Variable loadings: sklearn components_ shape (n_components, n_features)
    comp = pca.components_
    n_feat = comp.shape[1]
    names = np.array(feature_names)
    nv = min(plot_n_vars, n_feat)

    # Top |loading| on PC1 and PC2 (two panels)
    idx1 = np.argsort(-np.abs(comp[0]))[:nv]
    idx2 = np.argsort(-np.abs(comp[1]))[:nv]

    fig4, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12, 6))
    y_pos = np.arange(nv)[::-1]
    ax_a.barh(y_pos, comp[0, idx1], color="steelblue", alpha=0.85)
    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels([str(names[i])[:48] for i in idx1], fontsize=7)
    ax_a.set_xlabel("Loading (poids sur PC1)")
    ax_a.set_title(f"Top {nv} variables — contribution à PC1")

    ax_b.barh(y_pos, comp[1, idx2], color="darkseagreen", alpha=0.85)
    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels([str(names[i])[:48] for i in idx2], fontsize=7)
    ax_b.set_xlabel("Loading (poids sur PC2)")
    ax_b.set_title(f"Top {nv} variables — contribution à PC2")
    fig4.suptitle("Variables (transcrits) les plus liées à PC1 et PC2", y=1.02)
    fig4.tight_layout()
    load_path = prefix.parent / f"{prefix.name}_pca_loadings_pc1_pc2_top.png"
    fig4.savefig(load_path, dpi=150, bbox_inches="tight")
    plt.close(fig4)
    print(f"Wrote {load_path}")

    n_pc_avail = scores.shape[1]
    bn = min(biplot_n_vars, n_feat)
    for hue_mode, file_tag in hue_runs():
        if n_pc_avail >= 2:
            save_biplot(
                prefix,
                "pc1_pc2",
                file_tag,
                scores,
                comp,
                ratios,
                sample_index,
                names,
                0,
                1,
                bn,
                n_feat,
                hue_mode=hue_mode,
                markers_by_cellline=markers_by_cellline,
            )
        if n_pc_avail >= 3:
            save_biplot(
                prefix,
                "pc1_pc3",
                file_tag,
                scores,
                comp,
                ratios,
                sample_index,
                names,
                0,
                2,
                bn,
                n_feat,
                hue_mode=hue_mode,
                markers_by_cellline=markers_by_cellline,
            )
        if n_pc_avail >= 4:
            save_biplot(
                prefix,
                "pc1_pc4",
                file_tag,
                scores,
                comp,
                ratios,
                sample_index,
                names,
                0,
                3,
                bn,
                n_feat,
                hue_mode=hue_mode,
                markers_by_cellline=markers_by_cellline,
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="PCA on RSEM transcript counts.")
    parser.add_argument(
        "counts_tsv",
        nargs="?",
        default=Path(__file__).resolve().parent
        / "GSE266566_COUNTS.rsem_human.transcripts_NRP1_KD_BrCa_llynam.csv",
        type=Path,
        help="Tab-separated counts matrix (default: bundled GSE266566 file).",
    )
    parser.add_argument(
        "--n-components",
        type=int,
        default=None,
        help="Number of PCs (default: min(n_samples, n_features)).",
    )
    parser.add_argument(
        "--scale",
        action="store_true",
        help="Z-score each feature (transcript) before PCA (recommended if scales differ).",
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=None,
        help=(
            "If set, write ratios CSV, scores CSV, and PCA plots (PNG; requires matplotlib). "
            "Plots: scree, PC1–PC2, PC1–PC3, loadings, biplots PC1–PC2 / PC1–PC3 / PC1–PC4. "
            "Filenames: <prefix>_pca_*.csv / *.png"
        ),
    )
    parser.add_argument(
        "--plot-n-pcs",
        type=int,
        default=None,
        help="Number of PCs to show on scree plot (default: min(20, n_components)).",
    )
    parser.add_argument(
        "--plot-n-vars",
        type=int,
        default=25,
        help="Top N transcripts on the loading bar charts (default: 25).",
    )
    parser.add_argument(
        "--biplot-n-vars",
        type=int,
        default=None,
        help="Top N transcripts on biplots (default: min(plot_n_vars, 18) for readable labels).",
    )
    parser.add_argument(
        "--annotate-samples",
        action="store_true",
        help="Draw sample name labels on score scatters (can look busy).",
    )
    parser.add_argument(
        "--color-by",
        choices=["detail", "kd", "cellline", "both"],
        default="both",
        help=(
            "Couleur des points (scores + biplots): detail=lignée|condition; kd=KD (sh/si NRP1) "
            "vs non-KD (sh/si NT); cellline=lignée seule; both=génère detail ET kd (suffixe _kd sur les PNG kd)."
        ),
    )
    parser.add_argument(
        "--markers-by-cellline",
        action="store_true",
        help="Sur les vues codées en kd: distinguer les lignées (HS578T=o, MDAMB231=s, SUM159=^).",
    )
    args = parser.parse_args()

    path: Path = args.counts_tsv
    if not path.is_file():
        raise SystemExit(f"File not found: {path}")

    df = pd.read_csv(path, sep="\t")
    meta_cols = [
        "Ensembl.114.Transcript.ID",
        "Ensembl.114.Transcript.Name",
        "Ensembl.114.Transcript.Type",
    ]
    tid_col = "Ensembl.114.Transcript.ID"
    counts = df.drop(columns=[c for c in meta_cols if c in df.columns])
    counts = counts.apply(pd.to_numeric, errors="coerce")
    if tid_col in df.columns:
        counts = counts.copy()
        counts.index = df[tid_col].astype(str).values
    if counts.isna().any().any():
        n_na = int(counts.isna().sum().sum())
        raise SystemExit(f"Found {n_na} non-numeric/NaN values in count columns; fix upstream.")

    log_counts = np.log2(counts + 1.0)
    X = log_counts.T
    X = X.loc[:, X.var(axis=0) > 0]

    n_samples, n_features = X.shape
    print(f"Samples: {n_samples}, features after variance filter: {n_features}")

    feature_names = list(X.columns.astype(str))

    X_arr = X.to_numpy(dtype=np.float64)
    if args.scale:
        X_arr = StandardScaler(with_mean=True, with_std=True).fit_transform(X_arr)

    max_pc = max(min(n_samples - 1, n_features), 1)
    n_comp = args.n_components
    if n_comp is None:
        n_comp = max_pc
    elif n_comp > max_pc:
        print(f"Note: n_components={args.n_components} capped to {max_pc} (centered PCA).")
        n_comp = max_pc

    pca = PCA(
        n_components=n_comp,
        svd_solver="full" if min(n_samples, n_features) <= 500 else "auto",
    )
    scores = pca.fit_transform(X_arr)

    ratios = pca.explained_variance_ratio_
    cumsum = np.cumsum(ratios)

    print("\nExplained variance ratio per PC:")
    for i, (r, c) in enumerate(zip(ratios, cumsum, strict=True), start=1):
        print(f"  PC{i}: {r:.6f}  (cumulative: {c:.6f})")

    plot_k = args.plot_n_pcs if args.plot_n_pcs is not None else min(20, len(ratios))
    plot_nv = max(1, args.plot_n_vars)
    biplot_nv = args.biplot_n_vars if args.biplot_n_vars is not None else min(plot_nv, 18)

    if args.out_prefix is not None:
        prefix = args.out_prefix
        prefix.parent.mkdir(parents=True, exist_ok=True)
        ratios_df = pd.DataFrame(
            {
                "PC": [f"PC{i}" for i in range(1, len(ratios) + 1)],
                "explained_variance_ratio": ratios,
                "cumulative_explained_variance_ratio": cumsum,
            }
        )
        ratios_path = prefix.parent / f"{prefix.name}_pca_variance_ratios.csv"
        ratios_df.to_csv(ratios_path, index=False)
        scores_df = pd.DataFrame(scores, index=X.index, columns=[f"PC{i}" for i in range(1, scores.shape[1] + 1)])
        scores_path = prefix.parent / f"{prefix.name}_pca_scores.csv"
        scores_df.to_csv(scores_path)
        print(f"\nWrote {ratios_path}")
        print(f"Wrote {scores_path}")

        save_pca_plots(
            prefix,
            pca,
            ratios,
            scores,
            X.index,
            feature_names,
            plot_n_pcs=plot_k,
            plot_n_vars=plot_nv,
            biplot_n_vars=biplot_nv,
            annotate_samples=args.annotate_samples,
            color_by=args.color_by,
            markers_by_cellline=args.markers_by_cellline,
        )


if __name__ == "__main__":
    main()
