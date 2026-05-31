#!/usr/bin/env python3
# Copyright (c) Rishabh Gupta 2026
# This is a test code used for learning.

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def robustSigma(values):
    x = np.asarray(values)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return np.nan
    q16, q84 = np.percentile(x, [16, 84])
    return 0.5 * (q84 - q16)


def gaussianCoreStats(values, nsigma=3.0, niter=3):
    x = np.asarray(values)
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return np.nan, np.nan, len(x)

    mask = np.ones(len(x), dtype=bool)
    for _ in range(niter):
        xx = x[mask]
        if len(xx) < 5:
            break
        mu = np.mean(xx)
        sig = np.std(xx)
        if sig <= 0:
            break
        mask = np.abs(x - mu) < nsigma * sig

    xx = x[mask]
    if len(xx) < 5:
        return np.nan, np.nan, len(xx)

    return np.mean(xx), np.std(xx), len(xx)


def saveHist(data, bins, xlabel, ylabel, title, outpath):
    data = np.asarray(data)
    data = data[np.isfinite(data)]

    plt.figure(figsize=(6.6, 4.5))
    plt.hist(data, bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def saveScatter(x, y, xlabel, ylabel, title, outpath, alpha=0.35, s=8):
    x = np.asarray(x)
    y = np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y)

    plt.figure(figsize=(6.6, 4.5))
    plt.scatter(x[mask], y[mask], alpha=alpha, s=s)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def saveLine(x, y, xlabel, ylabel, title, outpath):
    plt.figure(figsize=(6.6, 4.5))
    plt.plot(x, y, marker="o")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def requireColumns(df, columns):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise RuntimeError(
            "Missing CSV columns:\n"
            + "\n".join(missing)
            + "\n\nCheck the CSV format."
        )


def thresholdScan(df, outdir):
    thresholds = np.arange(0, 61, 1)
    rows = []

    for threshold in thresholds:
        passed = (df["npe_left"] >= threshold) & (df["npe_right"] >= threshold)
        sample = df[passed].copy()

        if len(sample) >= 10:
            timingSigma = sample["dt_ps"].std()
            robustTimingSigma = robustSigma(sample["dt_ps"])
            positionSigma = sample["dz_um"].std()
            robustPositionSigma = robustSigma(sample["dz_um"])
        else:
            timingSigma = np.nan
            robustTimingSigma = np.nan
            positionSigma = np.nan
            robustPositionSigma = np.nan

        rows.append(
            {
                "threshold_npe": threshold,
                "accepted_hits": len(sample),
                "efficiency": len(sample) / len(df),
                "timing_sigma_ps": timingSigma,
                "timing_robust_sigma_ps": robustTimingSigma,
                "position_sigma_um": positionSigma,
                "position_robust_sigma_um": robustPositionSigma,
            }
        )

    scan = pd.DataFrame(rows)
    scan.to_csv(outdir / "threshold_scan.csv", index=False)

    saveLine(
        scan["threshold_npe"],
        scan["efficiency"],
        "Npe threshold",
        "Hit fraction",
        "Efficiency vs threshold",
        outdir / "efficiency_vs_threshold.png",
    )
    saveLine(
        scan["threshold_npe"],
        scan["timing_robust_sigma_ps"],
        "Npe threshold",
        "Timing [ps]",
        "Timing vs threshold",
        outdir / "timing_resolution_vs_threshold.png",
    )
    saveLine(
        scan["threshold_npe"],
        scan["position_robust_sigma_um"],
        "Npe threshold",
        r"Position [$\mu$m]",
        "Position vs threshold",
        outdir / "position_resolution_vs_threshold.png",
    )

    return scan


def resolutionVsNpe(valid, outdir):
    bins = np.array([0, 5, 10, 15, 20, 30, 40, 60, 80, 120, 180, 260, 400, 700])
    rows = []

    for low, high in zip(bins[:-1], bins[1:]):
        chunk = valid[(valid["npe_total"] >= low) & (valid["npe_total"] < high)]
        if len(chunk) < 20:
            continue
        rows.append(
            {
                "npe_low": low,
                "npe_high": high,
                "npe_center": 0.5 * (low + high),
                "count": len(chunk),
                "timing_sigma_ps": chunk["dt_ps"].std(),
                "timing_robust_sigma_ps": robustSigma(chunk["dt_ps"]),
                "position_sigma_um": chunk["dz_um"].std(),
                "position_robust_sigma_um": robustSigma(chunk["dz_um"]),
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(outdir / "resolution_vs_npe.csv", index=False)

    if len(out) > 0:
        saveLine(
            out["npe_center"],
            out["timing_robust_sigma_ps"],
            r"Total Npe $N_{\rm pe,L}+N_{\rm pe,R}$",
            "Timing [ps]",
            "Timing vs Npe",
            outdir / "timing_resolution_vs_npe.png",
        )
        saveLine(
            out["npe_center"],
            out["position_robust_sigma_um"],
            r"Total Npe $N_{\rm pe,L}+N_{\rm pe,R}$",
            r"Position [$\mu$m]",
            "Position vs Npe",
            outdir / "position_resolution_vs_npe.png",
        )

    return out


def occupancyPlots(df, outdir):
    layerCounts = df["layerID"].value_counts().sort_index()

    plt.figure(figsize=(6.4, 4.4))
    plt.bar(layerCounts.index.astype(str), layerCounts.values)
    plt.xlabel("Layer ID")
    plt.ylabel("Hit rows")
    plt.title("Layer occupancy")
    plt.tight_layout()
    plt.savefig(outdir / "layer_occupancy.png", dpi=300)
    plt.close()

    plt.figure(figsize=(7.4, 4.8))
    for layer in sorted(df["layerID"].unique()):
        sample = df[df["layerID"] == layer]
        counts = sample["fiberID"].value_counts().sort_index()
        plt.plot(counts.index, counts.values, marker="o", label=f"Layer {layer}")

    plt.xlabel("Fibre ID")
    plt.ylabel("Hit rows")
    plt.title("Fiber occupancy")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "fiber_occupancy_by_layer.png", dpi=300)
    plt.close()

    pivot = (
        df.pivot_table(
            index="layerID",
            columns="fiberID",
            values="eventID",
            aggfunc="count",
            fill_value=0,
        )
        .sort_index()
        .sort_index(axis=1)
    )

    plt.figure(figsize=(9.0, 3.8))
    plt.imshow(pivot.values, aspect="auto", origin="lower")
    plt.colorbar(label="Hit rows")
    plt.xlabel("Fibre ID")
    plt.ylabel("Layer ID")
    plt.title("Occupancy map")
    plt.tight_layout()
    plt.savefig(outdir / "occupancy_map.png", dpi=300)
    plt.close()

    pivot.to_csv(outdir / "occupancy_map.csv")


def clusterAnalysis(valid, outdir):
    rows = []

    for eventId, group in valid.groupby("eventID"):
        weights = group["npe_total"].to_numpy()
        if np.sum(weights) <= 0:
            continue

        tReco = np.average(group["t_reco_ns"], weights=weights)
        tTrue = np.average(group["t_true_ns"], weights=weights)
        zReco = np.average(group["z_reco_mm"], weights=weights)
        zTrue = np.average(group["z_mm"], weights=weights)

        rows.append(
            {
                "eventID": eventId,
                "n_hits": len(group),
                "n_layers": group["layerID"].nunique(),
                "npe_total": group["npe_total"].sum(),
                "t_cluster_reco_ns": tReco,
                "t_cluster_true_ns": tTrue,
                "dt_cluster_ps": (tReco - tTrue) * 1000.0,
                "z_cluster_reco_mm": zReco,
                "z_cluster_true_mm": zTrue,
                "dz_cluster_um": (zReco - zTrue) * 1000.0,
            }
        )

    clusters = pd.DataFrame(rows)
    clusters.to_csv(outdir / "cluster_summary.csv", index=False)

    if len(clusters) > 0:
        saveHist(
            clusters["dt_cluster_ps"],
            bins=80,
            xlabel=r"Cluster $t_{\rm reco}-t_{\rm true}$ [ps]",
            ylabel="Events",
            title="Cluster timing",
            outpath=outdir / "cluster_timing_residual.png",
        )
        saveHist(
            clusters["dz_cluster_um"],
            bins=80,
            xlabel=r"Cluster $z_{\rm reco}-z_{\rm true}$ [$\mu$m]",
            ylabel="Events",
            title="Cluster position",
            outpath=outdir / "cluster_position_residual.png",
        )
        saveHist(
            clusters["n_hits"],
            bins=np.arange(0.5, clusters["n_hits"].max() + 1.5, 1),
            xlabel="Number of valid fibre hits in event",
            ylabel="Events",
            title="Hits per event",
            outpath=outdir / "cluster_hit_multiplicity.png",
        )
        saveHist(
            clusters["n_layers"],
            bins=np.arange(0.5, 4.5, 1),
            xlabel="Number of layers with valid hits",
            ylabel="Events",
            title="Layers per event",
            outpath=outdir / "cluster_layer_multiplicity.png",
        )

    return clusters


def main():
    parser = argparse.ArgumentParser(description="Analyze SciFi CSV output.")
    parser.add_argument(
        "csv",
        help="Input CSV, e.g. scifi_hits.csv",
    )
    parser.add_argument(
        "--outdir",
        default="results",
        help="Output folder for plots and CSVs",
    )
    parser.add_argument(
        "--min-npe-side",
        type=float,
        default=3.0,
        help="Min Npe per side",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.csv)

    required = [
        "eventID",
        "layerID",
        "fiberID",
        "copyNo",
        "edep_keV",
        "evis_keV",
        "x_um",
        "y_um",
        "z_mm",
        "t_true_ns",
        "npe_left",
        "npe_right",
        "npe_total",
        "t_left_ns",
        "t_right_ns",
        "t_reco_ns",
        "z_reco_mm",
        "dt_ps",
        "dz_um",
        "valid",
    ]

    requireColumns(df, required)

    for column in required:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(subset=["eventID", "edep_keV", "npe_left", "npe_right"])

    validMask = (
        (df["valid"] == 1)
        & (df["npe_left"] >= args.min_npe_side)
        & (df["npe_right"] >= args.min_npe_side)
        & np.isfinite(df["dt_ps"])
        & np.isfinite(df["dz_um"])
        & (df["t_reco_ns"] > -100)
    )

    valid = df[validMask].copy()
    valid["npe_total"] = valid["npe_left"] + valid["npe_right"]

    print("\nSciFi analysis")
    print(f"Input file: {args.csv}")
    print(f"Total rows: {len(df)}")
    print(f"Valid rows: {len(valid)}")
    print(f"Valid fraction: {len(valid) / len(df):.6f}")

    if len(valid) == 0:
        raise RuntimeError(
            "No valid hits passed. Lower --min-npe-side."
        )

    saveHist(
        df["edep_keV"],
        bins=80,
        xlabel="Deposited energy [keV]",
        ylabel="Rows",
        title="Energy deposition",
        outpath=outdir / "energy_deposition.png",
    )

    saveHist(
        df["evis_keV"],
        bins=80,
        xlabel="Visible energy [keV]",
        ylabel="Rows",
        title="Visible energy",
        outpath=outdir / "visible_energy.png",
    )

    saveHist(
        valid["npe_left"],
        bins=80,
        xlabel=r"Left side photoelectrons $N_{\rm pe,L}$",
        ylabel="Valid rows",
        title="Left Npe",
        outpath=outdir / "npe_left.png",
    )

    saveHist(
        valid["npe_right"],
        bins=80,
        xlabel=r"Right side photoelectrons $N_{\rm pe,R}$",
        ylabel="Valid rows",
        title="Right Npe",
        outpath=outdir / "npe_right.png",
    )

    saveHist(
        valid["npe_total"],
        bins=80,
        xlabel=r"Total photoelectrons $N_{\rm pe,L}+N_{\rm pe,R}$",
        ylabel="Valid rows",
        title="Total Npe",
        outpath=outdir / "npe_total.png",
    )

    saveHist(
        valid["dt_ps"],
        bins=100,
        xlabel=r"$t_{\rm reco}-t_{\rm true}$ [ps]",
        ylabel="Valid rows",
        title="Timing residual",
        outpath=outdir / "timing_residual.png",
    )

    saveHist(
        valid["dz_um"],
        bins=100,
        xlabel=r"$z_{\rm reco}-z_{\rm true}$ [$\mu$m]",
        ylabel="Valid rows",
        title="Position residual",
        outpath=outdir / "position_residual.png",
    )

    saveScatter(
        valid["npe_total"],
        valid["dt_ps"],
        xlabel=r"Total photoelectrons $N_{\rm pe,L}+N_{\rm pe,R}$",
        ylabel=r"$t_{\rm reco}-t_{\rm true}$ [ps]",
        title="Timing vs Npe",
        outpath=outdir / "dt_vs_npe.png",
    )

    saveScatter(
        valid["z_mm"],
        valid["dt_ps"],
        xlabel="True hit position along fibre z [mm]",
        ylabel=r"$t_{\rm reco}-t_{\rm true}$ [ps]",
        title="Timing vs z",
        outpath=outdir / "dt_vs_z.png",
    )

    saveScatter(
        valid["z_mm"],
        valid["dz_um"],
        xlabel="True hit position along fibre z [mm]",
        ylabel=r"$z_{\rm reco}-z_{\rm true}$ [$\mu$m]",
        title="Position vs z",
        outpath=outdir / "dz_vs_z.png",
    )

    thresholdScanDf = thresholdScan(df, outdir)
    npeScanDf = resolutionVsNpe(valid, outdir)
    occupancyPlots(df, outdir)
    clusters = clusterAnalysis(valid, outdir)

    valid.to_csv(outdir / "valid_hits.csv", index=False)

    print("\nResults")
    print(f"Timing sigma:        {valid['dt_ps'].std():.2f} ps")
    print(f"Timing robust:       {robustSigma(valid['dt_ps']):.2f} ps")
    print(f"Position sigma:      {valid['dz_um'].std():.2f} um")
    print(f"Position robust:     {robustSigma(valid['dz_um']):.2f} um")
    print(f"Mean Npe:            {valid['npe_total'].mean():.2f}")
    print(f"Saved to:            {outdir}")


if __name__ == "__main__":
    main()