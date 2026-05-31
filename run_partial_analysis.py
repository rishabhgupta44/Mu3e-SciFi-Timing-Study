#!/usr/bin/env python3
# Copyright (c) Rishabh Gupta 2026
# This is a test code used for learning.
from pathlib import Path
import sys
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import analyze_scifi as a

NROWS = 5_000_000
CSV = "scifi_hits.csv"
OUTDIR = Path("results_5M")
OUTDIR.mkdir(parents=True, exist_ok=True)

print(f"Reading first {NROWS} rows from {CSV}...")

df = pd.read_csv(CSV, nrows=NROWS)

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

a.requireColumns(df, required)

for column in required:
    df[column] = pd.to_numeric(df[column], errors="coerce")

df = df.dropna(subset=["eventID", "edep_keV", "npe_left", "npe_right"])

minNpeSide = 3.0
validMask = (
    (df["valid"] == 1)
    & (df["npe_left"] >= minNpeSide)
    & (df["npe_right"] >= minNpeSide)
    & pd.notna(df["dt_ps"])
    & pd.notna(df["dz_um"])
    & (df["t_reco_ns"] > -100)
)

valid = df[validMask].copy()
valid["npe_total"] = valid["npe_left"] + valid["npe_right"]

print("Running analysis pipeline on the loaded rows...")

print(f"Total rows loaded: {len(df)}")
print(f"Valid rows: {len(valid)}")

try:
    a.saveHist(
        df["edep_keV"], bins=80, xlabel="Deposited energy [keV]", ylabel="Rows",
        title="Energy deposition", outpath=OUTDIR / "energy_deposition.png",
    )

    a.saveHist(
        df["evis_keV"], bins=80, xlabel="Visible energy [keV]", ylabel="Rows",
        title="Visible energy", outpath=OUTDIR / "visible_energy.png",
    )

    a.saveHist(
        valid["npe_left"], bins=80, xlabel=r"Left side photoelectrons $N_{\rm pe,L}$",
        ylabel="Valid rows", title="Left Npe",
        outpath=OUTDIR / "npe_left.png",
    )

    a.saveHist(
        valid["npe_right"], bins=80, xlabel=r"Right side photoelectrons $N_{\rm pe,R}$",
        ylabel="Valid rows", title="Right Npe",
        outpath=OUTDIR / "npe_right.png",
    )

    a.saveHist(
        valid["npe_total"], bins=80, xlabel=r"Total photoelectrons $N_{\rm pe,L}+N_{\rm pe,R}$",
        ylabel="Valid rows", title="Total Npe",
        outpath=OUTDIR / "npe_total.png",
    )

    a.saveHist(
        valid["dt_ps"], bins=100, xlabel=r"$t_{\rm reco}-t_{\rm true}$ [ps]",
        ylabel="Valid rows", title="Timing residual",
        outpath=OUTDIR / "timing_residual.png",
    )

    a.saveHist(
        valid["dz_um"], bins=100, xlabel=r"$z_{\rm reco}-z_{\rm true}$ [$\mu$m]",
        ylabel="Valid rows", title="Position residual",
        outpath=OUTDIR / "position_residual.png",
    )

    a.saveScatter(
        valid["npe_total"], valid["dt_ps"], xlabel=r"Total photoelectrons $N_{\rm pe,L}+N_{\rm pe,R}$",
        ylabel=r"$t_{\rm reco}-t_{\rm true}$ [ps]", title="Timing vs Npe",
        outpath=OUTDIR / "dt_vs_npe.png",
    )

    a.saveScatter(
        valid["z_mm"], valid["dt_ps"], xlabel="True hit position along fibre z [mm]",
        ylabel=r"$t_{\rm reco}-t_{\rm true}$ [ps]", title="Timing vs z",
        outpath=OUTDIR / "dt_vs_z.png",
    )

    a.saveScatter(
        valid["z_mm"], valid["dz_um"], xlabel="True hit position along fibre z [mm]",
        ylabel=r"$z_{\rm reco}-z_{\rm true}$ [$\mu$m]", title="Position vs z",
        outpath=OUTDIR / "dz_vs_z.png",
    )

    thresholdScanDf = a.thresholdScan(df, OUTDIR)
    npeScanDf = a.resolutionVsNpe(valid, OUTDIR)
    a.occupancyPlots(df, OUTDIR)
    clusters = a.clusterAnalysis(valid, OUTDIR)

    valid.to_csv(OUTDIR / "valid_hits.csv", index=False)

    print("Analysis done; results saved to:", OUTDIR)
except Exception as e:
    print("Analysis error:", e)
    raise
