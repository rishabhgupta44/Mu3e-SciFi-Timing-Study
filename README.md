# Mu3e SciFi Timing Study

**Reference:** This project takes its detector motivation from the Mu3e SciFi detector discussion in arXiv:2501.14692. This repository is my own attempt to learn the detector response and timing reconstruction problem using Geant4 and Python analysis.

This repository contains a Geant4 study of scintillating fibre timing response inspired by the Mu3e SciFi detector.

The aim is to understand how a charged particle passing through thin scintillating fibres produces a timing signal. The simulation uses a three layer fibre ribbon with 250 micrometre fibre scale geometry. It records energy deposited in the fibres and then applies a detector response model to estimate light production, photoelectrons at both SiPM ends, and reconstructed hit time.

## Detector response chain

The basic response chain is:

```text
charged particle
→ energy deposited in scintillating fibre
→ visible scintillation energy
→ scintillation photons
→ photoelectrons at left and right SiPM ends
→ left and right signal times
→ reconstructed hit time and position along the fibre
```

In short:

```text
e± → Edep → Evis → photons → photoelectrons → timing signals → reconstructed time
```

## What this project includes

* Three layer scintillating fibre ribbon geometry in Geant4
* 250 micrometre fibre scale geometry
* Electron or positron beam through the detector
* Energy deposition recorded in each fibre
* Light production and SiPM response model
* Light attenuation along the fibre
* Poisson fluctuations in photoelectron production
* Time walk and timing smearing model
* Timing reconstruction using signals from both fibre ends
* CSV output from Geant4
* Python scripts for detector response plots

## What the analysis studies

The Python analysis reads the Geant4 output and produces plots for:

* energy deposition
* visible energy after Birks type correction
* photoelectron yield
* hit timing residual
* timing resolution versus photoelectron yield
* efficiency versus threshold
* fibre and layer occupancy
* cluster level timing studies

## Build and run

You need Geant4 11 installed and available to CMake.

Build the project:

```bash
cmake -S . -B build
cmake --build build -j
```

Run the simulation:

```bash
./build/scifi_sim macros/run.mac
```

This creates:

```text
scifi_hits.csv
```

## Python analysis

Install the Python requirements:

```bash
python3 -m pip install -r analysis/requirements.txt
```

Run the analysis:

```bash
python3 analysis/analyze_scifi.py scifi_hits.csv --outdir results
```

The plots and summary files will be saved in:

```text
results/
```

## Visualization

To open the detector visualization:

```bash
./build/scifi_sim
```

To run a live event display:

```bash
./build/scifi_sim macros/live.mac
```

To save a detector image:

```bash
./build/scifi_sim macros/save_geometry.mac
```

To make a short setup movie:

```bash
./scripts/make_setup_movie.sh
```

## Note

It is my attempt to learn how a scintillating fibre timing detector works, starting from energy deposition in Geant4 and ending with timing reconstruction and detector response plots in Python.

Future improvements may include optical photon tracking, fibre cladding, optical surfaces, SiPM dark counts, cross talk, electronics response, and accidental hit mixing.
