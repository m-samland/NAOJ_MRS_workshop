# NAOJ MIRI MRS Workshop

Materials for a one-day hands-on workshop on **JWST MIRI MRS data analysis**, taught by Matthias Samland at NAOJ (National Astronomical Observatory of Japan).

The workshop walks participants through the practical end-to-end workflow for working with JWST Mid-Infrared Instrument Medium-Resolution Spectrometer (MIRI MRS) integral-field data — from environment setup, through running the calibration pipeline, to scientific analysis of molecular lines and extended emission.

A basic working knowledge of Python is assumed; no prior JWST-specific experience is required.

---

## Before the workshop: install everything

Setup is a **two-step process**. Please complete both before the in-person session — installs (especially the JWST pipeline and CARTA) can take a while and are easier to debug at home than in a lecture room.

### Step 1 — Set up a scientific Python environment

Pick the guide for your operating system:

- **macOS / Linux:** [installation_instructions/macos_linux.md](installation_instructions/macos_linux.md)
- **Windows:** [installation_instructions/windows.md](installation_instructions/windows.md) (routes through WSL2; several workshop tools are Unix-only)

These install Miniforge (conda + mamba), create a baseline `astro` Python environment, and set up VS Code and git.

### Step 2 — Install the workshop packages

Once the OS-level environment is in place, follow:

- [installation_instructions/workshop_packages.md](installation_instructions/workshop_packages.md)

This installs **CARTA** (interactive spectral-cube viewer), creates an **`mrs_analysis`** mamba environment with the JWST pipeline and the analysis dependencies used during the workshop, and creates a separate **`iSLAT`** environment for the interactive molecular-line analysis GUI.

### Step 3 — Download the data

Downloading the data we will use as an example here: https://datashare.mpcdf.mpg.de/s/5W8PB89Yea7BWjY

If you're interested in a particular target that you want to analyze after we are done with the base tutorial best think about what you would like to look at in advance! Check the MAST archive if there is any MIRI/IFU data for it.

--

## Reference material

- [useful_tools.md](useful_tools.md) — a curated overview of tools, archives, and packages relevant to MIRI MRS work, including pointers to the official JWST documentation (JDox), MAST/ESA Sky archive access, CARTA, iSLAT, NIST ASD, and q3dfit. Useful both during the workshop and as a take-home reference afterwards.

---

## What you'll do during the hands-on session

The hands-on portion walks through:

1. **Data acquisition** — find MRS observations of SY Cha on MAST, check the proposal PDF for context, and download the reduced `s3d` (cube) and `x1d` (1D spectrum) products into a single folder.
2. **Inspect the cube** — open it in CARTA and look for obvious features: molecular bands, atomic lines, silicate emission/absorption.
3. **Continuum subtraction** — run the Temmink et al. (2024) continuum subtraction on the 1D spectrum and verify the fit visually; this step is required before line searches.
4. **Molecular-line analysis with iSLAT** — reformat the continuum-subtracted spectrum for iSLAT and overplot LTE slab models for different molecular species and temperatures to identify what's there.
5. *(Optional)* **Extended line emission** — use the `extended_structures` notebook to map spatially extended atomic and molecular emission in the cube.

Bring a target you care about if you have one — otherwise SY Cha works well as a default.

---

## Workshop logistics

- **Format:** one day, in person at NAOJ.
- **Prerequisites for participants:** working laptop, environment set up per the installation guides above, >20 GB free disk space for example data.

Hands-on materials (pipeline notebook, extended-emission analysis notebook, scripts) will be added to this repository closer to the workshop date.
