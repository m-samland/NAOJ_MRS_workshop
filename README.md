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
- **Prerequisites for participants:** working laptop, environment set up per the installation guides above, >20 GB free disk space for example data. If you want to try to reduce raw data from scratch, this needs about 30 GB for a target.

Hands-on materials (pipeline notebook, extended-emission analysis notebook, scripts) will be added to this repository closer to the workshop date.

# Hands-on session: MIRI Medium Resolution Spectrograph

The hands-on session is built around a single linear workflow that takes you from a calibrated archive product all the way to a molecular-line identification on a real protoplanetary disk spectrum. The default target is **SY Cha**, but the steps are identical for other targets. As a start do all steps for SY Cha to be on familiar ground, if time allows start doing the same for your favourite target.

### Step 1. Data acquisition (download Stage-3 products from MAST)

- Go to the MAST portal and search for MIRI MRS observations of your target. Use **SY Cha**. In the advanced search of MAST, filter for MIRI/IFU and the target coordinates, then narrow it down further.
- Read the proposal PDF for context — observing mode, dither pattern, exposure times, and the science goals of the original program. It can be acquired by clicking on "..." and then "Show detail" for a data product in MAST. Look for Proposal ID. Then google "JWST proposal id XXX" or go directly to: https://www.stsci.edu/jwst-program-info/program/?program=XXX
- Download the **Level-3 reduced science products** (`*_s3d.fits` cubes and `*_x1d.fits` 1D extractions). These are the calibration-pipeline outputs hosted by STScI; you do **not** need to run the pipeline yourself. There should be 12 x1d and s3d files, one for each MRS band.
- Collect everything for your target into a single folder (the workshop layout is `data/<target>/stage3_mast/`). The example data shipped in `data/SYCha/stage3_mast/` follows exactly this convention — you can mirror it. If MAST is slow, you can download via the link provided on this github (see README.md)

### Step 2. Inspect the data in CARTA

- Open the `*_s3d.fits` cubes in CARTA to confirm how the data looks: source centred in roughly in the FOV, no obvious artefacts. Go to Widgets -> Profiles and select "Spectral Profiler" to get a spectrum when you hover over the image.
- Look for science features already visible by eye: molecular line forests, atomic lines (e.g. [Ne II], [Ar II]), the broad silicate feature around 10 µm, extended emission around the central source.
- This is also a good moment to get a feel for the wavelength coverage of each MRS sub-band (1A → 4C) and how they tile in wavelength.

### Step 3. Stitch the x1d spectra into a single 1D spectrum

- Run [`MRSReduction/combine_x1d_spectra.py`](../MRSReduction/combine_x1d_spectra.py) pointing it at the folder of `*_x1d.fits` files you downloaded.
- The script rescales the 12 sub-bands so they agree in their wavelength overlap regions and writes a stitched, non-overlapping spectrum (1A → 4C) as `<source>_spectrum_full.csv` and `<source>_spectrum_stitched.csv`, plus a diagnostic PNG.
- Inspect the diagnostic plot — if sub-bands disagree strongly, that usually points to a calibration issue or a problem with the extraction aperture rather than a bug in the stitching.

### Step 4. Continuum-subtract the stitched spectrum

- Run [`MRSReduction/continuum_subtraction.py`](../MRSReduction/continuum_subtraction.py) (an implementation of the Temmink et al. 2024 algorithm) on the stitched spectrum from Step 3.
- This is a required preparation step before any molecular-line work in iSLAT — slab models are fit to the line-only spectrum, not the continuum-dominated one.
- Look at the diagnostic plots and confirm the continuum is sensibly traced: not over-fitting line forests, not under-fitting the silicate feature. If something looks wrong, adjust the parameters and rerun.

### Step 5. Molecular-line analysis with iSLAT

- One of the outputs of the stitched spectrums has "for_i"Rename the columns of the continuum-subtracted CSV so they match the iSLAT input format.
- Open the spectrum in iSLAT (separate `iSLAT` env — it's a standalone GUI).
- Overplot LTE slab models for different molecular species (H₂O, CO, CO₂, HCN, C₂H₂, OH, …) at a range of temperatures and column densities to identify which species are present and to get a first feel for their excitation.

### (Optional) Step 6. Search for extended line emission

- Open [`extended_emission/extended_structures.ipynb`](../extended_emission/extended_structures.ipynb).
- The notebook uses the curated `linelist_combined.csv` to scan the Level-3 cubes for spatially extended emission in known atomic and molecular transitions (jets, outflows, photoevaporative winds, …) — complementary to the spatially-integrated spectrum from Steps 3–5.

### (If there is time) Step 7. Repeat on a new, "unknown" object

- Pick another MIRI MRS target from MAST that you have **not** looked at yet, ideally a different object class (e.g. a debris disk vs. a Class II disk, or an evolved star).
- Repeat Steps 1 → 5 on the new data, this time without the workshop hand-holding. The point is to confirm the workflow is yours, not just SY Cha's, and to surface the questions that come up when you can no longer compare against a known answer.

### (Advanced, if disk + bandwidth allow) Step 8. Run the calibration pipeline from raw data

- Download the **raw Stage-0 uncalibrated products** for your target from MAST instead of the Level-3 ones (note: this can be tens of GB and takes a while over a workshop wifi connection — be realistic about whether to attempt it live).
- Run [`MRSReduction/simplified_pipeline_script.py`](../MRSReduction/simplified_pipeline_script.py) (or step through `JWPipeNB-MIRI-MRS.ipynb`) to reduce the raw data yourself through Stages 1–3.
- Compare your self-reduced `*_x1d.fits` outputs against the MAST Level-3 products. This is the easiest way to internalise what each pipeline stage actually does, and how robust (or not) the defaults are for your science case.