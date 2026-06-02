# Useful Tools for JWST MIRI MRS

A curated list of tools, references, and packages that are useful when working with **JWST Mid-Infrared Instrument Medium-Resolution Spectrometer (MIRI MRS)** data — and adjacent IFU/spectral-cube datasets. Most are also applicable to NIRSpec IFU.

---

## Reference documentation

### JWST Documentation (JDox)

<https://jwst-docs.stsci.edu/>

The official STScI documentation hub for all JWST instruments. The MIRI MRS pages cover observing modes, the data reduction pipeline, calibration status, known issues, and recommended strategies for cube building and spectral extraction. This should be the first stop for any instrument-specific question.

### JWST pipeline notebooks

<https://github.com/spacetelescope/jwst-pipeline-notebooks>

Official Jupyter notebooks from STScI that walk through running the `jwst` calibration pipeline end-to-end for each instrument and mode. The MIRI MRS notebooks demonstrate Stage 1 (detector ramp fitting), Stage 2 (per-exposure calibration), and Stage 3 (cube building and 1D spectral extraction), and are the most up-to-date practical reference for pipeline usage.

---

## Data acquisition

All JWST data — raw uncalibrated, intermediate, and pipeline-processed — is distributed through STScI's **MAST** archive. The portals below are for interactive browsing; the Python tools further down are for scripted, filtered, and bulk downloads.

### MAST Portal

<https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html>

The **Mikulski Archive for Space Telescopes**, hosted by STScI, is the primary archive for all JWST data (alongside HST, Kepler, TESS, and others). The web portal lets you search by program ID, target, coordinates, or filters; preview observations on a sky map; and download individual files or staged bundles. Best when you know roughly what you want and want to inspect a handful of observations interactively.

### ESA Sky

<https://sky.esa.int/>

A discovery-first sky-visualization portal from ESA that overlays footprints from many missions (HST, JWST, XMM, Gaia, ALMA, …) onto a multiwavelength sky background. Useful when you want to see *what has been observed* in a given region of sky — across missions — rather than search MAST by program. JWST observations flow through to ESA Sky once they appear in MAST. From a hit you can deep-link back to MAST to download the actual files.

### jwst_mast_query (programmatic, JWST-specific)

<https://github.com/spacetelescope/jwst_mast_query>

STScI's purpose-built command-line and Python tool for filtering and downloading JWST observations from MAST. It understands JWST-specific concepts (proposal/program/observation IDs, calibration levels 1/2/3, product types, pipeline versions) and is the most efficient route for **bulk downloads, fetching raw `_uncal.fits` products for reprocessing, or pulling everything matching a given filter** without clicking through the portal one observation at a time.

### astroquery.mast (programmatic, general)

<https://astroquery.readthedocs.io/en/latest/mast/mast.html>

The community-maintained `astroquery` library exposes MAST via a Pythonic interface (`astroquery.mast.Observations`). More general-purpose than `jwst_mast_query` — better when you want to write custom filtering logic, cross-match with non-MAST catalogs in the same script, or integrate MAST queries into a larger pipeline.

> **MAST account.** Downloading proprietary data (within the program's 12-month exclusive-access period) requires a [MyST account and authentication token](https://auth.mast.stsci.edu/info). Public data does not require authentication. The same token is used by both `jwst_mast_query` and `astroquery.mast`.

---

## Visualization

### CARTA

[CARTA](https://cartavis.org/) (Cube Analysis and Rendering Tool for Astronomy) is the standard interactive viewer for radio/sub-mm spectral cubes from ALMA, VLA, MeerKAT, and SKA pathfinders. Although originally developed for radio data, it is also an excellent viewer for **JWST integral-field cubes from NIRSpec IFU and MIRI/MRS** — its handling of large 3D datasets, channel maps, and spectral profile extraction translates directly. See the [workshop packages](installation_instructions/workshop_packages.md#1-carta-spectral-cube-viewer) guide for setup instructions.

---

## Molecular line analysis

### iSLAT

<https://github.com/spexod/iSLAT>

The **Interactive Spectral-Line Analysis Tool**, developed by the SpExoDisks collaboration. iSLAT provides an interactive GUI for identifying and fitting molecular emission lines in mid-infrared spectra, with built-in LTE slab models for common species (H₂O, CO, OH, CO₂, HCN, C₂H₂, organics). It is particularly well suited to **MIRI MRS spectra of protoplanetary disks and other warm molecular environments**, where the dense forest of ro-vibrational lines makes by-eye identification impractical. Line lists are taken from HITRAN and complementary databases.

### DuCKLinG

<https://github.com/tillkaeufer/DuCKLinG>

**DuCKLinG** (Dust Continuum Kit with Line emission from Gas) fits the dust continuum and molecular slab-model line emission *simultaneously* in a single forward model, rather than continuum-subtracting first and fitting slab models afterwards. This is particularly relevant for **MIRI MRS spectra of protoplanetary disks** where the silicate feature and the molecular line forest overlap in wavelength and the inferred slab parameters can depend sensitively on how the continuum was drawn. A natural next step beyond the iSLAT-style "subtract continuum, then fit lines" workflow when you want statistically defensible parameter posteriors that account for continuum/line degeneracies.

### NIST Atomic Spectra Database

<https://physics.nist.gov/PhysRefData/ASD/lines_form.html>

The authoritative reference database for atomic line wavelengths, energy levels, transition probabilities, and oscillator strengths, maintained by the National Institute of Standards and Technology. Useful for identifying atomic and ionic features (e.g. [Ne II] 12.81 µm, [Ne III] 15.55 µm, [S III], [Fe II]) that appear alongside the molecular forest in MIRI MRS spectra of ionized regions, jets, and AGN.

---

## Extended emission and IFU decomposition

### q3dfit

<https://q3dfit.readthedocs.io/>

A Python package for fitting JWST and ground-based IFU data, descended from the IDL `IFSFIT` code originally developed for quasar host-galaxy decomposition. q3dfit simultaneously models point-source (e.g. AGN/quasar) emission and the underlying spatially extended host-galaxy continuum and emission lines on a spaxel-by-spaxel basis. It is designed for **NIRSpec IFU and MIRI MRS cubes of AGN, quasars, and other systems with a bright unresolved component superimposed on extended emission**, and handles PSF subtraction, multi-component line fitting, and continuum modeling in a unified framework.

---

## Contributing

If you know of another tool that belongs here (line identification, cube manipulation, PSF modeling, spectral fitting, etc.), feel free to add it with a short description following the format above.
