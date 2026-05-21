# Workshop packages

This guide installs the science packages used during the **NAOJ JWST MIRI-MRS hands-on session**. It assumes you have already followed the OS-specific environment setup in [macos_linux.md](macos_linux.md) or [windows.md](windows.md) — i.e. you have a working `mamba`/`conda` install and the `astro` Python environment.

Each package is installed independently; you can skip any you do not need.

| Package / env | What it provides | Used by |
| --- | --- | --- |
| [CARTA](#1-carta-spectral-cube-viewer) | Standalone interactive viewer for 3D spectral cubes | All cube inspection (MRS, NIRSpec IFU, ALMA) |
| [`mrs_analysis` env](#2-mrs_analysis-environment-jwst-pipeline--analysis) | JWST pipeline + extended-emission / continuum-fitting deps | `MRSReduction/`, `extended_emission/`, `continuum_subtraction/` |
| [`iSLAT` env](#3-islat-environment-molecular-line-analysis) | Interactive molecular line identification and slab-model fitting | MIRI MRS molecular-line analysis sessions |

The workshop uses **two separate mamba environments** so the JWST pipeline's dependency set does not collide with iSLAT's tighter version pins — see the note in [§3](#3-islat-environment-molecular-line-analysis) for the rationale.

---

## 1. CARTA (spectral cube viewer)

[CARTA](https://cartavis.org/) (Cube Analysis and Rendering Tool for Astronomy) is the standard interactive viewer for radio/sub-mm spectral cubes from ALMA, VLA, MeerKAT, and SKA pathfinders. Although originally developed for radio data, it is also an excellent viewer for **JWST integral-field cubes from NIRSpec IFU and MIRI/MRS** — its handling of large 3D datasets, channel maps, and spectral profile extraction translates directly.

CARTA is a standalone desktop application — it is **not** installed into your Python environment.

### CARTA on macOS

Supported on macOS 14 (Sonoma) and 15 (Sequoia), Intel and Apple Silicon. The cleanest install is via Homebrew:

```bash
brew install --cask cartavis/tap/carta
```

If you don't use Homebrew, download the DMG directly:

- [Apple Silicon DMG](https://github.com/CARTAvis/carta/releases/latest/download/CARTA-arm64.dmg)
- [Intel DMG](https://github.com/CARTAvis/carta/releases/latest/download/CARTA-x64.dmg)

Open the DMG and drag CARTA to `/Applications`.

### CARTA on Linux

**Ubuntu 22.04 / 24.04** — install from the project's PPA:

```bash
sudo add-apt-repository ppa:cartavis-team/carta
sudo apt update
sudo apt install carta
```

**RHEL / Rocky / AlmaLinux 8 / 9** — install via Fedora Copr:

```bash
sudo dnf install 'dnf-command(copr)'
sudo dnf copr enable cartavis/carta
sudo dnf install epel-release carta
```

**Other distros** — use the official AppImage:

- [AppImage x86_64](https://github.com/CARTAvis/carta/releases/latest/download/carta.AppImage.x86_64.tgz)
- [AppImage aarch64](https://github.com/CARTAvis/carta/releases/latest/download/carta.AppImage.aarch64.tgz)

Untar, `chmod +x` the AppImage, and run it directly.

### CARTA on Windows (via WSL)

CARTA does not officially support Windows. The supported path is to install it inside WSL's Ubuntu environment via the project's PPA; WSLg will display the GUI in Windows just like any other Linux app.

Inside your Ubuntu (WSL) terminal:

```bash
sudo add-apt-repository ppa:cartavis-team/carta
sudo apt update
sudo apt install carta
```

The PPA officially supports Ubuntu 22.04 and 24.04 — the defaults shipped by `wsl --install` are fine.

If you prefer not to install via the PPA, the official Docker image is the next-best option (`docker run --rm -ti -p 3002:3002 -v $PWD:/images cartavis/carta:latest`); see the [CARTA download page](https://cartavis.org/#download) for details.

### Using CARTA

Launch with:

```bash
carta
# or open a cube directly:
carta /path/to/cube.fits
```

On Windows/WSL the GUI window opens directly on your Windows desktop via WSLg. The full download page with Docker and source-build options is at [cartavis.org](https://cartavis.org/#download).

---

## 2. `mrs_analysis` environment (JWST pipeline + analysis)

This environment holds the JWST calibration pipeline plus the extra packages needed by the workshop's reduction, continuum-subtraction, and extended-emission notebooks. Specifically:

- [`MRSReduction/JWPipeNB-MIRI-MRS.ipynb`](../MRSReduction/JWPipeNB-MIRI-MRS.ipynb) and [`MRSReduction/simplified_pipeline_script.py`](../MRSReduction/simplified_pipeline_script.py) need the **`jwst`** pipeline (which itself brings in `numpy`, `scipy`, `astropy`, `matplotlib`, `asdf`, `gwcs`, `stdatamodels`, `crds`, `photutils`), plus `astroquery` (MAST downloads) and a Jupyter kernel.
- [`continuum_subtraction/ContinuumSubtraction.py`](../continuum_subtraction/ContinuumSubtraction.py) adds **`pybaselines`** and **`spectres`** on top of the standard numpy/scipy/astropy/pandas/matplotlib stack.
- [`extended_emission/extended_structures_cleaner_v2.4.ipynb`](../extended_emission/extended_structures_cleaner_v2.4.ipynb) adds **`vip_hci`**, **`seaborn`**, **`tqdm`**, and **`packaging`**.

All of these are mutually compatible: the JWST pipeline (≥ 1.18) supports Python 3.11–3.14 with loose `numpy>=1.25`, `scipy>=1.14.1`, `astropy>=6.1` constraints, and every additional dependency has wheels for Python 3.12 and supports current numpy 2.x.

### Create the env

```bash
mamba create -n mrs_analysis python=3.12
mamba activate mrs_analysis
pip install -U pip
```

> **Why Python 3.12 and not 3.13?** Some downstream packages (notably `vip_hci`) lag a release behind, so 3.12 is the safest sweet spot for a workshop where everyone must end up with the same working env. Bump this once `vip_hci` and friends officially ship for 3.13.

### Install the JWST pipeline

```bash
pip install jwst
```

This pulls in the bulk of the scientific stack — `numpy`, `scipy`, `astropy`, `matplotlib`, `asdf`, `gwcs`, `stdatamodels`, `crds`, `photutils`, and friends — so they should not be installed separately first.

### Install the additional analysis packages

```bash
pip install astroquery pybaselines spectres vip-hci seaborn tqdm packaging \
            ipykernel jupyterlab
```

`pip` will skip anything `jwst` already pulled in. `astroquery` and `photutils` may both already be present; that's fine.

> **Why `jupyterlab` rather than the `jupyter` / `notebook` meta-package?** JupyterLab is the actively developed interface from Project Jupyter; the classic single-document Notebook 7 is now built on JupyterLab's own components, so installing `jupyterlab` gets you the modern multi-pane UI and full notebook compatibility in one package. If you intend to run all notebooks from VS Code, you can drop `jupyterlab` from the line above — only `ipykernel` is required to make the env show up in VS Code's kernel picker.

### Register the env as a Jupyter kernel

So that VS Code and JupyterLab can find this environment from the kernel picker:

```bash
python -m ipykernel install --user --name mrs_analysis --display-name "Python (mrs_analysis)"
```

### Configure CRDS

The JWST pipeline needs to download reference files from STScI's CRDS server. Tell it where to cache them and which server to talk to by adding the following two lines to the rc file of whichever shell you use (`~/.zshrc` on macOS, `~/.bashrc` on most Linux / WSL):

```bash
export CRDS_PATH=$HOME/crds_cache
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
```

Then either `source` the rc file or open a new terminal. The first pipeline run will populate the cache; subsequent runs reuse it.

> **Disk usage.** The CRDS cache grows to several GB once you've run the MRS pipeline end-to-end. Point `CRDS_PATH` somewhere with room (an external drive is fine — just keep it on a filesystem with normal POSIX semantics, not a cloud-synced folder).

---

## 3. `iSLAT` environment (molecular line analysis)

[iSLAT](https://github.com/spexod/iSLAT) — the **Interactive Spectral-Line Analysis Tool**, developed by the SpExoDisks collaboration — provides an interactive GUI for identifying and fitting molecular emission lines in mid-infrared spectra, with built-in LTE slab models for common species (H₂O, CO, OH, CO₂, HCN, C₂H₂, organics). It is the workhorse tool for the molecular-line analysis session of the workshop.

iSLAT pulls in a relatively narrow set of dependencies (`numpy`, `scipy`, `astropy`, `pandas`, `lmfit`, `astroquery`, `matplotlib`, `tk`) but is **kept in its own environment, separate from `mrs_analysis`**, because:

1. iSLAT pins fairly recent minimum versions (e.g. `numpy>=2.3.5`, `scipy>=1.16.3`, `pandas>=2.3.3`, `matplotlib>=3.10.7`) that may not stay in lockstep with whatever `jwst` resolves to.
2. iSLAT installs from its own `requirements.txt`, so a single `pip install -r` against a shared env could downgrade or upgrade packages that the JWST pipeline depends on.
3. iSLAT is launched as a standalone GUI; nothing in the workshop notebooks imports from it, so there is no integration benefit to sharing an env.

### Create a dedicated `iSLAT` environment

```bash
mamba create -n iSLAT python=3.13
mamba activate iSLAT
```

Install `tk` from conda-forge (it ships with Python on most systems but is occasionally missing or broken in fresh conda envs, and iSLAT's GUI depends on it):

```bash
mamba install -c conda-forge tk
```

> **Windows / WSL note.** Run all of the steps below inside your WSL Ubuntu terminal, not in native Windows. The iSLAT GUI will display on your Windows desktop via WSLg, the same way CARTA does.

### Clone and install iSLAT from source

Pick a directory where you'd like to keep iSLAT (e.g. `~/code/` or `~/projects/`), then:

```bash
git clone https://github.com/spexod/iSLAT
cd iSLAT
pip install -r requirements.txt
```

This installs the Python dependencies into the active `iSLAT` env. Optionally, install iSLAT itself as a package (recommended — gives you the `iSLAT` launcher on your PATH):

```bash
pip install -e .
```

### Launch iSLAT

If you ran the optional `pip install -e .` step:

```bash
mamba activate iSLAT
iSLAT
```

Otherwise, launch via the script in the repo:

```bash
mamba activate iSLAT
cd /path/to/iSLAT
python iSLAT-launch.py
```

### Updating iSLAT

iSLAT is actively developed — pull the latest version periodically:

```bash
cd /path/to/iSLAT
git pull
pip install -r requirements.txt   # in case dependencies changed
```

---

## Troubleshooting

**`pip install jwst` fails on a dependency build.** Make sure you used `python=3.12` (not 3.13 yet) and that you ran `pip install -U pip` first. If the error mentions a missing C compiler, install build tools per the OS guide ([macOS/Linux](macos_linux.md#1-install-build-tools) or the WSL section of [windows.md](windows.md#2-install-build-tools-inside-wsl)).

**Pipeline run fails with `CrdsError` / "no reference file found".** The `CRDS_PATH` and `CRDS_SERVER_URL` variables are not set in the shell that launched the notebook. Re-check your `~/.zshrc` / `~/.bashrc`, then restart the kernel (in Jupyter / VS Code: *Kernel → Restart*).

**VS Code's kernel picker doesn't show `Python (mrs_analysis)`.** Re-run `python -m ipykernel install --user --name mrs_analysis --display-name "Python (mrs_analysis)"` with the env active, then reload VS Code (`Cmd/Ctrl+Shift+P` → *Developer: Reload Window*).

**`mamba activate iSLAT` followed by `iSLAT` says "command not found".** You skipped the `pip install -e .` step. Either run it from inside the iSLAT repo (with the env active), or launch via `python iSLAT-launch.py` from the repo directory.

**iSLAT GUI opens but throws a `tkinter` / `tk` error.** Install `tk` explicitly into the env: `mamba install -c conda-forge tk`. On WSL, also make sure WSLg is working (test by launching CARTA or any other GUI app).

**CARTA on WSL fails to open a window.** Run `wsl --update` from PowerShell to refresh the WSL kernel, then reopen the WSL terminal. WSLg requires a recent Windows 10 or Windows 11.
