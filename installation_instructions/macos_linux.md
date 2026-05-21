# Setup for macOS and Linux

For scientific Python work in astronomy, macOS and Linux are both first-class platforms — most astronomy packages target Unix-like environments, so almost everything works out of the box.

The setup below is essentially the same on macOS and Linux; differences are called out per section. Apple Silicon (M1/M2/M3/M4) Macs and ARM Linux machines are fully supported.

> **For the NAOJ MRS workshop:** this document covers the general-purpose scientific Python environment. Once that is in place, continue to [workshop_packages.md](workshop_packages.md) to install the workshop-specific tools (CARTA, iSLAT, …).

---

## 1. Install build tools

Some scientific Python packages need a C/C++/Fortran compiler at install time.

### macOS

Install the Xcode Command Line Tools (this also gives you `git`):

```bash
xcode-select --install
```

A graphical prompt will appear; accept it. This is a one-time, multi-GB install. You do *not* need the full Xcode app — only the Command Line Tools.

For Fortran-extension packages (less common, but some MCMC and instrument codes use them), install `gfortran` via Homebrew or via conda-forge later:

```bash
# via Homebrew (if you use it)
brew install gfortran

# or — preferred — get it from conda-forge inside your env (see section 4)
```

### Linux

On Debian/Ubuntu:

```bash
sudo apt update
sudo apt install build-essential gfortran git
```

On Fedora/RHEL/Rocky:

```bash
sudo dnf groupinstall "Development Tools"
sudo dnf install gcc-gfortran git
```

On Arch:

```bash
sudo pacman -S base-devel gcc-fortran git
```

---

## 2. Install Miniforge (conda + mamba)

We recommend **Miniforge** rather than the regular Anaconda/Miniconda distributions, because:

- It defaults to the community-maintained **conda-forge** channel, which has the largest and best-maintained selection of scientific Python packages.
- It avoids Anaconda Inc.'s commercial-use restrictions on the `defaults` channel that affect many institutions.
- It ships with **mamba**, a much faster solver alongside `conda` (this is also what the JWST pipeline recommends).

> **One install, two commands.** Miniforge installs *both* `conda` and `mamba`. They share the same environments, the same config (`~/.condarc`), and the same shell integration — `mamba` is not a separate ecosystem. Use whichever feels natural; `mamba` is just faster.

### Pick the right installer for your machine

| Platform | Installer |
| --- | --- |
| macOS, Apple Silicon (M1/M2/M3/M4) | `Miniforge3-MacOSX-arm64.sh` |
| macOS, Intel | `Miniforge3-MacOSX-x86_64.sh` |
| Linux, x86_64 (most desktops/servers) | `Miniforge3-Linux-x86_64.sh` |
| Linux, ARM64 (Raspberry Pi 4+, ARM servers) | `Miniforge3-Linux-aarch64.sh` |

> **On Apple Silicon, make sure you grab the `arm64` installer.** Using the `x86_64` installer would force every package to run through Rosetta translation — much slower and a frequent source of "this works for everyone else but not me" bugs.

To check your Mac's architecture if unsure: `uname -m` returns `arm64` on Apple Silicon and `x86_64` on Intel.

### Download and install

```bash
# Replace the filename below with the one from the table above.
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh
```

If `wget` isn't installed on macOS, use `curl -L -O <url>` instead, or install it via `brew install wget`.

Accept the license, accept the default install location (`~/miniforge3`), and answer **yes** when asked whether to run `conda init`. Then close and reopen the terminal so the shell picks up the new configuration.

Verify the install:

```bash
conda --version
mamba --version
```

> **Shell note.** `conda init` writes to the rc file of whichever shell you're using. macOS uses **zsh** by default (since Catalina), so the hook lands in `~/.zshrc`. Most Linux distributions use **bash** by default, so it lands in `~/.bashrc`. If you switch shells later, re-run `conda init <shell>` to add the hook to the new rc file.

---

## 3. Create a Python environment

Pick a recent Python (3.13 is a good default in 2026):

```bash
mamba create -n astro python=3.13
mamba activate astro
mamba install numpy scipy matplotlib astropy jupyter
```

`mamba install` resolves environments in seconds rather than the minutes `conda install` typically takes. `mamba` and `conda` are interchangeable in current Miniforge — `mamba activate`, `mamba install`, and `mamba create` all work, and you can freely mix them with `conda` commands against the same environments.

---

## 4. Configure git

```bash
git config --global user.name "Your Name"
git config --global user.email "you@example.com"
```

You do **not** need `core.autocrlf` settings on macOS or Linux — both use LF line endings natively, which is what git wants.

If you want a sensible default branch name for new repos:

```bash
git config --global init.defaultBranch main
```

---

## 5. Install VS Code

### VS Code on macOS

Either download the `.zip` from [code.visualstudio.com](https://code.visualstudio.com/) and drag VS Code to `/Applications`, or install via Homebrew:

```bash
brew install --cask visual-studio-code
```

After installing, open the Command Palette (`Cmd+Shift+P`) and run `Shell Command: Install 'code' command in PATH` so you can launch VS Code from a terminal with `code .`.

### VS Code on Linux

Download the `.deb` or `.rpm` from [code.visualstudio.com](https://code.visualstudio.com/) and install it via your package manager, e.g.:

```bash
sudo apt install ./code_*.deb
# or
sudo dnf install ./code-*.rpm
```

Snap and Flatpak builds also exist; the official Microsoft `.deb`/`.rpm` is usually the smoothest path.

### Configure Python in VS Code

Install the official **Python** extension from the Extensions panel. Then, with a `.py` or `.ipynb` file open:

- For `.py` files: command palette → `Python: Select Interpreter` → choose your `astro` env.
- For notebooks: top-right kernel picker → choose your `astro` env.

Run `code .` from your project directory to launch VS Code on that project.

---

## 6. Where to put your files

There's no special advice here — put project code wherever fits your normal workflow. A common convention is to keep code repos under `~/code/` or `~/projects/` and data elsewhere (often on an external/network drive).

Avoid putting active git repos directly inside cloud-synced folders (Dropbox, iCloud Drive, OneDrive). The sync clients aggressively touch and re-touch files, which can confuse git, conflict with file locks, and silently corrupt `.git/` internals. If you must, sync the *exports* (notebooks, plots, papers) rather than the working tree.

---

## Troubleshooting

**`conda activate` or `mamba activate` doesn't work in a new shell.** Either you skipped `conda init` during install, or you didn't restart the shell after. The shell hook is shared between conda and mamba — installing or repairing it is always done via `conda init`. Run `~/miniforge3/bin/conda init <your-shell>` (e.g. `zsh` on macOS, `bash` on most Linux) and reopen the terminal. If you're not sure which shell you have, run `echo $SHELL`.

**Packages install but fail to import with "Symbol not found" or "image not found" errors on Apple Silicon.** You're probably running an x86_64 build of conda on an arm64 Mac (or vice versa). Check with `conda info | grep platform` — it should report `osx-arm64` on Apple Silicon. If it doesn't match `uname -m`, the cleanest fix is to remove `~/miniforge3` and reinstall using the correct installer from the table in section 2.

**`xcode-select --install` fails or hangs on macOS.** Sometimes the install state gets wedged. Reset with `sudo rm -rf /Library/Developer/CommandLineTools` and re-run `xcode-select --install`.

**A package needs a Fortran compiler.** Install it from conda-forge inside your env (works on both macOS and Linux, no admin needed):

```bash
mamba install -c conda-forge gfortran_osx-arm64   # Apple Silicon macOS
mamba install -c conda-forge gfortran_osx-64      # Intel macOS
mamba install -c conda-forge gfortran_linux-64    # Linux x86_64
```

**GUI plots from `matplotlib` open but freeze, or you get a backend error.** Make sure a GUI backend is installed in your env (`mamba install pyqt`) and try `matplotlib.use("QtAgg")` at the top of your script. Inside Jupyter, `%matplotlib inline` (static) or `%matplotlib widget` (interactive) avoid the backend issue entirely.
