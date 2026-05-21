# Setup for Windows

> **Attribution.** This guide is adapted from the [Code/Astro workshop's setup instructions](https://github.com/semaphoreP/codeastro/tree/main), organized by Jason Wang and collaborators. The original codeastro repo remains a great general reference for getting started with astronomy coding (project structure, testing, packaging) — recommended reading once your environment is up.

For scientific Python work in astronomy, we strongly recommend developing inside **WSL2** (Windows Subsystem for Linux) rather than native Windows. Many astronomy packages do not ship Windows wheels or rely on Unix-only tooling, and the developer experience inside WSL is essentially identical to native Linux.

If you really need to stay on native Windows, see the [Native Windows notes](#native-windows-notes) at the bottom — but expect to hit incompatibilities sooner or later.

> **For the NAOJ MRS workshop:** this document covers the general-purpose scientific Python environment under WSL. Once that is in place, continue to [workshop_packages.md](workshop_packages.md) to install the workshop-specific tools (CARTA, iSLAT, …).

---

## 1. Install WSL2

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

This installs WSL2 and the default Ubuntu distribution in one step. Restart your computer when prompted. On first launch of Ubuntu, you will be asked to create a Linux username and password — pick anything; this is independent of your Windows account.

Reference: [Microsoft's WSL install guide](https://learn.microsoft.com/windows/wsl/install).

> **GUI apps just work.** On Windows 11 (and Windows 10 with up-to-date WSL), GUI applications launched from WSL — `matplotlib` windows, Jupyter, VS Code, etc. — display natively via WSLg. No X server (XMing, VcXsrv) and no `DISPLAY` variable are needed.

### Windows Terminal

Windows Terminal ships by default with Windows 11. On Windows 10, install it from the Microsoft Store (see [Microsoft's install guide](https://learn.microsoft.com/windows/terminal/install)). After installation, the dropdown next to the new-tab button lets you open an Ubuntu (WSL) session.

---

## 2. Install build tools inside WSL

Open an Ubuntu terminal and run:

```bash
sudo apt update
sudo apt install build-essential git
```

`build-essential` provides `gcc`, `make`, and the C headers that many scientific Python packages (e.g. ones with Cython or C extensions) need at install time.

---

## 3. Install Miniforge (conda + mamba)

We recommend **Miniforge** rather than the regular Anaconda/Miniconda distributions, because:

- It defaults to the community-maintained **conda-forge** channel, which has the largest and best-maintained selection of scientific Python packages.
- It avoids Anaconda Inc.'s commercial-use restrictions on the `defaults` channel that affect many institutions.
- It ships with **mamba**, a much faster solver alongside `conda` (this is also what the JWST pipeline recommends).

> **One install, two commands.** Miniforge installs *both* `conda` and `mamba`. They share the same environments, the same config (`~/.condarc`), and the same shell integration — `mamba` is not a separate ecosystem. Use whichever feels natural; `mamba` is just faster.

Download and run the installer inside WSL:

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

> **ARM Windows (Surface Pro X, Snapdragon laptops):** use `Miniforge3-Linux-aarch64.sh` instead.

Accept the license, accept the default install location, and answer **yes** when asked whether to run `conda init`. Then close and reopen the terminal so the shell picks up the new configuration.

Verify the install:

```bash
conda --version
mamba --version
```

---

## 4. Create a Python environment

Pick a recent Python (3.12 is a good default in 2026):

```bash
mamba create -n astro python=3.12
mamba activate astro
mamba install numpy scipy matplotlib astropy jupyter
```

`mamba install` resolves environments in seconds rather than the minutes `conda install` typically takes. As noted above, `mamba` and `conda` are interchangeable in current Miniforge — `mamba activate`, `mamba install`, and `mamba create` all work, and you can freely mix them with `conda` commands against the same environments.

---

## 5. Configure git

Inside WSL, set a sensible line-ending policy:

```bash
git config --global core.autocrlf input
git config --global user.name "Your Name"
git config --global user.email "you@example.com"
```

`core.autocrlf=input` tells git to leave LF line endings alone on checkout but normalize CRLF to LF on commit — the right setting for Linux/WSL. **Do not** use `true` here; that is the Windows-side setting and will mangle files in a WSL working tree.

---

## 6. Where to put your files

**Keep your project files inside the WSL filesystem** (i.e. under `~/`), not under `/mnt/c/...`.

Cross-filesystem I/O between WSL and the Windows drive is dramatically slower — often 10× or more for git operations, `pip install`, and anything that touches many small files. Microsoft's own guidance is now to keep code in the Linux filesystem.

You can still access your WSL files easily from Windows:

- From a Windows Explorer address bar, type `\\wsl$\Ubuntu\home\<your-wsl-username>` (or `\\wsl.localhost\Ubuntu\...`).
- From a WSL terminal, run `explorer.exe .` to open the current directory in Windows Explorer.

If you genuinely need a file in both worlds (e.g. an observing log on OneDrive), you can still reach `/mnt/c/Users/<you>/...` from WSL — just don't run repos or installs from there.

---

## 7. Install VS Code and the WSL extension

1. Download and install VS Code on **Windows** from [code.visualstudio.com](https://code.visualstudio.com/). (You do not install VS Code inside WSL; the Windows install talks to a small server it deploys into WSL automatically.)
2. From the Extensions panel, install the official **WSL** extension (formerly "Remote - WSL") by Microsoft.
3. From your WSL terminal, `cd` into a project and run:

   ```bash
   code .
   ```

   The first time, this auto-installs the VS Code server inside WSL and opens a window connected to it. You'll see a green status indicator in the bottom-left reading `WSL: Ubuntu`.

4. Inside that WSL-connected VS Code window, install the **Python** extension. Extensions installed in WSL are separate from your Windows extensions — this is correct and expected.

5. When opening a notebook or `.py` file, use the interpreter picker (bottom-right or command palette `Python: Select Interpreter`) to choose your `astro` environment.

---

## Native Windows notes

If you choose to use native Windows Python rather than WSL:

- Install Miniforge for Windows from the [Miniforge releases page](https://github.com/conda-forge/miniforge/releases).
- Most pure-Python and well-packaged scientific libraries (`numpy`, `scipy`, `astropy`, `matplotlib`, `astroquery`, `photutils`) work fine on Windows.
- Expect trouble with: radio-astronomy tools (CASA, casatools), packages with Fortran extensions, anything that assumes a POSIX shell, and many MCMC/HPC-oriented packages.
- For git, on the **Windows side** the correct line-ending setting is `git config --global core.autocrlf true` (the opposite of the WSL recommendation above).

---

## Troubleshooting

**`wsl --install` says the command is not recognized.** Your Windows version is too old. Update Windows, or follow the [manual install instructions](https://learn.microsoft.com/windows/wsl/install-manual).

**`conda activate` or `mamba activate` doesn't work in a new shell.** Either you skipped `conda init` during install, or you didn't restart the shell after. The shell hook is shared between conda and mamba — installing or repairing it is always done via `conda init`. Run `~/miniforge3/bin/conda init bash` and reopen the terminal.

**`pip install` or `git clone` is unbearably slow.** You're probably working under `/mnt/c/...`. Move the repo into `~/` and try again — see [section 6](#6-where-to-put-your-files).

**GUI plots don't appear.** Check that you're on Windows 11 or a recent Windows 10 with WSLg. Run `wsl --update` from PowerShell to refresh the WSL kernel.
