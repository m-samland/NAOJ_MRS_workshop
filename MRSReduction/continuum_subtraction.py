## Continuum subtraction for JWST MIRI spectra - created by: Milou Temmink (adapted by Matthias Samland)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from astropy.io import fits

import pybaselines as pb
from spectres import spectres

import os
import pickle
import itertools

from scipy.interpolate import interp1d as I1D
from scipy.signal import savgol_filter as sf
from astropy.stats import sigma_clipped_stats as scs
import warnings

## --- Helper utilities
def _fill_endpoint_nans(arr):
    """Replace leading/trailing NaNs in ``arr`` with the nearest finite value.

    ``spectres`` fills new wavelengths that fall outside the source range with
    NaN, which then propagates through downstream cubic interpolations
    (``scipy.interpolate.interp1d``) and turns the whole baseline into NaN.
    Snap those endpoint NaNs to the nearest finite value before returning the
    baseline so the rest of the pipeline sees a fully finite array.
    """
    arr = np.asarray(arr, dtype=float).copy()
    finite = np.isfinite(arr)
    if not finite.any():
        return arr
    first = int(np.argmax(finite))
    last  = len(arr) - 1 - int(np.argmax(finite[::-1]))
    if first > 0:
        arr[:first] = arr[first]
    if last < len(arr) - 1:
        arr[last+1:] = arr[last]
    return arr


def _adaptive_window_length(desired:int, n:int, polyorder:int=3) -> int:
    """Return a valid Savitzky-Golay window_length (odd, > polyorder, <= n).
    Falls back to the largest permissible odd length. If n too small, return -1.
    """
    if n <= polyorder+1:  # cannot filter meaningfully
        return -1
    wl = min(desired, n)
    if wl % 2 == 0:
        wl -= 1
    if wl <= polyorder:
        wl = polyorder + 2  # make it the minimal odd > polyorder
        if wl % 2 == 0:
            wl += 1
        if wl > n:
            return -1
    return wl


## === INSTRUCTIONS/NOTES; DO NOT SKIP
## - This code relies on the "pybaselines" package (https://pybaselines.readthedocs.io/en/latest/index.html), you can easily install this package with: "pip install pybaselines"
## - This particular version requires the .fits files and stitches them together (using the red part of the following subband). You can provide your own stitched spectrum. BUT NOTE
## 	1) Please unsure that you, in the regions of overlapping subbands, only keep the contribution of one. This can either be the blue or red side.
## 	2) If you use your own stitched spectrum, please cut off the longest wavelengths (where the spectra become incredibly noise). The version below cuts the spectrum of at ~27.8 micron.
##	   If you do not cut off the spectrum, this method will fail.

## - On a side node; the continuum subtracted method is saved using 'pickle'. Feel free to change this to your preferred way of saving.
##   If you are using pickle and cannot open the saved results, feel free to reach out to me!

## - The method is (partially) described in my CO paper (Temmink et al. 2024). If you use this method, a reference is appreciated!
##   In addition, this version is slightly different compared to what was originally put in said paper. If you plan to publish a paper including this method, please reach out and I can give you some comments on what has changed.

## - The code may not work with older versions of Scipy; an error may occur with the savgol-filter.
##   To ensure that the code works, make sure that your environment has Scipy v1.10.0 installed!

## Any questions/suggestions/comments are always welcome!
## ===

## === ADDITIONAL NOTE:
## - THIS CODE CONTAINS A LINE AT THE END WHICH CHANGES THE ESTIMATED NOISE;
##   THE NOISE WILL NOW BE DERIVED BY ASSUMING A S/N RATIO OF 300 ON THE CONTINUUM.


## === USER-PARAMETERS OVERVIEW (what to tune and their impact)
## Data ingestion and stitching
## - Stitch (bool): True stitches JWST subbands from FITS; False loads a single ASCII spectrum.
## - StitchSide ('Blue'|'Red'): which side to keep in overlap regions when stitching; affects flux continuity at subband joins.
## - DataFile (str): path to .dat/.txt when Stitch=False; expects wavelength in col 0 and flux in col 1.
## - CutWav (float): wavelength cutoff for ch4-long; trimming noisy long-λ tail stabilizes baseline fits.
##
## Outlier detection and smoothing
## - WL (int): desired Savitzky–Golay window length for initial smoothing while detecting positive outliers.
##   The code adapts WL downward if segments are short and warns when reduced significantly.
##   Larger WL → smoother trend (less sensitive to broad features); smaller WL → follows structure more closely.
## - Thresholds (fixed in code):
##     • Positive outliers masked at +2σ relative to smoothed spectrum (removes spikes/emission lines).
##     • Negative outliers flagged at −3σ when forming StitchedMask (excludes deep dips from baseline).
##
## Global baseline fit (pybaselines.spline.irsqr)
## - NKS (int): knot-spacing divisor; NK = ceil(len(StitchedWav)/NKS). Smaller NKS → more knots → more flexible baseline; larger NKS → smoother baseline.
## - Quant (float): quantile for robust spline; lower values bias toward the lower envelope. Defaults 0.2 globally; 0.5 in 8–12 μm silicate region.
## - SD (int): spline_degree; default 3.
## - DO (int): diff_order; default 3.
## - MI (int): max_iter; default 9999 for robust convergence.
## - Lam (float): smoothness parameter; higher → smoother baseline.
## - Tol (float): tolerance; convergence threshold.
##
## Region-specific controls
## - SKS_SF (bool): use smaller knot spacing in 8.25–11.25 μm (NK ≈ len/25) to better follow the silicate feature envelope.
## - IP_SF (bool): interpolate across ±DInt around SF edges to smooth transitions; IP_Kind selects 'linear' or 'cubic'.
## - SKS_Ch4 (bool): apply smaller knot spacing for λ ≳ 17.7/17.98 μm (depending on StitchSide) to fit broad long-λ structure.
## - IP_Ch4 (bool): small interpolation across the Ch4 boundary to avoid discontinuities; uses IP_Kind.
## - IP_Kind ('linear'|'cubic'): interpolation method; 'cubic' may produce NaNs on ill-conditioned segments—fallback to 'linear' if needed.
## - IP_Mols (bool) & MolRegs ([[λ1,λ2], ...]): optionally override baseline in broad molecular bands by interpolation; use cautiously and verify visually.
##
## Plotting/output
## - dY (float): vertical padding in plots.
## - SavePlot (str): path prefix for saved PNGs.
## - Source (str): used for output directories and possible per-source defaults (e.g., NKS for specific targets).
##
## Notes
## - WL must be odd and > polyorder; the code enforces valid values and will issue a warning when reducing WL significantly for short segments, including the affected wavelength range.
## - When resampling (spectres), values outside the source range are filled (default 0); consider switching to fill=np.nan and ignoring NaNs if edge effects appear.


## === Functions:
def Baseline(Wavelength = [],
             Flux       = [],
             Mask       = [],
             NK         = 100,
             Quant      = 0.05,
             SD         = 3,
             DO         = 3,
             MI         = 100,
             Lam        = 1e2,
             Tol        = 1e-6):
    ## Determining the baseline of JWST spectra using a consistent way:
    ## https://pybaselines.readthedocs.io/en/latest/
    Wavelength = np.asarray(Wavelength)
    Flux       = np.asarray(Flux)
    Mask       = np.asarray(Mask, dtype=bool)
    if len(Wavelength) != len(Flux) or len(Wavelength) != len(Mask):
        raise ValueError("Wavelength, Flux, and Mask must have the same length")
    if not np.any(Mask):
        # If everything got masked, fall back to using all points
        Mask = np.ones_like(Wavelength, dtype=bool)

    Baseline = pb.spline.irsqr(x_data        = Wavelength[Mask],
                               data          = Flux[Mask],
                               num_knots     = NK,
                               quantile      = Quant,
                               spline_degree = SD,
                               diff_order    = DO,
                               max_iter      = MI,
                               lam           = Lam,
                               tol           = Tol)[0]
    Baseline = spectres(Wavelength, Wavelength[Mask], Baseline)
    return _fill_endpoint_nans(Baseline)

def PureBaseline(Wavelength = [],
                 MaskedWav  = [],
                 MaskedFlux = [],
                 NK         = 100,
                 Quant      = 0.05,
                 SD         = 3,
                 DO         = 3,
                 MI         = 100,
                 Lam        = 1e2,
                 Tol        = 1e-6):
    Baseline = pb.spline.irsqr(x_data        = MaskedWav,
                               data          = MaskedFlux,
                               num_knots     = NK,
                               quantile      = Quant,
                               spline_degree = SD,
                               diff_order    = DO,
                               max_iter      = MI,
                               lam           = Lam,
                               tol           = Tol)[0]
    Baseline = spectres(Wavelength, MaskedWav, Baseline)
    return _fill_endpoint_nans(Baseline)

def PlotSpectrum(Wav      = [],
                 Flux     = [],
                 Baseline = [],
                 Mask     = [],
                 dY       = 0.25,
                 SavePlot = ''):
    
    XL = [np.nanmin(Wav), np.nanmax(Wav)]
    YL = [np.nanmin(Flux), np.nanmax(Flux)+dY]
    
    if len(Baseline) > 0:
        # Link x-axes for interactive zoom/pan
        fig, ax = plt.subplots(2, 1, sharex=True, figsize=(10,7))

        # Top panel: original spectrum and baseline
        ax[0].step(Wav, Flux, color='k', where='mid', zorder=5)
        ax[0].plot(Wav, Baseline, color='firebrick', zorder=5)
        ax[0].scatter(Wav[~Mask], Flux[~Mask], color='red', marker='x', s=50)
        ax[0].set(xlim=XL, ylim=YL)
        ax[0].tick_params(labelbottom=False)  # hide top x tick labels without altering formatter

        # Bottom panel: continuum-subtracted
        ax[1].step(Wav, Flux-Baseline, color='k', where='mid', zorder=5)
        ax[1].axhline(y=0, color='firebrick', ls='dashed')
        ax[1].scatter(Wav[~Mask], (Flux-Baseline)[~Mask], color='red', marker='x', s=50)
        YL = [np.nanmin(Flux-Baseline), np.nanmax(Flux-Baseline)+dY]
        ax[1].set(xlabel=r'Wavelength [$\mu$m]', ylabel='Flux [Jy]', xlim=XL, ylim=YL)
        ax[1].tick_params(labelbottom=True)
    else:
        fig, ax = plt.subplots(figsize=(10,4))
        ax.step(Wav, Flux, color='k', where='mid')
        ax.set(xlabel=r'Wavelength [$\mu$m]', ylabel='Flux [Jy]', xlim=XL, ylim=YL)
    
    # Support both single and two-panel layouts
    ax_list = ax if isinstance(ax, (list, tuple, np.ndarray)) else [ax]
    for a in ax_list:
        a.xaxis.set_minor_locator(AutoMinorLocator(5))
        a.yaxis.set_minor_locator(AutoMinorLocator(5))
        

    ## --- Plot the blue and red ends over the overlapping wavelength regions:
    Reds  = [5.74, 6.63, 7.65, 8.77, 10.13, 11.70, 13.47, 15.57, 17.98, 20.95, 24.48]
    Blues = [5.66, 6.53, 7.51, 8.67, 10.02, 11.55, 13.34, 15.41, 17.70, 20.69, 24.19]

    if isinstance(ax, (list, tuple, np.ndarray)):
        for i in range(len(Reds)):
            ax[0].axvline(x=Reds[i], color='red', lw=2, alpha=0.6, ls='dashed')
            ax[1].axvline(x=Reds[i], color='red', lw=2, alpha=0.6, ls='dashed')
            ax[0].axvline(x=Blues[i], color='royalblue', lw=2, alpha=0.6, ls='dashed')
            ax[1].axvline(x=Blues[i], color='royalblue', lw=2, alpha=0.6, ls='dashed')
            
    
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(f'{SavePlot}.png', dpi=250)
    plt.show()
    plt.close()
    


## === Data handling
Source = 'SYCha'
if not os.path.exists(f'./data/{Source}/stage3/'):
    os.mkdir(f'./data/{Source}/stage3/')
    pass
OutDir = f'./data/{Source}/stage3/ContSub/' ## Folder where the continuum subtracted spectrum will be placed.
if not os.path.exists(OutDir):
    os.mkdir(OutDir)
    pass

## Load in spectrum, two options
## First option is to supply a file (.dat/.txt that contains an already stitched spectrum)
## Second option is to supply all datafiles of the various subbands, those will then be stitched together.
Stitch     = False ## If True, it assumes the second option, if False it assumes the Second option.
StitchSide = 'Red' ## Choose what side is used in the overlapping regions: "Blue" for the blue wavelengthside of the following subband, "Red" for the redside of the earlier subband.
DataFile   = f'./data/{Source}/stage3/{Source}_spectrum_stitched.csv' ## CSV produced by combine_x1d_spectra.py (used when Stitch=False).
## When Stitch=False: which column of the stitched CSV to fit the baseline on.
## 'rf_flux' uses the residual-fringe-corrected flux (cleaner continuum); 'flux' uses the plain extracted flux.
## In either case the resulting baseline is subtracted from BOTH columns, so the output contains continuum-subtracted versions of each.
## Falls back to 'flux' automatically if 'rf_flux' is missing or entirely NaN (e.g. MAST products without RFC columns).
FitColumn = 'rf_flux'
FluxCol = 'RF_FLUX'  ## Used by the Stitch=True branch (FITS column name); independent from FitColumn above.
# FluxCol = 'FLUX'
FileNameRFC = True

if FileNameRFC:
    FileNameAdd = '_rfc1d'

CutWav = 27.5 ## Cut-off wavelength for subband 4C. 
## Including the longest wavelengths (>27.8 or >27.5 for some sources) may cause the continuum subtraction to fail.

## Use a smaller knot spacing at the silicate feature (8-12 micron); i.e. one spline every 25 data points.
## This may improve the continuum subtraction for broader features seen at the longer wavelengths.
## Additionally, you can select the option to use a small interpolation to have a smoother transition between
## the initial baseline fit (Ch1-3) and the one for Ch4.
SKS_SF = True
IP_SF  = True

## Use a smaller knot spacing in Channel 4? I.e. one spline not every 25 datapoints.
## This may improve the continuum subtraction for broader features seen at the longer wavelengths.
## Additionally, you can select the option to use a small interpolation to have a smoother transition between
## the initial baseline fit (Ch1-3) and the one for Ch4.
SKS_Ch4 = True
IP_Ch4  = True

## Interpolation method (i.e. 'linear' or 'cubic') to use, see scipy for a better description.
## Using 'cubic' might lead to nans! If that error occurs, please set the interpolation to 'linear'
IP_Kind = 'cubic'

## For broad molecular features, you can select regions (starting line 382) where to adjust the baseline using an interpolation.
## I.e. for a small pseudocontinuum of C2H2 (and/or HCN), you can select the region (recommended: [[13.45, 14.2]]) where you want to use a cubic interpolation
## The code is set to iterate over various regions, [[13.45, 14.2], [xx, yy]], etc. (See MolRegs parameter below)
## By setting 'IP_Mols = True' you will apply this interpolation.
## NOTE, it may not give the desired results for every spectrum. Double check whether you trust the interpolation when you use it!
IP_Mols = False
MolRegs = [[5.1, 5.148], [5.865, 5.934], [6.28, 6.33], [6.75, 6.8], [7.42, 7.54], [13.45, 14.2], [14.88, 15.01]]

## Window length for the savgol filter. Masking over the entire spectrum used WL = 100.
WL = 100

## Defaults for the optional second-flux-column outputs. Overridden in the
## Stitch=False branch when the input CSV carries both ``flux`` and ``rf_flux``.
StitchedFluxAlt = None
AltColumn = None

if Stitch:
    ## === STITCHING FROM FITS SUBBANDS ===
    ## Reads JWST MIRI subbands (Ch 1–4; short/medium/long), detects and removes positive outliers per subband,
    ## flags strong negative deviations for baseline masking, trims long-λ tail in ch4-long,
    ## and concatenates into StitchedWav/Flux/Err/Mask. Overlap handling is controlled by StitchSide (Blue/Red).
    ObsDir = f'./data/{Source}/stage3/'  ## Folder that contains the per-subband x1d FITS files.
    StitchedWav     = np.array([])
    StitchedFlux    = np.array([])
    StitchedFluxErr = np.array([])
    StitchedMask    = np.array([], dtype=bool)
    Chans, Ranges   = [1, 2, 3, 4], ['short', 'medium', 'long']
    if StitchSide == 'Blue':
        for i, (Chan, Range) in enumerate(itertools.product(Chans, Ranges)):
            print(f'Detection outliers in {Chan}-{Range}')
            #Level3_ifua_rfc1d_psf_ch4-medium_x1d.fits
            with fits.open(f'{ObsDir}Level3_ifua{FileNameAdd}_psf_ch{Chan}-{Range}_x1d.fits') as hdul:
                Data = hdul[1].data
            Spectrum = {'Wavelength':Data['WAVELENGTH'], 'Flux':Data[FluxCol], 'Uncertainty':Data['FLUX_ERROR']}

            if Chan == 4 and Range == 'long':
                WavMask                = (Spectrum['Wavelength'] <= CutWav)
                Spectrum['Wavelength'] = Spectrum['Wavelength'][WavMask]
                Spectrum['Flux']       = Spectrum['Flux'][WavMask]
                pass
            
            W, F        = Spectrum['Wavelength'], Spectrum['Flux']
            OutlierMask = np.ones_like(W, dtype=bool)
            Outliers    = 1
            iter_count = 0
            reduced_warned = False
            ## --- Iterative positive-outlier clipping (> +2σ) using SG smoothing ---
            while Outliers > 0:
                iter_count += 1
                wl_eff = _adaptive_window_length(WL, len(F), polyorder=3)
                if wl_eff == -1:
                    # Not enough points left to continue filtering
                    break
                if (not reduced_warned) and wl_eff < 0.75*WL:
                    wmin, wmax = float(np.nanmin(W)), float(np.nanmax(W))
                    warnings.warn(f"Savitzky-Golay window reduced from {WL} to {wl_eff} due to segment length (≈{len(F)} pts) in {wmin:.2f}–{wmax:.2f} μm.")
                    reduced_warned = True
                Filtered = sf(F, window_length=wl_eff, polyorder=3)
                STD      = scs(F-Filtered, sigma=3)[2]
                Mask     = (F-Filtered) > 2*STD   ## 2 Sigma lines are masked.

                Outliers = len(F[Mask])
                print(Outliers)
                if Outliers > 0:
                    W, F = W[~Mask], F[~Mask]
                if iter_count > 200:  # safety valve
                    print('Stopping iterations (Blue) after 200 loops to avoid oscillation.')
                    break
            ## Final smoothing and negative-outlier mask (< −3σ) for this subband
            wl_eff_final = _adaptive_window_length(WL, len(F), polyorder=3)
            if wl_eff_final == -1:
                Filtered = F.copy()
            else:
                if wl_eff_final < 0.75*WL:
                    wmin, wmax = float(np.nanmin(W)), float(np.nanmax(W))
                    warnings.warn(f"Final smoothing window reduced from {WL} to {wl_eff_final} (≈{len(F)} pts) in {wmin:.2f}–{wmax:.2f} μm.")
                Filtered    = sf(F, window_length=wl_eff_final, polyorder=3)
            STD         = scs(F-Filtered, sigma=3)[2]  ## STD over the masked Savgol-filter.
            FilterFull  = spectres(Spectrum['Wavelength'], W, Filtered)
            Mask        = (Spectrum['Flux']-FilterFull) < -3*STD
            
            OutlierMask[Mask] = False
            
            ## Append current subband into stitched arrays (keeping Blue side in overlaps)
            if i == 0:
                StitchedWav     = np.append(StitchedWav, Spectrum['Wavelength'])
                StitchedFlux    = np.append(StitchedFlux, Spectrum['Flux'])
                StitchedFluxErr = np.append(StitchedFluxErr, Spectrum['Uncertainty'])
                StitchedMask    = np.append(StitchedMask, OutlierMask)
                pass
            else:
                Mask            = (StitchedWav <= np.nanmin(Spectrum['Wavelength']))
                StitchedWav     = np.append(StitchedWav[Mask], Spectrum['Wavelength'])
                StitchedFlux    = np.append(StitchedFlux[Mask], Spectrum['Flux'])
                StitchedFluxErr = np.append(StitchedFluxErr[Mask], Spectrum['Uncertainty'])
                StitchedMask    = np.append(StitchedMask[Mask], OutlierMask)
                pass
            pass
        pass
    elif StitchSide == 'Red':
        for i, (Chan, Range) in enumerate(itertools.product(Chans, Ranges)):
            print(f'Detection outliers in {Chan}-{Range}')
            #Level3_ifua_rfc1d_psf_ch4-medium_x1d.fits
            with fits.open(f'{ObsDir}Level3_ifua{FileNameAdd}_psf_ch{Chan}-{Range}_x1d.fits') as hdul:
                Data = hdul[1].data
            Spectrum = {'Wavelength':Data['WAVELENGTH'], 'Flux':Data[FluxCol], 'Uncertainty':Data['FLUX_ERROR']}
            
            if Chan == 4 and Range == 'long':
                WavMask                 = (Spectrum['Wavelength'] <= CutWav)
                Spectrum['Wavelength']  = Spectrum['Wavelength'][WavMask]
                Spectrum['Flux']        = Spectrum['Flux'][WavMask]
                Spectrum['Uncertainty'] = Spectrum['Uncertainty'][WavMask]
                pass
            
            W, F        = Spectrum['Wavelength'], Spectrum['Flux']
            OutlierMask = np.ones_like(W, dtype=bool)
            Outliers    = 1
            iter_count = 0
            reduced_warned = False
            ## --- Iterative positive-outlier clipping (> +2σ) using SG smoothing ---
            while Outliers > 0:
                iter_count += 1
                wl_eff = _adaptive_window_length(WL, len(F), polyorder=3)
                if wl_eff == -1:
                    break
                if (not reduced_warned) and wl_eff < 0.75*WL:
                    wmin, wmax = float(np.nanmin(W)), float(np.nanmax(W))
                    warnings.warn(f"Savitzky-Golay window reduced from {WL} to {wl_eff} due to segment length (≈{len(F)} pts) in {wmin:.2f}–{wmax:.2f} μm.")
                    reduced_warned = True
                Filtered = sf(F, window_length=wl_eff, polyorder=3)
                STD      = scs(F-Filtered, sigma=3)[2]
                Mask     = (F-Filtered) > 2*STD  ## 2-sigma positive deviations
                Outliers = len(F[Mask])
                print(Outliers)
                if Outliers > 0:
                    W, F = W[~Mask], F[~Mask]
                if iter_count > 200:
                    print('Stopping iterations (Red) after 200 loops to avoid oscillation.')
                    break
            ## Final smoothing and negative-outlier mask (< −3σ) for this subband
            wl_eff_final = _adaptive_window_length(WL, len(F), polyorder=3)
            if wl_eff_final == -1:
                Filtered = F.copy()
            else:
                if wl_eff_final < 0.75*WL:
                    wmin, wmax = float(np.nanmin(W)), float(np.nanmax(W))
                    warnings.warn(f"Final smoothing window reduced from {WL} to {wl_eff_final} (≈{len(F)} pts) in {wmin:.2f}–{wmax:.2f} μm.")
                Filtered    = sf(F, window_length=wl_eff_final, polyorder=3)
            STD         = scs(F-Filtered, sigma=3)[2]  ## STD over the masked Savgol-filter.
            FilterFull  = spectres(Spectrum['Wavelength'], W, Filtered)
            Mask        = (Spectrum['Flux']-FilterFull) < -3*STD
            
            OutlierMask[Mask] = False

            ## Append current subband into stitched arrays (keeping Red side in overlaps)
            if i == 0:
                StitchedWav     = np.append(StitchedWav, Spectrum['Wavelength'])
                StitchedFlux    = np.append(StitchedFlux, Spectrum['Flux'])
                StitchedFluxErr = np.append(StitchedFluxErr, Spectrum['Uncertainty'])
                StitchedMask    = np.append(StitchedMask, OutlierMask)
                pass
            else:
                Mask            = (Spectrum['Wavelength'] >= np.nanmax(StitchedWav))
                StitchedWav     = np.append(StitchedWav, Spectrum['Wavelength'][Mask])
                StitchedFlux    = np.append(StitchedFlux, Spectrum['Flux'][Mask])
                StitchedFluxErr = np.append(StitchedFluxErr, Spectrum['Uncertainty'][Mask])
                StitchedMask    = np.append(StitchedMask, OutlierMask[Mask])
                pass
            pass
        pass
    pass
else:
    ## === SINGLE STITCHED CSV SPECTRUM ===
    ## Load a pre-stitched spectrum written by ``combine_x1d_spectra.py``: a CSV
    ## with a header row including columns ``wavelength``, ``flux``, and
    ## (optionally) ``flux_error``. Detect positive outliers over the full band
    ## and build the negative-outlier mask used to anchor the baseline fit.
    ## The same SG-based clipping and final smoothing logic applies here as for
    ## the per-subband path.
    ## NOTE, as opposed to loading in the spectrum band-by-band, the outliers
    ## will be inferred over the entire spectrum.
    Spectrum = pd.read_csv(DataFile)
    print(Spectrum)

    StitchedWav = Spectrum['wavelength'].to_numpy()

    ## Load both flux columns when available so the baseline (fit on FitColumn)
    ## can be subtracted from BOTH afterwards.
    HasFlux   = 'flux' in Spectrum.columns
    HasRFFlux = 'rf_flux' in Spectrum.columns
    PlainFlux = Spectrum['flux'].to_numpy()    if HasFlux   else np.full(len(StitchedWav), np.nan)
    RFFlux    = Spectrum['rf_flux'].to_numpy() if HasRFFlux else np.full(len(StitchedWav), np.nan)

    ## Auto-fallback when rf_flux is unavailable (e.g. MAST products without RFC).
    if FitColumn == 'rf_flux' and not (HasRFFlux and np.any(np.isfinite(RFFlux))):
        print("FitColumn='rf_flux' requested but the rf_flux column is missing/all-NaN; falling back to 'flux'.")
        FitColumn = 'flux'

    ## ``StitchedFlux`` is the column the rest of the script fits on; ``StitchedFluxAlt`` is the other column
    ## (the same baseline is later subtracted from it as well).
    if FitColumn == 'rf_flux':
        StitchedFlux,    StitchedFluxAlt = RFFlux,    PlainFlux
        AltColumn = 'flux'
    else:
        StitchedFlux,    StitchedFluxAlt = PlainFlux, RFFlux
        AltColumn = 'rf_flux'

    if 'flux_error' in Spectrum.columns:
        StitchedFluxErr = Spectrum['flux_error'].to_numpy()
    else:
        StitchedFluxErr = np.full_like(StitchedFlux, np.nan, dtype=float)

    WavMask         = (StitchedWav <= CutWav)
    StitchedWav     = StitchedWav[WavMask]
    StitchedFlux    = StitchedFlux[WavMask]
    StitchedFluxAlt = StitchedFluxAlt[WavMask]
    StitchedFluxErr = StitchedFluxErr[WavMask]

    W, F           = StitchedWav.copy(), StitchedFlux.copy()
    Outliers       = 1
    iter_count     = 0
    reduced_warned = False
    while Outliers > 0:
        iter_count += 1
        wl_eff = _adaptive_window_length(WL, len(F), polyorder=3)
        if wl_eff == -1:
            break
        if (not reduced_warned) and wl_eff < 0.75*WL:
            wmin, wmax = float(np.nanmin(W)), float(np.nanmax(W))
            warnings.warn(f"Savitzky-Golay window reduced from {WL} to {wl_eff} due to segment length (≈{len(F)} pts) in {wmin:.2f}–{wmax:.2f} μm.")
            reduced_warned = True
        Filtered = sf(F, window_length=wl_eff, polyorder=3)
        STD      = scs(F-Filtered, sigma=3)[2]
        Mask     = (F-Filtered) > 2*STD   ## 2-sigma positive deviations
        Outliers = len(F[Mask])
        if Outliers > 0:
            W, F = W[~Mask], F[~Mask]
        if iter_count > 200:
            print('Stopping iterations (Single spectrum) after 200 loops to avoid oscillation.')
            break
    wl_eff_final = _adaptive_window_length(WL, len(F), polyorder=3)
    if wl_eff_final == -1:
        Filtered = F.copy()
    else:
        if wl_eff_final < 0.75*WL:
            wmin, wmax = float(np.nanmin(W)), float(np.nanmax(W))
            warnings.warn(f"Final smoothing window reduced from {WL} to {wl_eff_final} (≈{len(F)} pts) in {wmin:.2f}–{wmax:.2f} μm.")
        Filtered = sf(F, window_length=wl_eff_final, polyorder=3)
    STD         = scs(F-Filtered, sigma=3)[2]  ## STD over the masked Savgol-filter.
    FilterFull  = spectres(StitchedWav, W, Filtered)
    ## Mask convention: True = keep (used as the include-mask in Baseline()).
    StitchedMask = (StitchedFlux - FilterFull) >= -3*STD
    StitchedMask = StitchedMask.astype(bool)

## === Estimate the baseline:
## ACTUAL GLOBAL BASELINE FIT HAPPENS HERE
## - Change GLOBAL baseline behavior here (applies across all channels):
##   • NKS (sets NK via NK = ceil(len(StitchedWav)/NKS))
##   • Quant (robust quantile, lower → lower envelope)
##   • SD, DO, MI, Lam, Tol (spline controls)
StitchedMask = StitchedMask.astype(bool)
NKS = 75
if Source in ['CXTau', 'DNTau']:
    NKS = 25
    pass

BL = Baseline(Wavelength = StitchedWav,
              Flux       = StitchedFlux,  ## The interpolated flux is being used for obtaining the continuum; if MaskFeature = True
              Mask       = StitchedMask,
              NK         = int(np.ceil(len(StitchedWav)/NKS)),  ## CHANGE HERE (GLOBAL NK via NKS)
              Quant      = 0.2,                                 ## CHANGE HERE (GLOBAL Quant)
              SD         = 3,
              DO         = 3,
              MI         = 9999,
              Lam        = 1e2,
              Tol        = 1e-9)
              
StitchedBL = BL

## === Adjusting the basline for broader molecular features (i.e. C2H2):
if IP_Mols:
    dMask   = 0.05 ## The additional bit that will be ignored when determining the interpolation.
    for i, Reg in enumerate(MolRegs):
        WavMask = (StitchedWav >= Reg[0]-dMask) & (StitchedWav <= Reg[1]+dMask)
        AdjWav  = StitchedWav[~WavMask].copy()
        AdjBL   = StitchedBL[~WavMask].copy()
    
        Interpolation       = I1D(x    = AdjWav,
                                  y    = AdjBL,
                                  kind = 'cubic')
        NewBL               = Interpolation(StitchedWav[WavMask])
        StitchedBL[WavMask] = NewBL
        pass
    pass

## === Silicate feature: 8-12 micron, smaller knot spacing and a different quantile value (0.5).
## REGIONAL BASELINE OVERRIDE (SHORT λ: 8.25–11.25 μm)
## - Tune short-region baseline by changing NK and Quant below.
if SKS_SF:
    MinWav, MaxWav = 8.25, 11.25
    DInt           = 0.05
    
    WavMask1 = (StitchedWav >= MinWav) & (StitchedWav <= MaxWav)
    WavMask2 = (StitchedWav[StitchedMask] >= MinWav) & (StitchedWav[StitchedMask] <= MaxWav)
    BLSF     = PureBaseline(Wavelength = StitchedWav[WavMask1],
                            MaskedWav  = StitchedWav[StitchedMask][WavMask2],
                            MaskedFlux = StitchedFlux[StitchedMask][WavMask2],
                            NK         = int(np.ceil(len(StitchedWav[WavMask1])/25)),  ## CHANGE HERE (SHORT-REGION NK)
                            Quant      = 0.5,                                          ## CHANGE HERE (SHORT-REGION Quant)
                            SD         = 3,
                            DO         = 3,
                            MI         = 9999,
                            Lam        = 1e2,
                            Tol        = 1e-9)

    BL[WavMask1] = BLSF

    if IP_SF:
        ## Add a small interpolation to the mask edges --> smoother transition.
        IPMask        = (StitchedWav >= MinWav-DInt) & (StitchedWav <= MinWav+DInt)
        Interpolation = I1D(x=StitchedWav[~IPMask], y=BL[~IPMask], kind=IP_Kind)
        BL            = Interpolation(StitchedWav)

        IPMask        = (StitchedWav >= MaxWav-DInt) & (StitchedWav <= MaxWav+DInt)
        Interpolation = I1D(x=StitchedWav[~IPMask], y=BL[~IPMask], kind=IP_Kind)
        BL            = Interpolation(StitchedWav)
        pass
        
    StitchedBL = BL
    pass

## === Channel 4: smaller knot spacing
## REGIONAL BASELINE OVERRIDE (LONG λ: Channel 4)
## - Starts ~17.7 μm (StitchSide='Blue') or ~17.98 μm (StitchSide='Red').
## - Tune long-region baseline by changing NK and Quant below.
if SKS_Ch4:
    # Initialize masks and guard variables to avoid unbound usage
    WavMask1 = np.zeros_like(StitchedWav, dtype=bool)
    # For WavMask2 we need masked array size; will reset below accordingly
    WavMask2 = np.zeros_like(StitchedWav[StitchedMask], dtype=bool)
    IPMask = None
    if Stitch:
        if StitchSide == 'Blue':
            WavMask1 = (StitchedWav >= 17.70)
            WavMask2 = (StitchedWav[StitchedMask] >= 17.70)
            if IP_Ch4:
                IPMask = (StitchedWav >= 17.65) & (StitchedWav <= 17.75)
                pass
            pass
        elif StitchSide == 'Red':
            WavMask1 = (StitchedWav >= 17.98)
            WavMask2 = (StitchedWav[StitchedMask] >= 17.98)
            if IP_Ch4:
                IPMask = (StitchedWav >= 17.93) & (StitchedWav <= 18.03)
                pass
            pass
        pass
    else:
        WavMask1 = (StitchedWav >= 17.70)  ## Change 17.70 into 17.98 if you stitched using the red-side of the previous subbands during stitches.
        WavMask2 = (StitchedWav[StitchedMask] >= 17.70)
        
        if IP_Ch4:
            IPMask = (StitchedWav >= 17.65) & (StitchedWav <= 17.75)
            pass
        pass

    if np.any(WavMask1) and np.any(WavMask2):
        BL4ABC = PureBaseline(
            Wavelength = StitchedWav[WavMask1],
            MaskedWav  = StitchedWav[StitchedMask][WavMask2],
            MaskedFlux = StitchedFlux[StitchedMask][WavMask2],
            NK         = int(np.ceil(len(StitchedWav[WavMask1])/20)),  ## CHANGE HERE (LONG-REGION NK)
            Quant      = 0.2,                                           ## CHANGE HERE (LONG-REGION Quant)
            SD         = 3,
            DO         = 3,
            MI         = 9999,
            Lam        = 1e2,
            Tol        = 1e-9
        )
        BL[WavMask1] = BL4ABC

        if IP_Ch4 and (IPMask is not None) and np.any(~IPMask):
            Interpolation = I1D(x=StitchedWav[~IPMask], y=BL[~IPMask], kind=IP_Kind)
            BL            = Interpolation(StitchedWav)
            pass
        
    StitchedBL = BL
    pass
    
StitchedFluxCS = StitchedFlux - StitchedBL
## When a second flux column is present (Stitch=False CSV with both ``flux`` and ``rf_flux``),
## apply the same baseline to it so the output carries continuum-subtracted versions of both.
StitchedFluxAltCS = (StitchedFluxAlt - StitchedBL) if StitchedFluxAlt is not None else None

## === Plot results
PlotSpectrum(Wav      = StitchedWav,
             Flux     = StitchedFlux,  ## Using the original stitched spectrum in the plotting.
             Baseline = StitchedBL,
             Mask     = StitchedMask,
             dY       = 0.025,
             SavePlot = f'{OutDir}FullSpectrum_CS')

## === Save the results:
## Canonical column names regardless of which one drove the baseline fit:
##   Flux     / CSFlux    — plain extracted flux and its continuum-subtracted version
##   RFFlux   / CSRFFlux  — residual-fringe-corrected flux and its continuum-subtracted version
## ``FitColumn`` tags which column was actually fit; the same baseline is in both CS columns.
ContSubData = {'wave': StitchedWav,
               'FluxErr': StitchedFluxErr,
               'Mask': StitchedMask,
               'Baseline': StitchedBL,
               'StitchSide': StitchSide,
               'FitColumn': FitColumn if not Stitch else FluxCol.lower()}

if not Stitch:
    ## Stitch=False: both columns are loaded → emit all four canonical keys.
    if FitColumn == 'rf_flux':
        ContSubData['RFFlux']   = StitchedFlux
        ContSubData['CSRFFlux'] = StitchedFluxCS
        ContSubData['Flux']     = StitchedFluxAlt
        ContSubData['CSFlux']   = StitchedFluxAltCS
    else:
        ContSubData['Flux']     = StitchedFlux
        ContSubData['CSFlux']   = StitchedFluxCS
        ContSubData['RFFlux']   = StitchedFluxAlt
        ContSubData['CSRFFlux'] = StitchedFluxAltCS
else:
    ## Stitch=True: only one FITS column is loaded → emit just that pair.
    if FluxCol == 'RF_FLUX':
        ContSubData['RFFlux']   = StitchedFlux
        ContSubData['CSRFFlux'] = StitchedFluxCS
    else:
        ContSubData['Flux']     = StitchedFlux
        ContSubData['CSFlux']   = StitchedFluxCS

pickle.dump(ContSubData, open(f'{OutDir}FullSpectrum_CS.p', 'wb'))

import pandas as pd
a = pd.DataFrame(ContSubData)
a.to_csv(f'{OutDir}FullSpectrum_CS.csv', index=False)

## Companion file for iSLAT: two columns only (``wave``, ``flux``), where
## ``flux`` is the continuum-subtracted version of the column that drove the
## baseline fit (CSRFFlux when FitColumn='rf_flux', CSFlux otherwise).
## ``StitchedFluxCS`` is, by construction, that same array in both branches.
pd.DataFrame({
    'wave': StitchedWav,
    'flux': StitchedFluxCS,
}).to_csv(f'{OutDir}FullSpectrum_CS_for_iSLAT.csv', index=False)

