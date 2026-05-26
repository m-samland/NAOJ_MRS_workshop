#!/usr/bin/env python
"""Post-process MIRI MRS 1D spectra: rescale and stitch the 12 sub-band x1d files.

Takes a folder containing per-band ``*_x1d.fits`` files for a single target
(either pipeline output or MAST-downloaded Level-3 products), reads them into a
combined table, rescales each band so adjacent bands match in their overlap
regions, and writes:

  - ``<source>_spectrum_full.csv``      all 12 bands, rescaled, with overlaps
  - ``<source>_spectrum_stitched.csv``  single non-overlapping spectrum (1A -> 4C)

The script does not require the JWST pipeline at runtime; it only needs
``numpy``, ``pandas``, ``astropy``, and (optionally) ``stdatamodels`` for the
``WAVELENGTH`` column.

Examples
--------
Command line (from the workshop root, with the ``mrs_analysis`` env active)::

    # Write CSVs alongside the input FITS files
    python MRSReduction/combine_x1d_spectra.py data/SYCha/stage3 --source SYCha

    # Write CSVs to a separate output directory
    python MRSReduction/combine_x1d_spectra.py data/SYCha/stage3 \\
        --source SYCha \\
        --output-dir data/SYCha/spectra

From Python or a notebook::

    from MRSReduction.combine_x1d_spectra import post_process

    spectrum_full, spectrum_stitched = post_process(
        input_dir='data/SYCha/stage3',
        source='SYCha',
        output_dir='data/SYCha/spectra',
        wl_upper=27.5,
    )

    import matplotlib.pyplot as plt
    plt.plot(spectrum_stitched['wavelength'], spectrum_stitched['flux'])
    plt.xlabel('Wavelength [um]'); plt.ylabel('Flux [Jy]')
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

try:
    from stdatamodels.jwst import datamodels
    _HAS_DATAMODELS = True
except ImportError:
    _HAS_DATAMODELS = False


BANDS = ['1A', '1B', '1C', '2A', '2B', '2C',
         '3A', '3B', '3C', '4A', '4B', '4C']

_BAND_LETTER = {'short': 'A', 'medium': 'B', 'long': 'C'}


def convert_to_native_endian(data):
    """Copy a FITS column array into native-endian numpy form (pandas-safe)."""
    native_data = np.empty(len(data))
    native_data[:] = data
    return native_data


def _parse_channel_band(basename: str):
    """Return (channel:int|None, band_name:str|None, band_letter:str|None)."""
    channel = next((int(c) for c in '1234' if f'ch{c}' in basename), None)
    band_name = next((b for b in _BAND_LETTER if b in basename), None)
    band_letter = _BAND_LETTER[band_name] if band_name is not None else None
    return channel, band_name, band_letter


def x1d_files_to_table(x1d_fnames):
    """Read a list of MIRI MRS x1d FITS files into a single DataFrame.

    Channel and band (A/B/C) are parsed from the filename; ``BAND_ID`` is
    e.g. ``'1A'``, ``'2B'``, etc.
    """
    spectra = []
    for fname in x1d_fnames:
        basename = os.path.splitext(os.path.basename(fname))[0]
        channel, band_name, band_letter = _parse_channel_band(basename)
        if channel is None or band_letter is None:
            raise ValueError(
                f"Could not parse channel/band from filename {basename!r}. "
                "Expected tokens like 'ch1' / 'short' in the name."
            )

        # Wavelength: prefer datamodels (also gives canonical channel), fall back to fits.
        # datamodels.open can fail for many reasons on partial/non-pipeline files;
        # any failure is fine because the fits-based read below covers it.
        wl = None
        if _HAS_DATAMODELS:
            try:
                dm = datamodels.open(fname)  # type: ignore[name-defined]
                wl = np.array(dm.spec[0].spec_table['WAVELENGTH'])
                channel = int(dm.meta.instrument.channel)
            except Exception:  # noqa: BLE001  -- datamodels surface area is wide
                wl = None

        with fits.open(fname) as hdul:
            hdul.verify('ignore')
            table = hdul[1].data
            if wl is None:
                wl = table['WAVELENGTH']

            # RF_FLUX / RF_BACKGROUND are written by the residual-fringe-corrected
            # extract_1d step. Many MAST-downloaded x1d products do not carry them;
            # fall back to NaN columns so downstream code keeps the same schema.
            def _maybe_column(name):
                try:
                    return convert_to_native_endian(table[name])
                except KeyError:
                    nan_col = np.empty(len(table['FLUX']))
                    nan_col[:] = np.nan
                    return nan_col

            spectrum_data = {
                'wavelength': convert_to_native_endian(wl),
                'flux': convert_to_native_endian(table['FLUX']),
                'flux_error': convert_to_native_endian(table['FLUX_ERROR']),
                'background': convert_to_native_endian(table['BACKGROUND']),
                'bkgd_error': convert_to_native_endian(table['BKGD_ERROR']),
                'rf_flux': _maybe_column('RF_FLUX'),
                'rf_background': _maybe_column('RF_BACKGROUND'),
            }

        spectrum = pd.DataFrame(spectrum_data)
        spectrum['CHANNEL'] = channel
        spectrum['BAND'] = band_name
        spectrum['BAND_NUM'] = band_letter
        spectrum['BAND_ID'] = f'{channel}{band_letter}'
        spectrum['wavelength_index'] = spectrum.index
        spectra.append(spectrum)

    return pd.concat(spectra, ignore_index=True)


def _create_line_segments(spectrum, bands, flux_column):
    return [
        (
            spectrum.loc[spectrum['BAND_ID'] == b, 'wavelength'].values,
            spectrum.loc[spectrum['BAND_ID'] == b, flux_column].values,
        )
        for b in bands
    ]


def _rescale_line_segments(line_segments):
    """Walk band-to-band, scaling each band to match the previous one in their
    overlap region. Returns the list of per-band scaling factors (first = 1)."""
    scaling_factors = [1.0]
    for i in range(len(line_segments) - 1):
        x1, y1 = line_segments[i]
        x2, y2 = line_segments[i + 1]
        mean_y1 = np.nanmean(y1[(x1 >= x2[0]) & (x1 <= x2[-1])])
        mean_y2 = np.nanmean(y2[(x2 >= x1[0]) & (x2 <= x1[-1])])
        scaling_factor = mean_y1 / mean_y2
        scaling_factors.append(scaling_factor)
        line_segments[i + 1] = (x2, y2 * scaling_factor)
    return line_segments, scaling_factors


def scale_bands_to_match(spectrum_input, bands=None):
    """Rescale each band so adjacent bands match in their overlap region.

    Applies independent factors to (``flux``, ``flux_error``) and to ``rf_flux``.
    """
    if bands is None:
        bands = BANDS
    spectrum = spectrum_input.copy()
    _, sf_flux = _rescale_line_segments(_create_line_segments(spectrum, bands, 'flux'))
    _, sf_rf = _rescale_line_segments(_create_line_segments(spectrum, bands, 'rf_flux'))

    for idx, band in enumerate(bands):
        mask = spectrum['BAND_ID'] == band
        if np.isfinite(sf_flux[idx]):
            spectrum.loc[mask, 'flux'] *= sf_flux[idx]
            spectrum.loc[mask, 'flux_error'] *= sf_flux[idx]
        if np.isfinite(sf_rf[idx]):
            spectrum.loc[mask, 'rf_flux'] *= sf_rf[idx]

    return spectrum, {'flux': sf_flux, 'rf_flux': sf_rf}


def stitch_bands(spectrum, wl_upper=27.5, bands=None):
    """Concatenate per-band spectra (in order ``bands``) into a single non-
    overlapping spectrum, keeping both ``flux`` and ``rf_flux`` from the first
    band that covers each wavelength. Drops wavelengths above ``wl_upper``
    (default 27.5 um, where MRS sensitivity collapses)."""
    if bands is None:
        bands = BANDS
    data = spectrum[spectrum['wavelength'] <= wl_upper].reset_index(drop=True)

    sw, sf, srf, se, sb = (np.array([]) for _ in range(5))
    for band in bands:
        sub = data[data['BAND_ID'] == band].reset_index(drop=True)
        if sub.empty:
            continue
        if sw.size == 0:
            sw  = np.append(sw,  sub['wavelength'])
            sf  = np.append(sf,  sub['flux'])
            srf = np.append(srf, sub['rf_flux'])
            se  = np.append(se,  sub['flux_error'])
            sb  = np.append(sb,  sub['BAND_ID'])
        else:
            mask = np.nanmax(sw) <= sub['wavelength']
            sw  = np.append(sw,  sub['wavelength'][mask])
            sf  = np.append(sf,  sub['flux'][mask])
            srf = np.append(srf, sub['rf_flux'][mask])
            se  = np.append(se,  sub['flux_error'][mask])
            sb  = np.append(sb,  sub['BAND_ID'][mask])

    return pd.DataFrame({
        'wavelength': sw,
        'flux': sf,
        'rf_flux': srf,
        'flux_error': se,
        'BAND_ID': sb,
    })


def find_x1d_files(input_dir, pattern='*_x1d.fits'):
    files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not files:
        raise FileNotFoundError(
            f"No files matching {pattern!r} in {input_dir!r}."
        )
    return files


def plot_stitched_spectrum(spectrum, source, output_path, flux_unit='Jy'):
    """Plot a stitched MRS spectrum and save as PNG.

    Single thin line on a linear-wavelength / log-flux axis so the full ~5-28 um
    dynamic range is legible (warm continuum + faint long-wavelength features).
    """
    # Local import keeps the FITS-only code paths usable without matplotlib.
    import matplotlib.pyplot as plt

    wl = spectrum['wavelength'].to_numpy()
    # Prefer the residual-fringe-corrected flux for the diagnostic plot; fall
    # back to the plain flux when rf_flux is missing or entirely NaN (e.g. for
    # MAST products that don't ship the RF columns).
    if 'rf_flux' in spectrum.columns and np.any(np.isfinite(spectrum['rf_flux'])):
        fl = spectrum['rf_flux'].to_numpy()
    else:
        fl = spectrum['flux'].to_numpy()

    fig, ax = plt.subplots(figsize=(9, 3.5), constrained_layout=True)
    ax.plot(wl, fl, color='black', linewidth=0.6)

    ax.set_xlabel(r'Wavelength [$\mu$m]')
    ax.set_ylabel(f'Flux density [{flux_unit}]')
    ax.set_title(f'{source} – MIRI MRS stitched spectrum')

    finite = np.isfinite(fl) & (fl > 0)
    if finite.any():
        ax.set_yscale('log')
        ymin = np.nanpercentile(fl[finite], 1)
        ymax = np.nanmax(fl[finite])
        ax.set_ylim(ymin * 0.5, ymax * 2)
    ax.set_xlim(np.nanmin(wl), np.nanmax(wl))

    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.grid(which='major', alpha=0.25, linewidth=0.5)

    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def post_process(input_dir, source, output_dir=None, wl_upper=27.5,
                 pattern='*_x1d.fits'):
    """End-to-end post-processing: read x1d files, rescale, stitch, write CSVs + PNG."""
    output_dir = Path(output_dir or input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    x1d_files = find_x1d_files(input_dir, pattern=pattern)
    print(f"Found {len(x1d_files)} x1d files in {input_dir}")
    for f in x1d_files:
        print(f"  {os.path.basename(f)}")

    spectrum_original = x1d_files_to_table(x1d_files)
    spectrum_matched, scaling = scale_bands_to_match(spectrum_original)
    spectrum_combined = stitch_bands(spectrum_matched, wl_upper=wl_upper)

    print("Flux scaling factors:   ", scaling['flux'])
    print("RF-flux scaling factors:", scaling['rf_flux'])

    full_csv = output_dir / f'{source}_spectrum_full.csv'
    stitched_csv = output_dir / f'{source}_spectrum_stitched.csv'
    stitched_png = output_dir / f'{source}_spectrum_stitched.png'
    spectrum_matched.to_csv(full_csv, index=False)
    spectrum_combined.to_csv(stitched_csv, index=False)
    plot_stitched_spectrum(spectrum_combined, source, stitched_png)
    print(f"Wrote {full_csv}")
    print(f"Wrote {stitched_csv}")
    print(f"Wrote {stitched_png}")

    return spectrum_matched, spectrum_combined


def parse_args(argv=None):
    """Parse command-line arguments."""
    summary = (__doc__ or '').splitlines()[0] if __doc__ else None
    p = argparse.ArgumentParser(description=summary)
    p.add_argument('input_dir',
                   help='Folder containing per-band x1d FITS files for one target.')
    p.add_argument('--source', required=True,
                   help='Source name; used as the CSV filename prefix.')
    p.add_argument('--output-dir', default=None,
                   help='Output directory for the CSVs (default: same as input_dir).')
    p.add_argument('--pattern', default='*_x1d.fits',
                   help='Glob pattern for x1d files inside input_dir '
                        '(default: %(default)s).')
    p.add_argument('--wl-upper', type=float, default=27.5,
                   help='Drop wavelengths above this value, in um '
                        '(default: %(default)s).')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    post_process(
        input_dir=args.input_dir,
        source=args.source,
        output_dir=args.output_dir,
        wl_upper=args.wl_upper,
        pattern=args.pattern,
    )


if __name__ == '__main__':
    main()
