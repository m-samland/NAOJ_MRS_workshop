#!/usr/bin/env python
# coding: utf-8

import jwst
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.associations import asn_from_list as afl
from jwst.pipeline import Detector1Pipeline, Spec2Pipeline, Spec3Pipeline
from jwst import __version__ as vjwst
import os
import multiprocessing
import numpy as np
import glob
import time
import logging
import sys
from datetime import datetime
from os.path import join, isdir, isfile, basename
from packaging import version
import pandas as pd
from astropy.io import fits
from stdatamodels.jwst import datamodels
import matplotlib.pyplot as plt

# Set path to cache of JWST pipeline
crds_cache_path = './crds_cache'

# Configure CRDS environment to fix connection issues
# import pathlib
# ROOT_DIR = pathlib.Path('~/').expanduser()
# crds_cache_dir = ROOT_DIR / "crds_cache"
# crds_cache_dir.mkdir(exist_ok=True)
# os.environ["CRDS_PATH"] = str(crds_cache_dir)
# os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# set the jwst pipeline crds context (leave blank to use latest)
# either leave blank or refer to a pmap file version >= 1094 (e.g. os.environ["CRDS_CONTEXT"] = "jwst_1190.pmap")
os.environ["CRDS_CONTEXT"] = ""


# Multiprocessing implementation to replace missing runmany function

def runmany(max_processes, func, file_list, *args, **kwargs):
    """
    Parallel execution function to replace the missing runmany from minds.mp

    Parameters:
    max_processes: maximum number of processes to use
    func: function to apply to each file
    file_list: list of files to process
    *args, **kwargs: additional arguments to pass to func
    """
    if max_processes == 1:
        # Sequential execution
        for file in file_list:
            try:
                func(file, *args, **kwargs)
            except Exception as e:
                print(f"Error processing {file}: {str(e)}")
                if enable_logging:
                    logger = logging.getLogger('jwst_reduction_master')
                    logger.error("Error processing %s: %s", file, str(e))
                continue
    else:
        # For now, use sequential execution to avoid multiprocessing pickle issues
        # TODO: Implement proper multiprocessing with worker functions
        print("Note: Running sequentially instead of parallel due to multiprocessing constraints")
        for file in file_list:
            try:
                func(file, *args, **kwargs)
            except Exception as e:
                print(f"Error processing {file}: {str(e)}")
                if enable_logging:
                    logger = logging.getLogger('jwst_reduction_master')
                    logger.error("Error processing %s: %s", file, str(e))
                continue


# directories
source = 'SYCha'
work_dir = f'../data/{source}/' # Write the full path here! Not relative path!
# Where the uncalibrated FITS files (for the science observation) are located.     Note: you can run the pipeline starting from rate files, by placing them in a work_dir+'stage1_004/' directory - in which case the definition of input_dir is not used.
input_dir = work_dir+'/stage0/'
# Where subdirectories corresponding to other processing stages will be created.
output_dir = work_dir+'/pipeline_standard/'
# Where the uncalibrated FITS files (for the background observation) are located.  If no background observation, leave this blank.
input_bgdir = work_dir+'/BG/stage0/'
# Where the output background observations should go.                              If no background observation, leave this blank
output_bgdir = work_dir+'/BG/'
# source name, only used in final output file names.

# Where final outputs will be located (cubes & spectra)
final_outdir = os.path.join(output_dir, 'final_outputs_bgsub/')

# suffix you'd like to add to output directories to better identify your reductions
ver = '_bgsub'

# Whether or not to run a given pipeline stage, and overwrite its products?
do_det1 = True  # Det1 Pipeline Stage & associated steps
do_spec2 = True  # Spec2 Pipeline Stage & associated steps
do_spec3 = True  # Spec3 Pipeline Stage & associated steps

# Simplified overwrite controls for the three main pipeline stages
# Overwrite Detector1Pipeline outputs? False will skip if outputs already exist
overwrite_det1 = True
# Overwrite Spec2Pipeline outputs? False will skip if outputs already exist
overwrite_spec2 = True
# Overwrite Spec3Pipeline outputs? False will skip if outputs already exist
overwrite_spec3 = True

# Multiprocessing:
usage = None  # {None, 'quarter', 'half', 'all'}. Takes precedence over maxp.
if usage is None:
    maxp = 4      # [if 'usage' is None] Set the maximum number of CPUs to use.

# Background subtraction - using annulus method only
bg_observation = True  # whether dedicated background observations or other reference observations are available to estimate the BKG. If True, these should be provided through the 'input_bgdir' parameter. If False, the pipeline will assume a dithering strategy was followed and will estimate the BKG from the dithered images themselves.

# Point source optimized reference files - param below are first used in Sec 4.1
# PSFF parameter removed - deprecated functionality that was never used

# Spectrum extraction + spike filtering - params below are first used in Sec 5.3
skip_spectral_leak = False  # whether to correct for spectral leakage
# Extract1d parameters
# whether [True] to use the 'autocen' feature of the pipeline, or [False] the centroid inferred in a weighted average frame for each band, where weights are proportional to the integrated flux (i.e., robust to faint sources)
autocen = True
# Float setting the size of the aperture in FWHM.
apsize = 2
# Whether to subtract the background in the spectrum extraction step
subtract_background = True
# Whether to apply the aperture correction in the spectrum extraction step
apply_apcorr = True
# Whether to overwrite target classification in header. This will influence the behaviour of default extract1d routine at the moment of extracting the spectrum: for an extended source, all spaxels are integrated, while for a point source, only spaxels in a small aperture are used.
overwrite_target_classification = True
# {'POINT','EXTENDED'} The new source type, if target classification is requested to change.
new_sourcetype = 'POINT'

# These parameters are for use of Danny Gasman's point source fringe correction code
# version of fringe to use, 01.05.00 is currently the only valid option
# ONLY USABLE WITH TARGET ACQUISITION
# Section from Gasman+25
do_G25 = False
# Give list of coordinates if known (omits centroiding, useful for fainter targets)
# Format for known coordinates is list in same order as list of files, nested per channel on the detector, e.g.:
# alpha_list = [[0.1,0.4], [0.3,0.2] ...]
# Put None if not known, otherwise alpha (along-slice) coordinate in arcsec
alpha_list = None
# Put None if not known, otherwise beta (across-slice) coordinate in arcsec
beta_list = None

if do_G25:
    # This should be skipped for the tutorial
    # Imports for Danny Gasman's fringe correction code
    # https://github.com/DannyG20/MIRI-MRS-Library-Pipeline
    # Works only for point sources with target acquisition
    import core.linearity_coeff_ref_files as lin_ref
    from core.distortion import d2cMapping
    from core.funcs import point_source_centroiding
    import core.flux_cal as phot
    G25_version = '01.05.00'
    fringedir = './references/FRINGE/'
    photomdir = './references/PHOTOM/'
    distdir = './references/DISTORTION/'
    distVers = 'flt8'
    linVers = 'v2'

# =======================================================================================
# Import packages for multiprocessing.  These won't be used on the online demo, but can be
# very useful for local data processing unless/until they get integrated natively into
# the cube building code.

num_cores = multiprocessing.cpu_count()
print('number of cores', num_cores)
if usage == 'quarter':
    maxp = num_cores // 4 or 1
elif usage == 'half':
    maxp = num_cores // 2 or 1
elif usage == 'all':
    maxp = num_cores

print('We will use '+str(maxp)+' CPUs in this run')

# Bad pixel correction/flagging
# Note: do_bpc1 parameter removed as it was deprecated and could only be False
# {'fit_profile', 'mingrad'} Algorithm used in PipelineStep pixel_replace.
pixel_replace_algo = 'mingrad'

# Detector1 pipeline parameters
# Whether to use more aggressive Det1 parameters to flag cosmic rays
use_agg_det1_params = False

# Spec3 options (params below first used in Sec 5.1)
skip_outlier_det_s3 = False  # Whether to skip outlier detection step in stage 3
# {'drizzle', 'emsm', 'msm'}. Method to combine different dithers in cube building.
dith_combi_method = 'drizzle'

# Run adaptive trace model code if desired to mitigate resampling noise, new feature in JWST pipeline, off by default.
# For extracting extended line emission or working on extended structures, this can be switched on to see if it works better, it should!
skip_adaptive_trace_model = True

# Logging configuration
enable_logging = True            # Master switch for logging
# Console log level: DEBUG, INFO, WARNING, ERROR, CRITICAL
log_level = 'WARNING'
log_to_console = True            # Also display logs in console
detailed_logging = True          # If True, file logs at DEBUG; else at log_level


# Detector1 pipeline outputs will go here
det1_dir = join(output_dir, 'stage1{}/'.format(ver))
# join(output_dir, 'stage1{}/'.format('_008_2024-02-08')) #Location of rate file inputs (if they already exist from a previous version)
det1_dir_ori = None
# Spec2 pipeline outputs will go here
spec2_dir = join(output_dir, 'stage2{}/'.format(ver))
# Spec3 pipeline outputs will go here
spec3_dir = join(output_dir, 'stage3{}/'.format(ver))
# Final Spec3 pipeline outputs will go here (final cubes and spectra)
fin3_dir = join(final_outdir, 'stage3{}/'.format(ver))
figs_dir = join(final_outdir, 'figures/')  # whether the figures will be save

# Logging directory and configuration
logs_dir = join(output_dir, 'logs/')  # Directory for all log files
timestamp = datetime.now().strftime('%Y_%m_%d_%H%M%S')
log_suffix = f"{timestamp}{ver}"

# **Do not touch anything below in this script - except if you know what you're doing**

# We need to check that the desired output directories exist, and if not create them
if not os.path.exists(det1_dir):
    os.makedirs(det1_dir)
if det1_dir_ori is None:
    det1_dir_ori = det1_dir
if not os.path.exists(spec2_dir):
    os.makedirs(spec2_dir)
if not os.path.exists(spec3_dir):
    os.makedirs(spec3_dir)
if not os.path.exists(fin3_dir):
    os.makedirs(fin3_dir)
if not os.path.exists(figs_dir):
    os.makedirs(figs_dir)
if enable_logging and not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

# Don't change these
firstframe_skip = False
lastframe_skip = False
rscd_skip = False

# list of authorized aperture sizes (in FWHM)
list_aps = [1.5, 2.0, 2.5, 3.0, 3.5]

# Import the base JWST package
vjwst = jwst.__version__
print("JWST pipeline version: ", vjwst)

if version.parse(vjwst) < version.parse("1.13.0"):
    raise ValueError(
        "Please update JWST package to a version larger or equal than 1.13.0")

# check the pmap version is > 1094, when not default
if os.environ["CRDS_CONTEXT"]:
    if int(os.environ["CRDS_CONTEXT"][-9:-5]) < 1095:
        raise ValueError("Please use a pmap version larger or equal than 1095")

# JWST pipelines (encompassing many steps)
# JWST pipeline utilities
# Function definitions
# Logging utility functions


def setup_master_logger():
    """Configure root logging per JWST docs and return an app logger.

    - Root logger captures everything at DEBUG
    - File handler writes time-stamped DEBUG (or chosen) logs to a file
    - Stream handler prints concise messages to stdout at the selected level
    """
    if not enable_logging:
        return logging.getLogger('dummy')

    # Configure root logger
    root_log = logging.getLogger()
    root_log.setLevel(logging.DEBUG)  # allow all messages; handlers filter
    root_log.handlers.clear()

    # File handler (captures JWST + app logs)
    file_level = logging.DEBUG if detailed_logging else getattr(
        logging, log_level.upper())
    file_path = join(logs_dir, f'pipeline_{log_suffix}.log')
    fh = logging.FileHandler(file_path)
    fh.setLevel(file_level)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    root_log.addHandler(fh)

    # Console handler (message-only)
    if log_to_console:
        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setLevel(getattr(logging, log_level.upper()))
        ch.setFormatter(logging.Formatter("%(message)s"))
        root_log.addHandler(ch)

    # Return a named logger for this script; messages will propagate to root
    return logging.getLogger('jwst_reduction_master')

# Define a function to check existing outputs (useful for overwrite parameter)


def check_output_exist(input_paths, output_dir, input_ending='uncal.fits', output_ending='rate.fits'):
    exists = []
    output_paths = [os.path.join(output_dir, os.path.basename(
        x).replace(input_ending, output_ending)) for x in input_paths]
    for out_path in output_paths:
        exists.append(isfile(out_path))
    exists = np.array(exists)
    return exists, output_paths


# First we'll define a function that will call the detector1 pipeline with our desired set of parameters
def rundet1(filename, outdir, use_agg_det1_params=False):
    print(filename)

    # skip if file exists and no overwrite is desired:
    expected_output = join(outdir, basename(
        filename).replace('_uncal.fits', '_rate.fits'))
    print("Expected output:", expected_output)

    # Use the recommended .call() syntax for pipeline execution
    # This eliminates the need for pipeline instantiation and configuration workarounds
    try:
        Detector1Pipeline.call(
            filename,
            output_dir=outdir,
            save_results=True,
            steps={'ramp_fit': {'algorithm': 'OLS_C'}},
        )
    except Exception as e:
        error_msg = f"Error processing {filename} in Detector1Pipeline: {str(e)}"
        print(f"ERROR: {error_msg}")
        if enable_logging:
            logger = logging.getLogger('jwst_reduction_master')
            logger.error(error_msg)
        raise


def rundet1_G25(list_of_files, input_dir, output_dir, rate_dir, lin_file,
                ramp_save=False):
    for i in range(len(list_of_files)):

        file = list_of_files[i]
        print(list_of_files[i])

        filesplit = file.split('_')
        ratefile = filesplit[0] + '_' + filesplit[1] + '_' + \
            filesplit[2] + '_' + filesplit[3] + '_rate.fits'
        num = int((filesplit[2])[-1])

        rate = fits.open(rate_dir+ratefile)

        if lin_file == 'interp':
            lin_ref.gen_custom_linearity(rate, distdir, crds_cache_path,
                                         num=num, file_ver=linVers, dist_ver=distVers)
            print('file generated')
            linearity_file = './references/LINEARITY/custom_linearity_ref_{}_{}.fits'.format(
                num, linVers)
        elif lin_file == 'grid':
            linearity_file = lin_ref.find_nearest_grid(rate, distdir, './references/LINEARITY/',
                                                       num=num, dist_ver=distVers)
            linearity_file = './references/LINEARITY/'+linearity_file

        rate.close()

        pipe = Detector1Pipeline()

        pipe.save_calibrated_ramp = ramp_save

        pipe.linearity.skip = False
        pipe.group_scale.skip = True
        pipe.dq_init.skip = False
        pipe.superbias.skip = False
        pipe.gain_scale.skip = True
#             pipe.emicorr.skip = False
        pipe.saturation.skip = False
        pipe.ipc.skip = True
        pipe.firstframe.skip = True      # <-- Changed
        pipe.lastframe.skip = False
        pipe.reset.skip = False
        pipe.saturation.n_pix_grow_sat = 0
        pipe.charge_migration.skip = True
        pipe.linearity.override_linearity = linearity_file
        pipe.rscd.skip = True      # <-- Changed
        pipe.dark_current.skip = False
        pipe.refpix.skip = False

        pipe.save_results = True
        pipe.output_dir = output_dir

        result = pipe.run(input_dir+file)


def runspec2(filename, outdir, pixel_replace_algo='mingrad',
             overwrite_target_classification=True, new_sourcetype='POINT', apsize=3.0):
    # Use the recommended .call() syntax for pipeline execution
    # This eliminates the need for pipeline instantiation and configuration workarounds

    # Set up the steps configuration for Spec2Pipeline
    steps_config = {
        'flat_field': {'skip': False},
        'bkg_subtract': {'skip': True},
        'pixel_replace': {
            'skip': False,
            'algorithm': pixel_replace_algo
        },
        'cube_build': {'skip': False},
        # Changed to get the 1d files! For BG sub.
        'extract_1d': {
            'skip': False,
            'subtract_background': False,
            'ifu_autocen': True,
            'ifu_rfcorr': False,
            'apply_apcorr': True,
            'ifu_rscale': apsize
        }
    }

    # Add conditional target classification override
    if overwrite_target_classification:
        steps_config['extract_1d']['ifu_set_srctype'] = new_sourcetype

    try:
        Spec2Pipeline.call(
            filename,
            output_dir=outdir,
            save_results=True,
            steps=steps_config,
        )
    except Exception as e:
        error_msg = f"Error processing {filename} in Spec2Pipeline: {str(e)}"
        print(f"ERROR: {error_msg}")
        if enable_logging:
            logger = logging.getLogger('jwst_reduction_master')
            logger.error(error_msg)
        raise


def runspec2_G25(scifile, input_dir, output_dir):
    d2cMaps = {}
    for b in ['1A', '1B', '1C', '2A', '2B', '2C', '3A', '3B', '3C', '4A']:
        d2cMaps[b] = d2cMapping(
            b, './references/DISTORTION/', slice_transmission='10pc', fileversion=distVers)

    fringe_file = {}
    alpha = {}
    beta = {}

    for i in range(len(scifile)):
        alpha[i] = {}
        beta[i] = {}
        fringe_file[i] = {}
        hdu = fits.open(input_dir + scifile[i])

        detector = hdu[0].header['DETECTOR']
        subbandl = hdu[0].header['BAND']
        dithdir = hdu[0].header['DITHDIRC']
        n = (scifile[i].split('_')[2])[-1]

        if subbandl == 'SHORT':
            subband = 'A'
        elif subbandl == 'MEDIUM':
            subband = 'B'
        else:
            subband = 'C'

        if detector == 'MIRIFUSHORT':
            band = ['1{}'.format(subband), '2{}'.format(subband)]
        else:
            band = ['3{}'.format(subband), '4{}'.format(subband)]

        data = hdu['SCI'].data
        hdu.close()

        for b in band:
            print(b, n)

            if b[0] in ['1', '2', '3', '4']:

                if b not in ['4B', '4C']:
                    if alpha_list is None:
                        try:
                            alpha[i][b], beta[i][b] = lin_ref.return_centroids(
                                data, b, d2cMaps[b])
                        except:
                            alpha[i][b], beta[i][b] = None, None
                            print('Centroiding failed, using central grid point...')
                    else:
                        print('Using predefined centroid...')
                        if b[0] in ['1', '3']:
                            ch_ind = 0
                        else:
                            ch_ind = 1

                        alpha[i][b] = alpha_list[i][ch_ind]
                        beta[i][b] = beta_list[i][ch_ind]

                    if b[0] in ['1', '2', '3']:
                        fringe_file[i][b], alpha[i][b], beta[i][b] = lin_ref.find_nearest_grid_fringe(
                            input_dir+scifile[i], alpha[i][b], beta[i][b], b, fringedir, G25_version)

    # Do fringe correction
    for i in range(len(scifile)):
        print(scifile[i])

        hdu = fits.open(input_dir + scifile[i])

        detector = hdu[0].header['DETECTOR']
        subbandl = hdu[0].header['BAND']
        dithdir = hdu[0].header['DITHDIRC']
        n = (scifile[i].split('_')[2])[-1]

        if subbandl == 'SHORT':
            subband = 'A'
        elif subbandl == 'MEDIUM':
            subband = 'B'
        else:
            subband = 'C'

        if detector == 'MIRIFUSHORT':
            band = ['1{}'.format(subband), '2{}'.format(subband)]
        else:
            band = ['3{}'.format(subband), '4{}'.format(subband)]

        data = hdu['SCI'].data
        err = hdu['ERR'].data

        hdu.close()

        if detector == 'MIRIFUSHORT':
            for b in band:
                fringe_hdu = fits.open(fringedir + fringe_file[i][b])
                fringe_flat = fringe_hdu[1].data
                err_fringe = fringe_hdu[2].data
                fringe_hdu.close()

                if b[0] == '1':
                    err[:, :512] = np.abs(data[:, :512]/fringe_flat[:, :512]) * np.sqrt(
                        (err[:, :512]/data[:, :512])**2 + (err_fringe[:, :512]/fringe_flat[:, :512])**2)
                    data[:, :512] = data[:, :512]/fringe_flat[:, :512]
                else:
                    err[:, 512:] = np.abs(data[:, 512:]/fringe_flat[:, 512:]) * np.sqrt(
                        (err[:, 512:]/data[:, 512:])**2 + (err_fringe[:, 512:]/fringe_flat[:, 512:])**2)
                    data[:, 512:] = data[:, 512:]/fringe_flat[:, 512:]
        else:
            fringe_hdu = fits.open(fringedir + fringe_file[i][band[0]])
            fringe_flat = fringe_hdu[1].data
            err_fringe = fringe_hdu[2].data
            fringe_hdu.close()

            err = np.abs(data/fringe_flat) * np.sqrt((err/data)
                                                     ** 2 + (err_fringe/fringe_flat)**2)
            data = data/fringe_flat

        hdu = fits.open(input_dir + scifile[i])
        hdu['SCI'].data = data
        hdu['ERR'].data = err

        hdu.writeto(output_dir + scifile[i], overwrite=True)

        hdu.close()

    # Apply the spectrophotometric calibration to the data
    for i in range(len(scifile)):
        file = scifile[i]

        dist_dir = './references/DISTORTION/'

        hdu = fits.open(input_dir+file)
        ifu = hdu[0].header['DETECTOR']
        subbandl = hdu[0].header['BAND']
        dithdirc = hdu[0].header['DITHDIRC']

        pipe = Spec2Pipeline()
        pipe.bkg_subtract.skip = True
        pipe.master_background_mos.skip = True

        pipe.fringe.skip = True
        pipe.residual_fringe.skip = True

        # generate pixel flat and find best photom file
        dithnum = (file.split('_')[2])[-1]
        if ifu == 'MIRIFULONG':
            channel = 'long'
        elif ifu == 'MIRIFUSHORT':
            channel = 'short'
            b = '1'

        if subbandl == 'SHORT':
            subband = 'A'
        elif subbandl == 'MEDIUM':
            subband = 'B'
        else:
            subband = 'C'

        phot.gen_custom_pixel_flat(input_dir+file, dist_dir, crds_cache_path,
                                   alpha[i], beta[i], dithnum, dist_ver=distVers, version=G25_version)
        pipe.flat_field.user_supplied_flat = './references/PIXEL/custom_pixel_ref_{}_{}.fits'.format(
            channel, dithnum)
        pipe.flat_field.skip = False
        pipe.straylight.skip = False

        if ifu == 'MIRIFULONG':
            pipe.photom.override_photom = photomdir + \
                '{}_{}_PS_PHOTOM_{}_{}.fits'.format(
                    ifu, subbandl, dithdirc, G25_version)
        else:
            pipe.photom.override_photom = photomdir + \
                phot.gen_photom(input_dir+file, dist_dir, alpha[i][b+subband],
                                beta[i][b+subband], dithnum, dist_ver=distVers, version=G25_version)

        pipe.photom.skip = False
        pipe.photom.mrs_time_correction = True
        pipe.cube_build.skip = True
        pipe.extract_1d.skip = True
        try:
            pipe.badpix_selfcal.skip = True
            pipe.nsclean.skip = True
        except:
            pass
        pipe.pixel_replace.skip = True

        pipe.output_dir = output_dir
        pipe.save_results = True
        result = pipe.run(output_dir + file)


def writel3asn(scifiles, asnfile, prodname, bg_files=None):
    # Define the basic association of science files
    asn = afl.asn_from_list(
        scifiles, rule=DMS_Level3_Base, product_name=prodname)

    if bg_files is not None:   # Add background files to the association
        nbg = len(bg_files)
        for ii in range(0, nbg):
            asn['products'][0]['members'].append(
                {'expname': bg_files[ii], 'exptype': 'background'})

    # Write the association to a json file
    _, serialized = asn.dump()
    with open(asnfile, 'w', encoding='utf-8') as outfile:
        outfile.write(serialized)


def runspec3(filename, outname, output_dir, coord_sys='skyalign', overwrite=True,
             pixel_replace_algo='mingrad', dith_combi_method='drizzle',
             overwrite_target_classification=True, new_sourcetype='POINT', apsize=3.0,
             subtract_background=True, skip_spectral_leak=False, apply_apcorr=True,
             skip_adaptive_trace_model=True):

    outlist = glob.glob(join(output_dir, '{}*_s3d.fits'.format(outname)))
    if (len(outlist) >= 12) and not overwrite:
        return None

    # Use the recommended .call() syntax for pipeline execution
    # This eliminates the need for pipeline instantiation and configuration workarounds

    # Set up the steps configuration for Spec3Pipeline
    steps_config = {
        'cube_build': {
            'output_type': 'band',
            'weighting': dith_combi_method,
            'coord_system': coord_sys
        },
        'pixel_replace': {
            'skip': False,
            'algorithm': pixel_replace_algo
        },
        'master_background': {'skip': False},
        'outlier_detection': {'skip': False},
        'spectral_leak': {'skip': skip_spectral_leak},
        'adaptive_trace_model': {'skip': skip_adaptive_trace_model},
        'extract_1d': {
            'skip': False,
            'subtract_background': subtract_background,
            'apply_apcorr': apply_apcorr,
            'ifu_autocen': True,
            'ifu_rfcorr': True,
            'ifu_rscale': apsize
        }
    }

    # Add conditional target classification override
    if overwrite_target_classification:
        steps_config['extract_1d']['ifu_set_srctype'] = new_sourcetype

    try:
        Spec3Pipeline.call(
            filename,
            output_dir=output_dir,
            steps=steps_config,
        )
    except Exception as e:
        error_msg = f"Error processing {filename} in Spec3Pipeline: {str(e)}"
        print(f"ERROR: {error_msg}")
        if enable_logging:
            logger = logging.getLogger('jwst_reduction_master')
            logger.error(error_msg)
        raise


def convert_to_native_endian(data):
    """Convert big-endian FITS data to numpy native format for pandas compatibility."""
    native_data = np.empty(len(data))
    native_data[:] = data
    return native_data


def x1d_files_to_table(x1d_fnames):
    spectra = []
    for fname_x1d in x1d_fnames:
        basename, _ = os.path.splitext(os.path.basename(fname_x1d))
        print(basename)
        # Initialize default values
        channel = None
        band = None
        band_num = None
        # Parse channel from filename
        if 'ch1' in basename:
            channel = 1
        elif 'ch2' in basename:
            channel = 2
        elif 'ch3' in basename:
            channel = 3
        elif 'ch4' in basename:
            channel = 4
        # Parse band from filename
        if 'long' in basename:
            band = 'long'
            band_num = 'C'
        elif 'medium' in basename:
            band = 'medium'
            band_num = 'B'
        elif 'short' in basename:
            band = 'short'
            band_num = 'A'
        # Try to read wavelength using datamodels first, then fall back to fits
        try:
            dm = datamodels.open(fname_x1d)
            wl = np.array(dm.spec[0].spec_table['WAVELENGTH'])
            channel = dm.meta.instrument.channel
        except Exception:
            hdu = fits.open(fname_x1d)
            wl = hdu[1].data['WAVELENGTH']
        # Read the spectrum data from FITS file
        hdul = fits.open(fname_x1d)
        hdul.verify('ignore')
        table = hdul[1].data
        # Convert all arrays to native endian format and create spectrum data dictionary
        spectrum_data = {
            'wavelength': convert_to_native_endian(wl),
            'flux': convert_to_native_endian(table['FLUX']),
            'flux_error': convert_to_native_endian(table['FLUX_ERROR']),
            'background': convert_to_native_endian(table['BACKGROUND']),
            'bkgd_error': convert_to_native_endian(table['BKGD_ERROR']),
            'rf_flux': convert_to_native_endian(table['RF_FLUX']),
            'rf_background': convert_to_native_endian(table['RF_BACKGROUND'])
        }
        hdul.close()
        # Create DataFrame from the spectrum data
        spectrum = pd.DataFrame(spectrum_data)
        # Add metadata columns
        spectrum['CHANNEL'] = channel
        spectrum['BAND'] = band
        spectrum['BAND_NUM'] = band_num
        band_id = f'{channel}{band_num}' if channel and band_num else None
        spectrum['BAND_ID'] = band_id
        spectrum['wavelength_index'] = spectrum.index
        spectra.append(spectrum)
    spectrum_df = pd.concat(spectra, ignore_index=True)
    return spectrum_df


def make_file_list(folder, prefix="Level3_psf_", suffix='_s3d.fits'):
    file_paths = []
    for channel in ['1', '2', '3', '4']:
        for band in ['short', 'medium', 'long']:
            file_path = os.path.join(
                folder, f"{prefix}ch{channel}-{band}{suffix}")
            file_paths.append(file_path)
    return file_paths


def rescale_line_segments(line_segments):
    # Sort the line segments by their start x value
    # line_segments.sort(key=lambda seg: seg[0][0])
    # Iterate over the line segments
    scaling_factors = [1]
    for i in range(len(line_segments) - 1):
        x1, y1 = line_segments[i]
        x2, y2 = line_segments[i + 1]
        # Check if the line segments overlap
        # print(x1[-1], x2[-1])
        # overlap = np.maximum(0, np.minimum(x1[-1], x2[-1]) - np.maximum(x1[0], x2[0]))
        # if overlap > 0:
        # print(i)
        # Calculate the mean y value in the overlapping region for the first line segment
        mean_y1 = np.nanmean(y1[(x1 >= x2[0]) & (x1 <= x2[-1])])
        # Calculate the mean y value in the overlapping region for the second line segment
        mean_y2 = np.nanmean(y2[(x2 >= x1[0]) & (x2 <= x1[-1])])
        # print(y2)
        # Rescale the y values of the second line segment such that the mean y value in the overlapping region is the same as for the first line segment
        scaling_factor = mean_y1 / mean_y2
        scaling_factors.append(scaling_factor)
        y2 = y2 * mean_y1 / mean_y2
        print(scaling_factor)
        # Update the second line segment
        line_segments[i + 1] = (x2, y2)
    return line_segments, scaling_factors


def create_line_segments(spectrum, bands, flux_column):
    """Create line segments for a given flux column."""
    line_segments = []
    for band in bands:
        line_segment = [
            spectrum.loc[spectrum['BAND_ID'] == band, 'wavelength'].values,
            spectrum.loc[spectrum['BAND_ID'] == band, flux_column].values,
        ]
        line_segments.append(line_segment)
    return line_segments


def scale_bands_to_match(spectrum_input, bands=['1A', '1B', '1C', '2A', '2B', '2C', '3A', '3B', '3C', '4A', '4B', '4C']):
    spectrum = spectrum_input.copy()
    # Create line segments for both flux types
    line_segments_flux = create_line_segments(spectrum, bands, 'flux')
    line_segments_rf_flux = create_line_segments(spectrum, bands, 'rf_flux')
    # Compute scaling factors for both flux types
    _, scaling_factors_flux = rescale_line_segments(line_segments_flux)
    _, scaling_factors_rf_flux = rescale_line_segments(line_segments_rf_flux)
    print("Flux scaling factors:", scaling_factors_flux)
    print("RF Flux scaling factors:", scaling_factors_rf_flux)
    # Apply scaling factors to respective columns
    for idx, band in enumerate(bands):
        if np.isfinite(scaling_factors_flux[idx]):
            spectrum.loc[spectrum['BAND_ID'] == band,
                         'flux'] *= scaling_factors_flux[idx]
            spectrum.loc[spectrum['BAND_ID'] == band,
                         'flux_error'] *= scaling_factors_flux[idx]
        if np.isfinite(scaling_factors_rf_flux[idx]):
            spectrum.loc[spectrum['BAND_ID'] == band,
                         'rf_flux'] *= scaling_factors_rf_flux[idx]
    return spectrum, {'flux': scaling_factors_flux, 'rf_flux': scaling_factors_rf_flux}


# def shift_by_radial_velocity(wavelength: Quantity, radial_velocity: Quantity):
#     if radial_velocity is None:
#         print("No radial velocity provided.")
#         return None
#     shifted_wavelength = (
#         wavelength + radial_velocity / const.c * wavelength
#     )
# #     shifted_wavelength = shifted_wavelength.to(self.spectrum[wavelength_column].unit)
#     print("Applying wavelength shift.")
#     return shifted_wavelength


def stitch_bands(spectrum, wl_upper=27.5):
    """Takes a spectrum data frame with BAND ID and stitches the spectra of different bands together. """
    Data = spectrum.copy()
    Data = Data[Data['wavelength'] <= wl_upper].reset_index(drop=True)
    SB1A, SB1B, SB1C = Data[Data['BAND_ID'] == '1A'].reset_index(drop=True), Data[Data['BAND_ID'] == '1B'].reset_index(
        drop=True), Data[Data['BAND_ID'] == '1C'].reset_index(drop=True)
    SB2A, SB2B, SB2C = Data[Data['BAND_ID'] == '2A'].reset_index(drop=True), Data[Data['BAND_ID'] == '2B'].reset_index(
        drop=True), Data[Data['BAND_ID'] == '2C'].reset_index(drop=True)
    SB3A, SB3B, SB3C = Data[Data['BAND_ID'] == '3A'].reset_index(drop=True), Data[Data['BAND_ID'] == '3B'].reset_index(
        drop=True), Data[Data['BAND_ID'] == '3C'].reset_index(drop=True)
    SB4A, SB4B, SB4C = Data[Data['BAND_ID'] == '4A'].reset_index(drop=True), Data[Data['BAND_ID'] == '4B'].reset_index(
        drop=True), Data[Data['BAND_ID'] == '4C'].reset_index(drop=True)
    Bands = [SB1A, SB1B, SB1C, SB2A, SB2B, SB2C,
             SB3A, SB3B, SB3C, SB4A, SB4B, SB4C]
    StitchedWav = np.array([])
    StitchedFlux = np.array([])
    StitchedError = np.array([])
    StitchedBand = np.array([])
    for i, Band in enumerate(Bands):
        if i == 0:
            StitchedWav = np.append(StitchedWav, Band['wavelength'])
            StitchedFlux = np.append(StitchedFlux, Band['rf_flux'])
            StitchedError = np.append(StitchedError, Band['flux_error'])
            StitchedBand = np.append(StitchedBand, Band['BAND_ID'])
        else:
            # Mask = (StitchedWav <= np.nanmin(Band['wavelength']))
            # StitchedWav = np.append(StitchedWav[Mask], Band['wavelength'])
            # StitchedFlux = np.append(StitchedFlux[Mask], Band['rf_flux'])
            # StitchedError = np.append(StitchedError[Mask], Band['rf_flux_error'])
            Mask = (np.nanmax(StitchedWav) <= Band['wavelength'])
            StitchedWav = np.append(StitchedWav, Band['wavelength'][Mask])
            StitchedFlux = np.append(StitchedFlux, Band['rf_flux'][Mask])
            StitchedError = np.append(
                StitchedError, Band['flux_error'][Mask])
            StitchedBand = np.append(StitchedBand, Band['BAND_ID'][Mask])
    stitched_spectrum = {'wavelength': StitchedWav, 'flux': StitchedFlux,
                         'flux_error': StitchedError, 'BAND_ID': StitchedBand}
    stitched_spectrum = pd.DataFrame(stitched_spectrum)
    return stitched_spectrum


def main():
    # Start timing
    time0 = time.perf_counter()

    # Set up logging
    master_logger = setup_master_logger()

    if enable_logging:
        master_logger.info("Starting JWST MIRI MRS reduction for %s", source)
        master_logger.info("Working directory: %s", work_dir)
        master_logger.info("Output directory: %s", output_dir)
        master_logger.info("JWST pipeline version: %s", vjwst)
        master_logger.info("Using %d CPU cores", maxp)

    # Output subdirectories to keep background data products organized
    if bg_observation:
        # Detector1 pipeline outputs will go here
        det1_bgdir = join(output_bgdir, 'stage1{}/'.format(ver))
        # Detector1 pipeline outputs will go here
        spec2_bgdir = join(output_bgdir, 'stage2{}/'.format(ver))
        print(det1_bgdir, spec2_bgdir, output_bgdir)
        # We need to check that the desired output directories exist, and if not create them
        if not os.path.exists(det1_bgdir):
            os.makedirs(det1_bgdir)
        if not os.path.exists(spec2_bgdir):
            os.makedirs(spec2_bgdir)
    else:
        det1_bgdir = None
        spec2_bgdir = None

    # Let's look for input files of the form *uncal.fits from the science observation
    if do_det1:
        # Simplified: always use uncal.fits files directly
        input_ending = 'uncal.fits'
        sstring = os.path.join(input_dir, 'jw*mirifu*uncal.fits')
        lvl1b_files = sorted(glob.glob(sstring))
        print('Found ' + str(len(lvl1b_files)) +
              ' science input files to process')

        if not overwrite_det1:
            exists, output_paths = check_output_exist(input_paths=lvl1b_files, output_dir=det1_dir_ori,
                                                      input_ending=input_ending, output_ending='rate.fits')
            # limit file list
            lvl1b_files = [lvl1b_files[i]
                           for i in range(len(lvl1b_files)) if not exists[i]]

        if len(lvl1b_files) > 0:
            if enable_logging:
                master_logger.info(
                    "Starting Detector1 pipeline processing for %d files", len(lvl1b_files))

            if maxp > 1:
                # Run the pipeline on these input files by a simple loop over our pipeline function
                runmany(maxp, rundet1, lvl1b_files, outdir=det1_dir,
                        use_agg_det1_params=use_agg_det1_params)
            else:
                for file in lvl1b_files:
                    rundet1(file, outdir=det1_dir,
                            use_agg_det1_params=use_agg_det1_params)

            if enable_logging:
                master_logger.info("Completed Detector1 pipeline processing")
    else:
        print('Skipping Detector1 processing')

    if do_G25 and do_det1:
        uncal_files = [os.path.basename(x) for x in glob.glob(
            f'{input_dir}*mirifushort_uncal.fits')]
        print(uncal_files)
        rundet1_G25(uncal_files, input_dir, det1_dir, det1_dir,
                    lin_file='grid', ramp_save=False)

    # Now let's look for input files of the form *uncal.fits from the background observation
    if do_det1 and bg_observation:
        # Simplified: always use uncal.fits files directly
        input_ending = 'uncal.fits'
        sstring = os.path.join(input_bgdir, 'jw*mirifu*'+input_ending)
        lvl1b_files = sorted(glob.glob(sstring))
        print('Found ' + str(len(lvl1b_files)) +
              ' background input files to process')

        if not overwrite_det1:
            exists, output_paths = check_output_exist(input_paths=lvl1b_files, output_dir=det1_bgdir,
                                                      input_ending=input_ending, output_ending='rate.fits')
            # limit file list
            lvl1b_files = [lvl1b_files[i]
                           for i in range(len(lvl1b_files)) if not exists[i]]

        if len(lvl1b_files) > 0:
            if maxp > 1:
                # Run the pipeline on these input files by a simple loop over our pipeline function
                runmany(maxp, rundet1, lvl1b_files, outdir=det1_bgdir,
                        use_agg_det1_params=use_agg_det1_params)
            else:
                for file in lvl1b_files:
                    rundet1(file, outdir=det1_bgdir)
    else:
        print('Skipping Detector1 processing for dedicated BKG observation')

    # Print out the time benchmark
    time1 = time.perf_counter()
    print(f"Runtime so far: {time1 - time0:0.4f} seconds")

    # For annulus background subtraction, we only use '_psf' suffix
    suffix = ['_psf']

    # Look for rate files from the Detector1 pipeline
    if do_spec2 and not do_G25:
        # DO PSF files
        psf_files = glob.glob(os.path.join(det1_dir, "jw*_rate.fits"))
        rate_file_list = [psf_files]
        # For annulus background method, we only process PSF files

        for rr, rate_files in enumerate(rate_file_list):
            # Simplified suffix handling - use rate files directly
            input_ending = 'rate.fits'
            output_ending = 'cal.fits'

            rate_files = np.array(rate_files)

            print('Found ' + str(len(rate_files)) +
                  ' input files to process for suffix: {}'.format(suffix[rr]))

            if not overwrite_spec2:
                exists, output_paths = check_output_exist(input_paths=rate_files, output_dir=spec2_dir,
                                                          input_ending=input_ending, output_ending=output_ending)
                # limit file list
                rate_files = [rate_files[i]
                              for i in range(len(rate_files)) if not exists[i]]

            if len(rate_files) > 0:
                if enable_logging:
                    master_logger.info("Starting Spec2 pipeline processing for %d files with suffix: %s", len(
                        rate_files), suffix[rr])

                if maxp > 1:
                    # Commented out: deprecated PSFF functionality
                    # runmany(maxp, runspec2_partial, rate_files)
                    runmany(maxp, runspec2, rate_files, outdir=spec2_dir,
                            pixel_replace_algo=pixel_replace_algo,
                            overwrite_target_classification=overwrite_target_classification,
                            new_sourcetype=new_sourcetype, apsize=apsize)
                else:
                    for file in rate_files:
                        runspec2(file, outdir=spec2_dir,
                                 pixel_replace_algo=pixel_replace_algo,
                                 overwrite_target_classification=overwrite_target_classification,
                                 new_sourcetype=new_sourcetype, apsize=apsize)

                if enable_logging:
                    master_logger.info(
                        "Completed Spec2 pipeline processing for suffix: %s", suffix[rr])
    else:
        print('Skipping Spec2 processing')

    if do_spec2 and bg_observation and not do_G25:
        # Look for background rate files from the Detector1 pipeline
        sstring = det1_bgdir + 'jw*mirifu*_rate.fits'
        bpcfiles = sorted(glob.glob(sstring))

        # Simplified suffix handling - use rate files as input
        input_ending = 'rate.fits'
        output_ending = 'cal.fits'

        bpcfiles = np.array(bpcfiles)
        print('Found ' + str(len(bpcfiles)) +
              ' input background files to process')

        if not overwrite_spec2:
            exists, output_paths = check_output_exist(input_paths=bpcfiles, output_dir=spec2_bgdir,
                                                      input_ending=input_ending, output_ending=output_ending)
            # limit file list
            bpcfiles = [bpcfiles[i]
                        for i in range(len(bpcfiles)) if not exists[i]]

        if len(bpcfiles) > 0:
            if maxp > 1:
                #             # Commented out: deprecated PSFF functionality
                #             runmany(maxp, runspec2_partial, bpcfiles)
                runmany(maxp, runspec2, bpcfiles, outdir=spec2_bgdir,
                        pixel_replace_algo=pixel_replace_algo,
                        overwrite_target_classification=overwrite_target_classification,
                        new_sourcetype=new_sourcetype, apsize=apsize)
            else:
                for file in bpcfiles:
                    runspec2(file, outdir=spec2_bgdir,
                             pixel_replace_algo=pixel_replace_algo,
                             overwrite_target_classification=overwrite_target_classification,
                             new_sourcetype=new_sourcetype, apsize=apsize)
    else:
        print('Skipping Spec2 processing for background files')

    if do_G25 and do_spec2:
        rate_files = [os.path.basename(x)
                      for x in glob.glob(f'{det1_dir}*_rate.fits')]
        runspec2_G25(rate_files, input_dir=det1_dir, output_dir=spec2_dir)

    # For annulus background method, we only use '_psf' suffix
    suffix = ['_psf']

    # Simplified suffix handling - use _cal.fits directly
    suff = ['{}_cal'.format(suffix[0])]

    # if needed, consider 2 lists for PSF and BKG, respectively
    if do_spec3:
        l3asn_dir = os.path.join(spec3_dir, "l3asn_files")
        if not isdir(l3asn_dir):
            os.makedirs(l3asn_dir)
        asnfiles = []
        asnfiles_ifua = []
        outnames = []
        outnames_ifua = []
        for s in range(len(suff)):
            sstring = os.path.join(spec2_dir, 'jw*_cal.fits')
            calfiles = np.array(sorted([os.path.abspath(f) for f in glob.glob(sstring)]))

            # Since do_bpc1 was deprecated and could only be False, we exclude '_bpc_' files
            calfiles = [calfiles[i] for i in range(
                len(calfiles)) if '_bpc_' not in calfiles[i]]

            print('Found ' + str(len(calfiles)) + ' science files to process')

            outname = 'Level3{}'.format(suffix[s])
            outname_ifua = 'Level3_ifua{}'.format(suffix[s])
            asnfile = os.path.join(l3asn_dir, 'l3asn{}.json'.format(suffix[s]))
            asnfile_ifua = os.path.join(
                l3asn_dir, 'l3asn_ifua{}.json'.format(suffix[s]))
            asnfiles.append(asnfile)
            asnfiles_ifua.append(asnfile_ifua)
            outnames.append(outname)
            outnames_ifua.append(outname_ifua)
            if bg_observation:
                # These have to be the x1d files!!
                bg_files = [os.path.abspath(f) for f in sorted(glob.glob(spec2_bgdir+'jw*x1d.fits'))]
                writel3asn(calfiles, asnfile, outname, bg_files=bg_files)
                writel3asn(calfiles, asnfile_ifua,
                           outname_ifua, bg_files=bg_files)
            else:
                writel3asn(calfiles, asnfile, outname)
                writel3asn(calfiles, asnfile_ifua, outname_ifua)

    # now let's run spec3
    if do_spec3:
        if enable_logging:
            master_logger.info("Starting Spec3 pipeline processing")

        for aa, asnfile in enumerate(asnfiles):
            runspec3(
                asnfile,
                outname=outnames[aa],
                output_dir=spec3_dir,
                overwrite=overwrite_spec3,
                coord_sys='skyalign',
                pixel_replace_algo=pixel_replace_algo,
                dith_combi_method=dith_combi_method,
                overwrite_target_classification=overwrite_target_classification,
                new_sourcetype=new_sourcetype,
                apsize=apsize,
                subtract_background=subtract_background,
                skip_spectral_leak=skip_spectral_leak,
                apply_apcorr=apply_apcorr,
                skip_adaptive_trace_model=skip_adaptive_trace_model,
            )
            runspec3(
                asnfiles_ifua[aa],
                outname=outnames_ifua[aa],
                output_dir=spec3_dir,
                coord_sys='ifualign',
                pixel_replace_algo=pixel_replace_algo,
                dith_combi_method=dith_combi_method,
                overwrite_target_classification=overwrite_target_classification,
                new_sourcetype=new_sourcetype,
                apsize=apsize,
                subtract_background=subtract_background,
                skip_spectral_leak=skip_spectral_leak,
                apply_apcorr=apply_apcorr,
                skip_adaptive_trace_model=skip_adaptive_trace_model,
            )  # ifualign

        if enable_logging:
            master_logger.info("Completed Spec3 pipeline processing")

        # Print out the time benchmark
        time1 = time.perf_counter()
        print(f"Runtime so far: {time1 - time0:0.4f} seconds")

        if enable_logging:
            master_logger.info("Total runtime: %.2f seconds", time1 - time0)
            master_logger.info(
                "JWST MIRI MRS reduction completed successfully")

    # Last steps stitching together the spectra
    x1d_files = glob.glob(f'{output_dir}stage3{ver}/Level3_ifua*_x1d.fits')
    spectrum_original = x1d_files_to_table(x1d_files)
    spectrum_matched, flux = scale_bands_to_match(spectrum_original)
    spectrum_combined = stitch_bands(spectrum_matched, wl_upper=27.5)

    spectrum_matched.to_csv(f'{final_outdir}/{source}_spectrum_full.txt')
    spectrum_combined.to_csv(f'{final_outdir}/{source}_spectrum_stitched.txt')

    # Final time
    time1 = time.perf_counter()
    print(f'Elapsed time: {(time1-time0)/60:.2f}min')


if __name__ == "__main__":
    main()
