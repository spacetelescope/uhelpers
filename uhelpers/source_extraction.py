"""Functions to call external source extraction packages, e.g. source extractor, SEP

Authors
-------

    Johannes Sahlmann

Use
---

"""

import glob
import os

import numpy as np


def system_command(sysStr):
    ret = os.system(sysStr)
    if (ret != 0):
        print(sysStr)
        raise RuntimeError('Error in system command:\n{}'.format(ret))

def run_source_extractor(file_list, source_extractor_config_dir, config_prefix=None, file_dir=None,
                         source_extractor_source_dir=None, weight=None, overwrite=0, use_psfex=0, fits_extension=None,
                         out_dir=None, prep_saturation_level=50000., saturation_level=60000., verbose=False):
    """System call to a source extractor executable.

    written 2014-11-17 JSA ESAC
    updated 2017-07-06 J. Sahlmann STScI

    Parameters
    ----------
    file_list
    source_extractor_config_dir : str
        directory containing the parameter and configuration files
    config_prefix
    file_dir
    source_extractor_source_dir : str
        directory containing the executables
    weight
    overwrite
    use_psfex
    fits_extension
    out_dir
    prep_saturation_level
    saturation_level
    verbose

    Returns
    -------

    """
    if config_prefix is None:
        config_prefix = 'default'

    sex_config_file = os.path.join(source_extractor_config_dir, '%s.sex_config' % config_prefix)
    sex_param_file = os.path.join(source_extractor_config_dir, '%s.sex_param' % config_prefix)

    # convolution filter for source detection in main source extraction
    sex_filter_file = glob.glob(os.path.join(source_extractor_config_dir, '%s*.sex_conv' % config_prefix))[0]

    # convolution filter for source detection in source extraction for psfex
    prep_filter_file = sex_filter_file

    if fits_extension is not None:
        extension_string = '[%d]' % fits_extension
    else:
        extension_string = ''

    if file_dir is None:
        file_dir = os.path.dirname(file_list[0])
        file_list = [os.path.basename(f) for f in file_list]

    original_weight = weight

    current_dir = os.getcwd()

    if out_dir is not None:
        out_dir = out_dir
    else:
        out_dir = file_dir

    for i, file in enumerate(file_list):
        print('')
        print('='*50)
        print('Working on %s: ' % file)

        ofile = os.path.join(out_dir, os.path.basename(file))
        name_cat = ofile.replace('.fits', '.cat')
        name_prep_cat = ofile.replace('.fits', '_prep.cat')
        name_psf = ofile.replace('.fits', '_prep.psf')
        name_psfxml = ofile.replace('.fits', '_psfex.xml')
        name_prepsfxml = ofile.replace('.fits', '_prepsfex.xml')
        name_sexxml = ofile.replace('.fits', '_sex.xml')

        if source_extractor_source_dir is not None:
            os.chdir(source_extractor_source_dir)
        else:
            source_extractor_source_dir = ''

        if use_psfex:

            prepare_for_psfex_config_file = os.path.join(source_extractor_config_dir, '%s_prepare_for_psfex.sex_config' % config_prefix)
            prepare_for_psfex_param_file = os.path.join(source_extractor_config_dir, '%s_prepare_for_psfex.sex_param' % config_prefix)

            psfex_config_file = os.path.join(source_extractor_config_dir, '%s.psfex_config' % config_prefix)

            if (not os.path.isfile(name_psf)) or (overwrite):

                sysString = '%s/sex %s%s -c %s -CATALOG_NAME %s -SATUR_LEVEL %3.1f -XML_NAME %s -WEIGHT_TYPE NONE -PARAMETERS_NAME %s -FILTER_NAME %s' % (
                source_extractor_source_dir, os.path.join(file_dir, file), extension_string,
                prepare_for_psfex_config_file,
                name_prep_cat,
                prep_saturation_level, name_prepsfxml,
                prepare_for_psfex_param_file, prep_filter_file)

                if verbose:
                    print('\nexecuting first sex command {}\n\n'.format(sysString))
                system_command(sysString)

                checkimage_file_names = np.array([os.path.join(out_dir, '%s_' % s + os.path.basename(file)) for s in
                                                  ['psf_prototype', 'psf_samples', 'psf_residuals', 'psf_snapshots']])

                # Compute the PSF with PSFex
                sysString = '%s/psfex %s -c %s -XML_NAME %s -CHECKIMAGE_NAME %s,%s,%s,%s' % (
                source_extractor_source_dir,
                name_prep_cat,
                psfex_config_file, name_psfxml, checkimage_file_names[0],
                checkimage_file_names[1], checkimage_file_names[2], checkimage_file_names[3])
                system_command(sysString)
                if verbose:
                    print('\nexecuting psfex command {}\n\n'.format(sysString))

                # # Extract final astrometry and photometry
                if weight is not None:
                    sysString = '%s/sex %s -c %sdefaultpsf.sex -CATALOG_NAME %s -PSF_NAME %s -XML_NAME %s -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE %s' % (
                    source_extractor_source_dir, file, source_extractor_source_dir, name_cat, name_psf, name_sexxml, weight)
                else:
                    sysString = '%s/sex %s%s -c %s -CATALOG_NAME %s -PSF_NAME %s -XML_NAME %s -WEIGHT_TYPE NONE -PARAMETERS_NAME %s -FILTER_NAME %s -SATUR_LEVEL %3.1f' % (
                    source_extractor_source_dir, os.path.join(file_dir, file), extension_string, sex_config_file,
                    name_cat, name_psf, name_sexxml, sex_param_file, sex_filter_file, saturation_level)
                    os.remove(name_prep_cat)

                if verbose:
                    print('\nexecuting second sex command {}\n\n'.format(sysString))

                sexout = system_command(sysString)

        else:
            if (not os.path.isfile(name_cat)) or (overwrite):
                sysString = '%s/sex %s%s -c %s -CATALOG_NAME %s -XML_NAME %s -WEIGHT_TYPE NONE -FILTER_NAME %s' % (
                source_extractor_source_dir, os.path.join(file_dir, file), extension_string, sex_config_file, name_cat,
                name_sexxml, sex_filter_file)

                if verbose:
                    print(sysString)
                system_command(sysString)

        weight = original_weight

    os.chdir(current_dir)
    return
