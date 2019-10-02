"""Helpers for interfacing astronomical archives.

Authors
-------

    Johannes Sahlmann

Use
---


"""

import netrc
import os
import warnings

import astropy.units as u
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import ascii
from astroquery.simbad import Simbad
import numpy as np

import pygacs.public.publicAccessTools as pgp
import pygacs.authen.manip as pga


CASJOBS_BASE_URL = 'http://mastweb.stsci.edu/gcasjobs/services/jobs.asmx'


def get_exoplanet_orbit_database(data_dir, overwrite=False, keep_only_most_massive_planet=True,
                                 verbose=True):
    """Return an astropy table with the exoplanet orbit database (exoplanets.org).

    Parameters
    ----------
    data_dir : str
        output directory
    overwrite : bool
    keep_only_most_massive_planet : bool
        if system contains several planets, return only the most massive one

    Returns
    -------
    eod_table : astropy.table
        The orbit database

    """
    eod_file = os.path.join(data_dir, 'exoplanet_orbit_database_table.csv')
    if not os.path.isfile(eod_file) or overwrite:

        from astroquery.exoplanet_orbit_database import ExoplanetOrbitDatabase

        eod_table = ExoplanetOrbitDatabase.get_table()

        gaia_dr2_source_ids = []
        for simbad_id in eod_table['SIMBADNAME']:
            value = 0
            result_table = None
            if not np.ma.is_masked(simbad_id) and simbad_id not in ['Qatar-1', 'Qatar-2', 'TYC', 'KIC', 'CoRoT 12', 'CoRoT-Exo 12']:
                if verbose:
                    print('Working on {}'.format(simbad_id))
                result_table = Simbad.query_objectids(simbad_id)
                try:
                    value = [np.int(s.split('Gaia DR2')[1]) for s in list(result_table['ID']) if
                     'Gaia DR2' in s][0]
                except IndexError:
                    pass

                if '114762' in simbad_id:
                    value = 3937211745904553600

            gaia_dr2_source_ids.append(value)

        eod_table['gaia_dr2_source_id'] = gaia_dr2_source_ids
        remove_index = np.where(eod_table['gaia_dr2_source_id'] == 0)[0]
        eod_table.remove_rows(remove_index)

        if keep_only_most_massive_planet:
            # keep only most massive planet
            for source_id in eod_table['gaia_dr2_source_id']:
                index = np.where(eod_table['gaia_dr2_source_id'] == source_id)[0]
                if len(index) != 1:
                    remove_index = None
                    remove_index = index[np.where(eod_table['MASS'][index] < np.max(eod_table['MASS'][index]))[0]]
                    eod_table.remove_rows(remove_index)

            # remove planets without mass estimate
            for source_id in eod_table['gaia_dr2_source_id']:
                index = np.where(eod_table['gaia_dr2_source_id'] == source_id)[0]
                if len(index) != 1:
                    remove_index = None
                    if len(np.unique(eod_table['MASS'][index])) == 1:
                        remove_index = index
                    else:
                        remove_index = index[np.where(eod_table['MASS'][index] == 0)[0]]
                    eod_table.remove_rows(remove_index)

        eod_table.write(eod_file, format='ascii.fixed_width', delimiter=',', bookend=False, overwrite=True)
    else:
        eod_table = Table.read(eod_file, format='ascii.fixed_width', delimiter=',')

    return eod_table


def return_vizier_catalog_as_table(out_dir, catalog_string, catalog_name=None, max_row_number=1e6,
                                   overwrite=False):
    """QUERY VIZIER CATALOGUE using vizquery command line tools and return astropy table.

    Parameters
    ----------
    out_dir
    catalog_string
    catalog_name
    max_row_number
    overwrite

    Returns
    -------

    """
    queryfile = os.path.join(out_dir, '{}_query.txt'.format(catalog_name))
    fits_file = queryfile.replace('.txt', '.fits')
    if not(os.path.isfile(queryfile)) | overwrite:
        query_string = 'vizquery -mime=csv -site=http://cdsarc.u-strasbg.fr -source=%s -out.all -out.max=%d > %s' % (catalog_string, max_row_number, queryfile)
        os.system(query_string)

    if not (os.path.isfile(fits_file)) | overwrite:
        t = Table.read(queryfile, format='ascii.basic', delimiter=';', data_start=3, guess=False, header_start=0)
        t.write(fits_file, format='fits', overwrite=overwrite)
    else:
        t = Table.read(fits_file, format='fits')
    return t


def return_gacs_query_as_table(query_string, output_file_seed, overwrite=False, verbose=True):
    """Run query on Gaia archive (GACS) and return result as astropy table.

    Parameters
    ----------
    query_string : str
        ADQL query
    output_file_seed : str
        naming seed for results
    overwrite : bool
    verbose : bool

    Returns
    -------
    d : astropy.table.Table


    """
    if (overwrite) | (not os.path.isfile(output_file_seed+'.vot')):
        pgp.retrieveQueryResult(query_string, output_file_seed+'.vot')
        d = Table.read(output_file_seed+'.vot', format='votable')

        problematic_columns = [c.name for c in d.columns.values() if c.dtype=='O'] # dtype = Object causes problems for fits files
        for colname in problematic_columns:
            if colname in d.colnames:
                tmp = np.array(d[colname]).astype(np.str)
                d.remove_column(colname)
                d[colname] = tmp
        d.write(output_file_seed+'.fits', overwrite=True)
    else:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', AstropyWarning)
            d = Table.read(output_file_seed+'.fits')
            # d = Table.read(output_file_seed + '.vot', format='votable')
    if verbose:
        print('Retrieved {} sources from Gaia catalog in this field'.format(len(d)))
    return d


def gacs_list_query(username, password, out_dir, input_table, input_table_name, gacs_table_name, id_str_gacs_table,
                    id_str_input_table, verbose=False, dr_string='gaiadr2', overwrite=True, outputFileName=None):
    """Upload file that contains list of target identifiers to Gaia archive (GACS),
    and crossmatch it with given GACS table using the identifier.

    requires GACS account (username+password)

    input_table has to be in votable format


    Parameters
    ----------
    username
    password
    out_dir
    input_table
    input_table_name
    gacs_table_name
    id_str_gacs_table
    id_str_input_table
    verbose
    dr_string
    overwrite
    outputFileName

    Returns
    -------
    astropy table with query result

    TODO
    ----
        check that requested gacs table field exists before executing query

    """
    xmlFileName = os.path.join(out_dir, 'gacsTableProperties.xml')
    tableProps = pga.GacsTableProperties(username, password, xmlFileName)
    if verbose:
        tableProps.printSchemaNames()

    # get and print available table in a specific schema
    publicTableNames = tableProps.getTableNames('public', verbose=0)

    # the user's schema name
    userSchemaName = 'user_%s' % username

    # check whether user table or xmatch table already exist in GACS. if yes delete them
    if userSchemaName in tableProps.schemaNames:
        userTableNames = tableProps.getTableNames(userSchemaName, verbose=1)

        if (input_table_name not in userTableNames):
            #   upload user table
            if verbose:
                print("*********************  Uploading Table")
            pga.authenticatedGacsCommand(username, password, pga.str_uploadTable(input_table, input_table_name))

        elif (overwrite):
            if verbose:
                print("*********************  Deleting Table")
            pga.authenticatedGacsCommand(username, password, pga.str_deleteTable(input_table_name))
            if verbose:
                print("*********************  Uploading Table")
            pga.authenticatedGacsCommand(username, password, pga.str_uploadTable(input_table, input_table_name))

    # construct query
    queryString = '''
    SELECT * FROM
    (select * FROM %s.%s) AS t
    INNER JOIN
    (select * FROM user_%s.%s) AS t2
    ON (t.%s = t2.%s)    
    ; ''' % (dr_string, gacs_table_name, username, input_table_name, id_str_gacs_table, id_str_input_table)

    if outputFileName is None:
        outputFileName = os.path.join(out_dir, '%s_%s_matched_on_%s.vot' % (
        gacs_table_name, input_table_name, id_str_gacs_table))

    if overwrite | (not os.path.isfile(outputFileName)):
        pga.authenticatedQuery(username, password, queryString, outputFileName, retrieve=True)

    T = Table.read(outputFileName, format='votable')
    return T


def execute_casjobs_query(userid, password, query, table_name, out_file,
                          overwrite_casjobs_query=False, verbose=True, download_result=True,
                          context='TwoMassNew'):
    """Execute query on casjobs, return status, and download data if query has finished.

    casjobs status codes (from http://casjobs.sdss.org/casjobs/services/jobs.asmx?op=GetJobStatus)
         0 = ready
         1 = started
         2 = canceling
         3 = cancelled
         4 = failed
         5 = finished

    Parameters
    ----------
    userid
    password
    query
    table_name
    out_file
    overwrite_casjobs_query
    verbose
    download_result
    context

    Returns
    -------
    job_is_ready : int
        Flag 0/1

    """
    casjobs_status_codes = np.array(['ready','started','canceling','cancelled','failed','finished'])

    job_is_ready = 1
    submit_query = 0

    if verbose:
        print(query)

    #         see https://galex.stsci.edu/casjobs/download/casjobs.config.x
    jobs = casjobs.CasJobs(userid=userid, password=password,
                           base_url=CASJOBS_BASE_URL)

    # get information on all previous jobs
    job_info = jobs.job_info()

    # reverse order (latest first)
    job_info = job_info[::-1]

    job_index = find_latest_query(query, job_info, verbose=verbose)
    if job_index is not None:
        print('Status is %s=%s' % (job_info[job_index]['Status'], casjobs_status_codes[np.int(job_info[job_index]['Status'])]))

    if job_index is None:
        # query has never been submitted
        submit_query = 1
    else:
        job_status = job_info[job_index]['Status']
        # if job still in the works, do nothing
        if job_status in ['0', '1']:
            if verbose:
                print('Query has not finished')
            job_is_ready = 0
        elif (job_info[job_index]['Status'] in ['5']) and (not overwrite_casjobs_query):
            # query has finished, let's download the data
            if download_result:
                if verbose:
                    print('Query has finished, let\'s download the data')
                try:
                    job_id = jobs.request_and_get_output(table_name, 'FITS', out_file)
                except Exception as e:
                    raise Exception('casjobs table download failed with message: {}\nPlease check that table {} is not empty.'.format(e, table_name))
            else:
                print('Query has finished and is ready for download')
        elif job_info[job_index]['Status'] in ['2', '3', '4']:
            submit_query = 1

    if overwrite_casjobs_query:
        submit_query = 1

    if submit_query == 1:
        try:
            # if table already existed, drop it
            jobs.drop_table(table_name)
        except:
            pass
        # submit job
        job_id = jobs.submit(query, context=context)
        if verbose:
            print('Submitted job')
        job_is_ready = 0

    return job_is_ready


def find_latest_query(query, job_info, verbose=True):
    """Return the most recent index of casjobs query in job_info if it exists, otherwise return None

    Parameters
    ----------
    query
    job_info

    Returns
    -------

    """
    job_index = None

    # look whether query has already been submitted as job
    index_list = []
    for jj in range(len(job_info)):
        if job_info[jj]['Query'].strip() == query.strip():
            index_list.append(jj)

    numer_of_matching_queries = len(index_list)
    if numer_of_matching_queries != 0:
        if verbose:
            print('Query has already been submitted {} times'.format(len(index_list)))
        if numer_of_matching_queries == 1:
            job_index = index_list[0]
        else:
            time_submit = np.zeros(numer_of_matching_queries)
            for j,index in enumerate(index_list):
                time_submit[j] = Time(job_info[index]['TimeSubmit'], format='isot').mjd
            latest_index = np.argmax(time_submit)
            job_index = index_list[latest_index]
            if verbose:
                print('Most recent query was submitted on {}'.format(job_info[job_index]['TimeSubmit']))

    return job_index


def inspect_casjobs_query(userid, password, query, verbose=True):
    """Inspect the status of a query on casjobs server and return the index of the most recent job

    Parameters
    ----------
    userid
    password
    query

    Returns
    -------
    job_index : int
        index of the most recent job in


    """
    jobs = casjobs.CasJobs(userid=userid, password=password,
                           base_url=CASJOBS_BASE_URL)

    # get information on all previous jobs
    job_info = jobs.job_info()

    # reverse order (latest first)
    job_info = job_info[::-1]

    job_index = find_latest_query(query, job_info, verbose=verbose)

    return job_index


def get_gaia_sources(RA_center_deg, DE_center_deg, search_radius, out_dir, overwrite=False,
                     tag='', catalog_name='gaiadr1.gaia_source', verbose=False):
    """Retrieve Gaia sources from Gaia archive within circular region.

    Parameters
    ----------
    RA_center_deg
    DE_center_deg
    search_radius
    out_dir
    overwrite
    tag
    catalog_name
    verbose

    Returns
    -------

    """
    query_string = '''SELECT * FROM {} WHERE 1 = CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {}, {}, {}))'''.format(
        catalog_name, RA_center_deg, DE_center_deg, search_radius.to(u.deg).value)


    outputFileSeed = os.path.join(out_dir, '%s_%s_gaia_query_result_area%2.3f' % (catalog_name, tag, (search_radius.to(u.deg).value*2)**2))
    if verbose:
        print(query_string)
    gaia_source = return_gacs_query_as_table(query_string, outputFileSeed, overwrite=overwrite)

    gaia_source_ra = np.array(gaia_source['ra'])
    gaia_source_de = np.array(gaia_source['dec'])
    gaia_cat = SkyCoord(ra=gaia_source_ra * u.degree, dec=gaia_source_de * u.degree)

    return gaia_source, gaia_cat

