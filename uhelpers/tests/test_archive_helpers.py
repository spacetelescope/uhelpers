#!/usr/bin/env python
"""Tests for the jwcf hawki module.

Authors
-------
    Johannes Sahlmann

"""
import netrc
import os

from astropy.table import Table
import pytest

from ..archive_helpers import get_exoplanet_orbit_database, gacs_list_query

local_dir = os.path.dirname(os.path.abspath(__file__))

ON_TRAVIS = os.environ.get('TRAVIS') == 'true'

@pytest.mark.skipif(ON_TRAVIS, reason='timeout issue.')
def test_eod():
    """Test the access to the exoplanet orbit database."""
    catalog = get_exoplanet_orbit_database(local_dir, verbose=False)

    assert len(catalog) > 100

@pytest.mark.skipif(ON_TRAVIS, reason='Requires access to .netrc file.')
def test_gacs_list_query():
    # print('test gacs list query')
    # Define which host in the .netrc file to use
    HOST = 'http://gea.esac.esa.int'
    # Read from the .netrc file in your home directory
    secrets = netrc.netrc()
    username, account, password = secrets.authenticators(HOST)

    out_dir = os.path.dirname(__file__)
    T = Table()
    id_str_input_table = 'ID_HIP'
    T[id_str_input_table] = [1, 2, 3, 4, 5, 6, 7]

    gacs_table_name = 'tgas_source'
    id_str_gacs_table = 'hip'

    input_table_name = 'hip_star_list'
    input_table = os.path.join(out_dir, 'hip_star_list.vot')
    T[[id_str_input_table]].write(input_table, format='votable', overwrite=1)

    T_out = gacs_list_query(username, password, out_dir, input_table, input_table_name, gacs_table_name,
                            id_str_gacs_table, id_str_input_table)
    T_out.pprint()
