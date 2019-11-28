"""Helper functions and classes of general use in astronomy

Authors
-------

    Johannes Sahlmann

Use
---

"""

import os
from astropy.table import Table
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
import numpy as np


class AstroSource(object):
    """Class for astronomical sources.

    Notes
    -----
        from astroquery.simbad import Simbad
        Simbad.list_votable_fields()
    """

    def __init__(self, identifier):
        self.identifier = identifier

    def __repr__(self):
        return '<AstroSource object id={0}>'.format(self.identifier)

    def set_simbad_fields(self, simbad_id, out_dir, overwrite=False,
                          votable_fields=('ra(d)', 'dec(d)', 'pmdec', 'pmra', 'parallax', 'sptype', 'ids')):
        """Retrieve source information from Simbad.

        Parameters
        ----------
        simbad_id
        out_dir
        overwrite
        votable_fields

        Returns
        -------

        """

        out_file = os.path.join(out_dir, '{}_simbad_parameters.txt'.format(self.identifier))
        if (not (os.path.isfile(out_file))) | overwrite:
            if os.path.isdir(out_dir) is False:
                os.makedirs(out_dir)
            simb = Simbad()
            if votable_fields is not None:
                simb.add_votable_fields(*votable_fields)
            pt = simb.query_object(simbad_id)
            pt.write(out_file, format='ascii.basic', delimiter=',')
        else:
            pt = Table.read(out_file, format='ascii.basic', delimiter=',')

        for c in pt.colnames:
            try:
                setattr(self, c, pt[c][0])
            except:
                setattr(self, c, None)


    def set_gaia_source_id(self):
        """Add the Gaia source_id attribute."""

        if hasattr(self, 'IDS') is False:
            raise RuntimeError('RUN set_simbad_fields() first!')
        try:
            self.simbad_identifiers = self.IDS.decode().split('|')
        except AttributeError:
            self.simbad_identifiers = self.IDS.split('|')

        self.gaia_dr2_id = None
        self.gaia_dr1_id = None
        try:
            self.gaia_dr2_id = [np.int(id.split(' ')[2].replace("'", "").replace("\"", "")) for id in self.simbad_identifiers if 'Gaia DR2 ' in id][0]
        except IndexError:
            try:
                self.gaia_dr1_id = [id for id in self.simbad_identifiers if 'Gaia DR1 ' in id][0]
            except IndexError:
                print('No Gaia identifier listed in Simbad!')


    def add_gaia(self):
        """Retrieve and add Gaia DR2 or DR1 parameters."""

        if hasattr(self, 'IDS') is False:
            raise RuntimeError('RUN set_simbad_fields() first!')
        try:
            self.simbad_identifiers = self.IDS.decode().split('|')
        except AttributeError:
            self.simbad_identifiers = self.IDS.split('|')

        try:
            gaia_dr2_id = [id for id in self.simbad_identifiers if 'Gaia DR2 ' in id][0].replace("'", "")
            # print(gaia_dr2_id)
            gacs_query = "SELECT * FROM gaiadr2.gaia_source WHERE source_id={}".format(gaia_dr2_id.split(' ')[-1])
            job = Gaia.launch_job_async(gacs_query)
            self.gaiadr2_table = job.get_results()
        except IndexError:
            try:
                gaia_dr1_id = [id for id in self.simbad_identifiers if 'Gaia DR1 ' in id][0]
                # return_gacs_query_as_table(query_string, output_file_seed, overwrite=False, verbose=True):
                gacs_query = "SELECT * FROM gaiadr1.gaia_source WHERE source_id={}".format(gaia_dr1_id.split(' ')[-1])
                job = Gaia.launch_job_async(gacs_query)
                self.gaiadr1_table = job.get_results()
            except IndexError:
                print('No Gaia identifier listed in Simbad!')

