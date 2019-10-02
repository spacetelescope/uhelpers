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
                          votable_fields=('ra(d)', 'dec(d)', 'pmdec', 'pmra', 'parallax', 'sptype')):

        out_file = os.path.join(out_dir, '{}_simbad_parameters.txt'.format(self.identifier))
        if (not (os.path.isfile(out_file)))| overwrite:
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
