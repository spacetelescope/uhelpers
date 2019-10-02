"""Helpers for astrometry tasks

Authors
-------

    Johannes Sahlmann

Use
---

"""

import scipy
import numpy as np

def get_hull(points):
    """Return the convex hull of a collection of points.

    Parameters
    ----------
    points

    Returns
    -------

    """
    hull = scipy.spatial.ConvexHull(points)
    return hull


def get_hull_points(points):
    """Return coordinates of a closed polygon aroind a collection of points.

    Parameters
    ----------
    points

    Returns
    -------

    """
    hull = get_hull(points)
    hull_points = points[hull.vertices,]
    hull_points_closed = np.vstack((hull_points, hull_points[0]))
    return hull_points_closed


def write_hull_to_ds9_region(points, region_file, coordinate_system='icrs'):
    """Get the hull around the region covered by points, write it to ds9 region file.

    See file parameters at http://ds9.si.edu/doc/ref/region.html#RegionProperties

    Parameters
    ----------
    points
    region_file
    coordinate_system

    Returns
    -------

    """
    hull_points = get_hull_points(points)
    with open(region_file, 'w') as f:
        # Region file format: DS9 version 4.1
        f.write('global color=yellow dashlist=8 3 width=2 ')
        f.write('select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('{}\n'.format(coordinate_system))
        f.write('polygon(')
        for j in range(len(hull_points)):
            f.write('{},{},'.format(hull_points[j, 0], hull_points[j, 1]))
        f.write(')\n')
