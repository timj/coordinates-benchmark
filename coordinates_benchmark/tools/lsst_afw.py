# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Coordinate conversions with the LSST AFW Python package.

https://github.com/lsst/afw
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from astropy.table import Table
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom

SUPPORTED_SYSTEMS = 'fk5 icrs galactic ecliptic'.split()


def transform_celestial(coords, systems):
    out = Table()
    out['lon'] = np.zeros(len(coords), dtype='float64')
    out['lat'] = np.zeros(len(coords), dtype='float64')

    for ii, (lon, lat) in enumerate(zip(coords['lon'], coords['lat'])):

        # Create relevant coordinate object
        if systems['in'] == 'icrs':
            coord = afwCoord.IcrsCoord(afwGeom.Point2D(lon, lat), afwGeom.degrees)
        elif systems['in'] == 'galactic':
            coord = afwCoord.GalacticCoord(afwGeom.Point2D(lon, lat), afwGeom.degrees)
        elif systems['in'] == 'ecliptic':
            coord = afwCoord.EclipticCoord(afwGeom.Point2D(lon, lat), afwGeom.degrees, 2000.0)
        elif systems['in'] == 'fk5':
            coord = afwCoord.Fk5Coord(afwGeom.Point2D(lon, lat), afwGeom.degrees, 2000.0)

        # Now convert to the output system
        if systems['out'] == 'fk5':
            outsys = afwCoord.FK5
        elif systems['out'] == 'icrs':
            outsys = afwCoord.ICRS
        elif systems['out'] == 'galactic':
            outsys = afwCoord.GALACTIC
        elif systems['out'] == 'ecliptic':
            outsys = afwCoord.ECLIPTIC

        out_coord = coord.convert(outsys)

        out[ii]['lon'] = out_coord.getLongitude().asDegrees()
        out[ii]['lat'] = out_coord.getLatitude().asDegrees()

    return out
