""" Routines to work with model grids.
"""
import os
import netCDF4 as nc
import numpy as np
from matplotlib import path

from analysispkg import pkg_grid

def expandf(glamf, gphif):
    # Expand the f points grid so the f points form complete boxes around the t points
    # This is needed because the coordinates file truncates the f points by one.
    NY, NX = glamf.shape[0], glamf.shape[1]
    glamfe = np.zeros([NY+1, NX+1])
    gphife = np.zeros([NY+1, NX+1])
    # Long
    glamfe[1:,1:] = glamf
    glamfe[0,1:] = glamf[0,:] - (glamf[1,:] - glamf[0,:])     # extrapolation
    glamfe[:,0] = glamfe[:,1] - (glamfe[:,2] - glamfe[:,1])   # extrapolation
    # Lat
    gphife[1:,1:] = gphif
    gphife[0,1:] = gphif[0,:] - (gphif[1,:] - gphif[0,:])     # extrapolation
    gphife[:,0] = gphife[:,1] - (gphife[:,2] - gphife[:,1])   # extrapolation
    return glamfe, gphife

# Build a polygon following the encompassing the domain
def makebox(glamfe,gphife,il,ir,jl,jr):
    # Builds a bounding polygon enclosing grid boxes il..ir by jl..jr
    n, L = 0, 2*(ir-il+1) + 2*(jr-jl-1) + 1
    px, py = np.zeros([L,1]), np.zeros([L,1])
    for i in range(il,ir+1):
        px[n] = glamfe[jl,i]
        py[n] = gphife[jl,i]
        n+=1
    for j in range(jl+1,jr+1):
        px[n] = glamfe[j,ir]
        py[n] = gphife[j,ir]
        n+=1
    for i in range(ir-1,il-1,-1):
        px[n] = glamfe[jr,i]
        py[n] = gphife[jr,i]
        n+=1
    for j in range(jr-1,jl-1,-1):
        px[n] = glamfe[j,il]
        py[n] = gphife[j,il]
        n+=1
    if n < L:
        print("n,L:",n,L)
    return np.hstack((px,py))

def domain_polygon(file_coord):
    """
    Returns a matplotlib Path object that traces the outer boundary of the
    model domain. The f-points are used here, and they are extrapolated by one
    via expandf(). The purpose of this polygon is for testing if points are
    inside the domain, eg:

        bbox = pkg_geo.domain_polygon('mesh_mask.nc')
        if bbox.contains_point((lon,lat)):
            # do something with it

    The input can be a coordinates file or a mesh mask file.

    Note: this is slightly approximate because we're not using a great circle path
    between each pair of f points. If we added a number of great circle waypoints on
    each segment it would be marginally improved. Unlikely to matter.
    """
    glamf, gphif = pkg_grid.fgrid(file_coord)
    glamfe, gphife = expandf(glamf, gphif)
    NY, NX = glamf.shape[0], glamf.shape[1]
    p = makebox(glamfe,gphife,0,NX,0,NY)
    poly = path.Path(p, closed=True)
    return poly

def search(glamfe,gphife,cx,cy):
    """ Conducts a binary search to find which grid box contains point cx,cy.
    """
    NY,NX=glamfe.shape[0], glamfe.shape[1]
    il, ir = 0, NX-1
    jl, jr = 0, NY-1
    # Conducts a binary search to find which grid box contains point cx,cy
    while (ir-il > 1) or (jr-jl > 1):
        # Do a refinement in i
        if ir-il>1:
            irr = il + (ir-il)//2
            p = makebox(glamfe,gphife,il,irr,jl,jr)
            poly = path.Path(p, closed=True)
            test = poly.contains_points([(cx,cy)])
            if test: ir = irr
            else: il = irr
        # Do a refinement in j
        if jr-jl > 1:
            jrr = jl + (jr-jl)//2
            p = makebox(glamfe,gphife,il,ir,jl,jrr)
            poly = path.Path(p, closed=True)
            test = poly.contains_points([(cx,cy)])
            if test: jr = jrr
            else: jl = jrr
    return il,jl


def magnetic_declination(lon, lat, time):
    """ magnetic declination

    Uses the NOAA National Geophysical Data Center, epoch 2015 data.
    Software used:
    https://github.com/cmweiss/geomag.git

    Parameters
    ----------
    lon,lat : int
    time : datetime.datetime object

    Returns
    -------
    int, magnetic declination, positive is East
    """
    import geomag
    gm = geomag.geomag.GeoMag(os.path.join(geomag.__path__[0], "WMM.COF"))
    mag = gm.GeoMag(lat, lon, time=time.date())
    return mag.dec
