""" Spatial interpolation routines
"""
import numpy as np
import netCDF4 as nc

from analysispkg import pkg_geo, pkg_data, pkg_grid

# TODO: split this into get_weights_2d and get_weights_3d
def get_weights(lonm,latm,lono,lato,tmask,zm=None,zo=0):
    """ Horizontal and vertical interpolation weights from model grid to obs location.

    Parameters
    ----------

    lonm,latm : array_like (nj,ni)
        Model grid coords.
    lono,lato : float
        Observation location.
    tmask : array_like (nk,nj,ni), dtype=int
        Model land mask.
    zm : array_like (nk,), optional, used for 3D
        Nominal vertical levels in the model, 'nav_lev'.
    zo : float, optional
        Observation depth.

    Returns
    -------
        Dict with k,j,i indices and (2,2,2) array of weights for interpolation.
    """
    if zm is None:
        zm = [0, 1]
    i1,j1 = pkg_geo.search(lonm,latm,lono,lato)
    i2=i1+1
    j2=j1+1

    k2 = np.digitize(zo,zm)
    if k2 == 0:
        k2 = 1
        k1 = 0
    elif k2 == len(zm) - 1:
        k2 = len(zm) - 1
        k1 = k2 - 1
    else:
        k1 = k2 - 1

    # distance from corners to data point
    d11 = pkg_data.haversine(lonm[j1,i1],latm[j1,i1], lono, lato)
    d12 = pkg_data.haversine(lonm[j1,i2],latm[j1,i2], lono, lato)
    d21 = pkg_data.haversine(lonm[j2,i1],latm[j2,i1], lono, lato)
    d22 = pkg_data.haversine(lonm[j2,i2],latm[j2,i2], lono, lato)
    d = np.array([[d11,d21],[d12,d22]])

    # find longest diagonal length and use it to form an inverse distance weighting length scale
    # this could be swapped for bilinear weights if someone wants to code it up
    da = pkg_data.haversine(lonm[j1,i1],latm[j1,i1], lonm[j2,i2],latm[j2,i2])
    db = pkg_data.haversine(lonm[j2,i1],latm[j2,i1], lonm[j1,i2],latm[j1,i2])
    dd = np.maximum(da,db)/2

    # horizontal weights
    wh = np.exp(-d/dd)
    wh /= wh.sum()

    # vertical linear interp weight
    z1 = 1 - (zo - zm[k1]) / (zm[k2] - zm[k1])
#    print(k1,k2,z1,1-z1)

    # extract corresponding tmask
    tm = tmask[k1:k2+1,j1:j2+1,i1:i2+1]

    # assemble weights
    w= np.zeros([2,2,2])

    # Only use lower level if there are data points
    if tm[1,...].any():
        levels = [0,1]
        w[0,...] = (  z1)*wh
        w[1,...] = (1-z1)*wh
    else:
        w[0,...] = wh
        w[1,...] = 0
        levels = [0]

    # Remove land values and renormalize
    for ii in levels:
        w0 = w[ii,...].sum()                # original sum of weights for this level
        w[ii, tm[ii,...]==0] = 0            # remove weights for land points
        wnew = w[ii,...].sum()
        if wnew != 0: # don't divide by zero for levels with all-zero weights
            w[ii,...] *=  w0 / wnew

    return {'k':[k1,k2],'j':[j1,j2],'i':[i1,i2], "w":w}


def get_weights_nearest(lono, lato, coordfile, bathyfile, max_boxes=6, extrap=True):
    """ Weights for use with "nearest" interpolation

    Parameters
    ----------

    lono,lato : float
        Observation location.
    coordfile, bathyfile : str
        Model coordinates and bathy file.
    max_boxes : int, optional
        Maximum number of boxes to search around the box containing the observation.
    extrap : bool, optional
        Include nodes within 20% of the longest grid dimension in the search. Default is True.

    Returns
    -------
        Dict with k,j,i indices and (1,1,1) array of weights for interpolation.
    """

    # Coordinate file
    glamt, gphit = pkg_grid.tgrid(coordfile)
    glamf, gphif = pkg_grid.fgrid(coordfile)
    glamfe, gphife = pkg_geo.expandf(glamf, gphif)

    # Bathymetry file
    with nc.Dataset(bathyfile, "r") as ncid:
        bathy = ncid.variables["Bathymetry"][:, :]
    lm = np.where(bathy > 0, 1, 0) #TODO Differs from tmask along outer margin

    i, j, dist = pkg_data.find_nearest_point_fast(
        lono, lato, glamt, gphit, lm, glamfe, gphife, max_boxes=max_boxes
    )
    # find_nearest_point_fast searches only within 6 nearest nodes,
    # use slow search if fast search fails
    if extrap and np.isnan(i):
        # search within 20% of the longest grid dimension
        # TODO What is the rationale behind 20%? Is find_nearest_point_fast with max_boxes param enough?
        radius = max(0.2*(glamt.max() - glamt.min()),
                     0.2*(gphit.max() - gphit.min()))
        i, j, dist = pkg_data.find_nearest_point(lono, lato, glamt, gphit, lm, radius)

    k = 0

    w = np.zeros([1, 1, 1]) #TODO: 1 seems to work, check the result, otherwise  w = np.zeros([2, 2, 2])
    if np.isnan(i):
        w[0, 0, 0] = 0
    else:
        w[0, 0, 0] = 1

    return {'k':[k,k],'j':[j,j],'i':[i,i], "w":w}


def interp(v,w):
    """ Apply spatial interpolation weights.

    This interpolates from the nearest 8 data points in 3D (or 4 points in 2D)
    that form a box around the observation location.

    Parameters
    ----------
    v : array_like, (nt,2,2,2) or (nt,2,2)
        Source values for the surrounding box along dimensions (t,z,y,x) or (t,y,x).
    w:  array_like, (2,2,2)
        Interpolation weights.

    Returns
    -------
    r : array_like (nt,)
        Spatially interpolated series.
    """
    if np.ndim(v)==2:
        r=np.einsum('ji,ji->',v,w[0,...]) # single 2D interp (grid angle)
    elif np.ndim(v)==3:
        r=np.einsum('tji,ji->t',v,w[0,...]) # only upper layer is used for 2D
    elif np.ndim(v)==4:
        r=np.einsum('tkji,kji->t',v,w)
    else:
        raise ValueError('Source array should have either 3 or 4 dimensions.')
    return r


def interp_weights_1d(zi, z, zmask=None, extrap_threshold=None):
    """ 1d linear interpolation indices and weights.

    Works on n-d arrays along 1st dimension.

    Parameters
    ----------
        zi : array_like of shape (n,d1,...)
            N-d array or scalar. Coordinate of `n` points to interpolate to.
        z : array_like of shape (m,d1,...)
            Points to interpolate from. Second and subsequent dimensions must
            match those of `zi`.
        zmask : array_like of shape (m,d1,...), optional
            Mask for `z` with 0's for invalid points.
        extrap_threshold : float, optional
            Extrapolate no further that the threshold (same units as `z` and `zi`).
            No threshold by default.

    Returns
    -------
        i1 : array_like of shape (n,d1,...)
        i2 : array_like of shape (n,d1,...)
            Interpolation indices, ravelled to use with (m,d1,...) arrays.
        w1 : array_like of shape (n,d1,...)
        w2 : array_like of shape (n,d1,...)
            Corresponding weights.

    Example
    -------
    To apply indices and weights:
        vi = np.take(v,i1)*w1 + np.take(v,i2)*w2 # where v has z.shape
    """
    # generalize for inputs of various shapes, including scalars
    scalar = np.isscalar(zi)
    if not scalar:
        sz = zi.shape
    zi, z = atleast_2d0(zi, z)

    n = zi.shape[0]
    i1, i2 = [np.zeros(zi.shape, dtype=np.int32) for i in range(2)]  # initialize
    w1, w2 = [np.full(zi.shape, np.nan) for i in range(2)]  # initialize

    # deal with completely masked nodes
    if zmask is not None:
        allmasked = np.all(zmask == 0, axis=0)
        nodes = np.where(~allmasked)[0]  # not masked
    else:
        nodes = range(zi.shape[1])  # all nodes

    for kl in nodes:
        if zmask is not None:
            zm = zmask[:, kl]
            i1[:, kl], i2[:, kl], w1[:, kl], w2[:, kl] = vinterp1d(zi[:, kl], z[zm, kl], extrap_threshold)
        else:
            i1[:, kl], i2[:, kl], w1[:, kl], w2[:, kl] = vinterp1d(zi[:, kl], z[:, kl], extrap_threshold)

    # ravel indices for subsequent indexing of arrays
    dim1 = zi.shape[1:]  # 2nd and subsequent dimensions
    dim1prod = np.prod(dim1)  # number of nodes
    # indices of all columns (ravelled for 2nd and subsequent dimensions) tiled n times
    ic = np.repeat(np.arange(dim1prod, dtype=int).reshape(dim1)[None, :], n, axis=0)
    i1 = np.ravel_multi_index((i1, ic), (z.shape[0], dim1prod))
    i2 = np.ravel_multi_index((i2, ic), (z.shape[0], dim1prod))

    if scalar:
        i1 = np.asscalar(i1)
        i2 = np.asscalar(i2)
    else:
        # keep dimensions of the input zi
        i1 = i1.reshape(sz)
        i2 = i2.reshape(sz)

    return i1, i2, w1, w2


def vinterp1d(gdepr, z, extrap_threshold):
    n = len(z)
    i2 = np.searchsorted(z, gdepr, side='right')
    ileft = i2==0  # dst layers shallower than 1st src layer
    irght = i2==n  # dst layers deeper than last src layer
    i2[ileft] = 1
    i2[irght] = n-1
    i1 = i2 - 1
    w1 = (z[i2] - gdepr) / (z[i2] - z[i1])
    w1[ileft] = 1  # this is nearest neighbour extrapolation to shallower layers
    w1[irght] = 0  # this is nearest neighbour extrapolation to deeper layers
    if extrap_threshold is not None:
        # drop points beyond extrap threshold
        invalid = np.logical_or(gdepr < z[ 0] - extrap_threshold,
                                gdepr > z[-1] + extrap_threshold)
        w1[invalid] = np.nan
    w2 = 1 - w1
    return i1,i2,w1,w2


def atleast_2d0(*arys):
    """ As numpy.atleast_2d, but populates 1st dimension first.
    View inputs as arrays with at least two dimensions.
    Code is from numpy.atleast_2d
    """
    res = []
    for ary in arys:
        ary = np.asanyarray(ary)
        if ary.ndim == 0:
            result = ary.reshape(1, 1)
        elif ary.ndim == 1:
            result = ary[:,np.newaxis] #mk: the only change
        else:
            result = ary
        res.append(result)
    if len(res) == 1:
        return res[0]
    else:
        return res
