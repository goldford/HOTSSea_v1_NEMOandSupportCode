import netCDF4 as nc
import numpy as np


def tgrid(file_coord):
    with nc.Dataset(file_coord, "r") as ncid:
        glamt = np.squeeze(ncid.variables["glamt"][..., :, :])
        gphit = np.squeeze(ncid.variables["gphit"][..., :, :])
    return glamt, gphit

def ugrid(file_coord):
    with nc.Dataset(file_coord, "r") as ncid:
        glamu = np.squeeze(ncid.variables["glamu"][..., :, :])
        gphiu = np.squeeze(ncid.variables["gphiu"][..., :, :])
    return glamu, gphiu

def vgrid(file_coord):
    with nc.Dataset(file_coord, "r") as ncid:
        glamv = np.squeeze(ncid.variables["glamv"][..., :, :])
        gphiv = np.squeeze(ncid.variables["gphiv"][..., :, :])
    return glamv, gphiv

def fgrid(file_coord):
    with nc.Dataset(file_coord, "r") as ncid:
        glamf = np.squeeze(ncid.variables["glamf"][..., :, :])
        gphif = np.squeeze(ncid.variables["gphif"][..., :, :])
    return glamf, gphif

def tmask(file_mesh):
    with nc.Dataset(file_mesh) as ncid:
        t_mask = ncid['tmask'][0,...]
    return t_mask

def umask(file_mesh):
    with nc.Dataset(file_mesh) as ncid:
        u_mask = ncid['umask'][0,...]
    return u_mask

def vmask(file_mesh):
    with nc.Dataset(file_mesh) as ncid:
        v_mask = ncid['vmask'][0,...]
    return v_mask

def nav_lev(file_mesh):
    with nc.Dataset(file_mesh) as ncid:
        zm = ncid['nav_lev'][:].filled()
    return zm

def bathymetry(file_bathy, maskfile=None):
    """
    Loads bathymetry from bathy_file. If provided, applies tmaskutil from maskfile and returns bathy as a masked array.
    """
    with nc.Dataset(file_bathy) as ncid:
        bathy = ncid.variables['Bathymetry'][...]

    if maskfile is not None:
        with nc.Dataset(maskfile) as ncm:
            mask = ncm['tmaskutil'][0, ...].filled()
        bathy = np.ma.masked_where(mask == 0, bathy)
    return bathy

def landmask(file_bathy):
    """
    Loads bathymetry from file_bathy, converts to a landmask where land=1 and ocean=NaN.
    """
    bathy = bathymetry(file_bathy)
    land = np.where(bathy > 0, np.nan, 1.0)
    return land

def tgrid_and_landmask(file_coord, file_bathy):
    lon, lat = tgrid(file_coord)
    land = landmask(file_bathy)
    return lon, lat, land

def grid_angle_t(file_coord):
    """ Compute angles between model grid lines and the North direction

    Sines and cosines of the angle between the north-south direction and the
    j-direction for T-grid.

    Computation done on the north stereographic polar plane.
    """

    glamt, gphit = tgrid(file_coord)
    glamv, gphiv = vgrid(file_coord)
    gsint, gcost = grid_angle(glamt, gphit, glamv, gphiv)
    return gsint, gcost


def grid_angle(glamt, gphit, glamv, gphiv):
    """ Compute angles between model grid lines and the North direction """
    # Code is adapted from nemo_tools.py by Christoph Renkl,
    # which is adapted from fortran code by
    #  !!   7.0  !  96-07  (O. Marti )  Original code
    #  !!   8.0  !  98-06  (G. Madec )
    #  !!   8.5  !  98-06  (G. Madec )  Free form, F90 + opt.
    #  !!   9.2  !  07-04  (S. Masson)  Add T, F points and bugfix in cos lateral boundary

    # initialize output arrays, same shape as grid dimensions
    gsint = np.zeros_like(glamt)
    gcost = np.zeros_like(glamt)

    # degree to radians conversion
    rad = np.pi/180.  # == np.deg2rad(1)

    # north pole direction & modulous (at t-point)
    zlam = rad*glamt[1:-1, 1:]
    zphi = rad*gphit[1:-1, 1:]
    tan_phi = np.tan( np.pi/4. - zphi/2. )
    zxnpt = 0. - 2. * np.cos(zlam) * tan_phi
    zynpt = 0. - 2. * np.sin(zlam) * tan_phi
    znnpt = zxnpt*zxnpt + zynpt*zynpt

    # j-direction: v-point segment direction (around t-point)
    zlam = rad*glamv[1:-1, 1:]
    zphi = rad*gphiv[1:-1, 1:]
    zlan = rad*glamv[:-2, 1:]
    zphh = rad*gphiv[:-2, 1:]
    tan_phi = np.tan( np.pi/4. - zphi/2. )
    tan_phh = np.tan( np.pi/4. - zphh/2. )
    zxvvt =  2. * np.cos(zlam) * tan_phi  -  2. * np.cos(zlan) * tan_phh
    zyvvt =  2. * np.sin(zlam) * tan_phi  -  2. * np.sin(zlan) * tan_phh

    znvvt = np.sqrt( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt ) )
    znvvt = np.maximum( znvvt, 1.e-14 )

    # cosine and sine using scalar and vectorial products
    gsint[1:-1, 1:] = ( zxnpt*zyvvt - zynpt*zxvvt ) / znvvt
    gcost[1:-1, 1:] = ( zxnpt*zxvvt + zynpt*zyvvt ) / znvvt

    # Lateral boundary conditions
    # copy second to first row
    gsint[0,:] = gsint[1,:]
    gcost[0,:] = gcost[1,:]
    # copy second last to last row
    gsint[-1,:] = gsint[-2,:]
    gcost[-1,:] = gcost[-2,:]
    # copy second to first column
    gsint[:,0] = gsint[:,1]
    gcost[:,0] = gcost[:,1]
    # copy second last to last column
    gsint[:,-1] = gsint[:,-2]
    gcost[:,-1] = gcost[:,-2]

    # in case mesh mask with land processor elimiation is supplied
    lpe_mask = ((glamt == 0) & (gphit == 0)) | ((glamv == 0) & (gphiv == 0))
    # expand mask one row down (see zlan,zphh calculations)
    lpe_mask[1:, :] = lpe_mask[1:, :] | lpe_mask[:-1, :]
    gsint[lpe_mask] = 0
    gcost[lpe_mask] = 1

    return gsint,gcost

def rotate(u, v, sina, cosa, direction):
    """
    Rotates u and v by angle a (provided as sina, cosa)
    If direction is "ij->en" we rotate from grid-aligned to east-north
    If direction is "en->ij" we rotate from east-north to grid-aligned

    Returns
    -------
    ur, vr  -- rotated velocities
    """
    if direction == 'ij->en':
        ur = u * cosa - v * sina
        vr = u * sina + v * cosa
    elif direction == 'en->ij':
        ur =   u * cosa + v * sina
        vr = - u * sina + v * cosa
    else:
        raise ValueError('Unknown rotation direction {}'.format(direction))
    return ur, vr