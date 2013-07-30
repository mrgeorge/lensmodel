#! env python

import numpy as np
import profiles

def readData(filename):
    """Read Rachel's data files.
    Returns:
        xshear - comoving Mpc/h
        yshear - delta sigma in h Msun / comoving pc **2
        errshear - error in delta sigma
    """
    xshear,yshear,errshear=np.loadtxt(filename,unpack=True)
    return (xshear,yshear,errshear)

def convertCosmo(xshear,yshear,errshear,redshift=0.1):
    """Convert from Rachel's comoving h=1 coordinates to physical h=0.7 coordinates.
    Returns:
        xshear - physical transverse distance kpc (/h=0.7)
        yshear - delta sigma in (h=0.7) Msun / pc**2
        errshear - error in delta sigma
    """
    h=profiles.cosmo.H0()/100.
    xOut=xshear / (1.+redshift) / h * 1.e3
    yOut=yshear * h * (1.+redshift)**2
    errOut=errshear * h * (1.+redshift)**2

    return (xOut,yOut,errOut)
