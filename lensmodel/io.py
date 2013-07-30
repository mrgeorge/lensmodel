#! env python

import numpy as np

def readData(filename):
    """Read Rachel's data files.
    Returns:
       xshear - comoving Mpc/h
       yshear - delta sigma in h Msun / comoving pc **2
       errshear - error in delta sigma
    """
    xshear,yshear,errshear=np.loadtxt(filename,unpack=True)
    return (xshear,yshear,errshear)
