#! env python

import numpy as np
import profiles
import fitsio

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

def parsToRec(pars,labels=np.array(["logMstar", "logRstar", "logMhalo", "conc", "innerSlopeGNFW", "nuDutton", "AGnedin", "wGnedin"])):
    dtype=[(label,float) for label in labels]
    rec=np.recarray(len(pars),dtype=dtype)
    for ii in range(len(labels)):
        rec[labels[ii]]=pars[:,ii]
    return rec

def chainToRec(chain,lnprob,labels=np.array(["logMstar", "logRstar", "logMhalo", "conc", "innerSlopeGNFW", "nuDutton", "AGnedin", "wGnedin"])):
    nGal=chain.shape[0]
    nPars=chain.shape[1]
    arr=np.zeros((nGal,nPars+1))
    arr[:,:-1]=chain
    arr[:,-1]=lnprob
    labels=np.append(labels,"lnprob")
    rec=parsToRec(arr,labels=labels)
    return rec

def recToPars(rec,labels=np.array(["logMstar", "logRstar", "logMhalo", "conc", "innerSlopeGNFW", "nuDutton", "AGnedin", "wGnedin"])):
    recLabels=rec.dtype.fields.keys() # note, this list is unordered since rec is a dict, so we need to use labels (which should be sorted to match the order of columns in pars array)
    pars=np.zeros((len(rec),len(labels)))
    for ii in range(len(labels)):
        pars[:,ii]=rec[labels[ii]]
    return pars

def writeRec(rec,filename,clobber=True,compress="GZIP"):
    fitsio.write(filename,rec,clobber=clobber,compress=compress)

def readRec(filename):
    rec=fitsio.read(filename)
    return rec
