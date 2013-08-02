#! env python

import numpy as np
import os
import sys

try:
    import lensmodel
except ImportError:
    path, filename = os.path.split(__file__)
    sys.path.append(os.path.abspath(os.path.join(path,"../..")))
    import lensmodel

# Input model pars
logMstars=11.
logRstars=0.5
logMhalo=12.0
conc=5.
innerSlopeGNFW=1.
nuDutton=0.
AGnedin=-1.
wGnedin=-1.
cenType="hernquist"
odType="critical"
delta=200.
redshift=0.1

od=lensmodel.profiles.overdensity(redshift=redshift,delta=delta,type=odType)
mstars=10.**logMstars
rstars=10.**logRstars
mhalo=10.**logMhalo

Rmin=40. # kpc - note this is the geometric bin center, so measurements would need to go to radii 0.5*bin width smaller
Rmax=1000. # kpc - Use only bins at R<1Mpc
dlog10R=0.17 # log bin spacing
nBins=np.floor(np.log10(Rmax/Rmin)/dlog10R)+1
Rkpc=Rmin*10.**(np.arange(nBins) * dlog10R)
rkpc=Rkpc

decimal=3 # almost_equal asserstions valid when difference is < 0.5 * 10**(-decimal)

def test_rho_nfw_vs_gnfw():
    rhoNFW=lensmodel.profiles.rhoNFW(rkpc, mhalo, conc, od)
    rhoGNFW=lensmodel.profiles.rhoGNFW(rkpc, mhalo, conc, 1., od)
    np.testing.assert_array_almost_equal(rhoNFW,rhoGNFW,decimal=decimal,err_msg="rhoNFW != rhoGNFW(beta=1)")
def test_sigma_nfw_vs_gnfw():
    sigmaNFW=lensmodel.profiles.sigmaNFW(rkpc, mhalo, conc, od)
    sigmaGNFW=lensmodel.profiles.sigmaGNFW(rkpc, mhalo, conc, 1., od)
    np.testing.assert_array_almost_equal(sigmaNFW,sigmaGNFW,decimal=decimal,err_msg="sigmaNFW != sigmaGNFW(beta=1)")
def test_deltaSigma_nfw_vs_gnfw():
    deltaSigmaNFW=lensmodel.profiles.deltaSigmaNFW(rkpc, mhalo, conc, od)
    deltaSigmaGNFW=lensmodel.profiles.deltaSigmaGNFW(rkpc, mhalo, conc, 1., od)
    np.testing.assert_array_almost_equal(deltaSigmaNFW,deltaSigmaGNFW,decimal=decimal,err_msg="deltaSigmaNFW != deltaSigmaGNFW(beta=1)")

def test_nfw_vs_contra():
    pass

def test_gnfw_vs_contra():
    pass

def test_blum_vs_contra():
    pass
