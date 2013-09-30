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

Rmin=1. # kpc - note this is the geometric bin center, so measurements would need to go to radii 0.5*bin width smaller
Rmax=1000. # kpc - Use only bins at R<1Mpc
nBins=50
dlog10R=np.log10(Rmax/Rmin)/(nBins-1)
Rkpc=Rmin*10.**(np.arange(nBins) * dlog10R)
rkpc=Rkpc


#Rmin=40. # kpc - note this is the geometric bin center, so measurements would need to go to radii 0.5*bin width smaller
#Rmax=1000. # kpc - Use only bins at R<1Mpc
#dlog10R=0.15 # log bin spacing
#nBins=np.floor(np.log10(Rmax/Rmin)/dlog10R)+1
#Rkpc=Rmin*10.**(np.arange(nBins) * dlog10R)
#rkpc=Rkpc
nradContra=200

decimal=2 # almost_equal assertions valid when abs(desired-actual) < 0.5 * 10**(-decimal)

def test_rho_nfw_vs_gnfw():
    rhoNFW=lensmodel.profiles.rhoNFW(rkpc, mhalo, conc, od)
    rhoGNFW=lensmodel.profiles.rhoGNFW(rkpc, mhalo, conc, 1., od)
    assert(~np.isnan(rhoNFW).any())
    np.testing.assert_array_almost_equal((rhoGNFW-rhoNFW)/rhoNFW,np.zeros(len(rkpc)),decimal=decimal,err_msg="rhoNFW != rhoGNFW(beta=1)")
def test_sigma_nfw_vs_gnfw():
    sigmaNFW=lensmodel.profiles.sigmaNFW(Rkpc, mhalo, conc, od)
    sigmaGNFW=lensmodel.profiles.sigmaGNFW(Rkpc, mhalo, conc, 1., od)
    assert(~np.isnan(sigmaNFW).any())
    np.testing.assert_array_almost_equal((sigmaGNFW-sigmaNFW)/sigmaNFW,np.zeros(len(Rkpc)),decimal=decimal,err_msg="sigmaNFW != sigmaGNFW(beta=1)")
def test_deltaSigma_nfw_vs_gnfw():
    deltaSigmaNFW=lensmodel.profiles.deltaSigmaNFW(Rkpc, mhalo, conc, od)
    deltaSigmaGNFW=lensmodel.profiles.deltaSigmaGNFW(Rkpc, mhalo, conc, 1., od)
    assert(~np.isnan(deltaSigmaNFW).any())
    np.testing.assert_array_almost_equal((deltaSigmaGNFW-deltaSigmaNFW)/deltaSigmaNFW,np.zeros(len(Rkpc)),decimal=decimal,err_msg="deltaSigmaNFW != deltaSigmaGNFW(beta=1)")

def test_rho_nfw_vs_contra():
    rhoNFW=lensmodel.profiles.rhoNFW(rkpc, mhalo, conc, od)
    rhoAC=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 2, 1., 0., -1., -1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(rhoNFW).any())
    np.testing.assert_array_almost_equal((rhoAC-rhoNFW)/rhoNFW,np.zeros(len(rkpc)),decimal=decimal,err_msg="rhoNFW != rhoAC(beta=1,nu=0)")
def test_sigma_nfw_vs_contra():
    sigmaNFW=lensmodel.profiles.sigmaNFW(Rkpc, mhalo, conc, od)
    sigmaAC=lensmodel.profiles.sigmaAC(Rkpc, mhalo, conc, od, 2, 1., 0., -1., -1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(sigmaNFW).any())
    np.testing.assert_array_almost_equal((sigmaAC-sigmaNFW)/sigmaNFW,np.zeros(len(Rkpc)),decimal=decimal,err_msg="sigmaNFW != sigmaAC(beta=1,nu=0)")
def test_deltaSigma_nfw_vs_contra():
    deltaSigmaNFW=lensmodel.profiles.deltaSigmaNFW(Rkpc, mhalo, conc, od)
    deltaSigmaAC=lensmodel.profiles.deltaSigmaAC(Rkpc, mhalo, conc, od, 2, 1., 0., -1., -1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(deltaSigmaNFW).any())
    np.testing.assert_array_almost_equal((deltaSigmaAC-deltaSigmaNFW)/deltaSigmaNFW,np.zeros(len(Rkpc)),decimal=decimal,err_msg="deltaSigmaNFW != deltaSigmaAC(beta=1,nu=0)")

def test_rho_gnfw_vs_contra():
    rhoGNFW=lensmodel.profiles.rhoGNFW(rkpc, mhalo, conc, 0.8, od)
    rhoAC=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 2, 0.8, 0., -1., -1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(rhoGNFW).any())
    np.testing.assert_array_almost_equal((rhoAC-rhoGNFW)/rhoGNFW,np.zeros(len(rkpc)),decimal=decimal,err_msg="rhoGNFW(beta=0.8) != rhoAC(beta=0.8,nu=0)")
def test_sigma_gnfw_vs_contra():
    sigmaGNFW=lensmodel.profiles.sigmaGNFW(Rkpc, mhalo, conc, 0.8, od)
    sigmaAC=lensmodel.profiles.sigmaAC(Rkpc, mhalo, conc, od, 2, 0.8, 0., -1., -1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(sigmaGNFW).any())
    np.testing.assert_array_almost_equal((sigmaAC-sigmaGNFW)/sigmaGNFW,np.zeros(len(Rkpc)),decimal=decimal,err_msg="sigmaGNFW(beta=0.8) != sigmaAC(beta=0.8,nu=0)")
def test_deltaSigma_gnfw_vs_contra():
    deltaSigmaGNFW=lensmodel.profiles.deltaSigmaGNFW(Rkpc, mhalo, conc, 0.8, od)
    deltaSigmaAC=lensmodel.profiles.deltaSigmaAC(Rkpc, mhalo, conc, od, 2, 0.8, 0., -1., -1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(deltaSigmaGNFW).any())
    np.testing.assert_array_almost_equal((deltaSigmaAC-deltaSigmaGNFW)/deltaSigmaGNFW,np.zeros(len(Rkpc)),decimal=decimal,err_msg="deltaSigmaGNFW(beta=0.8) != deltaSigmaAC(beta=0.8,nu=0)")

def test_rho_blum_vs_gnedin():
    rhoBlum=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 0, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    rhoGnedin=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 1, 1., 1., 1., 1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(rhoBlum).any())
    np.testing.assert_array_almost_equal((rhoGnedin-rhoBlum)/rhoBlum,np.zeros(len(rkpc)),decimal=decimal,err_msg="rhoBlum != rhoGnedin(A=1.,w=1.)")
def test_sigma_blum_vs_gnedin():
    sigmaBlum=lensmodel.profiles.sigmaAC(rkpc, mhalo, conc, od, 0, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    sigmaGnedin=lensmodel.profiles.sigmaAC(rkpc, mhalo, conc, od, 1, 1., 1., 1., 1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(sigmaBlum).any())
    np.testing.assert_array_almost_equal((sigmaGnedin-sigmaBlum)/sigmaBlum,np.zeros(len(rkpc)),decimal=decimal,err_msg="sigmaBlum != sigmaGnedin(A=1.,w=1.)")
def test_deltaSigma_blum_vs_gnedin():
    deltaSigmaBlum=lensmodel.profiles.deltaSigmaAC(rkpc, mhalo, conc, od, 0, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    deltaSigmaGnedin=lensmodel.profiles.deltaSigmaAC(rkpc, mhalo, conc, od, 1, 1., 1., 1., 1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(deltaSigmaBlum).any())
    np.testing.assert_array_almost_equal((deltaSigmaGnedin-deltaSigmaBlum)/deltaSigmaBlum,np.zeros(len(rkpc)),decimal=decimal,err_msg="deltaSigmaBlum != deltaSigmaGnedin(A=1.,w=1.)")

def test_rho_blum_vs_dutton():
    rhoBlum=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 0, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    rhoDutton=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 2, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(rhoBlum).any())
    np.testing.assert_array_almost_equal((rhoDutton-rhoBlum)/rhoBlum,np.zeros(len(rkpc)),decimal=decimal,err_msg="rhoBlum != rhoDutton(nu=1)")
def test_sigma_blum_vs_dutton():
    sigmaBlum=lensmodel.profiles.sigmaAC(rkpc, mhalo, conc, od, 0, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    sigmaDutton=lensmodel.profiles.sigmaAC(rkpc, mhalo, conc, od, 2, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(sigmaBlum).any())
    np.testing.assert_array_almost_equal((sigmaDutton-sigmaBlum)/sigmaBlum,np.zeros(len(rkpc)),decimal=decimal,err_msg="sigmaBlum != sigmaDutton(nu=1)")
def test_deltaSigma_blum_vs_dutton():
    deltaSigmaBlum=lensmodel.profiles.deltaSigmaAC(rkpc, mhalo, conc, od, 0, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    deltaSigmaDutton=lensmodel.profiles.deltaSigmaAC(rkpc, mhalo, conc, od, 2, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    assert(~np.isnan(deltaSigmaBlum).any())
    np.testing.assert_array_almost_equal((deltaSigmaDutton-deltaSigmaBlum)/deltaSigmaBlum,np.zeros(len(rkpc)),decimal=decimal,err_msg="deltaSigmaBlum != deltaSigmaDutton(nu=1)")


def test_mnfw():
    rhoNFW=lensmodel.profiles.rhoNFW(rkpc, mhalo, conc, od)
    rhalo=lensmodel.profiles.haloRadius(mhalo, od)
    massNFW=lensmodel.profiles.rhoToMenclosed(np.array([rhalo]), rkpc, rhoNFW)[0]
    np.testing.assert_almost_equal((massNFW-mhalo)/mhalo, 0., decimal=decimal,err_msg="massNFW(<rhalo) != mhalo")
    
def test_mblum():
    rhoBlum=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 0, 1., 1., -1., -1., mstars, rstars, nrad=nradContra)
    rhalo=lensmodel.profiles.haloRadius(mhalo, od)
    massBlum=lensmodel.profiles.rhoToMenclosed(np.array([rhalo]), rkpc, rhoBlum)[0]
    np.testing.assert_almost_equal((massBlum-mhalo)/mhalo, 0., decimal=decimal,err_msg="massBlum(<rhalo) != mhalo")
    
    #def test_rho_gnedin_vs_dutton():
    #    rhoGnedin=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 1, 1., 1., 0.85, 0.8, mstars, rstars)
    #    rhoDutton=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 2, 1., 0.8, -1., -1., mstars, rstars)
    #    np.testing.assert_array_almost_equal((rhoDutton-rhoGnedin)/rhoGnedin,np.zeros(len(rkpc)),decimal=decimal,err_msg="rhoGnedin(A=0.85,w=0.8) != rhoDutton(nu=0.8)")

