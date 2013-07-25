#! env python

import numpy as np
import matplotlib.pyplot as plt
import esutil
import collections
import scipy.interpolate
import scipy.integrate

cosmo=esutil.cosmology.Cosmo(h=0.7,omega_m=0.3,omega_l=0.7)
# distance convention will be physical kpc/h70
# r for 3d, R for 2d
# 200c for default overdensity


####
# Top Level (generic) Density and Mass Profiles
####
def rho(rkpc, mass, conc, redshift, haloType="nfw", cenType="ps", delta=200., odType="critical", beta=None, rhern=None):
    pass

def sigma(Rkpc):
    pass

def deltaSigma(Rkpc):
    pass

def massEnclosed(rkpc):
    pass

####
# Profile conversions via integration
####
def rhoToSigma(rkpc, rho):
    """Return surface mass density profile by interpolating and integrating rho."""
    logInterp=scipy.interpolate.UnivariateSpline(np.log10(rkpc),np.log10(rho),s=0) # cubic spline interpolation on log axes. Allows extrapolation.
    integrand=lambda theta,r: (r/np.cos(theta)**2) * 10.**logInterp(np.log10(r/np.cos(theta))) # convert linear coords for logInterp
    sigma=2.*np.array([scipy.integrate.quad(integrand,0,np.pi/2.,args=(rr))[0] for rr in rkpc])
    return sigma

def sigmaToDeltaSigma(Rkpc, sigma):
    """Return DS profile by interpolating and integrating sigma.
    Based on e.g. Wright & Brainerd 2000, Eq. 13 line 1
    """
    logInterp=scipy.interpolate.UnivariateSpline(np.log10(Rkpc),np.log10(sigma),s=0) # cubic spline interpolation on log axes. Allows extrapolation.
    integrand=lambda R: R*10.**logInterp(np.log10(R)) # convert linear coords for logInterp
    minRkpc=np.min([1.e-3,0.1*np.min(Rkpc)]) # minimum radius for integration
    sigmaInt=np.array([(2./RR**2) * scipy.integrate.quad(integrand,minRkpc,RR)[0] for RR in Rkpc])
    ds=sigmaInt-sigma
    return ds

####
# Parameter conversions
####
def hubble(redshift):
    """Return Hubble parameter at given redshift in km/s/Mpc."""
    Hz=cosmo.H0()*np.sqrt(cosmo.omega_m()*(1.+redshift)**3 + cosmo.omega_l())
    return Hz

def rhoCritical(redshift):
    """Return rhoCrit = 3H(z)**2 / [8*Pi*G] in Msun/kpc**3."""
    gravConst=4.302e-9 # (km/s)**2 Mpc / Msun
    rhoCrit=3.*hubble(redshift)**2 / (8.*np.pi*gravConst) * 1.e-9 # Msun / kpc**3
    return rhoCrit

def rhoMatter(redshift):
    """Return rhoM = Omega_m(z)*rho_crit(z) in Msun/kpc**3."""
    rhoM=rhoCritical(redshift)*cosmo.omega_m()*(1.+redshift)**3
    return rhoM

def overdensity(redshift, delta=200., type="critical"):
    """Return overdensity of given type in Msun/kpc**3."""
    if(delta=="vir"):
        # Bryan & Norman 1998, as written in footnote 5 of Behroozi 2010
        xx=1./(1.+cosmo.omega_l()/(cosmo.omega_m()*(1.+redshift)**3)) - 1.
        delta=(18.*np.pi**2 + 82.*xx -39*xx**2)/(1.+xx)

    if(type=="critical"):
        od=delta*rhoCritical(redshift)
    elif(type=="background"):
        od=delta*rhoMatter(redshift)
    else:
        raise ValueError(type)
    return od

def haloMass(rhalo, od):
    """Return halo mass = (4/3) pi rhalo**3 * overdensity using units given."""
    mhalo=(4./3)*np.pi*rhalo**3 * od
    return mhalo

def haloRadius(mhalo, od):
    """Return halo radius = [(3/4pi) mhalo/overdensity]**(1/3) using units given."""
    rhalo=((3./(4.*np.pi)) * mhalo / od)**(1./3.)
    return rhalo

def scaleRadius(rhalo, conc):
    """Return scale radius = rhalo/conc in same units as rhalo."""
    return rhalo/conc

####
# Specific Models
####

def rhoNFW(rkpc, mass, conc, od):
    """Return density of NFW profile at rkpc in Msun/kpc**3.
    Follows Wright & Brainerd 2000 Eqs. 1-2
    """
    xx=rkpc/scaleRadius(haloRadius(mass, od), conc)
    deltaC=(od/3.) * (conc**3)/(np.log(1.+conc) - conc/(1.+conc))
    rho=deltaC/(xx * (1.+xx)**2)
    return rho

def sigmaNFW(Rkpc, mass, conc, od):
    """Return surface mass density of NFW profile at Rkpc in Msun/kpc**2.
    Follows Wright & Brainerd 2000 Eq. 11
    Rkpc can be scalar or iterable
    """
    rs=scaleRadius(haloRadius(mass, od), conc)
    xx=Rkpc/rs
    low=(xx<1)
    mid=(xx==1)
    high=(xx>1)

    deltaC=(od/3.) * (conc**3)/(np.log(1.+conc) - conc/(1.+conc))

    if(isinstance(xx,collections.Iterable)):
        sigma=np.zeros_like(xx)
        sigma[low]=2.*rs*deltaC/(xx[low]**2-1.) * (1. - 2./np.sqrt(1.-xx[low]**2) * np.arctanh(np.sqrt((1.-xx[low])/(1.+xx[low]))))
        sigma[mid]=2.*rs*deltaC/3.
        sigma[high]=2.*rs*deltaC/(xx[high]**2-1.) * (1. - 2./np.sqrt(xx[high]**2-1.) * np.arctan(np.sqrt((xx[high]-1.)/(1.+xx[high]))))
    else: # Rkpc is a scalar
        if(low):
            sigma=2.*rs*deltaC/(xx**2-1.) * (1. - 2./np.sqrt(1.-xx**2) * np.arctanh(np.sqrt((1.-xx)/(1.+xx))))
        elif(mid):
            sigma=2.*rs*deltaC/3.
        elif(high):
            sigma=2.*rs*deltaC/(xx**2-1.) * (1. - 2./np.sqrt(xx**2-1.) * np.arctan(np.sqrt((xx-1.)/(1.+xx))))
        else:
            raise ValueError(xx)
    return sigma

def sigmaInteriorNFW(Rkpc, mass, conc, od):
    """Return mean surface mass density of NFW profile inside Rkpc in Msun/kpc**2.
    Follows Wright & Brainerd 2000 Eq. 13
    Rkpc can be scalar or iterable
    """
    rs=scaleRadius(haloRadius(mass, od), conc)
    xx=Rkpc/rs
    low=(xx<1)
    mid=(xx==1)
    high=(xx>1)

    deltaC=(od/3.) * (conc**3)/(np.log(1.+conc) - conc/(1.+conc))

    if(isinstance(xx,collections.Iterable)):
        sigmaInt=np.zeros_like(xx)
        sigmaInt[low]=(4./xx[low]**2)*rs*deltaC * (2./np.sqrt(1.-xx[low]**2) * np.arctanh(np.sqrt((1.-xx[low])/(1.+xx[low]))) + np.log(xx[low]/2.))
        sigmaInt[mid]=4.*rs*deltaC * (1. + np.log(1./2.))
        sigmaInt[high]=(4./xx[high]**2)*rs*deltaC * (2./np.sqrt(xx[high]**2-1.) * np.arctan(np.sqrt((xx[high]-1.)/(1.+xx[high]))) + np.log(xx[high]/2.))
    else: # Rkpc is a scalar
        if(low):
            sigmaInt=(4./xx**2)*rs*deltaC * (2./np.sqrt(1.-xx**2) * np.arctanh(np.sqrt((1.-xx)/(1.+xx))) + np.log(xx/2.))
        elif(mid):
            sigmaInt=4.*rs*deltaC * (1. + np.log(1./2.))
        elif(high):
            sigmaInt=(4./xx**2)*rs*deltaC * (2./np.sqrt(xx**2-1.) * np.arctan(np.sqrt((xx-1.)/(1.+xx))) + np.log(xx/2.))
        else:
            raise ValueError(xx)
    return sigmaInt
    
def deltaSigmaNFW(Rkpc, mass, conc, od):
    """Return DS for GNFW profile in Msun/kpc**2.
    Follows Wright & Brainerd 2000
    """
    return sigmaInteriorNFW(Rkpc, mass, conc, od) - sigmaNFW(Rkpc, mass, conc, od)

def rhoGNFW(rkpc, mass, conc, beta, od):
    """Return density of Generalized NFW profile at rkpc in Msun/kpc**3.
    Follows Eq. 1 of Wyithe, Turner, & Spergel 2001
    Allows generic overdensity (WTS assumes 200c)
    """
    rscale=scaleRadius(haloRadius(mass, od), conc)
    rho=od * rhoGNFW0(conc,beta) / ((rkpc/rscale)**beta * (1.+rkpc/rscale)**(3.-beta))
    return rho

def rhoGNFW0(conc,beta):
    """Return scale density of GNFW profile.
    Follows Eq. 4 of Wyithe, Turner, & Spergel 2001
    NOTE: Overdensity has been absorbed into rhoGNFW (Eq. 1), WTS assumes 200c
    """
    deltaC=conc**3 / (3.*gnfwF(conc,beta))
    return deltaC
    
def gnfwF(conc, beta):
    """Return F(conc) for GNFW profiles following Eq. 5 of Wyithe, Turner, & Spergel 2001."""
    ff,abserr=scipy.integrate.quad(gnfwFIntegrand,0,conc,args=(beta))
    return ff

def gnfwFIntegrand(xx, beta):
    """Compute integrand for Wyithe, Turner, & Spergel 2001 Eq. 5."""
    return xx**(2.-beta) * (1.+xx)**(beta-3.)

def sigmaGNFW(Rkpc, mass, conc, beta, od):
    """Return surface mass density for GNFW profile in Msun/kpc**2.
    Follows Eq. 8 of Wyithe, Turner, & Spergel 2001
    """
    Rscale=scaleRadius(haloRadius(mass, od), conc)
    xx=Rkpc/Rscale
    sigma=2. * od * rhoGNFW0(conc,beta) * Rscale * xx**(1.-beta) * gnfwF2(xx, beta)
    return sigma

def gnfwF2(xx, beta):
    """Compute integral used to define sigmaGNFW in Eq. 8 of Wyithe, Turner, & Spergel 2001."""
    ff,abserr=scipy.integrate.quad(gnfwF2Integrand,0,np.pi/2.,args=(xx,beta))
    return ff

def gnfwF2Integrand(theta, xx, beta):
    """Compute integrand for gnfwF2 used to get sigmaGNFW."""
    return np.sin(theta)*(np.sin(theta) + xx)**(beta-3.)

def deltaSigmaGNFW(Rkpc, mass, conc, beta, od):
    """Return DS for GNFW profile in Msun/kpc**2.
    Simply calls sigmaToDeltaSigma, lacking an analytic profile.
    """
    return sigmaToDeltaSigma(Rkpc, sigmaGNFW(Rkpc, mass, conc, beta, od))

def rhoHernquist(rkpc, mass, rhern):
    """Return density profile for Hernquist model in Msun/kpc**3.
    rhern is the scale radius 'a' used in Eq. 2 of Hernquist 1990
    a=rhalf/(1.+sqrt(2)) where rhalf is the 3d half-mass radius
    """
    rho=mass/(2.*np.pi) * (rhern/rkpc) * 1./(rkpc+a)**3
    return rho

def sigmaHernquist(Rkpc, mass, rhern):
    """Return the surface mass density for Hernquist model in Msun/kpc**2.
    Follows Eq. 32 of Hernquist 1990
    Assumes Upsilon=1 (i.e. consistently uses mass, not light)
    """
    ss=Rkpc/rhern
    sigma=mass/(2.*np.pi*rhern**2*(1.-ss**2)**2) * ((2.+ss**2)*hernX(ss) - 3.)

def hernX(ss):
    """Return X(s) from Hernquist 1990 Eqs. 33 and 34.
    ss can be a scalar or iterable type
    """
    low=((ss >= 0) & (ss <= 1))
    high=(ss>1)
    if(isinstance(ss,collections.Iterable)):
        XX=np.zero_like(ss)
        XX[low]=1./np.sqrt(1.-ss[low]**2) * np.log((1.+np.sqrt(1.-ss[low]**2))/ss[low])
        XX[high]=1./np.sqrt(ss[high]**2-1.) * np.arccos(1./ss[high])
    else: # ss is a scalar
        if(low):
            XX=1./np.sqrt(1.-ss**2) * np.log((1.+np.sqrt(1.-ss**2))/ss)
        elif(high):
            XX=1./np.sqrt(ss**2-1.) * np.arccos(1./ss)
    return XX

def deltaSigmaHernquist(Rkpc, mass, rhern):
    """Return DS for Hernquist profile in Msun/kpc**2.
    Simply calls sigmaToDeltaSigma, lacking an analytic profile.
    """
    return sigmaToDeltaSigma(Rkpc, sigmaHernquist(Rkpc, mass, rhern))

def deltaSigmaPS(Rkpc, mass):
    """Return DS for a point source in Msun/kpc**2."""
    return mass/(np.pi*Rkpc**2)

####
# Wrappers to CONTRA for contracted profiles
####
def rhoAC(rkpc):
    pass

####
# Scaling relations
####
def concentration(mass, redshift, type="Klypin200"):
    """Return concentration based on M-c relation specified by type.
    NOTE: Klypin redshift dependence is currently ignored (i.e. uses z=0)
    """
    if(type=="Klypin200"): # From Klypin 2011 for overdensity=200 (text on pg 7)
        conc=7.2*(mass / (1.e12/(cosmo.H0()/100.)))**-0.075
    elif(type=="KlypinVir"): # From Klypin 2011 for virial overdensity (Eq. 10)
        conc=9.60*(mass / (1.e12/(cosmo.H0()/100.)))**-0.075
    else:
        raise ValueError(type)
    return conc

def stellarMassToHaloMass(logSM, redshift):
    """Return halo mass in Msun given log(StellarMass) (in Msun), z from Behroozi 2010.
    Values drawn from Table 2
    NOTE: Behroozi uses overdensity=delta_vir*background for abundance matching
    """
    logSM00=10.72
    logSM0a=0.55
    logM10=12.35
    logM1a=0.28
    beta0=0.44
    betaa=0.18
    delta0=0.57
    deltaa=0.17
    gamma0=1.56
    gammaa=2.51

    a=1./(1.+redshift)

    logM1=logM10 + logM1a*(a-1.)
    logSM0=logSM00 + logSM0a*(a-1.)
    beta=beta0 + betaa*(a-1.)
    delta=delta0 + deltaa*(a-1.)
    gamma=gamma0 + gammaa*(a-1.)

    # Eq. 21
    logMhalo=logM1 + beta*(logSM - logSM0) + 10.**(delta*(logSM-logSM0)) / (1.+10.**(-gamma*(logSM-logSM0))) - 0.5
    return 10.**logMhalo

