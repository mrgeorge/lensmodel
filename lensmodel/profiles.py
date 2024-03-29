#! env python

import numpy as np
import collections
import scipy.interpolate
import scipy.integrate
import contra
import esutil

cosmo=esutil.cosmology.Cosmo(h=0.7,omega_m=0.3,omega_l=0.7)
# distance convention will be physical kpc/h70
# r for 3d, R for 2d
# 200c for default overdensity

extFactor=100. # used to define extfactor*rvir; contracted profiles are extended following NFW here
truncFactor=100 # used to define truncFactor*rvir; rhoToSigma integral is truncated here
Rkpcfine=np.logspace(0,3.5,num=100) # used for more extended / finer sampling when interpolating
rkpcfine=Rkpcfine

####
# Top Level (generic) Density and Mass Profiles
####
def rho(pars, rkpc, redshift=0.1, cenType="hernquist", delta=200., odType="critical", rhoStarsOut=None, rhoHaloOut=None):

    # Unpack model parameters
    logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin = pars
    mstars=10.**logMstars
    rstars=10.**logRstars
    mhalo=10.**logMhalo
    od=overdensity(redshift, delta=delta, type=odType)

    # Stellar term
    if(cenType=="ps"):
        rhoStars=rkpc * 0.
    elif(cenType=="hernquist"):
        rhoStars=rhoHernquist(rkpc, mstars, rstars)        
    else:
        raise ValueError(cenType)

    # Halo term
    if((innerSlopeGNFW == 1) & (nuDutton == 0) & (AGnedin <= 0) & (wGnedin <= 0)): # vanilla NFW
        rhoHalo=rhoNFW(rkpc, mhalo, conc, od)
    elif((innerSlopeGNFW != 1) & (nuDutton == 0) & (AGnedin <= 0) & (wGnedin <= 0)): # GNFW
        rhoHalo=rhoGNFW(rkpc, mhalo, conc, innerSlopeGNFW, od)
    elif((nuDutton == 0) & (AGnedin == 1) & (wGnedin == 1)): # Blumenthal AC
        MAC=0
        rhoHalo=rhoAC(rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars)
    elif((nuDutton == 0) & (AGnedin > 0) & (wGnedin > 0)): # Gnedin modified AC
        MAC=1
        rhoHalo=rhoAC(rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars)
    elif((nuDutton != 0.) & (AGnedin <= 0) & (wGnedin <= 0)): # Dutton modified AC
        MAC=2
        rhoHalo=rhoAC(rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars)
    else:
        raise ValueError(innerSlopeGNFW, nuDutton, AGnedin, wGnedin)

    rhoTot=rhoStars+rhoHalo

    # Save profiles to outputs if desired (probably won't work for non-iterable Rkpc)
    if(rhoStarsOut is not None):
       rhoStarsOut[:]=rhoStars
    if(rhoHaloOut is not None):
       rhoHaloOut[:]=rhoHalo

    return rhoTot

def sigma(pars, Rkpc, redshift=0.1, cenType="hernquist", delta=200., odType="critical", sigmaStarsOut=None, sigmaHaloOut=None):

    # Unpack model parameters
    logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin = pars
    mstars=10.**logMstars
    rstars=10.**logRstars
    mhalo=10.**logMhalo
    od=overdensity(redshift, delta=delta, type=odType)

    # Stellar term
    if(cenType=="ps"):
        sigmaStars=Rkpc * 0.
    elif(cenType=="hernquist"):
        sigmaStars=sigmaHernquist(Rkpc, mstars, rstars)
    else:
        raise ValueError(cenType)

    # Halo term
    if((innerSlopeGNFW == 1) & (nuDutton == 0) & (AGnedin <= 0) & (wGnedin <= 0)): # vanilla NFW
        sigmaHalo=sigmaNFW(Rkpc, mhalo, conc, od)
    elif((innerSlopeGNFW != 1) & (nuDutton == 0) & (AGnedin <= 0) & (wGnedin <= 0)): # GNFW
        sigmaHalo=sigmaGNFW(Rkpc, mhalo, conc, innerSlopeGNFW, od)
    elif((nuDutton == 0) & (AGnedin == 1) & (wGnedin == 1)): # Blumenthal AC
        MAC=0
        sigmaHalo=sigmaAC(Rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars)
    elif((nuDutton == 0) & (AGnedin > 0) & (wGnedin > 0)): # Gnedin modified AC
        MAC=1
        sigmaHalo=sigmaAC(Rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars)
    elif((nuDutton != 0.) & (AGnedin <= 0) & (wGnedin <= 0)): # Dutton modified AC
        MAC=2
        sigmaHalo=sigmaAC(Rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars)
    else:
        raise ValueError(innerSlopeGNFW, nuDutton, AGnedin, wGnedin)

    sigmaTot=sigmaStars+sigmaHalo

    # Save profiles to outputs if desired (probably won't work for non-iterable Rkpc)
    if(sigmaStarsOut is not None):
       sigmaStarsOut[:]=sigmaStars
    if(sigmaHaloOut is not None):
       sigmaHaloOut[:]=sigmaHalo
    
    return sigmaTot


def deltaSigma(pars, Rkpc, redshift=0.1, cenType="hernquist", delta=200., odType="critical", deltaSigmaStarsOut=None, deltaSigmaHaloOut=None):

    # Unpack model parameters
    logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin = pars
    mstars=10.**logMstars
    rstars=10.**logRstars
    mhalo=10.**logMhalo
    od=overdensity(redshift, delta=delta, type=odType)

    # Stellar term
    if(cenType=="ps"):
        deltaSigmaStars=Rkpc * 0.
    elif(cenType=="hernquist"):
        deltaSigmaStars=deltaSigmaHernquist(Rkpc, mstars, rstars)
    else:
        raise ValueError(cenType)

    # Halo term
    if((innerSlopeGNFW == 1) & (nuDutton == 0) & (AGnedin <= 0) & (wGnedin <= 0)): # vanilla NFW
        deltaSigmaHalo=deltaSigmaNFW(Rkpc, mhalo, conc, od)
    elif((innerSlopeGNFW != 1) & (nuDutton == 0) & (AGnedin <= 0) & (wGnedin <= 0)): # GNFW
        deltaSigmaHalo=deltaSigmaGNFW(Rkpc, mhalo, conc, innerSlopeGNFW, od)
    elif((nuDutton == 0) & (AGnedin == 1) & (wGnedin == 1)): # Blumenthal AC
        MAC=0
        deltaSigmaHalo=deltaSigmaAC(Rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars)
    elif((nuDutton == 0) & (AGnedin > 0) & (wGnedin > 0)): # Gnedin modified AC
        MAC=1
        deltaSigmaHalo=deltaSigmaAC(Rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars)
    elif((nuDutton != 0.) & (AGnedin <= 0) & (wGnedin <= 0)): # Dutton modified AC
        MAC=2
        deltaSigmaHalo=deltaSigmaAC(Rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars)
    else:
        raise ValueError(innerSlopeGNFW, nuDutton, AGnedin, wGnedin)

    deltaSigmaTot=deltaSigmaStars+deltaSigmaHalo

    # Save profiles to outputs if desired (probably won't work for non-iterable Rkpc)
    if(deltaSigmaStarsOut is not None):
       deltaSigmaStarsOut[:]=deltaSigmaStars
    if(deltaSigmaHaloOut is not None):
       deltaSigmaHaloOut[:]=deltaSigmaHalo

    return deltaSigmaTot


####
# Profile conversions via integration
####
def rhoToSigma(Rkpc_out, rkpc_in, rho, rtrunc=np.inf):
    """Return surface mass density profile by interpolating and integrating rho.
    Inputs:
        Rkpc_out - radius in kpc for DS to be sampled on
        rkpc_in - radius in kpc at which input rho is sampled
        rho - mass density in Msun/pc**3, same length as rkpc_in
        rtrunc - radius to integrate density profile along LOS (default=inf)
    Returns:
        sigma - surface mass density in Msun/pc**2, same length as Rkpc_out
    """
    logInterp=scipy.interpolate.UnivariateSpline(np.log10(rkpc_in),np.log10(rho),s=0) # cubic spline interpolation on log axes. Allows extrapolation.
    sigma=2.*np.array([scipy.integrate.quad(rhoTruncIntegrand,0,np.pi/2.,args=(RR,rtrunc,logInterp))[0] for RR in Rkpc_out]) * 1.e3
    return sigma

def rhoTruncIntegrand(theta, Rkpc, rtrunc, logInterp):
    rlos=Rkpc/np.cos(theta)
    low=(rlos <= rtrunc)
    high=(rlos > rtrunc)
    if(isinstance(rlos,collections.Iterable)):
        integrand=np.zeros(len(rlos))
        integrand[low]=(rlos[low]/np.cos(theta[low])) * 10.**logInterp(np.log10(rlos[low]))
        integrand[high]=0.
    else:
        if(low):
            integrand=(rlos/np.cos(theta)) * 10.**logInterp(np.log10(rlos))
        elif(high):
            integrand=0.
        else:
            raise ValueError(rlos)

    return integrand

def rhoToMenclosed(rkpc_out, rkpc_in, rho):
    logInterp=scipy.interpolate.UnivariateSpline(np.log10(rkpc_in),np.log10(rho),s=0) # cubic spline interpolation on log axes. Allows extrapolation.
    integrand=lambda r: r**2 * 10.**logInterp(np.log10(r)) * 1.e9
    massEnclosed=4.*np.pi*np.array([scipy.integrate.quad(integrand,0,rr)[0] for rr in rkpc_out])
    return massEnclosed
    
def sigmaToDeltaSigma(Rkpc_out, Rkpc_in, sigma):
    """Return DS profile by interpolating and integrating sigma.
    Based on e.g. Wright & Brainerd 2000, Eq. 13 line 1
    Inputs:
        Rkpc_out - radius in kpc for DS to be sampled on
        Rkpc_in - radius in kpc at which input sigma is sampled
        sigma - surface mass density in Msun/pc**2, same length as Rkpc_in
    Returns:
        ds - delta sigma in Msun/pc**2, same length as Rkpc_out
    """
    logInterp=scipy.interpolate.UnivariateSpline(np.log10(Rkpc_in),np.log10(sigma),s=0) # cubic spline interpolation on log axes. Allows extrapolation.
    integrand=lambda R: R*10.**logInterp(np.log10(R)) # convert linear coords for logInterp
    minRkpc=np.min([1.e-3,0.1*np.min(Rkpc_out)]) # minimum radius for integration
    sigmaInt=np.array([(2./RR**2) * scipy.integrate.quad(integrand,minRkpc,RR)[0] for RR in Rkpc_out])
    ds=sigmaInt-10.**logInterp(np.log10(Rkpc_out))
    return ds

####
# Parameter conversions
####
def hubble(redshift):
    """Return Hubble parameter at given redshift in km/s/Mpc."""
    Hz=cosmo.H0()*np.sqrt(cosmo.omega_m()*(1.+redshift)**3 + cosmo.omega_l())
    return Hz

def rhoCritical(redshift):
    """Return rhoCrit = 3H(z)**2 / [8*Pi*G] in Msun/pc**3."""
    gravConst=4.302e-9 # (km/s)**2 Mpc / Msun
    rhoCrit=3.*hubble(redshift)**2 / (8.*np.pi*gravConst) * 1.e-18 # Msun / pc**3
    return rhoCrit

def rhoMatter(redshift):
    """Return rhoM = Omega_m(z)*rho_crit(z) in Msun/pc**3."""
    rhoM=rhoCritical(redshift)*cosmo.omega_m()*(1.+redshift)**3
    return rhoM

def deltaVir(redshift):
    # Bryan & Norman 1998, as written in footnote 5 of Behroozi 2010
    xx=1./(1.+cosmo.omega_l()/(cosmo.omega_m()*(1.+redshift)**3)) - 1.
    delta=(18.*np.pi**2 + 82.*xx -39*xx**2)/(1.+xx)
    return delta
    

def overdensity(redshift, delta=200., type="critical"):
    """Return overdensity of given type in Msun/pc**3."""
    if(delta=="vir"):
        delta=deltaVir(redshift)

    if(type=="critical"):
        od=delta*rhoCritical(redshift)
    elif(type=="background"):
        od=delta*rhoMatter(redshift)
    else:
        raise ValueError(type)
    return od

def haloMass(rhalo, od):
    """Return halo mass = (4/3) pi rhalo**3 * overdensity in Msun [rhalo in kpc, od in Msun/pc**3]."""
    mhalo=(4./3)*np.pi*rhalo**3 * (od*1.e9)
    return mhalo

def haloRadius(mhalo, od):
    """Return halo radius = [(3/4pi) mhalo/overdensity]**(1/3) in kpc [mhalo in Msun, od in Msun/pc**3]."""
    rhalo=((3./(4.*np.pi)) * mhalo / (od*1.e9))**(1./3.)
    return rhalo

def scaleRadius(rhalo, conc):
    """Return scale radius = rhalo/conc in same units as rhalo."""
    return rhalo/conc

####
# Specific Models
####

def rhoNFW(rkpc, mass, conc, od):
    """Return density of NFW profile at rkpc in Msun/pc**3.
    Follows Wright & Brainerd 2000 Eqs. 1-2
    """
    xx=rkpc/scaleRadius(haloRadius(mass, od), conc)
    deltaC=(od/3.) * (conc**3)/(np.log(1.+conc) - conc/(1.+conc))
    rho=deltaC/(xx * (1.+xx)**2)
    return rho

def sigmaNFW(Rkpc, mass, conc, od):
    """Return surface mass density of NFW profile at Rkpc in Msun/pc**2.
    Follows Wright & Brainerd 2000 Eq. 11
    Rkpc can be scalar or iterable
    """
    rs=scaleRadius(haloRadius(mass, od), conc) # kpc
    xx=Rkpc/rs
    low=(xx<1)
    mid=(xx==1)
    high=(xx>1)

    deltaC=(od/3.) * (conc**3)/(np.log(1.+conc) - conc/(1.+conc)) # Msun/pc**3

    if(isinstance(xx,collections.Iterable)):
        sigma=np.zeros_like(xx)
        sigma[low]=2.*(rs*1.e3)*deltaC/(xx[low]**2-1.) * (1. - 2./np.sqrt(1.-xx[low]**2) * np.arctanh(np.sqrt((1.-xx[low])/(1.+xx[low]))))
        sigma[mid]=2.*(rs*1.e3)*deltaC/3.
        sigma[high]=2.*(rs*1.e3)*deltaC/(xx[high]**2-1.) * (1. - 2./np.sqrt(xx[high]**2-1.) * np.arctan(np.sqrt((xx[high]-1.)/(1.+xx[high]))))
    else: # Rkpc is a scalar
        if(low):
            sigma=2.*(rs*1.e3)*deltaC/(xx**2-1.) * (1. - 2./np.sqrt(1.-xx**2) * np.arctanh(np.sqrt((1.-xx)/(1.+xx))))
        elif(mid):
            sigma=2.*(rs*1.e3)*deltaC/3.
        elif(high):
            sigma=2.*(rs*1.e3)*deltaC/(xx**2-1.) * (1. - 2./np.sqrt(xx**2-1.) * np.arctan(np.sqrt((xx-1.)/(1.+xx))))
        else:
            raise ValueError(xx)
    return sigma

def sigmaInteriorNFW(Rkpc, mass, conc, od):
    """Return mean surface mass density of NFW profile inside Rkpc in Msun/pc**2.
    Follows Wright & Brainerd 2000 Eq. 13
    Rkpc can be scalar or iterable
    """
    rs=scaleRadius(haloRadius(mass, od), conc) # kpc
    xx=Rkpc/rs
    low=(xx<1)
    mid=(xx==1)
    high=(xx>1)

    deltaC=(od/3.) * (conc**3)/(np.log(1.+conc) - conc/(1.+conc)) # Msun/pc**3

    if(isinstance(xx,collections.Iterable)):
        sigmaInt=np.zeros_like(xx)
        sigmaInt[low]=(4./xx[low]**2)*(rs*1.e3)*deltaC * (2./np.sqrt(1.-xx[low]**2) * np.arctanh(np.sqrt((1.-xx[low])/(1.+xx[low]))) + np.log(xx[low]/2.))
        sigmaInt[mid]=4.*(rs*1.e3)*deltaC * (1. + np.log(1./2.))
        sigmaInt[high]=(4./xx[high]**2)*(rs*1.e3)*deltaC * (2./np.sqrt(xx[high]**2-1.) * np.arctan(np.sqrt((xx[high]-1.)/(1.+xx[high]))) + np.log(xx[high]/2.))
    else: # Rkpc is a scalar
        if(low):
            sigmaInt=(4./xx**2)*(rs*1.e3)*deltaC * (2./np.sqrt(1.-xx**2) * np.arctanh(np.sqrt((1.-xx)/(1.+xx))) + np.log(xx/2.))
        elif(mid):
            sigmaInt=4.*(rs*1.e3)*deltaC * (1. + np.log(1./2.))
        elif(high):
            sigmaInt=(4./xx**2)*(rs*1.e3)*deltaC * (2./np.sqrt(xx**2-1.) * np.arctan(np.sqrt((xx-1.)/(1.+xx))) + np.log(xx/2.))
        else:
            raise ValueError(xx)
    return sigmaInt
    
def deltaSigmaNFW(Rkpc, mass, conc, od):
    """Return DS for GNFW profile in Msun/pc**2.
    Follows Wright & Brainerd 2000
    """
    return sigmaInteriorNFW(Rkpc, mass, conc, od) - sigmaNFW(Rkpc, mass, conc, od)

def rhoGNFW(rkpc, mass, conc, beta, od):
    """Return density of Generalized NFW profile at rkpc in Msun/pc**3.
    Follows Eq. 1 of Wyithe, Turner, & Spergel 2001
    Allows generic overdensity (WTS assumes 200c)
    """
    rscale=scaleRadius(haloRadius(mass, od), conc) # kpc
    rho=od * rhoGNFW0(conc,beta) / ((rkpc/rscale)**beta * (1.+rkpc/rscale)**(3.-beta))
    return rho

def rhoGNFW0(conc,beta):
    """Return dimensionless scale density of GNFW profile.
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
    """Return surface mass density for GNFW profile in Msun/pc**2.
    Follows Eq. 8 of Wyithe, Turner, & Spergel 2001
    """
    Rscale=scaleRadius(haloRadius(mass, od), conc) # kpc
    xx=Rkpc/Rscale
    sigma=2. * od * rhoGNFW0(conc,beta) * (Rscale*1.e3) * xx**(1.-beta) * gnfwF2(xx, beta)
    return sigma

def gnfwF2(xx, beta):
    """Compute integral used to define sigmaGNFW in Eq. 8 of Wyithe, Turner, & Spergel 2001."""
    if(isinstance(xx,collections.Iterable)):
        ff=np.array([scipy.integrate.quad(gnfwF2Integrand,0,np.pi/2.,args=(xxi,beta))[0] for xxi in xx])
    else:
        ff,abserr=scipy.integrate.quad(gnfwF2Integrand,0,np.pi/2.,args=(xx,beta))
    return ff

def gnfwF2Integrand(theta, xx, beta):
    """Compute integrand for gnfwF2 used to get sigmaGNFW."""
    return np.sin(theta)*(np.sin(theta) + xx)**(beta-3.)

def deltaSigmaGNFW(Rkpc, mass, conc, beta, od):
    """Return DS for GNFW profile in Msun/pc**2.
    Simply calls sigmaToDeltaSigma, lacking an analytic profile.
    """
    return sigmaToDeltaSigma(Rkpc, Rkpcfine, sigmaGNFW(Rkpcfine, mass, conc, beta, od))

def rhoHernquist(rkpc, mass, rhern):
    """Return density profile for Hernquist model in Msun/pc**3.
    rhern is the scale radius 'a' used in Eq. 2 of Hernquist 1990
    a=rhalf/(1.+sqrt(2)) where rhalf is the 3d half-mass radius
    """
    rho=mass/(2.*np.pi) * (rhern/rkpc) * 1./(1.e3*(rkpc+rhern))**3
    return rho

def sigmaHernquist(Rkpc, mass, rhern):
    """Return the surface mass density for Hernquist model in Msun/pc**2.
    Follows Eq. 32 of Hernquist 1990
    Assumes Upsilon=1 (i.e. consistently uses mass, not light)
    """
    ss=Rkpc/rhern
    mid=(ss == 1.)
    if(isinstance(ss,collections.Iterable)):
        sigma=np.zeros_like(ss)
        sigma[~mid]=mass/(2.*np.pi*(1.e3*rhern)**2*(1.-ss[~mid]**2)**2) * ((2.+ss[~mid]**2)*hernX(ss[~mid]) - 3.)
        sigma[mid]=2.*mass/(15.*np.pi*(1.e3*rhern)**2)
    else:
        if(mid):
            sigma=2.*mass/(15.*np.pi*(1.e3*rhern)**2)
        else:
            sigma=mass/(2.*np.pi*(1.e3*rhern)**2*(1.-ss**2)**2) * ((2.+ss**2)*hernX(ss) - 3.)
    return sigma

def hernX(ss):
    """Return X(s) from Hernquist 1990 Eqs. 33 and 34.
    ss can be a scalar or iterable type
    """
    low=((ss >= 0) & (ss <= 1))
    mid=(ss==1.)
    high=(ss>1)
    if(isinstance(ss,collections.Iterable)):
        XX=np.zeros_like(ss)
        XX[low]=1./np.sqrt(1.-ss[low]**2) * np.log((1.+np.sqrt(1.-ss[low]**2))/ss[low])
        XX[mid]=1.
        XX[high]=1./np.sqrt(ss[high]**2-1.) * np.arccos(1./ss[high])
    else: # ss is a scalar
        if(low):
            XX=1./np.sqrt(1.-ss**2) * np.log((1.+np.sqrt(1.-ss**2))/ss)
        elif(mid):
            XX=1.
        elif(high):
            XX=1./np.sqrt(ss**2-1.) * np.arccos(1./ss)
    return XX

def deltaSigmaHernquist(Rkpc, mass, rhern):
    """Return DS for Hernquist profile in Msun/pc**2.
    Simply calls sigmaToDeltaSigma, lacking an analytic profile.
    """
    return sigmaToDeltaSigma(Rkpc, Rkpcfine, sigmaHernquist(Rkpcfine, mass, rhern))

def deltaSigmaPS(Rkpc, mass):
    """Return DS for a point source in Msun/pc**2."""
    return mass/(np.pi*(1.e3*Rkpc)**2)

def RdevTorHern(Rdev):
    """Return 3d scale length for Hernquist profile from 2d de Vaucouleurs radius
    Hernquist 1990, Eq. 38
    """
    rhern=Rdev/1.8153
    return rhern

####
# Wrappers to CONTRA for contracted profiles
####
def rhoAC(rkpc, mhalo, conc, od, MAC, isl, nuac, Aac, wac, mstars, rstars, nrad=200):
    """Return contracted profile at rkpc in Msun/pc**3.
    Calls Gnedin's CONTRA code and goes to GNFW beyond rhalo.
    """

    # Get contracted profile from Contra
    BAR=2 # 1 = exponential disk, 2 = Hernquist
    fb=mstars/(mhalo+mstars) # baryon/total mass fraction
    rhalo=haloRadius(mhalo, od)
    rb=rstars/rhalo
    acProfile=contra.contra(MAC, BAR, conc, fb, rb, isl, nuac, Aac, wac, nrad)

    # convert profile to physical units
    rcontra=acProfile[0]*rhalo # kpc
    rhocontra=acProfile[1]*(mhalo+mstars)/(1.e3*rhalo)**3 # Msun/pc**3 - note, contra's mass is defined such that Mtot = Mhalo+Mstars = 1. The rho returned is DM only.

    # extend AC profile using GNFW outside virial radius
    rExt,rhoExt=appendRhoGNFW(rcontra,rhocontra,rhalo,conc,isl,od,extFactor)
    
    # interpolate onto rkpc
    logInterp=scipy.interpolate.UnivariateSpline(np.log10(rExt),np.log10(rhoExt),s=0) # cubic spline interpolation on log axes. Allows extrapolation.
    rhoAC=10.**logInterp(np.log10(rkpc))

    return rhoAC


def appendRhoGNFW(rkpc,rho,rhalo,conc,isl,od,extFactor):
    """Append GNFW density profile from rhalo to extFactor*rhalo with log-spacing similar to input profile."""
    nOld=len(rkpc)
    nExt=nOld * np.log10(extFactor/1.1)/np.log10(rkpc[-1]/rkpc[0]) # how many points on the extension
    rExt=rhalo*1.1 * np.logspace(0,np.log10(extFactor),num=nExt)

    mhalo=haloMass(rhalo, od)
    if(isl==1.):
        rhoExt=rhoNFW(rExt, mhalo, conc, od)
    else:
        rhoExt=rhoGNFW(rExt, mhalo, conc, isl, od)

    rNew=np.append(rkpc,rExt)
    rhoNew=np.append(rho,rhoExt)
    return (rNew,rhoNew)

def sigmaAC(Rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars, nrad=200):
    """Return surface mass density for contracted profile in Msun/pc**2.
    Simply calls rhoToSigma, lacking an analytic profile.
    """
    return rhoToSigma(Rkpc, rkpcfine, rhoAC(rkpcfine, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars, nrad=nrad), rtrunc=truncFactor*haloRadius(mhalo,od))

def deltaSigmaAC(Rkpc, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars, nrad=200):
    """Return DS for contracted profile in Msun/pc**2.
    Simply calls sigmaToDeltaSigma, lacking an analytic profile.
    """
    return sigmaToDeltaSigma(Rkpc, Rkpcfine, sigmaAC(Rkpcfine, mhalo, conc, od, MAC, innerSlopeGNFW, nuDutton, AGnedin, wGnedin, mstars, rstars, nrad=nrad))


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
    # TO DO - check that h-dependence is right - Behroozi assumes h=0.7
    
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
    return logMhalo

def schechterSMF(logSM,phiStar,alpha,logMstar):
    return phiStar * 10.**((logSM - logMstar)* alpha) * np.exp(-10.**(logSM - logMstar))
    
def liwhiteSMF(logSM):

    logSM_overh2 = logSM + np.log10((cosmo.H0()/100.)**2)

    if((np.min(logSM_overh2 < 8.00)) | (np.max(logSM_overh2) > 12.00)):
        raise ValueError(np.min(logSM_overh2),np.max(logSM_overh2))
    
    low=((logSM_overh2 > 8.00) & (logSM_overh2 <= 9.33))
    mid=((logSM_overh2 > 9.33) & (logSM_overh2 <= 10.67))
    high=((logSM_overh2 > 10.67) & (logSM_overh2 <= 12.00))

    phiStar=[0.0146,0.0132,0.0044]
    alpha=[-1.13,-0.90,-1.99]
    logMstar=[9.61,10.37,10.71]
    
    if(isinstance(logSM,collections.Iterable)):
        phi_h3=np.zeros(len(logSM))
        phi_h3[low]=schechterSMF(logSM_overh2[low],phiStar[0],alpha[0],logMstar[0])
        phi_h3[mid]=schechterSMF(logSM_overh2[mid],phiStar[1],alpha[1],logMstar[1])
        phi_h3[high]=schechterSMF(logSM_overh2[high],phiStar[2],alpha[2],logMstar[2])
    else:
        if(low):
            phi_h3=schechterSMF(logSM_overh2,phiStar[0],alpha[0],logMstar[0])
        elif(mid):
            phi_h3=schechterSMF(logSM_overh2,phiStar[1],alpha[1],logMstar[1])
        elif(high):
            phi_h3=schechterSMF(logSM_overh2,phiStar[2],alpha[2],logMstar[2])
        
    phi=phi_h3 * (cosmo.H0()/100.)**3
    
    return phi

def hkC3(xx):
    """Hu & Kravtsov equation C3"""
    return xx**3 * (np.log(1.+1./xx) - 1./(1.+xx))

def hkC11(ff):
    """Hu & Kravtsov equation C11"""
    
    a1=0.5116
    a2=-0.4283
    a3=-3.13e-3
    a4=-3.52e-5

    pp=a2 + a3 * np.log(ff) + a4 * np.log(ff)**2
    
    xx=(a1 * ff**(2.*pp) + (3./4)**2)**(-0.5) + 2 * ff
    return xx
    
def convertHaloMass(logMvir,cvir,redshift,delta,odType):
    """Convert virial halo mass to another overdensity definition.

    Follows Appendix C of Hu & Kravtsov 2003
    Inputs:
        logMvir - log10 virial halo mass (assumes overdensity is deltaVir * rhoMatter)
        cvir - virial radius / scale radius
        redshift
        delta - desired overdensity of output halo mass
        odType - output mass overdensity relative to critical or background
    Return:
        logMhalo - log10 halo mass for specified overdensity
    """

    if(delta=="vir"):
        delta=deltaVir(redshift)
    
    deltaV=deltaVir(redshift)

    f1c=hkC3(1./cvir) # f(1/c) in C8 and C9

    if(odType=="background"):
        dHdV=delta/deltaV # DeltaH/DeltaV, both using background
    elif(odType=="critical"):
        dHdV=delta/(deltaV*cosmo.omega_m()) # DeltaH*rho_c/DeltaV*rho_m
    else:
        raise ValueError(odType)

    fh=dHdV*f1c # rhs arg in C9
    
    xx=hkC11(fh) # C9

    logMhalo=logMvir + np.log10(dHdV * (cvir*xx)**(-3)) # C10

    return logMhalo

def massSize(logSM,type="early"):
    """Mass-size relation from Shen et al 2003 (SDSS) with 2007 Erratum
    Returns logR (kpc) [z-band Sersic half-light radii]
    """

    # TO DO - test , this doesn't seem right, since half-light radii should be larger than deV but these values are pretty small
    
    if(type=="early"): # Eq 17
        bb=2.88e-6 # correction from 2007 Erratum
        aa=0.56
        logR=np.log10(bb) + aa * logSM
    elif(type=="late"): # Eq 18
        gamma=0.10
        alpha=0.14
        beta=0.39
        M0=3.98e10
        SM=10.**logSM
        Rkpc=gamma * SM**alpha * (1.+SM/M0)**(beta-alpha)
        logR=np.log10(Rkpc)

#    return logR
