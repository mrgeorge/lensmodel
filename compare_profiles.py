#! env python

import matplotlib.pyplot as plt
import lensmodel
import numpy as np

# Input model pars
redshift=0.11
mstars=7.e10
mhalo=0.7e13
conc=9.
rstars=4.3/(1.+redshift) # convert Schulz' comoving coords to physical
innerSlopeGNFW=1.
nuDuttonAC=1.
nuDuttonExp=-0.2
AGnedin=-1.
wGnedin=-1.
cenType="hernquist"
odType="background"
delta=200.

od=lensmodel.profiles.overdensity(redshift=redshift,delta=delta,type=odType)

Rmin=1. # kpc - note this is the geometric bin center, so measurements would need to go to radii 0.5*bin width smaller
Rmax=1000. # kpc - Use only bins at R<1Mpc
nBins=50
dlog10R=np.log10(Rmax/Rmin)/(nBins-1)
Rkpc=Rmin*10.**(np.arange(nBins) * dlog10R)
rkpc=Rkpc

nradAC=600

rhoNFW=lensmodel.profiles.rhoNFW(rkpc, mhalo, conc, od)
rhoAC=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 2, innerSlopeGNFW, nuDuttonAC, AGnedin, wGnedin, mstars, rstars, nrad=nradAC)
rhoExp=lensmodel.profiles.rhoAC(rkpc, mhalo, conc, od, 2, innerSlopeGNFW, nuDuttonExp, AGnedin, wGnedin, mstars, rstars, nrad=nradAC)

sigmaNFW=lensmodel.profiles.sigmaNFW(Rkpc, mhalo, conc, od)
sigmaAC=lensmodel.profiles.sigmaAC(Rkpc, mhalo, conc, od, 2, innerSlopeGNFW, nuDuttonAC, AGnedin, wGnedin, mstars, rstars, nrad=nradAC)
sigmaExp=lensmodel.profiles.sigmaAC(Rkpc, mhalo, conc, od, 2, innerSlopeGNFW, nuDuttonExp, AGnedin, wGnedin, mstars, rstars, nrad=nradAC)

deltaSigmaNFW=lensmodel.profiles.deltaSigmaNFW(Rkpc, mhalo, conc, od)
deltaSigmaAC=lensmodel.profiles.deltaSigmaAC(Rkpc, mhalo, conc, od, 2, innerSlopeGNFW, nuDuttonAC, AGnedin, wGnedin, mstars, rstars, nrad=nradAC)
deltaSigmaExp=lensmodel.profiles.deltaSigmaAC(Rkpc, mhalo, conc, od, 2, innerSlopeGNFW, nuDuttonExp, AGnedin, wGnedin, mstars, rstars, nrad=nradAC)

plt.clf()
figsize=(10,5)
logylim=(1.,1.e3)
linylim=(0.75,1.09)


plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes',linewidth=1.5)

fig,axarr=plt.subplots(2,3,figsize=figsize)
fig.subplots_adjust(hspace=0,wspace=0.5)

plt.sca(axarr[0,0])
lensmodel.plot.plotProfile(rkpc,rhoAC,xlim=(Rmin,Rmax),xscale="log",yscale="log",ylabel=r"$\rho~($M$_{\odot}~$pc$^{-3})$",label="AC",color="blue",ls="--")
lensmodel.plot.plotProfile(rkpc,rhoNFW,label="NFW",color="black",ls="-")
lensmodel.plot.plotProfile(rkpc,rhoExp,label="Expansion",color="green",ls=":")
plt.setp(axarr[0,0].get_xticklabels(),visible=False)
plt.legend(loc="lower left", frameon=False, prop={'size':12})

plt.sca(axarr[0,1])
lensmodel.plot.plotProfile(Rkpc,sigmaAC,xlim=(Rmin,Rmax),ylim=logylim,xscale="log",yscale="log",ylabel=r"$\Sigma~($M$_{\odot}~$pc$^{-2})$",color="blue",ls="--")
lensmodel.plot.plotProfile(Rkpc,sigmaNFW)
lensmodel.plot.plotProfile(Rkpc,sigmaExp,color="green",ls=":")
plt.setp(axarr[0,1].get_xticklabels(),visible=False)

plt.sca(axarr[0,2])
lensmodel.plot.plotProfile(Rkpc,deltaSigmaAC,xlim=(Rmin,Rmax),ylim=logylim,xscale="log",yscale="log",ylabel=r"$\Delta\Sigma~($M$_{\odot}~$pc$^{-2})$",color="blue",ls="--")
lensmodel.plot.plotProfile(Rkpc,deltaSigmaNFW)
lensmodel.plot.plotProfile(Rkpc,deltaSigmaExp,color="green",ls=":")
plt.setp(axarr[0,2].get_xticklabels(),visible=False)

plt.sca(axarr[1,0])
lensmodel.plot.plotProfile(rkpc,rhoAC/rhoNFW,xlim=(Rmin,Rmax),ylim=linylim,xscale="log",yscale="linear",ylabel=r"$\rho/\rho_{NFW}$",xlabel=r"$r$ (kpc)",color="blue",ls="--")
lensmodel.plot.plotProfile(rkpc,rhoExp/rhoNFW,yscale="linear",color="green",ls=":")

plt.sca(axarr[1,1])
lensmodel.plot.plotProfile(Rkpc,sigmaAC/sigmaNFW,xlim=(Rmin,Rmax),ylim=linylim,xscale="log",yscale="linear",ylabel=r"$\Sigma/\Sigma_{NFW}$",xlabel=r"$R$ (kpc)",color="blue",ls="--")
lensmodel.plot.plotProfile(rkpc,sigmaExp/sigmaNFW,yscale="linear",color="green",ls=":")

plt.sca(axarr[1,2])
lensmodel.plot.plotProfile(Rkpc,deltaSigmaAC/deltaSigmaNFW,xlim=(Rmin,Rmax),ylim=linylim,xscale="log",yscale="linear",ylabel=r"$\Delta\Sigma/\Delta\Sigma_{NFW}$",xlabel=r"$R$ (kpc)",color="blue",ls="--")
lensmodel.plot.plotProfile(rkpc,deltaSigmaExp/deltaSigmaNFW,yscale="linear",color="green",ls=":")


plt.savefig("../plots/compareAC.pdf")