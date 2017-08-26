import numpy as np
import scipy.integrate as itg
import py_cosmo_mad as csm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import brentq

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

RHOCRIT0=2.7744948E11
OM0=0.3
HH0=0.7
FB0=0.1666666
NS=0.96
S8=0.8
NRMAX=10.
LMMIN=8.
LMMAX=16.
DLM=0.1
NLM=64
B_CMASS=2.
Z_EFF=0.5

#Cosmological parameters outside the voide
pcs_BG0=csm.PcsPar()
pcs_BG0.unset_gsl_eh()
pcs_BG0.background_set(OM0,1-OM0,OM0*FB0,-1.,0.,HH0,2.7255)
dgrowth0=pcs_BG0.growth_factor(1.) #Normalization of the growth factor
#Factor needed to scale the void profile at z=Z_EFF to z=0
scale_delta=B_CMASS*pcs_BG0.growth_factor(1./(1+Z_EFF))/dgrowth0

def rDelta(m,z,Delta,cpar) :
    """Returns r_Delta
    """
    hn=(cpar['pcs']).hubble(1./(1+z))/(cpar['pcs']).hubble(1.)
    rhoc=RHOCRIT0*hn*hn
    return (3*m/(4*np.pi*Delta*rhoc))**0.333333333*(1+z)

def get_battaglia(m,z,cpar) :
    """Sets all parameters needed to compute the Battaglia et al. profile
    """
    fi=0.518
    fb=cpar['OB']/cpar['OM']
    ez2=((cpar['pcs']).hubble(1./(1+z))/(cpar['pcs']).hubble(1))**2
    h=cpar['h']
    mr=m/(1E14*h)
    p0=18.1*mr**0.154*(1+z)**(-0.758)
    delta=200
    r200=rDelta(m,z,delta,cpar)
    dic={'ups0':2.026212E-21*delta*m*ez2*h*fi*fb*p0*(1+z)/r200,
         'r200':r200,
         'xc':0.497*(mr**(-0.00865))*((1+z)**0.731),
         'beta':4.35*(mr**0.0393)*((1+z)**0.415)}
    return dic

def ups_battaglia(x,bp) :
    """Battaglia pressure profile in units of pressure * sigma_T/(m_e c^2)
    x = r/r200
    """
    xr=x/bp['xc']
    return bp['ups0']*(xr**(-0.3))*((1+xr)**(-bp['beta']))

def integrated_profile(bp) :
    """Volume integral of the Battaglia pressure profile
    """
    def integrand(x) :
        return x**2*ups_battaglia(x,bp)
    return 4*np.pi*(bp['r200'])**3*itg.quad(integrand,0,NRMAX)[0]

def get_zeff(om,ol,ok,hh,z0) :
    """ Returns equivalent redshift in effective universe
    by equating cosmic times
    """
    def t_intg(x,o_m,o_l,o_k,h_h) :
        return np.sqrt(x/(o_m+o_l*x**3+o_k*x))/h_h
    tBG=itg.quad(t_intg,0,1./(1+z0),args=(OM0,1-OM0,0.,HH0))[0]
    def fzero(zz) :
        return itg.quad(t_intg,0,1./(1+zz),args=(om,ol,ok,hh))[0]-tBG
    zeff=brentq(fzero,0.,10.)
    return zeff

def get_pz(om,ol,ob,hh,get_nm=False) :
    """Average electron pressure in the relevant redshift range
    for a set of cosmological parameters"""
    pcs=csm.PcsPar()
    pcs.unset_gsl_eh()
    pcs.background_set(om,ol,ob,-1.,0.,hh,2.7255)
    g_ratio=pcs.growth_factor(1.)/dgrowth0 #Growth ratio with respect to background universe
    pcs.set_linear_pk('EH',-5.,3.,0.01,NS,S8*g_ratio)
    pcs.set_mf_params(LMMIN,LMMAX,DLM)
    cospar={'OM':om,'OB':ob,'h':hh,'pcs':pcs}

    marr=np.logspace(LMMIN,LMMAX,NLM)
    zarr=np.array([Z_EFF])
    pzarr=np.zeros_like(zarr)
    for iz,z in enumerate(zarr) :
        zeff=get_zeff(om,ol,1-om-ol,hh,z) #Equivalent redshift in the effective Universe
        nmarr=np.zeros_like(marr)
        for im,m in enumerate(marr) :
            bp=get_battaglia(m,zeff,cospar)
            pmean=integrated_profile(bp)
            nmarr[im]=pmean*pcs.mass_function_logarithmic(m,zeff,'Tinker10_200')
        nmf=interp1d(np.log10(marr),nmarr,bounds_error=False,fill_value=0)
        pzarr[iz]=itg.quad(nmf,7.5,LMMAX+0.5)[0]
    if get_nm :
        return marr,nmarr,np.mean(pzarr)
    else :
        return np.mean(pzarr)

#Compute mean void radius
data_ngc_cut=np.genfromtxt('data/voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5.dat',
                           dtype='i8,f8,f8,f8,i8,f8,f8,f8,f8,f8,f8,f8',
                           names=['VoidID','RA','Dec','z','ng','Vol','Reff','n_dens',
                                  'delta_core','delta_void','prob','R_edge'])
data_sgc_cut=np.genfromtxt('data/voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5.dat',
                           dtype='i8,f8,f8,f8,i8,f8,f8,f8,f8,f8,f8,f8',
                           names=['VoidID','RA','Dec','z','ng','Vol','Reff','n_dens',
                                  'delta_core','delta_void','prob','R_edge'])
data_cut=np.hstack((data_ngc_cut,data_sgc_cut))
rv_mean=np.mean(data_cut['Reff']/(1+data_cut['z']))
print "Mean void scale : %.3lf Mpc/h"%rv_mean

#Read 3D void profile
dum,r,dum,dum,dum,d=np.loadtxt("data/profiles.txt",unpack=True)
rarrt=np.zeros(len(r)+1); rarrt[1:]=r
darrt=np.zeros(len(d)+1); darrt[0]=d[0]*0.9; darrt[1:]=d

#deltaf=interp1d(rarrt,darrt-1,bounds_error=False,fill_value=0,kind='cubic')
rarr=np.linspace(0.,5.,150)
deltaf=interp1d(rarrt,(darrt-1)/scale_delta,bounds_error=False,fill_value=0)

#Compute tBB for background cosmology
def tintg0(x) :
    return np.sqrt(x/(OM0+(1-OM0)*x**3))
tbb0=itg.quad(tintg0,0,1)[0]

def delta_fit(r) :
    """Fitting funtion to the void profile
    """
    delta_cmass=-0.886*np.exp(-(r/0.65)**2)+0.23*np.exp(-((r-0.93)/0.4)**2)
    return delta_cmass/scale_delta

def rhomean(r) :
    """Compute mean density within a radius r
    """
    if r==0 :
        return 1+delta_fit(r)
    else :
        def intg(s) :
            return s*s*(1+delta_fit(s))
        return itg.quad(intg,0,r)[0]/(r**3/3)

def tbb(eta2,rhom) :
    """tBB in an effective cosmology with rho_mean=rhom and eta^2=eta2
    """
    om=OM0*rhom/eta2
    ol=(1-OM0)/eta2
    ok=1-om-ol
    def tintg(x) :
        return np.sqrt(x/(eta2*(om+ol*x**3+ok*x)))
    return itg.quad(tintg,0,1)[0]-tbb0

def e2ofrhom(rhom) :
    """Compute eta^2 from a given rho_mean=rhom by matching tBB
    to the background cosmology"""
    return brentq(tbb,0.8,3.,args=(rhom))

rhomarr=np.array([rhomean(r) for r in rarr]) #Array of rho_mean(r)
eta2arr=np.array([e2ofrhom(rhom) for rhom in rhomarr]) #Array of eta^2(r)
omarr=OM0*rhomarr/eta2arr #Effective matter parameter
olarr=(1-OM0)/eta2arr #Effective cosmological constant
okarr=1-olarr-omarr #Effective curvature
harr=HH0*np.sqrt(eta2arr) #Effective expansion rate


omf=interp1d(rarr,omarr)
olf=interp1d(rarr,olarr)
hhf=interp1d(rarr,harr)
r_plot=[0.,0.3,0.6,1.0]
fmts=['r:','r-.','r--','r-']
plt.figure()
ax=plt.gca()
for r,fmt in zip(r_plot,fmts) :
    mm,nmm,pzz=get_pz(omf(r),olf(r),FB0*omf(r),hhf(r),get_nm=True)
    plt.plot(mm,nmm/pzz,fmt,lw=2,label='$r/r_v=%.2lf$'%r)
mm,nmm,pzz=get_pz(omf(4),olf(4),FB0*omf(4),hhf(4),get_nm=True)
ax.plot(mm,nmm/pzz,'k-',lw=1,label='$r/r_v=\\infty$')
ax.set_xlim([1E10,5E15])
ax.set_ylim([0,0.6])
ax.set_xlabel('$M\\,[h^{-1}\\,M_\\odot]$',fontsize=14)
ax.set_ylabel('$d\\,\\log\\langle P_e\\rangle/d\\log_{10}M$',fontsize=14)
ax.set_xscale('log')
plt.legend(loc='upper left',frameon=False,labelspacing=0.1,fontsize=14,ncol=2)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(12)
plt.savefig('doc/pe_mass.pdf',bbox_inches='tight')


parr=np.zeros_like(harr) #This will hold mean electron pressure(r)
for i in np.arange(len(parr)) :
    parr[i]=get_pz(omarr[i],olarr[i],FB0*omarr[i],harr[i])
parr*=harr/harr[-1] #Correct for h-inverse units
parr-=parr[-1] #Subtract background (which we don't measure)

#Integrate pressure profile along the line of sight
pf=interp1d(rarr,parr,bounds_error=False,fill_value=parr[-1])
def integ_project(rl,rt) :
    r=np.sqrt(rl**2+rt**2)
    return pf(r)
pparr=np.array([2*itg.quad(integ_project,0.,5.,args=(r,))[0] for r in rarr])*rv_mean

np.savetxt("data/data_y/void_univ_eff.txt",np.transpose([rarr,delta_fit(rarr),omarr,olarr,okarr,harr]))
np.savetxt("data/data_y/y_th_void.txt",np.transpose([rarr,pparr,rv_mean*parr]))

plt.figure(); plt.xlim([0,2]); plt.plot(rarrt,(darrt-1)/scale_delta+1); plt.plot(rarr,1+delta_fit(rarr));
plt.figure(); plt.xlim([0,5]); plt.plot(rarr,omarr,'r-'); plt.plot(rarr,olarr,'g-'); plt.plot(rarr,okarr,'y-')
plt.figure(); plt.xlim([0,2]); plt.plot(rarr,harr)
plt.figure(); plt.plot(rarr,rv_mean*parr,'r--'); plt.plot(rarr,pparr,'r-');
plt.xlim([0,2]); plt.ylim([-7E-8,7E-8])
plt.show()
