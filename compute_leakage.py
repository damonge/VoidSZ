import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pymaster as nmt
import os

nside=2048
n_rb=32
lmax=1000
lmin=50

def get_cells(fsky_mask) :
    fname_cls="data/cls_A%d.txt"%fsky_mask

    if os.path.isfile(fname_cls) :
        leff,cl_cc,cl_cy,cl_yy=np.loadtxt(fname_cls,unpack=True)
    else :
        mp_545=hp.read_map("data/HFI_SkyMap_545_2048_R2.02_full.fits",verbose=False)*1E6
        mp_y=hp.read_map("data/data_y/milca_ymaps.fits",verbose=False)
        msk=hp.read_map("data/data_y/mask_planck%d.fits"%fsky_mask,verbose=False)
        
        fsk=np.mean(msk)
        print "cc"
        cl_cc=hp.anafast(mp_545*msk,map2=mp_545*msk,iter=0)/fsk
        print "cy"
        cl_cy=hp.anafast(mp_y*msk,mp_545*msk,iter=0)/fsk
        print "yy"
        cl_yy=hp.anafast(mp_y*msk,mp_y*msk,iter=0)/fsk
        leff=np.arange(len(cl_yy))+0.

        np.savetxt(fname_cls,np.transpose([leff,cl_cc,cl_cy,cl_yy]))

    return leff,cl_cc,cl_cy,cl_yy

l_cc_mid,cc_mid=np.loadtxt("data/cl_cib_cib.txt",unpack=True); x_cc_mid=np.log10(l_cc_mid/100.)
cc_midf=interp1d(x_cc_mid,np.log10(cc_mid))
l_cc_midb,cc_midb=np.loadtxt("data/cl_cib_cib_c.txt",unpack=True); x_cc_midb=np.log10(l_cc_midb/100.)
cc_midbf=interp1d(x_cc_midb,np.log10(cc_midb))
l_cy_mid,cy_mid=np.loadtxt("data/cl_cib_tsz.txt",unpack=True); x_cy_mid=np.log10(l_cy_mid/100.)
cy_midf=interp1d(x_cy_mid,np.log10(cy_mid))
def cibcib(l,use_b=False) :
    """CIB auto power spectrum in units of Jy^2 sr^-1 (i.e. CIB is in units of Jy sr^-1)
    """
    x=np.log10(l/100.)

    if use_b :
        ids_lo=np.where(x<=x_cc_midb[0])[0]
        ids_mid=np.where((x>x_cc_midb[0]) & (x<x_cc_midb[-1]))[0]
        ids_hi=np.where(x>=x_cc_midb[-1])[0]
        
        ret=np.zeros_like(l)
        tilt_lo=np.log10(cc_midb[1]/cc_midb[0])/(x_cc_midb[1]-x_cc_midb[0])
        offs_lo=np.log10(cc_midb[0])
        tilt_hi=np.log10(cc_midb[-1]/cc_midb[-2])/(x_cc_midb[-1]-x_cc_midb[-2])
        offs_hi=np.log10(cc_midb[-1])
        ret[ids_lo]=offs_lo+tilt_lo*(x[ids_lo]-x_cc_midb[0])
        ret[ids_hi]=offs_hi+tilt_hi*(x[ids_hi]-x_cc_midb[-1])
        ret[ids_mid]=cc_midbf(x[ids_mid])
    else :
        ids_lo=np.where(x<=x_cc_mid[0])[0]
        ids_mid=np.where((x>x_cc_mid[0]) & (x<x_cc_mid[-1]))[0]
        ids_hi=np.where(x>=x_cc_mid[-1])[0]
        
        ret=np.zeros_like(l)
        tilt_lo=np.log10(cc_mid[1]/cc_mid[0])/(x_cc_mid[1]-x_cc_mid[0])
        offs_lo=np.log10(cc_mid[0])
        tilt_hi=np.log10(cc_mid[-1]/cc_mid[-2])/(x_cc_mid[-1]-x_cc_mid[-2])
        offs_hi=np.log10(cc_mid[-1])
        ret[ids_lo]=offs_lo+tilt_lo*(x[ids_lo]-x_cc_mid[0])
        ret[ids_hi]=offs_hi+tilt_hi*(x[ids_hi]-x_cc_mid[-1])
        ret[ids_mid]=cc_midf(x[ids_mid])
        ret+=3.

    return 10.**ret
def ciby(l) :
    """y-CIB power spectrum in units of Jy (i.e. CIB is in Jy sr^-1 and y is unitless)
    """

    x=np.log10(l/100.)

    ids_lo=np.where(x<=x_cy_mid[0])[0]
    ids_mid=np.where((x>x_cy_mid[0]) & (x<x_cy_mid[-1]))[0]
    ids_hi=np.where(x>=x_cy_mid[-1])[0]
    
    ret=np.zeros_like(l)
    tilt_lo=np.log10(cy_mid[1]/cy_mid[0])/(x_cy_mid[1]-x_cy_mid[0])
    offs_lo=np.log10(cy_mid[0])
    tilt_hi=np.log10(cy_mid[-1]/cy_mid[-2])/(x_cy_mid[-1]-x_cy_mid[-2])
    offs_hi=np.log10(cy_mid[-1])
    ret[ids_lo]=offs_lo+tilt_lo*(x[ids_lo]-x_cy_mid[0])
    ret[ids_hi]=offs_hi+tilt_hi*(x[ids_hi]-x_cy_mid[-1])
    ret[ids_mid]=cy_midf(x[ids_mid])

    return 1.6*10.**(ret-6.)


def rebin(a) :
    return np.mean(a.reshape([len(a)/n_rb,n_rb]),axis=1)

msk=hp.read_map("data/data_y/mask_planck80.fits",verbose=False)
fsky=np.mean(msk)

l80,cl_cc80,cl_cy80,cl_yy80=get_cells(80)
beam_y=np.exp(-0.5*l80*(l80+1)*(10./2.355*(np.pi/180./60.))**2)
beam_c=np.exp(-0.5*l80*(l80+1)*(5./2.355*(np.pi/180./60.))**2)
l80_b=rebin(l80); cl_cc80_b=rebin(cl_cc80); cl_cy80_b=rebin(cl_cy80); cl_yy80_b=rebin(cl_yy80)
cl_cy_th=rebin(ciby(l80)*beam_y*beam_c)
cl_cc_th=rebin(cibcib(l80)*beam_c*beam_c)
cl_gg_th=cl_cc80_b-cl_cc_th
ids=np.where((l80_b>lmin) & (l80_b<lmax))[0]
err_cy=np.sqrt((cl_yy80_b*cl_cc80_b+cl_cy80_b**2)/((2*l80_b+1.)*fsky*n_rb))

a_cc=np.sum(cl_cc_th[ids]*cl_cc_th[ids]/err_cy[ids]**2)
a_cg=np.sum(cl_cc_th[ids]*cl_gg_th[ids]/err_cy[ids]**2)
a_gg=np.sum(cl_gg_th[ids]*cl_gg_th[ids]/err_cy[ids]**2)
e_c=np.sum((cl_cy80_b[ids]-cl_cy_th[ids])*cl_cc_th[ids]/err_cy[ids]**2)
e_g=np.sum((cl_cy80_b[ids]-cl_cy_th[ids])*cl_gg_th[ids]/err_cy[ids]**2)
a_mat=np.array([[a_cc,a_cg],[a_cg,a_gg]])
e_vec=np.array([e_c,e_g])
covar=np.linalg.inv(a_mat)
sol=np.dot(covar,e_vec)*1E6
err=np.sqrt(np.diag(covar))*1E6
res=(cl_cy80_b[ids]-cl_cy_th[ids]-sol[0]*cl_cc_th[ids]*1E-6-sol[1]*cl_gg_th[ids]*1E-6)/err_cy[ids]
print "alpha_CIB = %lE +- %lE (MJy/sr)^-1"%(sol[0],err[0])
print "alpha_Gal = %lE +- %lE (MJy/sr)^-1"%(sol[1],err[1])
print "chi^2/dof = %lE"%(np.sum(res**2)/(len(res)-2))

f=open("data/data_y/y_contamination.txt","w")
stout ="#[0]-a_CIB [1]-error(a_CIB) [2]-a_Gal [3]-error(a_Gal)\n"
stout+=" %lE %lE %lE %lE\n"%(sol[0],err[0],sol[1],err[1])
f.write(stout)
f.close()

plt.figure()
plt.plot(l80_b[ids],(l80_b*(l80_b+1.)*cl_cc80_b)[ids]/(2*np.pi),'r-',label='545 x 545')
plt.plot(l80_b[ids],(l80_b*(l80_b+1.)*cl_cc_th)[ids]/(2*np.pi) ,'g-',label='CIB x CIB')
plt.plot(l80_b[ids],(l80_b*(l80_b+1.)*cl_gg_th)[ids]/(2*np.pi) ,'k-',label='Gal x Gal')
plt.loglog()

plt.figure()
plt.errorbar(l80_b[ids],cl_cy80_b[ids],yerr=err_cy[ids],color='r',label='545 x y')
plt.plot(l80_b[ids],cl_cy_th[ids],'g-',label='CIB x y')
plt.loglog()
plt.show()
