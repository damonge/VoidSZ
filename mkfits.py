import numpy as np
import pyfits as pf
import py_cosmo_mad as csm
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits

seed=1001
np.random.seed(seed)
r=hp.Rotator(coord=['C','G'])
nside=2048

print " Subsampling randoms"
data_n=(fits.open('data/random0_DR12v5_CMASS_North.fits'))[1].data
data_s=(fits.open('data/random0_DR12v5_CMASS_South.fits'))[1].data
data=np.hstack((data_n,data_s))
print len(data)
ids=np.random.choice(np.arange(len(data)),size=2200000,replace=False)
print len(ids)
tbhdu=pf.new_table([pf.Column(name='Z',format='D',array=(data['Z'])[ids]),
                    pf.Column(name='RA',format='D',array=(data['RA'])[ids]),
                    pf.Column(name='DEC',format='D',array=(data['DEC'])[ids])])
tbhdu.writeto("data/random0_DR12v5_CMASS_Subsample.fits",clobber=True)

print " Computing BOSS mask"
nside_lo=128
th_c=(90-data['DEC'][ids])*np.pi/180
phi_c=data['RA'][ids]*np.pi/180
th_g,phi_g=r(th_c,phi_c)
ipix=hp.ang2pix(nside_lo,th_g,phi_g)
ngals,b=np.histogram(ipix,range=[0,hp.nside2npix(nside_lo)],bins=hp.nside2npix(nside_lo))
msk_lo=np.zeros(hp.nside2npix(nside_lo)); msk_lo[ngals>0]=1;
msk_hi=hp.ud_grade(msk_lo,nside_out=nside)
hp.write_map("data/mask_DR12.fits",msk_hi)

print " Cosmological parameters"
RHOCRIT0=2.7744948E11
HH=0.7
pcs=csm.PcsPar()
pcs.background_set(0.3,0.7,0.05,-1.,0,HH,2.7255)

print " Making SZ point source mask"
def rDelta(m,zz,Delta) :
    """Returns r_Delta
    """
    hn=np.array([pcs.hubble(1./(1+z))/pcs.hubble(1.) for z in zz])
    rhoc=RHOCRIT0*hn*hn
    return (3*m/(4*np.pi*Delta*rhoc))**0.333333333*(1+zz)

#Read catalog and remove all clusters above z=0.43
data=(fits.open('data/data_y/HFI_PCCS_SZ-union_R2.08.fits'))[1].data
ids=np.where(data['REDSHIFT']>=0)[0]; data=data[ids]
idsz=np.where(data['REDSHIFT']<0.43)[0]; data=data[idsz]
#Compute their angular extent
r500=rDelta(data['MSZ']*HH*1E14,data['REDSHIFT'],500)
chi=np.array([pcs.radial_comoving_distance(1./(1+z)) for z in data['REDSHIFT']])
th500=r500/chi
#Compute angular positions for each cluster
theta=(90-data['GLAT'])*np.pi/180
phi=data['GLON']*np.pi/180
vx=np.sin(theta)*np.cos(phi)
vy=np.sin(theta)*np.sin(phi)
vz=np.cos(theta)
#Generate mask by cutting out a circle of radius
#3*theta_500 around each cluster
mask_sz=np.ones(hp.nside2npix(nside))
for i in np.arange(len(data)) :
    v=np.array([vx[i],vy[i],vz[i]])
    radius=3*th500[i]
    ip=hp.query_disc(nside,v,radius)
    mask_sz[ip]=0
hp.write_map("data/data_y/mask_sz.fits",mask_sz)

#Make mask
print " Mask"
mask_sz=hp.read_map("data/data_y/mask_sz.fits")
mask_gal_80=hp.read_map("data/data_y/HFI_Mask_GalPlane-apo0_2048_R2.00.fits",verbose=False,field=4)
mask_gal_60=hp.read_map("data/data_y/HFI_Mask_GalPlane-apo0_2048_R2.00.fits",verbose=False,field=2)
mask_gal_40=hp.read_map("data/data_y/HFI_Mask_GalPlane-apo0_2048_R2.00.fits",verbose=False,field=1)
mask_gal_20=hp.read_map("data/data_y/HFI_Mask_GalPlane-apo0_2048_R2.00.fits",verbose=False,field=0)
mask_p0=hp.read_map("data/data_y/LFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,hdu=1); print np.mean(mask_p0)
mask_p1=hp.read_map("data/data_y/LFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,hdu=2); print np.mean(mask_p1)
mask_p2=hp.read_map("data/data_y/LFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,hdu=3); print np.mean(mask_p2)
mask_pl=mask_p0*mask_p1*mask_p2
mask_p0=hp.read_map("data/data_y/HFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,field=0); print np.mean(mask_p0)
mask_p1=hp.read_map("data/data_y/HFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,field=1); print np.mean(mask_p1)
mask_p2=hp.read_map("data/data_y/HFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,field=2); print np.mean(mask_p2)
mask_p3=hp.read_map("data/data_y/HFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,field=3); print np.mean(mask_p3)
mask_ph=mask_p0*mask_p1*mask_p2*mask_p3
hp.write_map("data/data_y/mask_planck20.fits",mask_gal_20*mask_ph)
hp.write_map("data/data_y/mask_planck40.fits",mask_gal_40*mask_ph)
hp.write_map("data/data_y/mask_planck60.fits",mask_gal_60*mask_ph)
hp.write_map("data/data_y/mask_planck80.fits",mask_gal_80*mask_ph)
hp.write_map("data/data_y/mask_planck60L.fits",mask_gal_60*mask_ph*mask_pl)
hp.write_map("data/data_y/mask_planck80L.fits",mask_gal_80*mask_ph*mask_pl)
hp.write_map("data/data_y/mask_planck20S.fits",mask_gal_20*mask_ph*mask_sz)
hp.write_map("data/data_y/mask_planck40S.fits",mask_gal_40*mask_ph*mask_sz)
hp.write_map("data/data_y/mask_planck60S.fits",mask_gal_60*mask_ph*mask_sz)
hp.write_map("data/data_y/mask_planck80S.fits",mask_gal_80*mask_ph*mask_sz)
hp.write_map("data/data_y/mask_planck20LS.fits",mask_gal_20*mask_ph*mask_pl*mask_sz)
hp.write_map("data/data_y/mask_planck40LS.fits",mask_gal_40*mask_ph*mask_pl*mask_sz)
hp.write_map("data/data_y/mask_planck60LS.fits",mask_gal_60*mask_ph*mask_pl*mask_sz)
hp.write_map("data/data_y/mask_planck80LS.fits",mask_gal_80*mask_ph*mask_pl*mask_sz)

#Make null map
print " Null maps"
y_first,y_last=hp.read_map("data/data_y/milca_ymaps.fits",field=[1,2],verbose=False)
y_null=0.5*(y_first-y_last)
hp.write_map("data/data_y/y_null_milca.fits",y_null)
y_first,y_last=hp.read_map("data/data_y/nilc_ymaps.fits",field=[1,2],verbose=False)
y_null=0.5*(y_first-y_last)
hp.write_map("data/data_y/y_null_nilc.fits",y_null)

# [1] VoidID, [2] Center RA (deg), [3] Center Dec (deg), [4] Center Redshift, [5] Number of Galaxies, [6] Total Voronoi Volume (Mpc^3/h^3), [7] Effective Radius (Mpc/h), [8] Core Number Density (h^3/Mpc^3), [9] Core Density Contrast, [10] Void Density Contrast, [11] Probability [12] Distance from Center to Edge 

def mk_fits_data(fname_n,fname_s,fname_out) :
    data_ngc=np.genfromtxt(fname_n,
                           dtype='i8,f8,f8,f8,i8,f8,f8,f8,f8,f8,f8,f8',
                           names=['VoidID','RA','Dec','z','ng','Vol','Reff','n_dens',
                                  'delta_core','delta_void','prob','R_edge'])
    data_sgc=np.genfromtxt(fname_s,
                           dtype='i8,f8,f8,f8,i8,f8,f8,f8,f8,f8,f8,f8',
                           names=['VoidID','RA','Dec','z','ng','Vol','Reff','n_dens',
                                  'delta_core','delta_void','prob','R_edge'])
    data=np.hstack((data_sgc,data_ngc))
    #Rotate to galactic coordinates
    th_c=(90-data['Dec'])*np.pi/180.
    phi_c=data['RA']*np.pi/180.
    th_g,phi_g=r(th_c,phi_c)
    b=90-180*th_g/np.pi
    l=180*phi_g/np.pi
    #Distance
    chi=np.array([pcs.radial_comoving_distance(1./(1+z)) for z in data['z']])
    #Projected size
    theta_eff=data['Reff']/chi
    theta_edge=data['R_edge']/chi
    print len(data_ngc), len(data_sgc), len(data)

    tbhdu=pf.new_table([pf.Column(name='VoidID',format='K',array=data['VoidID']),
                        pf.Column(name='RA'    ,format='D',array=data['RA']),
                        pf.Column(name='Dec'   ,format='D',array=data['Dec']),
                        pf.Column(name='L'     ,format='D',array=l),
                        pf.Column(name='B'     ,format='D',array=b),
                        pf.Column(name='z'     ,format='D',array=data['z']),
                        pf.Column(name='ng'    ,format='K',array=data['ng']),
                        pf.Column(name='Vol'   ,format='D',array=data['Vol']),
                        pf.Column(name='Reff'  ,format='D',array=data['Reff']),
                        pf.Column(name='n_dens',format='D',array=data['n_dens']),
                        pf.Column(name='delta_core',format='D',array=data['delta_core']),
                        pf.Column(name='delta_void',format='D',array=data['delta_void']),
                        pf.Column(name='prob',format='D',array=data['prob']),
                        pf.Column(name='R_edge',format='D',array=data['R_edge']),
                        pf.Column(name='rdist',format='D',array=chi),
                        pf.Column(name='theta_eff',format='D',array=theta_eff),
                        pf.Column(name='theta_edge',format='D',array=theta_edge)])
    tbhdu.writeto(fname_out,clobber=True)

print " Data catalogs"
mk_fits_data('data/voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5.dat',
             'data/voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5.dat',
             'data/voids_BOSS_cmass_dr12v4_Om0.3_dt0.5.fits')
mk_fits_data('data/voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5_uncut.dat',
             'data/voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5_uncut.dat',
             'data/voids_BOSS_cmass_dr12v4_Om0.3_dt0.5_uncut.fits')
mk_fits_data('data/voids_BOSS_lowz_dr12v4_ngc_Om0.3_dt0.5.dat',
             'data/voids_BOSS_lowz_dr12v4_sgc_Om0.3_dt0.5.dat',
             'data/voids_BOSS_lowz_dr12v4_Om0.3_dt0.5.fits')
mk_fits_data('data/voids_BOSS_lowz_dr12v4_ngc_Om0.3_dt0.5_uncut.dat',
             'data/voids_BOSS_lowz_dr12v4_sgc_Om0.3_dt0.5_uncut.dat',
             'data/voids_BOSS_lowz_dr12v4_Om0.3_dt0.5_uncut.fits')

def mkfits_mock(imock) :
    fname_ngc="data/mocks/voids_QPM_a0.6452_dr12c_%04d_cmass_ngc_zspace_Om0.3_dt0.5.dat"%imock
    fname_sgc="data/mocks/voids_QPM_a0.6452_dr12c_%04d_cmass_sgc_zspace_Om0.3_dt0.5.dat"%imock
    fname_end="data/mocks/voids_QPM_a0.6452_dr12c_%04d_cmass_zspace_Om0.3_dt0.5.fits"%imock
    data_ngc_cut=np.genfromtxt(fname_ngc,
                               dtype='i8,f8,f8,f8,i8,f8,f8,f8,f8,f8,f8,f8',
                               names=['VoidID','RA','Dec','z','ng','Vol','Reff','n_dens',
                                      'delta_core','delta_void','prob','R_edge'])
    data_sgc_cut=np.genfromtxt(fname_sgc,
                               dtype='i8,f8,f8,f8,i8,f8,f8,f8,f8,f8,f8,f8',
                               names=['VoidID','RA','Dec','z','ng','Vol','Reff','n_dens',
                                      'delta_core','delta_void','prob','R_edge'])
    data_cut=np.hstack((data_sgc_cut,data_ngc_cut))
    #Rotate to galactic coordinates
    th_c=(90-data_cut['Dec'])*np.pi/180.
    phi_c=data_cut['RA']*np.pi/180.
    th_g,phi_g=r(th_c,phi_c)
    b_cut=90-180*th_g/np.pi
    l_cut=180*phi_g/np.pi
    #Distance
    chi_cut=np.array([pcs.radial_comoving_distance(1./(1+z)) for z in data_cut['z']])
    #Projected size
    theta_eff_cut=data_cut['Reff']/chi_cut
    theta_edge_cut=data_cut['R_edge']/chi_cut
    nkeep=int(774*(len(theta_edge_cut)+0.)/923.446)
    ids=np.random.choice(np.arange(len(theta_edge_cut)),size=nkeep,replace=False)
    if imock%100==0 :
        print len(ids)
    tbhdu=pf.new_table([pf.Column(name='VoidID',format='K',array=data_cut['VoidID'][ids]),
                        pf.Column(name='RA'    ,format='D',array=data_cut['RA'][ids]),
                        pf.Column(name='Dec'   ,format='D',array=data_cut['Dec'][ids]),
                        pf.Column(name='L'     ,format='D',array=l_cut[ids]),
                        pf.Column(name='B'     ,format='D',array=b_cut[ids]),
                        pf.Column(name='z'     ,format='D',array=data_cut['z'][ids]),
                        pf.Column(name='ng'    ,format='K',array=data_cut['ng'][ids]),
                        pf.Column(name='Vol'   ,format='D',array=data_cut['Vol'][ids]),
                        pf.Column(name='Reff'  ,format='D',array=data_cut['Reff'][ids]),
                        pf.Column(name='n_dens',format='D',array=data_cut['n_dens'][ids]),
                        pf.Column(name='delta_core',format='D',array=data_cut['delta_core'][ids]),
                        pf.Column(name='delta_void',format='D',array=data_cut['delta_void'][ids]),
                        pf.Column(name='prob',format='D',array=data_cut['prob'][ids]),
                        pf.Column(name='R_edge',format='D',array=data_cut['R_edge'][ids]),
                        pf.Column(name='rdist',format='D',array=chi_cut[ids]),
                        pf.Column(name='theta_eff',format='D',array=theta_eff_cut[ids]),
                        pf.Column(name='theta_edge',format='D',array=theta_edge_cut[ids])])
    tbhdu.writeto(fname_end,clobber=True)
    return len(data_cut)+0.

print " Mock catalogs"
ns=[]
for i in np.arange(1000)+1 :
    if i%100==0 :
        print i,
    n=mkfits_mock(i)
    ns.append(n)
print np.mean(ns)
