import numpy as np
import pyfits as pf
import py_cosmo_mad as csm
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits

print "Subsampling randoms"
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

print "Cosmological parameters"
pcs=csm.PcsPar()
pcs.background_set(0.30,0.70,0.05,-1.,0.,0.7,2.7255)
r=hp.Rotator(coord=['C','G'])

#Make mask
print "Mask"
mask_gal_80=hp.read_map("data/data_y/HFI_Mask_GalPlane-apo0_2048_R2.00.fits",verbose=False,field=0)
mask_gal_60=hp.read_map("data/data_y/HFI_Mask_GalPlane-apo0_2048_R2.00.fits",verbose=False,field=1)
mask_p0=hp.read_map("data/data_y/LFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,hdu=1); print np.mean(mask_p0)
mask_p1=hp.read_map("data/data_y/LFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,hdu=2); print np.mean(mask_p1)
mask_p2=hp.read_map("data/data_y/LFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,hdu=3); print np.mean(mask_p2)
mask_pl=mask_p0*mask_p1*mask_p2
mask_p0=hp.read_map("data/data_y/HFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,field=0); print np.mean(mask_p0)
mask_p1=hp.read_map("data/data_y/HFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,field=1); print np.mean(mask_p1)
mask_p2=hp.read_map("data/data_y/HFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,field=2); print np.mean(mask_p2)
mask_p3=hp.read_map("data/data_y/HFI_Mask_PointSrc_2048_R2.00.fits",verbose=False,field=3); print np.mean(mask_p3)
mask_ph=mask_p0*mask_p1*mask_p2*mask_p3
mask_b=hp.read_map("data/data_y/mask.fits",verbose=False);
hp.write_map("data/data_y/mask_planck60.fits",mask_gal_60*mask_ph)
hp.write_map("data/data_y/mask_planck80.fits",mask_gal_80*mask_ph)

#Make null map
print "Null map"
y_first,y_last=hp.read_map("data/data_y/milca_ymaps.fits",field=[1,2],verbose=False)
y_null=0.5*(y_first-y_last)
hp.write_map("data/data_y/y_null_milca.fits",y_null)

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

print "Data catalogs"
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

print "Mock catalogs"
ns=[]
for i in np.arange(1000)+1 :
    print i,
    n=mkfits_mock(i)
    ns.append(n)
print np.mean(ns)
