import py_cosmo_mad as csm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from collections import Counter
from astropy.io import fits
import pyfits as pf

pcs=csm.PcsPar()
pcs.set_verbosity(1)
pcs.background_set(0.3,0.7,0.05,-1.,0.,0.7,2.7255)

print "Reading galaxies"
data_gn=(fits.open('data/galaxy_DR12v5_CMASS_North.fits'))[1].data
data_gs=(fits.open('data/galaxy_DR12v5_CMASS_South.fits'))[1].data
data_g=np.hstack((data_gn,data_gs))
rgals=np.array([pcs.radial_comoving_distance(1./(1+z)) for z in data_g['Z']])
phigals=data_g['RA']*np.pi/180
thgals=(90-data_g['DEC'])*np.pi/180
xgals=rgals*np.sin(thgals)*np.cos(phigals)
ygals=rgals*np.sin(thgals)*np.sin(phigals)
zgals=rgals*np.cos(thgals)
ngals=len(zgals)

print "Read randoms"
data_r=(fits.open('data/random0_DR12v5_CMASS_Subsample.fits'))[1].data
rrans=np.array([pcs.radial_comoving_distance(1./(1+z)) for z in data_r['Z']])
phirans=data_r['RA']*np.pi/180
thrans=(90-data_r['DEC'])*np.pi/180
xrans=rrans*np.sin(thrans)*np.cos(phirans)
yrans=rrans*np.sin(thrans)*np.sin(phirans)
zrans=rrans*np.cos(thrans)
nrans=len(zrans)

print "Reading voids"
data_ngc_cut=np.genfromtxt('data/voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5.dat',
                           dtype='i8,f8,f8,f8,i8,f8,f8,f8,f8,f8,f8,f8',
                           names=['VoidID','RA','Dec','z','ng','Vol','Reff','n_dens',
                                  'delta_core','delta_void','prob','R_edge'])
data_sgc_cut=np.genfromtxt('data/voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5.dat',
                           dtype='i8,f8,f8,f8,i8,f8,f8,f8,f8,f8,f8,f8',
                           names=['VoidID','RA','Dec','z','ng','Vol','Reff','n_dens',
                                  'delta_core','delta_void','prob','R_edge'])
data_cut=np.hstack((data_ngc_cut,data_sgc_cut))
rvoids=np.array([pcs.radial_comoving_distance(1./(1+z)) for z in data_cut['z']])
phivoids=data_cut['RA']*np.pi/180
thvoids=(90-data_cut['Dec'])*np.pi/180
xvoids=rvoids*np.sin(thvoids)*np.cos(phivoids)
yvoids=rvoids*np.sin(thvoids)*np.sin(phivoids)
zvoids=rvoids*np.cos(thvoids)
nvoids=len(rvoids)

cols=['r','g','b','y','m','c']
nb=30
ngh=np.zeros(nb)
nrh=np.zeros(nb)
r0h=np.zeros(nb)
rfh=np.zeros(nb)
for i in np.arange(nvoids) :
    print i,nvoids
    xarr=(xgals-xvoids[i])/data_cut['Reff'][i]
    yarr=(ygals-yvoids[i])/data_cut['Reff'][i]
    zarr=(zgals-zvoids[i])/data_cut['Reff'][i]
    rarr=np.sqrt(xarr**2+yarr**2+zarr**2)
    h,b=np.histogram(rarr,bins=nb,range=[0,3])
    ngh+=h
    r0h+=b[:-1]
    rfh+=b[1:]
    xarr=(xrans-xvoids[i])/data_cut['Reff'][i]
    yarr=(yrans-yvoids[i])/data_cut['Reff'][i]
    zarr=(zrans-zvoids[i])/data_cut['Reff'][i]
    rarr=np.sqrt(xarr**2+yarr**2+zarr**2)
    h,b=np.histogram(rarr,bins=nb,range=[0,3])
    nrh+=h
ngh=(ngh+0.)/nvoids
nrh=(nrh+0.)/nvoids*(ngals+0.)/(nrans+0.)
r0h=(r0h+0.)/nvoids
rfh=(rfh+0.)/nvoids
rmh=0.5*(r0h+rfh)
densg=ngh/(4*np.pi*(rfh**3-r0h**3)/3)
densr=nrh/(4*np.pi*(rfh**3-r0h**3)/3)
plt.plot(rmh,densg/densr)
plt.show()

np.savetxt("data/profiles.txt",np.transpose([r0h,rmh,rfh,densg,densr,densg/densr]))
