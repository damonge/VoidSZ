#!/bin/bash

echo "Download y data"
#mkdir -p data
#cd data/
#wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_YSZ_R2.00.fits.tgz
#tar -xvf COM_CompMap_YSZ_R2.00.fits.tgz
#mv COM_CompMap_YSZ_R2.00.fits data_y
#rm COM_CompMap_YSZ_R2.00.fits.tgz
#cd data_y/
#rm MILCA_Csz_* nilc_weights_*
#wget irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/HFI_Mask_GalPlane-apo0_2048_R2.00.fits
#wget irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/LFI_Mask_PointSrc_2048_R2.00.fits
#wget irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/HFI_Mask_PointSrc_2048_R2.00.fits
#wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/lensing/COM_CompMap_Lensing_2048_R2.00.tar
#tar -xvf COM_CompMap_Lensing_2048_R2.00.tar
#cd data
#gunzip mask.fits.gz
#mv mask.fits ..
#cd ..
#rm -r data
#rm COM_CompMap_Lensing_2048_R2.00.tar
#cd ../
#wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_545_2048_R2.02_full.fits
#cd ../

echo "Download BOSS data"
#cd data
#wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5.dat
#wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5_uncut.dat
#wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5.dat
#wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5_uncut.dat
#wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_lowz_dr12v4_ngc_Om0.3_dt0.5.dat
#wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_lowz_dr12v4_ngc_Om0.3_dt0.5_uncut.dat
#wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_lowz_dr12v4_sgc_Om0.3_dt0.5.dat
#wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_lowz_dr12v4_sgc_Om0.3_dt0.5_uncut.dat
#wget lss.phy.vanderbilt.edu/voids/files/Mocks/QPM_mock_voids_rspace.tar.gz
#wget lss.phy.vanderbilt.edu/voids/files/Mocks/QPM_mock_voids_zspace.tar.gz
#tar -xvf QPM_mock_voids_zspace.tar.gz
#mkdir -p mocks
#mv voids_QPM_* mocks
#rm QPM_mock_voids_rspace.tar.gz QPM_mock_voids_zspace.tar.gz
#wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_CMASS_North.fits.gz
#wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_CMASS_South.fits.gz
#wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_LOWZ_North.fits.gz
#wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_LOWZ_South.fits.gz
#gunzip galaxy_DR12v5_CMASS_North.fits.gz
#gunzip galaxy_DR12v5_CMASS_South.fits.gz
#gunzip galaxy_DR12v5_LOWZ_North.fits.gz
#gunzip galaxy_DR12v5_LOWZ_South.fits.gz
#wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_CMASS_North.fits.gz
#wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_CMASS_South.fits.gz
#gunzip random0_DR12v5_CMASS_North.fits.gz
#gunzip random0_DR12v5_CMASS_South.fits.gz
#wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_LOWZ_North.fits.gz
#wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_LOWZ_South.fits.gz
#cd ..

echo "Reforming data"
#cd data
#python mkfits.py
#cd ../

echo "Compute y stacks"
#fname_y=data/data_y/milca_ymaps.fits
#fname_msk=data/data_y/mask.fits
#fname_data_voids=data/voids_BOSS_cmass_dr12v4_Om0.3_dt0.5.fits
#prefix_mocks=data/mocks/voids_QPM_a0.6452_dr12c
#suffix_mocks=cmass_zspace_Om0.3_dt0.5.fits
#prefix_out=output_voids_milca
#addqueue -q cmb -s -n 1x12 -m 1 /usr/local/shared/python/2.7.6-gcc/bin/python run_corr.py ${fname_y} ${fname_msk} ${fname_data_voids} ${prefix_mocks} ${suffix_mocks} ${prefix_out}

echo "Compute null stacks"
#fname_y=data/data_y/y_null_milca.fits
#fname_msk=data/data_y/mask.fits
#fname_data_voids=data/voids_BOSS_cmass_dr12v4_Om0.3_dt0.5.fits
#prefix_mocks=data/mocks/voids_QPM_a0.6452_dr12c
#suffix_mocks=cmass_zspace_Om0.3_dt0.5.fits
#prefix_out=output_voids_null
#addqueue -q cmb -s -n 1x12 -m 1 /usr/local/shared/python/2.7.6-gcc/bin/python run_corr.py ${fname_y} ${fname_msk} ${fname_data_voids} ${prefix_mocks} ${suffix_mocks} ${prefix_out}

echo "Compute 545 stacks"
#fname_y=data/HFI_SkyMap_545_2048_R2.02_full.fits
#fname_msk=data/data_y/mask.fits
#fname_data_voids=data/voids_BOSS_cmass_dr12v4_Om0.3_dt0.5.fits
#prefix_mocks=data/mocks/voids_QPM_a0.6452_dr12c
#suffix_mocks=cmass_zspace_Om0.3_dt0.5.fits
#prefix_out=output_voids_545
#addqueue -q cmb -s -n 1x12 -m 1 /usr/local/shared/python/2.7.6-gcc/bin/python run_corr.py ${fname_y} ${fname_msk} ${fname_data_voids} ${prefix_mocks} ${suffix_mocks} ${prefix_out}

echo "Compute NILC stacks"
#fname_y=data/data_y/nilc_ymaps.fits
#fname_msk=data/data_y/mask.fits
#fname_data_voids=data/voids_BOSS_cmass_dr12v4_Om0.3_dt0.5.fits
#prefix_mocks=data/mocks/voids_QPM_a0.6452_dr12c
#suffix_mocks=cmass_zspace_Om0.3_dt0.5.fits
#prefix_out=output_voids_nilc
#addqueue -q cmb -s -n 1x12 -m 1 /usr/local/shared/python/2.7.6-gcc/bin/python run_corr.py ${fname_y} ${fname_msk} ${fname_data_voids} ${prefix_mocks} ${suffix_mocks} ${prefix_out}

echo "Compute 3D void profile"
python voidprof.py

<<COMMENT
echo "Compute theoretical y stack"
python voidth.py
