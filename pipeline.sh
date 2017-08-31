#!/bin/bash

echo "Download y data"
mkdir -p data
cd data/
#Y maps
if [ ! -f data_y/milca_ymaps.fits ] ; then
    echo " Downloading Y maps"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_YSZ_R2.00.fits.tgz
    tar -xvf COM_CompMap_YSZ_R2.00.fits.tgz
    mv COM_CompMap_YSZ_R2.00.fits data_y
    rm COM_CompMap_YSZ_R2.00.fits.tgz
    cd data_y/
    rm MILCA_Csz_* nilc_weights_*
    cd ../
fi
cd data_y/
#Masks
if [ ! -f HFI_Mask_GalPlane-apo0_2048_R2.00.fits ] ; then
    echo " Downloading galatic mask" 
   wget irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/HFI_Mask_GalPlane-apo0_2048_R2.00.fits
fi
if [ ! -f LFI_Mask_PointSrc_2048_R2.00.fits ] ; then
    echo " Downloading LFI PS mask"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/LFI_Mask_PointSrc_2048_R2.00.fits
fi
if [ ! -f HFI_Mask_PointSrc_2048_R2.00.fits ] ; then
    echo " Downloading HFI PS mask"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/HFI_Mask_PointSrc_2048_R2.00.fits
fi
if [ ! -f HFI_PCCS_SZ-union_R2.08.fits ] ; then
    echo " Downloading HFI SZ catalog"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/HFI_PCCS_SZ-union_R2.08.fits.gz
    gunzip HFI_PCCS_SZ-union_R2.08.fits.gz
fi
cd ../
#545 map
if [ ! -f HFI_SkyMap_100_2048_R2.02_full.fits ] ; then
    echo " Downloading 100 map"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_100_2048_R2.02_full.fits
fi
if [ ! -f HFI_SkyMap_143_2048_R2.02_full.fits ] ; then
    echo " Downloading 143 map"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_143_2048_R2.02_full.fits
fi
if [ ! -f HFI_SkyMap_217_2048_R2.02_full.fits ] ; then
    echo " Downloading 217 map"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_217_2048_R2.02_full.fits
fi
if [ ! -f HFI_SkyMap_353_2048_R2.02_full.fits ] ; then
    echo " Downloading 353 map"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_353_2048_R2.02_full.fits
fi
if [ ! -f HFI_SkyMap_545_2048_R2.02_full.fits ] ; then
    echo " Downloading 545 map"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_545_2048_R2.02_full.fits
fi
if [ ! -f HFI_SkyMap_857_2048_R2.02_full.fits ] ; then
    echo " Downloading 857 map"
    wget irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_857_2048_R2.02_full.fits
fi
cd ../

echo "Download BOSS data"
cd data
echo " Void catalogs"
if [ ! -f voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5.dat ] ; then
    wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5.dat
fi
if [ ! -f voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5_uncut.dat ] ; then
    wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_cmass_dr12v4_ngc_Om0.3_dt0.5_uncut.dat
fi
if [ ! -f voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5.dat ] ; then
    wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5.dat
fi
if [ ! -f voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5_uncut.dat ] ; then
    wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_cmass_dr12v4_sgc_Om0.3_dt0.5_uncut.dat
fi
if [ ! -f voids_BOSS_lowz_dr12v4_ngc_Om0.3_dt0.5.dat ] ; then
    wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_lowz_dr12v4_ngc_Om0.3_dt0.5.dat
fi
if [ ! -f voids_BOSS_lowz_dr12v4_ngc_Om0.3_dt0.5_uncut.dat ] ; then
    wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_lowz_dr12v4_ngc_Om0.3_dt0.5_uncut.dat
fi
if [ ! -f voids_BOSS_lowz_dr12v4_sgc_Om0.3_dt0.5.dat ] ; then
    wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_lowz_dr12v4_sgc_Om0.3_dt0.5.dat
fi
if [ ! -f voids_BOSS_lowz_dr12v4_sgc_Om0.3_dt0.5_uncut.dat ] ; then
    wget http://lss.phy.vanderbilt.edu/voids/files/DR12/voids_BOSS_lowz_dr12v4_sgc_Om0.3_dt0.5_uncut.dat
fi
echo " Void mocks"
if [ ! -f mocks/voids_QPM_a0.6452_dr12c_0500_cmass_ngc_zspace_Om0.3_dt0.5.dat ] ; then
    wget lss.phy.vanderbilt.edu/voids/files/Mocks/QPM_mock_voids_rspace.tar.gz
    wget lss.phy.vanderbilt.edu/voids/files/Mocks/QPM_mock_voids_zspace.tar.gz
    tar -xvf QPM_mock_voids_zspace.tar.gz
    mkdir -p mocks
    mv voids_QPM_* mocks
    rm QPM_mock_voids_rspace.tar.gz QPM_mock_voids_zspace.tar.gz
fi
echo " Galaxy catalogs"
if [ ! -f galaxy_DR12v5_CMASS_North.fits ] ; then
    wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_CMASS_North.fits.gz
    gunzip galaxy_DR12v5_CMASS_North.fits.gz
fi
if [ ! -f galaxy_DR12v5_CMASS_South.fits ] ; then
    wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_CMASS_South.fits.gz
    gunzip galaxy_DR12v5_CMASS_South.fits.gz
fi
if [ ! -f galaxy_DR12v5_LOWZ_North.fits ] ; then
    wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_LOWZ_North.fits.gz
    gunzip galaxy_DR12v5_LOWZ_North.fits.gz
fi
if [ ! -f galaxy_DR12v5_LOWZ_South.fits ] ; then
    wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_LOWZ_South.fits.gz
    gunzip galaxy_DR12v5_LOWZ_South.fits.gz
fi
echo " Galaxy randoms"
if [ ! -f random0_DR12v5_CMASS_North.fits ] ; then
    wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_CMASS_North.fits.gz
    gunzip random0_DR12v5_CMASS_North.fits.gz
fi
if [ ! -f random0_DR12v5_CMASS_South.fits ] ; then
    wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_CMASS_South.fits.gz
    gunzip random0_DR12v5_CMASS_South.fits.gz
fi
cd ..

echo "Reforming data"
if [ ! -f data/data_y/mask_planck60LS.fits ] ; then
    python mkfits.py
fi

predir=pl60LS
fname_msk=data/data_y/mask_planck60LS.fits
fname_data_voids=data/voids_BOSS_cmass_dr12v4_Om0.3_dt0.5.fits
prefix_mocks=data/mocks/voids_QPM_a0.6452_dr12c
suffix_mocks=cmass_zspace_Om0.3_dt0.5.fits
echo "Compute stacks"
if [ ! -f output_${predir}_voids_milca_20/wth_mock_0632.txt ] ; then
    echo " MILCA stacks"
    fname_y=data/data_y/milca_ymaps.fits
    prefix_out=output_${predir}_voids_milca
    addqueue -q cmb -s -n 1x12 -m 1 /usr/local/shared/python/2.7.6-gcc/bin/python run_corr.py ${fname_y} ${fname_msk} ${fname_data_voids} ${prefix_mocks} ${suffix_mocks} ${prefix_out}
fi
if [ ! -f output_${predir}_voids_null_20/wth_mock_0632.txt ] ; then
    echo " NULL stacks"
    fname_y=data/data_y/y_null_milca.fits
    prefix_out=output_${predir}_voids_null
    addqueue -q cmb -s -n 1x12 -m 1 /usr/local/shared/python/2.7.6-gcc/bin/python run_corr.py ${fname_y} ${fname_msk} ${fname_data_voids} ${prefix_mocks} ${suffix_mocks} ${prefix_out}
fi
if [ ! -f output_${predir}_voids_nulln_20/wth_mock_0632.txt ] ; then
    echo " NULL NILC stacks"
    fname_y=data/data_y/y_null_nilc.fits
    prefix_out=output_${predir}_voids_nulln
    addqueue -q cmb -s -n 1x12 -m 1 /usr/local/shared/python/2.7.6-gcc/bin/python run_corr.py ${fname_y} ${fname_msk} ${fname_data_voids} ${prefix_mocks} ${suffix_mocks} ${prefix_out}
fi
for freq in 100 143 217 353 545 857
do
    if [ ! -f output_${predir}_voids_${freq}_20/wth_mock_0632.txt ] ; then
	echo " ${freq} stacks"
	fname_y=data/HFI_SkyMap_${freq}_2048_R2.02_full.fits
	prefix_out=output_${predir}_voids_${freq}
	addqueue -q cmb -s -n 1x12 -m 1 /usr/local/shared/python/2.7.6-gcc/bin/python run_corr.py ${fname_y} ${fname_msk} ${fname_data_voids} ${prefix_mocks} ${suffix_mocks} ${prefix_out}
    fi
done
if [ ! -f output_${predir}_voids_nilc_20/wth_mock_0632.txt ] ; then
    echo " NILC stacks"
    fname_y=data/data_y/nilc_ymaps.fits
    prefix_out=output_${predir}_voids_nilc
    addqueue -q cmb -s -n 1x12 -m 1 /usr/local/shared/python/2.7.6-gcc/bin/python run_corr.py ${fname_y} ${fname_msk} ${fname_data_voids} ${prefix_mocks} ${suffix_mocks} ${prefix_out}
fi

echo "Compute 3D void profile"
if [ ! -f data/profiles.txt ] ; then
    python voidprof.py
fi

echo "Compute CIB leakage"
if [ ! -f data/data_y/y_contamination.txt ] ; then
    python compute_leakage.py
fi

echo "Compute theory prediction"
if [ ! -f data/data_y/y_th_void.txt ] ; then
    python voidth.py
fi

echo "Analyze data"
python analysis.py 20
