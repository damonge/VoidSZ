#ifndef _COMMON_
#define _COMMON_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <omp.h>
#include <fitsio.h>
#include <chealpix.h>
/*
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <fftw3.h>
#include <sharp_almhelpers.h>
#include <sharp_geomhelpers.h>
#include <sharp.h>
*/

#define DTOR 0.017453292519943295 // x deg = x*DTORAD rad
#define MIN(a, b) (((a) < (b)) ? (a) : (b)) //Minimum of two numbers

void timer(int i);
void report_error(int level,char *fmt,...);
void *my_malloc(size_t size);
void *my_calloc(size_t nmemb,size_t size);
FILE *my_fopen(const char *path,const char *mode);
double *read_catalog(char *fname_cat,char *weight_name,char *cut_name,char *scale_name,long *ngal);
void compute_correlation(long ngal,double *pos,long nside,double *fld,double *msk,
			 double thmin,double thmax,int nth,int do_log,double *hf_th,double *hm_th);
void compute_correlation_scaled(long ngal,double *pos,long nside,double *fld,double *msk,
				double thmin,double thmax,int nth,int do_log,double *hf_th,double *hm_th);
void compute_correlation_scaled_2D(long ngal,double *pos,long nside,double *fld,double *msk,
				   double thmin,double thmax,int nth,int na,
				   int do_log,double *hf_th,double *hm_th);
void compute_correlation_scaled_2Db(long ngal,double *pos,long nside,double *fld,double *msk,
				    double thmin,double thmax,int nth,
				    int do_log,double *hf_th,double *hm_th);
void write_correlation(char *fname_out,double *hf_th,double *hm_th,
		       double thmin,double thmax,int nth,int na,int do_2d,int do_log);

long he_nside2npix(long nside);
void he_pix2vec_ring(long nside, long ipix, double *vec);
long he_ang2pix(long nside,double cth,double phi);
void he_write_healpix_map(double **tmap,int nfields,long nside,char *fname);
void he_get_file_params(char *fname,long *nside,int *nfields,int *isnest);
double *he_read_healpix_map(char *fname,long *nside,int nfield);
int he_ring_num(long nside,double z);
void he_query_strip(long nside,double theta1,double theta2,int *pixlist,long *npix_strip);
void he_ring2nest_inplace(double *map_in,long nside);
void he_nest2ring_inplace(double *map_in,long nside);
void he_in_ring(int nside,int iz,double phi0,double dphi,int *listir,int *nir);
void he_query_disc(int nside,double cth0,double phi,double radius,int *listtot,int *nlist,int inclusive);
void he_udgrade(double *map_in,long nside_in,double *map_out,long nside_out,int nest);
/*
long he_nalms(int lmax);
long he_indexlm(int l,int m,int lmax);
void he_alm2map(int nside,int lmax,int ntrans,int spin,double **maps,complex **alms);
void he_map2alm(int nside,int lmax,int ntrans,int spin,double **maps,complex **alms,int niter);
void he_alm2cl(complex **alms_1,complex **alms_2,int pol_1,int pol_2,double **cls,int lmax);
void he_anafast(double **maps_1,double **maps_2,int pol_1,int pol_2,double **cls,
		int nside,int lmax,int iter);
*/
/*
double *he_generate_beam_window(int lmax,double fwhm_amin);
void he_zero_alm(int lmax,complex *alm);
void he_alter_alm(int lmax,double fwhm_amin,complex *alm_in,complex *alm_out,double *window,int add_to_out);
void he_map_product(int nside,double *mp1,double *mp2,double *mp_out);
double he_map_dot(int nside,double *mp1,double *mp2);
complex **he_synalm(int nside,int nmaps,int lmax,double **cells,double **beam,int seed);
#ifdef _WITH_NEEDLET
#define HE_NBAND_NX 512
#define HE_NORM_FT 2.2522836206907617
#define HE_NL_INTPREC 1E-6
#define HE_NT_NSIDE_MIN 32
typedef struct {
  double b;
  double inv_b;
  gsl_spline *b_spline;
  gsl_interp_accel *b_intacc;
  int niter;
  int nside0;
  int jmax_min;
  int nj;
  int *nside_arr;
  int *lmax_arr;
  double **b_arr;
} he_needlet_params;
void he_nt_end(he_needlet_params *par);
he_needlet_params *he_nt_init(double b_nt,int nside0,int niter);
void he_free_needlet(he_needlet_params *par,int pol,double ***nt);
double ***he_alloc_needlet(he_needlet_params *par,int pol);
void he_nt_get_window(he_needlet_params *par,int j,double *b);
complex **he_needlet2map(he_needlet_params *par,double **map,double ***nt,
			  int return_alm,int pol,int input_TEB,int output_TEB);
complex **he_map2needlet(he_needlet_params *par,double **map,double ***nt,
			  int return_alm,int pol,int input_TEB,int output_TEB);

*/
#endif //_COMMON_
