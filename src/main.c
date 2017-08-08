#include "common.h"

int main(int argc,char **argv)
{
  char fname_field[256];
  char fname_mask[256];
  char fname_catalog[256];
  char fname_out[256];
  char field_weight[256];
  char field_cut[256];
  char field_scale[256]="none";
  double thmax,thmin=0;
  int nth,do_logbin=0,do_scale=0,type_2d=0;
  if((argc<11) || (argc>13)) {
    fprintf(stderr,"Usage: ./FieldXCorr fname_field fname_mask fname_catalog fname_out thmax nth"
	    " logbin thmin field_weight field_cut field_scale corrtype\n");
    exit(0);
  }

  timer(4);

  sprintf(fname_field  ,"%s",argv[1]);
  sprintf(fname_mask   ,"%s",argv[2]);
  sprintf(fname_catalog,"%s",argv[3]);
  sprintf(fname_out    ,"%s",argv[4]);
  thmax=atof(argv[5]);
  nth=atoi(argv[6]);
  do_logbin=atoi(argv[7]);
  thmin=atof(argv[8]);
  sprintf(field_weight,"%s",argv[9]);
  sprintf(field_cut   ,"%s",argv[10]);
  if(argc>11) {
    sprintf(field_scale,"%s",argv[11]);
    type_2d=atoi(argv[12]);
  }

  if(strcmp(field_scale,"none"))
    do_scale=1;

  printf("Reading field :\n");
  timer(0);
  long nside;
  double *field=he_read_healpix_map(fname_field,&nside,0);
  printf(" read field with Nside=  %ld\n",nside);
  timer(2);
  printf("\n");

  printf("Reading mask :\n");
  long nside_b;
  double *mask=he_read_healpix_map(fname_mask,&nside_b,0);
  if(nside_b!=nside)
    report_error(1,"Mask and field must have the same resolution\n");
  printf(" read mask with Nside=  %ld\n",nside);
  timer(2);
  printf("\n");

  long ip;
  double mean_field=0;
  for(ip=0;ip<he_nside2npix(nside);ip++) {
    field[ip]*=mask[ip];
    mean_field+=field[ip];
  }
  mean_field/=he_nside2npix(nside);
  double mean_mask=0;
  for(ip=0;ip<he_nside2npix(nside);ip++) {
    mean_mask+=mask[ip];
  }
  mean_mask/=he_nside2npix(nside);
  mean_field/=mean_mask;

  printf("Reading catalog :\n");
  long npart;
  double *pos=read_catalog(fname_catalog,field_weight,field_cut,field_scale,&npart);
  printf(" read positions for %ld objects\n",npart);
  timer(2);
  printf("\n");

  printf("Computing correlation function\n");
  double *hf_th,*hm_th;
  if(do_scale) {
    if(type_2d==0) {
      hf_th=my_calloc(nth,sizeof(double));
      hm_th=my_calloc(nth,sizeof(double));
      compute_correlation_scaled(npart,pos,nside,field,mask,thmin,thmax,nth,do_logbin,hf_th,hm_th);
    }
    else if(type_2d==1) {
      hf_th=my_calloc(16*nth,sizeof(double));
      hm_th=my_calloc(16*nth,sizeof(double));
      compute_correlation_scaled_2D(npart,pos,nside,field,mask,thmin,thmax,nth,16,do_logbin,hf_th,hm_th);
    }
    else {
      hf_th=my_calloc(4*nth*nth,sizeof(double));
      hm_th=my_calloc(4*nth*nth,sizeof(double));
      compute_correlation_scaled_2Db(npart,pos,nside,field,mask,thmin,thmax,nth,do_logbin,hf_th,hm_th);
    }
  }
  else {
    hf_th=my_calloc(nth,sizeof(double));
    hm_th=my_calloc(nth,sizeof(double));
    compute_correlation(npart,pos,nside,field,mask,thmin,thmax,nth,do_logbin,hf_th,hm_th);
  }
  timer(2);
  printf("\n");

  printf("Writing output\n");
  write_correlation(fname_out,hf_th,hm_th,thmin,thmax,nth,16,type_2d,do_logbin);
  timer(1);
  printf("\n");

  printf("All done!\n");
  free(hf_th);
  free(hm_th);
  free(pos);
  free(field);

  timer(5);

  return 0;
}
