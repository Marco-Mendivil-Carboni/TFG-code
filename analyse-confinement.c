//Libraries

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Alex_pol.h"
#include <glob.h>
#include <string.h>

//Defines

#define M_PI 3.14159265358979323846

#define hnb 128

//Structs

struct sim_params
{
  float T;
  int N;
  float R;
};

struct rand_var_props
{
  double m_1;
  double m_2;
  double var;
  double sem;
  int n_term;
  int term_end;
  int n_blocks;
};

typedef struct sim_params sim_params;
typedef struct rand_var_props rand_var_props;

//Functions

void read_parameters( sim_params *sp, FILE *f)
{
  if( fscanf(f,"T\t%f\n",&(sp->T))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"N\t%d\n",&(sp->N))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"R\t%f\n",&(sp->R))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( (sp->T)<__FLT_MIN__){ printf("T must be positive.\n"); exit(-1);}
  if( (sp->N)<__FLT_MIN__){ printf("N must be positive.\n"); exit(-1);}
  if( (sp->R)<__FLT_MIN__){ printf("R must be positive.\n"); exit(-1);}
}

int read_trajectory_positions( int N, float *r, float *t, FILE *f)
{
  //header
  int magickvalue;
  if( fread(&magickvalue,sizeof(int),1,f)==0){ return 0;}
  if( magickvalue!=1993){ printf("Wrong format.\n"); exit(-1);}
  char trrversion[13];
  int len_s_a;
  int len_s_b;
  fread(&len_s_a,sizeof(int),1,f);
  fread(&len_s_b,sizeof(int),1,f);
  fread(trrversion,sizeof(char),sizeof(trrversion)-1,f);
  int zero;
  for( int i = 0; i<7; i++)
  {
    fread(&zero,sizeof(int),1,f);
  }
  int x_size;
  fread(&x_size,sizeof(int),1,f);
  if( x_size!=3*N*sizeof(float)){ printf("Error reading trajectory.\n"); exit(-1);}
  int v_size;
  fread(&v_size,sizeof(int),1,f);
  if( v_size!=0){ printf("Error reading trajectory.\n"); exit(-1);}
  int f_size;
  fread(&f_size,sizeof(int),1,f);
  if( f_size!=0){ printf("Error reading trajectory.\n"); exit(-1);}
  int natoms;
  fread(&natoms,sizeof(int),1,f);
  if( natoms!=N){ printf("Error reading trajectory.\n"); exit(-1);}
  int step;
  fread(&step,sizeof(int),1,f);
  float time;
  fread(&zero,sizeof(int),1,f);
  fread(&time,sizeof(float),1,f);
  fread(&zero,sizeof(int),1,f);
  *t = time;
  //coordinates
  fread(r,sizeof(float),3*N,f);
  return 1;
}

void calc_cm( int N, float *r, double *r_cm)
{
  r_cm[0] = r_cm[1] = r_cm[2] = 0.0;
  for( int i_p = 0; i_p<N; i_p++)
  {
    r_cm[0] += r[3*i_p+0];
    r_cm[1] += r[3*i_p+1];
    r_cm[2] += r[3*i_p+2];
  }
  r_cm[0] /= N;
  r_cm[1] /= N;
  r_cm[2] /= N;
}

void calc_Rg2( int N, float *r, double *r_cm, double *Rg2)
{
  *Rg2 = 0.0;
  for( int i_p = 0; i_p<N; i_p++)
  {
    *Rg2 += (r[3*i_p+0]-r_cm[0])*(r[3*i_p+0]-r_cm[0]);
    *Rg2 += (r[3*i_p+1]-r_cm[1])*(r[3*i_p+1]-r_cm[1]);
    *Rg2 += (r[3*i_p+2]-r_cm[2])*(r[3*i_p+2]-r_cm[2]);
  }
  *Rg2 = *Rg2/N;
}

void block_data( int *n_data, double *x)
{
  int i;
  for( i = 0; (2*i)<(*n_data-1); i++)
  {
    x[i] = 0.5*(x[2*i]+x[2*i+1]);
  }
  *n_data = i;
}

void calc_moments( int n_data, double *x, double *m_1, double *m_2)
{
  *m_1 = 0.0;
  *m_2 = 0.0;
  for( int i = 0; i<n_data; i++)
  {
    *m_1 += x[i];
    *m_2 += x[i]*x[i];
  }
  *m_1 /= n_data;
  *m_2 /= n_data;
}

int calc_sig_mean( int n_data, double *x, double *sem)
{
  int n_data_b = n_data;
  double m_1, m_2, var_m_1, uplim_var_m_1;
  calc_moments(n_data_b,x,&m_1,&m_2);
  var_m_1 = (m_2-m_1*m_1)/(n_data_b-1.0);
  uplim_var_m_1 = var_m_1*(1.0+sqrt(2.0/(n_data_b-1.0)));
  while( n_data_b>3)
  {
    block_data(&n_data_b,x);
    calc_moments(n_data_b,x,&m_1,&m_2);
    var_m_1 = (m_2-m_1*m_1)/(n_data_b-1.0);
    if( var_m_1>uplim_var_m_1)
    {
      uplim_var_m_1 = var_m_1*(1.0+sqrt(2.0/(n_data_b-1.0)));
    }
    else
    {
      *sem = sqrt(var_m_1);
      return n_data_b;
    }
  }
  *sem = sqrt(var_m_1);
  return n_data_b; 
}

int estimate_term( int n_data, double *x, int *opt_n_term)
{
  int i_term;
  int n_term;
  double m_1, m_2, smer;
  double smer_min = INFINITY;
  for( i_term = 0; i_term<50; i_term++)
  {
    n_term = i_term*n_data/100;
    calc_moments(n_data-n_term,&x[n_term],&m_1,&m_2);
    smer = (m_2-m_1*m_1)/(n_data-n_term);
    if( smer<smer_min)
    {
      *opt_n_term = n_term;
      smer_min = smer;
    }
  }
  if( *opt_n_term==n_term)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

void calc_statistics( int n_data, double *x, rand_var_props *x_p)
{
  (*x_p).m_1 = 0.0;
  (*x_p).m_2 = 0.0;
  for( int i = 0; i<n_data; i++)
  {
    (*x_p).m_1 += x[i];
    (*x_p).m_2 += x[i]*x[i];
  }
  (*x_p).m_1 /= n_data;
  (*x_p).m_2 /= n_data;
  (*x_p).var = ((*x_p).m_2-(*x_p).m_1*(*x_p).m_1);
  (*x_p).var *= n_data/(n_data-1.0);
  (*x_p).sem = sqrt((*x_p).var/n_data);
}

void update_rpdfcs( int N, float R, float *r, double *rpdfcs)
{
  int i_h;
  float d_r;
  for( int i_p = 0; i_p<N; i_p++)
  {
    d_r = 0.0;
    for( int i_c = 0; i_c<3; i_c++)
    {
      d_r += r[3*i_p+i_c]*r[3*i_p+i_c];
    }
    d_r = sqrt(d_r);
    i_h = hnb*d_r/R;
    if(i_h>=hnb){ printf("Error calculating histogram.\n"); exit(-1);}
    rpdfcs[i_h] += hnb/(N*R);
  }
}

void update_rpdfcm( int N, float R, float *r, double *rpdfcm)
{
  int i_h;
  float d_r;
  double r_cm[3];
  calc_cm(N,r,r_cm);
  for( int i_p = 0; i_p<N; i_p++)
  {
    d_r = 0.0;
    for( int i_c = 0; i_c<3; i_c++)
    {
      d_r += (r[3*i_p+i_c]-r_cm[i_c])*(r[3*i_p+i_c]-r_cm[i_c]);
    }
    d_r = sqrt(d_r);
    i_h = hnb*d_r/R;
    if(i_h>=2*hnb){ printf("Error calculating histogram.\n"); exit(-1);}
    rpdfcm[i_h] += hnb/(N*R);
  }
}

void update_cm_rpd( int N, float R, float *r, double *cm_rpd)
{
  int i_h;
  float d_r;
  double r_cm[3];
  calc_cm(N,r,r_cm);
  d_r = 0.0;
  for( int i_c = 0; i_c<3; i_c++)
  {
    d_r += r_cm[i_c]*r_cm[i_c];
  }
  d_r = sqrt(d_r);
  i_h = hnb*d_r/R;
  if(i_h>=hnb){ printf("Error calculating histogram.\n"); exit(-1);}
  cm_rpd[i_h] += hnb/R;
}

void update_oop( int N, float *r, double *oop)
{
  float b[3];
  double cos;
  double d_r, d_b;
  for( int i_p = 0; i_p<(N-1); i_p++)
  {
    d_r = d_b = cos = 0.0;
    for( int i_c = 0; i_c<3; i_c++)
    {
      b[i_c] = r[3*(i_p+1)+i_c]-r[3*i_p+i_c];
      d_r += r[3*i_p+i_c]*r[3*i_p+i_c];
      d_b += b[i_c]*b[i_c];
      cos += b[i_c]*r[3*i_p+i_c];
    }
    d_r = sqrt(d_r);
    d_b = sqrt(d_b);
    cos /= d_b*d_r;
    *oop += 0.5*(3.0*cos*cos-1.0)/(N-1.0);
  }
}

int main( int argc, char const *argv[])
{
  if( argc!=2)
  {
    if( argc<2){ printf("You forgot the input.\n");}
    else{ printf("Too many arguments.\n");}
    exit(-1);
  }
  if( sizeof(argv[1])>128){ printf("Directory name too long.\n"); exit(-1);}

  char sim_dir[128];
  snprintf(sim_dir,sizeof(sim_dir),"%s",argv[1]);

  FILE *file_i1;
  FILE *file_o1;
  FILE *file_o2;

  char filename[256];

  //Simulation parameters and variables

  sim_params sp;

  snprintf(filename,sizeof(filename),"%s/adjustable-parameters.dat",sim_dir);
  file_i1 = fopen(filename,"rt");
  if( file_i1==NULL){ printf("Error opening parameters file.\n"); exit(-1);}
  read_parameters(&sp,file_i1);
  fclose(file_i1);

  float *r;
  r = (float*)malloc(3*sp.N*sizeof(float));

  float t;

  //Individual analysis

  int sim_num = 0;

  glob_t prev_sim_files;
  snprintf(filename,sizeof(filename),"%s/trajectory-positions-*",sim_dir);
  if( glob(filename,0,NULL,&prev_sim_files)==0)
  {
    sim_num = prev_sim_files.gl_pathc;
  }
  globfree(&prev_sim_files);

  double **rpdfcs = (double**)malloc(sim_num*sizeof(double*));
  double **rpdfcm = (double**)malloc(sim_num*sizeof(double*));
  double **cm_rpd = (double**)malloc(sim_num*sizeof(double*));

  double *oop = (double*)malloc(sim_num*sizeof(double));

  snprintf(filename,sizeof(filename),"%s/analysis-statistics.dat",sim_dir);
  file_o2 = fopen(filename,"wt");
  if( file_o2==NULL){ printf("Error opening file.\n"); exit(-1);}

  fprintf(file_o2,"#%14s%15s%15s","sim_idx","n_term","n_blocks");
  fprintf(file_o2,"%15s%15s%15s%15s\n","<Rg2>","SEM(Rg2)","Var(Rg2)","knot");

  for( int sim_idx = 0; sim_idx<sim_num; sim_idx++)
  {
    //Simulation output reading

    double r_cm[3];
    double Rg2;
    int n_data = 0;

    snprintf(filename,sizeof(filename),"%s/Rg2-time-series-%03d.dat",sim_dir,sim_idx);
    file_o1 = fopen(filename,"wt");
    if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}

    snprintf(filename,sizeof(filename),"%s/trajectory-positions-%03d.trr",sim_dir,sim_idx);
    file_i1 = fopen(filename,"rb");
    if( file_i1==NULL){ printf("Error opening file.\n"); exit(-1);}

    while( read_trajectory_positions(sp.N,r,&t,file_i1)) 
    {
      calc_cm(sp.N,r,r_cm);
      calc_Rg2(sp.N,r,r_cm,&Rg2);
      fprintf(file_o1,"%15.6f%15.6lf\n",t,Rg2);
      n_data++;
    }

    if( n_data==0){ printf("No data.\n"); exit(-1);}

    fclose(file_i1);

    fclose(file_o1);

    //Time series reading

    double *Rg2_ts;
    Rg2_ts = (double*)malloc(n_data*sizeof(double));

    snprintf(filename,sizeof(filename),"%s/Rg2-time-series-%03d.dat",sim_dir,sim_idx);
    file_i1 = fopen(filename,"rt");
    if( file_i1==NULL){ printf("Error opening file.\n"); exit(-1);}

    for( int i_m = 0; i_m<n_data; i_m++)
    {
      fscanf(file_i1,"%15f%15lf\n",&t,&Rg2_ts[i_m]);
    }

    fclose(file_i1);

    //Analysis of the gyration radius

    rand_var_props Rg2_p;

    Rg2_p.term_end = estimate_term(n_data,Rg2_ts,&Rg2_p.n_term);

    calc_moments(n_data-Rg2_p.n_term,&Rg2_ts[Rg2_p.n_term],&Rg2_p.m_1,&Rg2_p.m_2);
    Rg2_p.var = (Rg2_p.m_2-Rg2_p.m_1*Rg2_p.m_1);
    Rg2_p.var *= (n_data-Rg2_p.n_term)/(n_data-Rg2_p.n_term-1.0);

    Rg2_p.n_blocks = calc_sig_mean(n_data-Rg2_p.n_term,&Rg2_ts[Rg2_p.n_term],&Rg2_p.sem);

    fprintf(file_o2,"%15d%15d%15d",sim_idx,Rg2_p.n_term,Rg2_p.n_blocks);
    fprintf(file_o2,"%15.6lf%15.6lf%15.6lf",Rg2_p.m_1,Rg2_p.sem,Rg2_p.var);

    free(Rg2_ts);

    //Polymer distribution analysis

    rpdfcs[sim_idx] = (double*)malloc(hnb*sizeof(double));
    rpdfcm[sim_idx] = (double*)malloc(2*hnb*sizeof(double));
    cm_rpd[sim_idx] = (double*)malloc(hnb*sizeof(double));

    for( int i_h = 0; i_h<hnb; i_h++)
    {
      rpdfcs[sim_idx][i_h] = 0.0;
    }
    for( int i_h = 0; i_h<2*hnb; i_h++)
    {
      rpdfcm[sim_idx][i_h] = 0.0;
    }
    for( int i_h = 0; i_h<hnb; i_h++)
    {
      cm_rpd[sim_idx][i_h] = 0.0;
    }

    oop[sim_idx] = 0.0;

    snprintf(filename,sizeof(filename),"%s/trajectory-positions-%03d.trr",sim_dir,sim_idx);
    file_i1 = fopen(filename,"rb");
    if( file_i1==NULL){ printf("Error opening file.\n"); exit(-1);}

    for( int i_m = 0; i_m<Rg2_p.n_term; i_m++)
    {
      read_trajectory_positions(sp.N,r,&t,file_i1);
    }
    for( int i_m = Rg2_p.n_term; i_m<n_data; i_m++)
    {
      read_trajectory_positions(sp.N,r,&t,file_i1);
      update_rpdfcs(sp.N,sp.R,r,rpdfcs[sim_idx]);
      update_rpdfcm(sp.N,sp.R,r,rpdfcm[sim_idx]);
      update_cm_rpd(sp.N,sp.R,r,cm_rpd[sim_idx]);
      update_oop(sp.N,r,&oop[sim_idx]);
    }

    for( int i_h = 0; i_h<hnb; i_h++)
    {
      rpdfcs[sim_idx][i_h] /= n_data-Rg2_p.n_term;
    }
    for( int i_h = 0; i_h<2*hnb; i_h++)
    {
      rpdfcm[sim_idx][i_h] /= n_data-Rg2_p.n_term;
    }
    for( int i_h = 0; i_h<hnb; i_h++)
    {
      cm_rpd[sim_idx][i_h] /= n_data-Rg2_p.n_term;
    }

    oop[sim_idx] /= n_data-Rg2_p.n_term;

    char knot_name[128];
    calc_invariant(sp.N,r,knot_name);
    fprintf(file_o2,"%15s\n",knot_name);

    fclose(file_i1);
  }

  fclose(file_o2);

  free(r);

  //Global analysis

  double **rpdfcs_t = (double**)malloc(hnb*sizeof(double*));
  for( int i_h = 0; i_h<hnb; i_h++)
  {
    rpdfcs_t[i_h] = (double*)malloc(sim_num*sizeof(double));
  }

  for( int i_h = 0; i_h<hnb; i_h++)
  {
    for( int sim_idx = 0; sim_idx<sim_num; sim_idx++)
    {
      rpdfcs_t[i_h][sim_idx] = rpdfcs[sim_idx][i_h];
    }
  }

  rand_var_props *rpdfcs_t_p = (rand_var_props*)malloc(hnb*sizeof(rand_var_props));

  for( int i_h = 0; i_h<hnb; i_h++)
  {
    calc_statistics(sim_num,rpdfcs_t[i_h],&rpdfcs_t_p[i_h]);
  }

  double **rpdfcm_t = (double**)malloc(2*hnb*sizeof(double*));
  for( int i_h = 0; i_h<2*hnb; i_h++)
  {
    rpdfcm_t[i_h] = (double*)malloc(sim_num*sizeof(double));
  }

  for( int i_h = 0; i_h<2*hnb; i_h++)
  {
    for( int sim_idx = 0; sim_idx<sim_num; sim_idx++)
    {
      rpdfcm_t[i_h][sim_idx] = rpdfcm[sim_idx][i_h];
    }
  }

  rand_var_props *rpdfcm_t_p = (rand_var_props*)malloc(2*hnb*sizeof(rand_var_props));

  for( int i_h = 0; i_h<2*hnb; i_h++)
  {
    calc_statistics(sim_num,rpdfcm_t[i_h],&rpdfcm_t_p[i_h]);
  }

  double **cm_rpd_t = (double**)malloc(hnb*sizeof(double*));
  for( int i_h = 0; i_h<hnb; i_h++)
  {
    cm_rpd_t[i_h] = (double*)malloc(sim_num*sizeof(double));
  }

  for( int i_h = 0; i_h<hnb; i_h++)
  {
    for( int sim_idx = 0; sim_idx<sim_num; sim_idx++)
    {
      cm_rpd_t[i_h][sim_idx] = cm_rpd[sim_idx][i_h];
    }
  }

  rand_var_props *cm_rpd_t_p = (rand_var_props*)malloc(hnb*sizeof(rand_var_props));

  for( int i_h = 0; i_h<hnb; i_h++)
  {
    calc_statistics(sim_num,cm_rpd_t[i_h],&cm_rpd_t_p[i_h]);
  }

  rand_var_props oop_p;

  calc_statistics(sim_num,oop,&oop_p);

  //Results

  snprintf(filename,sizeof(filename),"%s/analysis-rpd-histogram.dat",sim_dir);
  file_o1 = fopen(filename,"wt");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}

  for( int i_h = 0; i_h<hnb; i_h++)
  {
    fprintf(file_o1,"%15.6lf%15.6lf%15.6lf\n",(i_h+0.5)*(sp.R/hnb),rpdfcs_t_p[i_h].m_1,rpdfcs_t_p[i_h].sem);
  }

  fprintf(file_o1,"\n\n");

  for( int i_h = 0; i_h<2*hnb; i_h++)
  {
    fprintf(file_o1,"%15.6lf%15.6lf%15.6lf\n",(i_h+0.5)*(sp.R/hnb),rpdfcm_t_p[i_h].m_1,rpdfcm_t_p[i_h].sem);
  }

  fprintf(file_o1,"\n\n");

  for( int i_h = 0; i_h<hnb; i_h++)
  {
    fprintf(file_o1,"%15.6lf%15.6lf%15.6lf\n",(i_h+0.5)*(sp.R/hnb),cm_rpd_t_p[i_h].m_1,cm_rpd_t_p[i_h].sem);
  }

  fclose(file_o1);

  snprintf(filename,sizeof(filename),"%s/analysis-order-parameter.dat",sim_dir);
  file_o1 = fopen(filename,"wt");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}

  fprintf(file_o1,"%15.6lf%15.6lf\n",oop_p.m_1,oop_p.sem);

  fclose(file_o1);

  printf("Analysed %d simulations from %s.\n",sim_num,sim_dir);

  return 0;
}
