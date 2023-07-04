//Libraries

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Alex_pol.h"
#include <glob.h>
#include <string.h>

//Defines

#define M_PI 3.14159265358979323846

#define n_b 2
#define bpm 8

//Structs

struct sim_params
{
  float T;
  int N;
  float V;
};

struct rand_var_props
{
  double m_1;
  double m_2;
  double var;
  double sem;
};

typedef struct sim_params sim_params;
typedef struct rand_var_props rand_var_props;

//Functions

void read_parameters( sim_params *sp, FILE *f)
{
  if( fscanf(f,"T\t%f\n",&(sp->T))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"N\t%d\n",&(sp->N))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"V\t%f\n",&(sp->V))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( (sp->T)<__FLT_MIN__){ printf("T must be positive.\n"); exit(-1);}
  if( (sp->N)<__FLT_MIN__){ printf("N must be positive.\n"); exit(-1);}
  if( (sp->V)<__FLT_MIN__){ printf("V must be positive.\n"); exit(-1);}
}

void read_initial_configuration( int N, float *r, int *sim_try, FILE *f)
{
  int N_f, i_p_f, error = 0; char name[5]; float box_f;
  if( fscanf(f,"sim_try=%d, t=0.0\n",sim_try)!=1){ error = 1;}
  if( (fscanf(f,"%5d\n",&N_f)!=1)||(N_f!=N)){ error = 1;}
  for( int i_p = 0; i_p<N; i_p++)
  {
    if( fscanf(f,"%5d%5s%5s%5d",&i_p_f,name,name,&i_p_f)!=4){ error = 1;}
    if( fscanf(f,"%8f%8f%8f\n",&r[3*i_p+0],&r[3*i_p+1],&r[3*i_p+2])!=3){ error = 1;}
  }
  if( fscanf(f,"%10f%10f%10f\n",&box_f,&box_f,&box_f)!=3){ error = 1;}
  if( error){ printf("Error reading initial configuration.\n"); exit(-1);}
}

void read_pore_borders_times( int N, float *t_b[n_b], FILE *f)
{
  int i_bp_f, error = 0;
  for( int i_bp = 0; i_bp<bpm*(N-1); i_bp++)
  {
    if( (fscanf(f,"%5d",&i_bp_f)!=1)){ error = 1;}
    for( int i_b = 0; i_b<n_b; i_b++)
    {
      if( fscanf(f,"%10f",&t_b[i_b][i_bp])!=1){ error = 1;}
    }
    if( fscanf(f,"\n")!=0){ error = 1;}
  }
  if( error){ printf("Error reading pore borders times.\n"); exit(-1);}
}

void calc_statistics( int n_data, float *x, rand_var_props *x_p)
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

  //Simulation output reading

  int sim_num = 0;
  int sim_n_f = 0;

  char *knot_list;
  knot_list = (char*)malloc(22*sizeof(char));

  glob_t prev_sim_files;
  snprintf(filename,sizeof(filename),"%s/pore-borders-times-*",sim_dir);
  if( glob(filename,0,NULL,&prev_sim_files)==0)
  {
    sim_num = prev_sim_files.gl_pathc;
  }
  globfree(&prev_sim_files);

  float ***t_b = (float***)malloc(sim_num*sizeof(float**));
  float **t_w = (float**)malloc(bpm*(sp.N-1)*sizeof(float*));
  for( int i_bp = 0; i_bp<bpm*(sp.N-1); i_bp++)
  {
    t_w[i_bp] = (float*)malloc(sim_num*sizeof(float));
  }
  float *tau = (float*)malloc(sim_num*sizeof(float));

  for( int sim_idx = 0; sim_idx<sim_num; sim_idx++)
  {
    int sim_try;
    char knot_name[16];

    t_b[sim_idx] = (float**)malloc(n_b*sizeof(float*));
    for( int i_b = 0; i_b<n_b; i_b++)
    {
      t_b[sim_idx][i_b] = (float*)malloc(bpm*(sp.N-1)*sizeof(float));
    }

    snprintf(filename,sizeof(filename),"%s/initial-configuration-%03d.gro",sim_dir,sim_idx);
    file_i1 = fopen(filename,"rt");
    if( file_i1==NULL){ printf("Error opening initial configuration file.\n"); exit(-1);}
    read_initial_configuration(sp.N,r,&sim_try,file_i1);
    fclose(file_i1);

    calc_invariant(sp.N,r,knot_name);
    if( knot_name[0]!='0')
    {
      snprintf(knot_list+strlen(knot_list),22,"%03d %-16s\n",sim_idx,knot_name);
      knot_list = (char*)realloc(knot_list,(strlen(knot_list)+22)*sizeof(char));
    }

    sim_n_f += sim_try;

    snprintf(filename,sizeof(filename),"%s/pore-borders-times-%03d.dat",sim_dir,sim_idx);
    file_i1 = fopen(filename,"rt");
    if( file_i1==NULL){ printf("Error opening pore borders times file.\n"); exit(-1);}
    read_pore_borders_times(sp.N,t_b[sim_idx],file_i1);
    fclose(file_i1);

    for( int i_bp = 0; i_bp<bpm*(sp.N-1); i_bp++)
    {
      t_w[i_bp][sim_idx] = t_b[sim_idx][0][i_bp]-t_b[sim_idx][1][i_bp];
    }

    tau[sim_idx] = t_b[sim_idx][0][bpm*(sp.N-1)-1];
  }

  //Statistics calculation

  rand_var_props tau_p;

  calc_statistics(sim_num,tau,&tau_p);

  rand_var_props *t_w_p;
  t_w_p = (rand_var_props*)malloc(bpm*(sp.N-1)*sizeof(rand_var_props));

  for( int i_bp = 0; i_bp<bpm*(sp.N-1); i_bp++)
  {
    calc_statistics(sim_num,t_w[i_bp],&t_w_p[i_bp]);
  }

  //Results

  snprintf(filename,sizeof(filename),"%s/analysis-statistics.dat",sim_dir);
  file_o1 = fopen(filename,"wt");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}

  fprintf(file_o1,"#Analysis of %d translocations from %s.\n",sim_num,sim_dir);
  fprintf(file_o1,"#%9s%10s%10s\n","T","N","V");
  fprintf(file_o1,"%10.4f%10d%10.4f\n",sp.T,sp.N,sp.V);

  fprintf(file_o1,"\n");
  fprintf(file_o1,"#%14s%15s%15s\n","<tau>","SEM(tau)","Var(tau)");
  fprintf(file_o1,"%15.6lf%15.6lf%15.6lf\n",tau_p.m_1,tau_p.sem,tau_p.var);

  fprintf(file_o1,"\n");
  fprintf(file_o1,"#SR=%5.1f%%\n",sim_num*100.0/(sim_n_f+sim_num));

  if( strlen(knot_list)>0)
  {
    fprintf(file_o1,"\n");
    fprintf(file_o1,"knots:\n");
    fprintf(file_o1,"%s",knot_list);
  }

  fclose(file_o1);

  snprintf(filename,sizeof(filename),"%s/analysis-waiting-times.dat",sim_dir);
  file_o1 = fopen(filename,"wt");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}

  for( int i_bp = 0; i_bp<bpm*(sp.N-1); i_bp++)
  {
    fprintf(file_o1,"%05d%15.6lf%15.6lf",i_bp,t_w_p[i_bp].m_1,t_w_p[i_bp].sem);
    fprintf(file_o1,"%15.6lf\n",t_w_p[i_bp].var);
  }

  fclose(file_o1);

  printf("Analysed %d translocations from %s.\n",sim_num,sim_dir);

  return 0;
}
