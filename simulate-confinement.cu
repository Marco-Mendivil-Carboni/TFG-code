//Libraries

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <curand_kernel.h>
#include <glob.h>

//Defines

#define xi  1.00
#define l_0 0.591
#define q_e 0.0521
#define k_e 550.0
#define k_b 6.66
#define k_B 0.00112
#define r_c 1.1224620

#define dt  1.0/2048.0
#define n_s 10*2048
#define n_m 1000

#define rct 100

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
typedef struct curandStatePhilox4_32_10 PRNGstate;
typedef struct rand_var_props rand_var_props;

//Functions

inline void cuda_check( cudaError_t result)
{
  if( result!=cudaSuccess)
  {
    fprintf(stderr,"CUDA Error: %s\n",cudaGetErrorString(result));
    exit(-1);
  }
}

void read_parameters( sim_params *sp, FILE *f)
{
  if( fscanf(f,"T\t%f\n",&(sp->T))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"N\t%d\n",&(sp->N))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"R\t%f\n",&(sp->R))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( (sp->T)<__FLT_MIN__){ printf("T must be positive.\n"); exit(-1);}
  if( (sp->N)<__FLT_MIN__){ printf("N must be positive.\n"); exit(-1);}
  if( (sp->R)<__FLT_MIN__){ printf("R must be positive.\n"); exit(-1);}
}

void generate_initial_configuration( int N, float T, float *r)
{
  curandGenerator_t gen;
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(NULL));
  float random, theta, varphi, bondlen, bondangle;
  float dir_old[3], dir_new[3], perpdir[3], perpdirnorm;
  float beta = 1.0/(k_B*T);
  //positions
  curandGenerateUniform(gen,&random,1); theta = acos(1.0-2.0*random); 
  curandGenerateUniform(gen,&random,1); varphi = 2.0*M_PI*random;
  dir_old[0]=sin(theta)*cos(varphi);
  dir_old[1]=sin(theta)*sin(varphi);
  dir_old[2]=cos(theta);
  r[0]=r[1]=r[2]=0.0;
  for( int i_p = 1; i_p<N; i_p++)
  {
    curandGenerateUniform(gen,&random,1); theta = acos(1.0-2.0*random);
    curandGenerateUniform(gen,&random,1); varphi = 2.0*M_PI*random;
    perpdir[0] = dir_old[1]*cos(theta)-dir_old[2]*sin(theta)*sin(varphi);
    perpdir[1] = dir_old[2]*sin(theta)*cos(varphi)-dir_old[0]*cos(theta);
    perpdir[2] = dir_old[0]*sin(theta)*sin(varphi)-dir_old[1]*sin(theta)*cos(varphi);
    perpdirnorm = sqrt(perpdir[0]*perpdir[0]+perpdir[1]*perpdir[1]+perpdir[2]*perpdir[2]);
    perpdir[0] /= perpdirnorm; perpdir[1] /= perpdirnorm; perpdir[2] /= perpdirnorm;
    curandGenerateUniform(gen,&random,1);
    bondangle = acos(1.0+log(1.0-(1.0-exp(-2.0*(k_b*beta)))*random)/(k_b*beta));
    dir_new[0] = dir_old[0]*cos(bondangle)+perpdir[0]*sin(bondangle);
    dir_new[1] = dir_old[1]*cos(bondangle)+perpdir[1]*sin(bondangle);
    dir_new[2] = dir_old[2]*cos(bondangle)+perpdir[2]*sin(bondangle);
    curandGenerateUniform(gen,&random,1);
    bondlen = l_0+sqrt(2.0/(k_e*beta))*erfinv(2.0*random-1.0);
    int accept = 1;
    r[3*i_p+0] = bondlen*dir_new[0]+r[3*(i_p-1)+0]; if( !isfinite(r[3*i_p+0])){ accept = 0;}
    r[3*i_p+1] = bondlen*dir_new[1]+r[3*(i_p-1)+1]; if( !isfinite(r[3*i_p+1])){ accept = 0;}
    r[3*i_p+2] = bondlen*dir_new[2]+r[3*(i_p-1)+2]; if( !isfinite(r[3*i_p+2])){ accept = 0;}
    for( int j_p = 0; j_p<(i_p-1); j_p++)
    {
      float dist = 0.0;
      dist += (r[3*j_p+0]-r[3*i_p+0])*(r[3*j_p+0]-r[3*i_p+0]);
      dist += (r[3*j_p+1]-r[3*i_p+1])*(r[3*j_p+1]-r[3*i_p+1]);
      dist += (r[3*j_p+2]-r[3*i_p+2])*(r[3*j_p+2]-r[3*i_p+2]);
      dist = sqrt(dist);
      if( dist<r_c){ accept = 0; break;}
    }
    if( accept)
    {
      dir_old[0] = dir_new[0];
      dir_old[1] = dir_new[1];
      dir_old[2] = dir_new[2];
    }
    else{ i_p--;}
  }
  curandDestroyGenerator(gen);
}

void write_initial_configuration( int N, float *r, FILE *f)
{
  fprintf(f,"Confinement simulation, t=0.0\n");
  fprintf(f,"%5d\n",N);
  for( int i_p = 0; i_p<N; i_p++)
  {
    fprintf(f,"%5d%-5s%5s%5d",i_p+1,"X","X",i_p+1);
    fprintf(f,"%8.3f%8.3f%8.3f\n",r[3*i_p+0],r[3*i_p+1],r[3*i_p+2]);
  }
  fprintf(f,"%10.5f%10.5f%10.5f\n",0.0,0.0,0.0);
}

void write_trajectory_positions( int N, float *r, float t, int i_f, FILE *f)
{
  //header
  int magickvalue = 1993;
  fwrite(&magickvalue,sizeof(int),1,f);
  char trrversion[] = "GMX_trn_file";
  int len_s_a = sizeof(trrversion);
  int len_s_b = sizeof(trrversion)-1;
  fwrite(&len_s_a,sizeof(int),1,f);
  fwrite(&len_s_b,sizeof(int),1,f);
  fwrite(trrversion,sizeof(char),sizeof(trrversion)-1,f);
  int zero = 0;
  for( int i = 0; i<7; i++)
  {
    fwrite(&zero,sizeof(int),1,f);
  }
  int x_size = 3*N*sizeof(float);
  fwrite(&x_size,sizeof(int),1,f);
  int v_size = 0;
  fwrite(&v_size,sizeof(int),1,f);
  int f_size = 0;
  fwrite(&f_size,sizeof(int),1,f);
  int natoms = N;
  fwrite(&natoms,sizeof(int),1,f);
  int step = i_f;
  fwrite(&step,sizeof(int),1,f);
  float time = t;
  fwrite(&zero,sizeof(int),1,f);
  fwrite(&time,sizeof(float),1,f);
  fwrite(&zero,sizeof(int),1,f);
  //coordinates
  fwrite(r,sizeof(float),3*N,f);
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
  for( i_term = 0; i_term<4; i_term++)
  {
    n_term = i_term*n_data/16;
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

//Kernels

__global__
void setup_PRNG( int seed, PRNGstate *state)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(seed, i_p, 0, &state[i_p]);
}

__global__
void call_PRNG( float c_rn, float *nrn, PRNGstate *state)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  for( int i_c = 0; i_c<3; i_c++)
  {
    nrn[3*i_p+i_c] = c_rn*curand_normal(&state[i_p]);
  }
}

__global__
void calc_extern_f( int N, float *f_c, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      f[3*i_p+i_c] = f_c[3*i_p+i_c];
    }
  }
}

__global__
void calc_sphere_f( int N, float R, float *r, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    float k_LJ;
    float d_r = 0.0;
    for( int i_c = 0; i_c<3; i_c++)
    {
      d_r += r[3*i_p+i_c]*r[3*i_p+i_c];
    }
    d_r = sqrt(d_r);
    float dsp = (R-d_r);
    float d2 = dsp*dsp;
    if( d2<r_c*r_c)
    {
      k_LJ = 48.0/(d2*d2*d2*d2*d2*d2*d2)-24.0/(d2*d2*d2*d2);
      for( int i_c = 0; i_c<3; i_c++)
      {
        f[3*i_p+i_c] += k_LJ*(-dsp/d_r)*r[3*i_p+i_c];
      }
    }
  }
}

__global__ 
void calc_bonds( int N, float *r, float *b, float *invlen)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N-1)
  {
    invlen[i_p+2] = 0.0;
    for( int i_c = 0; i_c<3; i_c++)
    {
      b[3*(i_p+2)+i_c] = r[3*(i_p+1)+i_c]-r[3*i_p+i_c];
      invlen[i_p+2] += b[3*(i_p+2)+i_c]*b[3*(i_p+2)+i_c];
    }
    invlen[i_p+2] = 1.0/sqrt(invlen[i_p+2]);
  }
}

__global__
void calc_cosines( int N, float *b, float *invlen, float *cosine)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N-2)
  {
    cosine[i_p+3] = 0.0;
    for( int i_c = 0; i_c<3; i_c++)
    {
      cosine[i_p+3] += b[3*(i_p+3)+i_c]*b[3*(i_p+2)+i_c];
    }
    cosine[i_p+3] *= invlen[i_p+3]*invlen[i_p+2];
  }
}

__global__
void calc_intern_f( int N, float *b, float *invlen, float *cosine, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      f[3*i_p+i_c] += k_e*(1.0-l_0*invlen[i_p+1])*(-b[3*(i_p+1)+i_c]);

      f[3*i_p+i_c] += k_e*(1.0-l_0*invlen[i_p+2])*(+b[3*(i_p+2)+i_c]);

      f[3*i_p+i_c] += k_b*(+b[3*(i_p+0)+i_c])*invlen[i_p+1]*invlen[i_p+0];

      f[3*i_p+i_c] += k_b*(-cosine[i_p+2]-cosine[i_p+1])*b[3*(i_p+1)+i_c]*invlen[i_p+1]*invlen[i_p+1];

      f[3*i_p+i_c] += k_b*(+b[3*(i_p+2)+i_c]-b[3*(i_p+1)+i_c])*invlen[i_p+2]*invlen[i_p+1];

      f[3*i_p+i_c] += k_b*(+cosine[i_p+3]+cosine[i_p+2])*b[3*(i_p+2)+i_c]*invlen[i_p+2]*invlen[i_p+2];

      f[3*i_p+i_c] += k_b*(-b[3*(i_p+3)+i_c])*invlen[i_p+3]*invlen[i_p+2];
    }
  }
}

__device__
int calc_LJ_f( int N, float *r, int i_p, int j_p, float *f)
{
  float k_LJ;
  float d2 = 0.0;
  for( int i_c = 0; i_c<3; i_c++)
  {
    d2 += (r[3*i_p+i_c]-r[3*j_p+i_c])*(r[3*i_p+i_c]-r[3*j_p+i_c]);
  }
  if( d2<r_c*r_c)
  {
    k_LJ = 48.0/(d2*d2*d2*d2*d2*d2*d2)-24.0/(d2*d2*d2*d2);
    for( int i_c = 0; i_c<3; i_c++)
    {
      f[3*i_p+i_c] += k_LJ*(r[3*i_p+i_c]-r[3*j_p+i_c]);
    }
    return 0;
  }
  else
  {
    return ((sqrt(d2)-r_c)/(1.25*l_0));
  }
}

__global__
void calc_exclvol_f( int N, float *r, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    int skip;
    for( int j_p = i_p-3; j_p>=0; j_p -= 1+skip)
    {
      skip = calc_LJ_f(N,r,i_p,j_p,f);
    }
    for( int j_p = i_p+3; j_p<N; j_p += 1+skip)
    {
      skip = calc_LJ_f(N,r,i_p,j_p,f);
    }
  }
}

__global__
void RK_stage_1( int N, float *r_1, float *r_2, float *f_2, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      r_1[3*i_p+i_c] = r_2[3*i_p+i_c]+f_2[3*i_p+i_c]*dt/xi+nrn[3*i_p+i_c]/xi;
    }
  }
}

__global__
void RK_stage_2( int N, float *r_2, float *f_1, float *f_2, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      r_2[3*i_p+i_c] = r_2[3*i_p+i_c]+0.5*(f_1[3*i_p+i_c]+f_2[3*i_p+i_c])*dt/xi+nrn[3*i_p+i_c]/xi;
    }
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

  char filename[256];

  //Simulation parameters and variables

  sim_params sp;

  snprintf(filename,sizeof(filename),"%s/adjustable-parameters.dat",sim_dir);
  file_i1 = fopen(filename,"rt");
  if( file_i1==NULL){ printf("Error opening parameters file.\n"); exit(-1);}
  read_parameters(&sp,file_i1);
  fclose(file_i1);

  size_t threads_block = 256;
  size_t n_blocks = (sp.N+threads_block-1)/threads_block;
  size_t n_threads = n_blocks*threads_block;

  float c_rn = sqrt(2.0*xi*k_B*sp.T*dt);

  float *r_2;
  float *r_1;

  float *f_2;
  float *f_1;

  float *nrn;
  PRNGstate *state;

  float *b;
  float *invlen;
  float *cosine;

  float *f_c;

  float R;

  double r_cm[3];
  double *Rg2 = NULL;
  rand_var_props Rg2_p;

  //Memory allocation

  cuda_check( cudaMallocManaged( &r_2, 3*sp.N*sizeof(float)));
  cuda_check( cudaMallocManaged( &r_1, 3*sp.N*sizeof(float)));

  cuda_check( cudaMallocManaged( &f_2, 3*sp.N*sizeof(float)));
  cuda_check( cudaMallocManaged( &f_1, 3*sp.N*sizeof(float)));

  cuda_check( cudaMallocManaged( &nrn, 3*n_threads*sizeof(float)));
  cuda_check( cudaMallocManaged( &state, n_threads*sizeof(PRNGstate)));

  cuda_check( cudaMallocManaged( &b, 3*(sp.N+3)*sizeof(float)));
  cuda_check( cudaMallocManaged( &invlen, (sp.N+3)*sizeof(float)));
  cuda_check( cudaMallocManaged( &cosine, (sp.N+4)*sizeof(float)));

  cuda_check( cudaMallocManaged( &f_c, 3*sp.N*sizeof(float)));

  //Exceptions for the polymer ends

  b[0] = b[1] = b[2] = 0.0;
  b[3] = b[4] = b[5] = 0.0;
  b[3*(sp.N+2)+0] = b[3*(sp.N+2)+1] = b[3*(sp.N+2)+2] = 0.0;
  b[3*(sp.N+1)+0] = b[3*(sp.N+1)+1] = b[3*(sp.N+1)+2] = 0.0;

  invlen[0] = invlen[1] = 0.0;
  invlen[sp.N+2] = invlen[sp.N+1] = 0.0;

  cosine[0] = cosine[1] = cosine[2] = 0.0;
  cosine[sp.N+3] = cosine[sp.N+2] = cosine[sp.N+1] = 0.0;

  //Constant force

  for( int i_p = 0; i_p<sp.N; i_p++)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      f_c[3*i_p+i_c] = 0.0;
    }
  }

  //PRNG initialization

  setup_PRNG<<<n_blocks,threads_block>>>(time(NULL),state);
  cuda_check( cudaDeviceSynchronize());

  //Simulation

  int sim_idx = 0;
  int sim_n_b = 0;
  int sim_end = 0;

  glob_t prev_sim_files;
  snprintf(filename,sizeof(filename),"%s/initial-configuration-*",sim_dir);
  if( glob(filename,0,NULL,&prev_sim_files)==0)
  {
    sim_idx = prev_sim_files.gl_pathc;
  }
  globfree(&prev_sim_files);

  printf("%s: T=%05.1f N=%04d R=%05.1f sim_idx=%03d\n",sim_dir,sp.T,sp.N,sp.R,sim_idx);

  float t = 0.0;

  R = sp.N*l_0;

  generate_initial_configuration(sp.N,sp.T,r_2);

  snprintf(filename,sizeof(filename),"%s/initial-configuration-%03d.gro",sim_dir,sim_idx);
  file_o1 = fopen(filename,"wt");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}
  write_initial_configuration(sp.N,r_2,file_o1);
  fclose(file_o1);

  snprintf(filename,sizeof(filename),"%s/trajectory-positions-%03d.trr",sim_dir,sim_idx);
  file_o1 = fopen(filename,"wb");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}

  int i_f = 0;

  while( !sim_end)
  {
    sim_n_b += 1;

    Rg2 = (double*)realloc(Rg2,sim_n_b*n_m*sizeof(double));

    for( int i_m = 0; i_m<n_m; i_m++)
    {
      printf("sim_n_b=%04d Progress: %5.1lf%%\r",sim_n_b,(100.0*i_m)/(1.0*n_m));
      fflush(stdout);

      for( int i_s = 0; i_s<n_s; i_s++)
      {
        call_PRNG<<<n_blocks,threads_block>>>(c_rn,nrn,state);

        calc_extern_f<<<n_blocks,threads_block>>>(sp.N,f_c,f_2);
        calc_sphere_f<<<n_blocks,threads_block>>>(sp.N,R,r_2,f_2);
        calc_bonds<<<n_blocks,threads_block>>>(sp.N,r_2,b,invlen);
        calc_cosines<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine);
        calc_intern_f<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine,f_2);
        calc_exclvol_f<<<n_blocks,threads_block>>>(sp.N,r_2,f_2);

        RK_stage_1<<<n_blocks,threads_block>>>(sp.N,r_1,r_2,f_2,nrn);

        calc_extern_f<<<n_blocks,threads_block>>>(sp.N,f_c,f_1);
        calc_sphere_f<<<n_blocks,threads_block>>>(sp.N,R,r_1,f_1);
        calc_bonds<<<n_blocks,threads_block>>>(sp.N,r_1,b,invlen);
        calc_cosines<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine);
        calc_intern_f<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine,f_1);
        calc_exclvol_f<<<n_blocks,threads_block>>>(sp.N,r_1,f_1);

        RK_stage_2<<<n_blocks,threads_block>>>(sp.N,r_2,f_1,f_2,nrn);

        if( R>sp.R){R += -(sp.N*l_0/rct)*dt;}
      }

      cuda_check( cudaDeviceSynchronize());

      t += n_s*dt;

      i_f += 1;

      write_trajectory_positions(sp.N,r_2,t,i_f,file_o1);

      calc_cm(sp.N,r_2,r_cm);
      calc_Rg2(sp.N,r_2,r_cm,&Rg2[(sim_n_b-1)*n_m+i_m]);
    }

    Rg2_p.term_end = estimate_term(sim_n_b*n_m,Rg2,&Rg2_p.n_term);

    if( !Rg2_p.term_end)
    {
      sim_end = 0;
    }
    else
    {
      sim_end = 1;
    }
  }

  fclose(file_o1);

  //Memory deallocation

  cudaFree(r_2);
  cudaFree(r_1);

  cudaFree(f_2);
  cudaFree(f_1);

  cudaFree(nrn);
  cudaFree(state);

  cudaFree(b);
  cudaFree(invlen);
  cudaFree(cosine);

  cudaFree(f_c);

  return 0;
}
