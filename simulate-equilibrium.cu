//Libraries

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <curand_kernel.h>

//Defines

#define BROWNIAN

#define m   1.00
#define xi  1.00
#define l_0 1.00
#define k_e 100.0
#define k_b 0.0
#define k_B 0.00336
#define sig 1.5
#define eps 3.0
#define r_c pow(2,1.0/6.0)*sig

#ifdef BROWNIAN
#define dt  1.0/2048.0
#else
#define dt  0.001
#endif

//Structs

struct sim_params
{
  float T;
  float eta;
  int N;
  float f_pull;
  char IC_type;
  float meas_freq;
  float sim_time;
};

typedef struct sim_params sim_params;
typedef struct curandStatePhilox4_32_10 PRNGstate;

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
  if( fscanf(f,"eta\t%f\n",&(sp->eta))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"N\t%d\n",&(sp->N))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"f_pull\t%f\n",&(sp->f_pull))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"IC_type\t%c\n",&(sp->IC_type))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"meas_freq\t%f\n",&(sp->meas_freq))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
  if( fscanf(f,"sim_time\t%f\n",&(sp->sim_time))!=1){ printf("Error reading parameters file.\n"); exit(-1);}
}

void write_gro_frame( int N, float *r, FILE *f)
{
  fprintf(f,"MD simulation of a polymer, t = 0.0\n");
  fprintf(f,"%5d\n",N);
  for( int i_p = 0; i_p<N; i_p++)
  {
    fprintf(f,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",i_p+1,"X","X",i_p+1,r[3*i_p+0],r[3*i_p+1],r[3*i_p+2]);
  }
  fprintf(f,"%10.5f%10.5f%10.5f\n",0.0,0.0,0.0);
}

void write_trr_frame( int N, float *r, float t, int i_f, FILE *f)
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

void write_checkpoint( int N, float *r, float *v, float t, int n_threads, PRNGstate *state, int f_idx, FILE *f)
{
  int natoms = N;
  fwrite(&natoms,sizeof(int),1,f);
  int time = t;
  fwrite(&time,sizeof(float),1,f);
  int index = f_idx;
  fwrite(&index,sizeof(int),1,f);
  fwrite(r,sizeof(float),3*N,f);
  fwrite(v,sizeof(float),3*N,f);
  fwrite(state,sizeof(PRNGstate),n_threads,f);
}

void read_checkpoint( int N, float *r, float *v, float *t, int n_threads, PRNGstate *state, int *f_idx, FILE *f)
{
  int natoms;
  fread(&natoms,sizeof(int),1,f);
  if( natoms!=N){ printf("Error reading checkpoint.\n"); exit(-1);}
  int time;
  fread(&time,sizeof(float),1,f);
  int index;
  fread(&index,sizeof(int),1,f);
  fread(r,sizeof(float),3*N,f);
  fread(v,sizeof(float),3*N,f);
  fread(state,sizeof(PRNGstate),n_threads,f);
  *t=time;
  *f_idx=index+1;
}

void set_random_IC( int N, float T, float *r, float *v)
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
    if( k_b<__FLT_MIN__){ bondangle = acos(1.0-2.0*random);}
    else{ bondangle = acos(1.0+log(1.0-(1.0-exp(-2.0*(k_b*beta)))*random)/(k_b*beta));}
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
  //velocities
  for( int i_p=0; i_p<N; i_p++)
  {
    for( int i_c=0; i_c<3; i_c++)
    {
      curandGenerateUniform(gen,&random,1);
      v[3*i_p+i_c] = sqrt(2.0/(m*beta))*erfinv(2.0*random-1.0);
    }
  }
  curandDestroyGenerator(gen);
}

void set_linear_IC( int N, float *r, float *v)
{
  r[0]=r[1]=r[2]=0.0;
  for( int i_p=1; i_p<N; i_p++)
  {
    r[3*i_p+0] = l_0+r[3*(i_p-1)+0];
    r[3*i_p+1] = r[3*(i_p-1)+1];
    r[3*i_p+2] = r[3*(i_p-1)+2];
  }
  for( int i_p=0; i_p<N; i_p++)
  {
    for( int i_c=0; i_c<3; i_c++)
    {
      v[3*i_p+i_c] = 0.0;
    }
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
void calc_extern_f( int N, float eta, float *v, float *f_c, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      f[3*i_p+i_c] = f_c[3*i_p+i_c];

      #ifndef BROWNIAN
      f[3*i_p+i_c] += -eta*v[3*i_p+i_c];
      #endif
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
    k_LJ = 4.0*eps*(12.0*pow(sig,12)/(d2*d2*d2*d2*d2*d2*d2)-6.0*pow(sig,6)/(d2*d2*d2*d2));
    for( int i_c = 0; i_c<3; i_c++)
    {
      f[3*i_p+i_c] += k_LJ*(r[3*i_p+i_c]-r[3*j_p+i_c]);
    }
    return 0;
  }
  else
  {
    return ((sqrt(d2)-r_c)/(1.5*l_0));
  }
}

__global__
void calc_exclvol_f( int N, float *r, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    int skip;
    for( int j_p = i_p-2; j_p>=0; j_p -= 1+skip)
    {
      skip = calc_LJ_f(N,r,i_p,j_p,f);
    }
    for( int j_p = i_p+2; j_p<N; j_p += 1+skip)
    {
      skip = calc_LJ_f(N,r,i_p,j_p,f);
    }
  }
}

__global__
void RK_stage_1( int N, float *r_1, float *r_2, float *v_1, float *v_2, float *f_2, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      #ifdef BROWNIAN
      r_1[3*i_p+i_c] = r_2[3*i_p+i_c]+f_2[3*i_p+i_c]*dt/xi+nrn[3*i_p+i_c]/xi;
      #else
      r_1[3*i_p+i_c] = r_2[3*i_p+i_c]+v_2[3*i_p+i_c]*dt;

      v_1[3*i_p+i_c] = v_2[3*i_p+i_c]+f_2[3*i_p+i_c]*dt/m+nrn[3*i_p+i_c]/m;
      #endif
    }
  }
}

__global__
void RK_stage_2( int N, float *r_2, float *v_1, float *v_2, float *f_1, float *f_2, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      #ifdef BROWNIAN
      r_2[3*i_p+i_c] = r_2[3*i_p+i_c]+0.5*(f_1[3*i_p+i_c]+f_2[3*i_p+i_c])*dt/xi+nrn[3*i_p+i_c]/xi;
      #else
      r_2[3*i_p+i_c] = r_2[3*i_p+i_c]+0.5*(v_1[3*i_p+i_c]+v_2[3*i_p+i_c])*dt;

      v_2[3*i_p+i_c] = v_2[3*i_p+i_c]+0.5*(f_1[3*i_p+i_c]+f_2[3*i_p+i_c])*dt/m+nrn[3*i_p+i_c]/m;
      #endif
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

  snprintf(filename,sizeof(filename),"%s/parameters.dat",sim_dir);
  file_i1 = fopen(filename,"rt");
  if( file_i1==NULL){ printf("Error opening parameters file.\n"); exit(-1);}
  read_parameters(&sp,file_i1);
  fclose(file_i1);

  size_t threads_block = 256;
  size_t n_blocks = (sp.N+threads_block-1)/threads_block;
  size_t n_threads = n_blocks*threads_block;

  #ifdef BROWNIAN
  float c_rn = sqrt(2.0*xi*k_B*sp.T*dt);
  #else
  float c_rn = sqrt(2.0*sp.eta*k_B*sp.T*dt);
  #endif

  float *r_2;
  float *r_1;

  float *v_2;
  float *v_1;

  float *f_2;
  float *f_1;

  float *nrn;
  PRNGstate *state;

  float *b;
  float *invlen;
  float *cosine;

  float *f_c;

  //Memory allocation

  cuda_check( cudaMallocManaged( &r_2, 3*sp.N*sizeof(float)));
  cuda_check( cudaMallocManaged( &r_1, 3*sp.N*sizeof(float)));

  cuda_check( cudaMallocManaged( &v_2, 3*sp.N*sizeof(float)));
  cuda_check( cudaMallocManaged( &v_1, 3*sp.N*sizeof(float)));

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
  f_c[3*0+0] = -sp.f_pull;
  f_c[3*(sp.N-1)+0] = sp.f_pull;

  //PRNG initialization

  setup_PRNG<<<n_blocks,threads_block>>>(time(NULL),state);
  cuda_check( cudaDeviceSynchronize());

  //Initial condition

  int f_idx = 0;

  float t = 0.0;

  snprintf(filename,sizeof(filename),"%s/checkpoint.bin",sim_dir);
  file_i1 = fopen(filename,"rb");

  if( file_i1==NULL)
  {
    if( sp.IC_type=='f')
    {
      snprintf(filename,sizeof(filename),"%s/initial-condition.bin",sim_dir);
      file_i1 = fopen(filename,"rb");
      if( file_i1==NULL){ printf("Error opening file.\n"); exit(-1);}
      read_checkpoint(sp.N,r_2,v_2,&t,n_threads,state,&f_idx,file_i1);
      fclose(file_i1);
      f_idx = 0;
      t = 0.0;
    }
    else if( sp.IC_type=='r')
    {
      set_random_IC(sp.N,sp.T,r_2,v_2);
    }
    else if( sp.IC_type=='l')
    {
      set_linear_IC(sp.N,r_2,v_2);
    }
    else
    {
      printf("Unknown IC.\n");
      exit(-1);
    }
  }
  else
  {
    read_checkpoint(sp.N,r_1,v_1,&t,n_threads,state,&f_idx,file_i1);
    fclose(file_i1);
  }

  //Simulation

  printf("%s parameters: T=%05.1f eta=%06.3f N=%04d f_pull=%06.3f\n",sim_dir,sp.T,sp.eta,sp.N,sp.f_pull);
  printf("IC_type=%c meas_freq=%06.3f sim_time=%08.1f f_idx=%d\n",sp.IC_type,sp.meas_freq,sp.sim_time,f_idx);

  if( f_idx==0)
  {
    snprintf(filename,sizeof(filename),"%s/initial-condition.gro",sim_dir);
    file_o1 = fopen(filename,"wt");
    if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}
    write_gro_frame(sp.N,r_2,file_o1);
    fclose(file_o1);
  }

  snprintf(filename,sizeof(filename),"%s/trajectory-file-%d.trr",sim_dir,f_idx);
  file_o1 = fopen(filename,"wb");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}

  int i_s;
  int n_steps = round((1.0/sp.meas_freq)/(dt));
  int i_f;
  int n_frames = round(sp.sim_time*sp.meas_freq);

  for( i_f = 0; i_f < n_frames; i_f++)
  {
    printf("Progress:%05.1lf%%\r",(100.0*i_f)/(1.0*n_frames));
    fflush(stdout);

    for( i_s = 0; i_s < n_steps; i_s++)
    {
      call_PRNG<<<n_blocks,threads_block>>>(c_rn,nrn,state);

      calc_extern_f<<<n_blocks,threads_block>>>(sp.N,sp.eta,v_2,f_c,f_2);
      calc_bonds<<<n_blocks,threads_block>>>(sp.N,r_2,b,invlen);
      calc_cosines<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine);
      calc_intern_f<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine,f_2);
      // calc_exclvol_f<<<n_blocks,threads_block>>>(sp.N,r_2,f_2);

      RK_stage_1<<<n_blocks,threads_block>>>(sp.N,r_1,r_2,v_1,v_2,f_2,nrn);

      calc_extern_f<<<n_blocks,threads_block>>>(sp.N,sp.eta,v_1,f_c,f_1);
      calc_bonds<<<n_blocks,threads_block>>>(sp.N,r_1,b,invlen);
      calc_cosines<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine);
      calc_intern_f<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine,f_1);
      // calc_exclvol_f<<<n_blocks,threads_block>>>(sp.N,r_1,f_1);

      RK_stage_2<<<n_blocks,threads_block>>>(sp.N,r_2,v_1,v_2,f_1,f_2,nrn);

      cuda_check( cudaDeviceSynchronize());
    }

    t += n_steps*dt;

    write_trr_frame(sp.N,r_2,t,i_f,file_o1);
  }

  fclose(file_o1);

  snprintf(filename,sizeof(filename),"%s/checkpoint.bin",sim_dir);
  file_o1 = fopen(filename,"wb");
  if( file_o1==NULL){ printf("Error opening file.\n"); exit(-1);}
  write_checkpoint(sp.N,r_2,v_2,t,n_threads,state,f_idx,file_o1);
  fclose(file_o1);

  //Memory deallocation

  cudaFree(r_2);
  cudaFree(r_1);

  cudaFree(v_2);
  cudaFree(v_1);

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
