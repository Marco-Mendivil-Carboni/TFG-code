//Libraries

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <curand_kernel.h>

//Defines

#define m   1.00
#define xi  1.00
#define l_0 0.591
#define q_e 0.0521
#define k_e 550.0
#define k_b 6.66
#define k_B 0.00112
#define r_c 1.1224620

#define dt  1.0/2048.0

#define r_p 1.15
#define l_p 1.19

#define T   298.0
#define V   250.0

#define n_s 2048
#define sim_num 32

//Structs

struct rand_var_props
{
  double m_1;
  double m_2;
  double var;
  double sem;
};

typedef struct rand_var_props rand_var_props;
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

void generate_initial_configuration( int N, float *r, float *v)
{
  curandGenerator_t gen;
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,rand());
  float random, theta, varphi, bondlen, bondangle;
  float dir_old[3], dir_new[3], perpdir[3], perpdirnorm;
  float beta = 1.0/(k_B*T);
  //positions
  int n_failures = 0;
  r[0] = r[1] = r[2] = 0.0;
  dir_old[0] = 1.0; dir_old[1] = 0.0; dir_old[2] = 0.0;
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
    if( r[3*i_p+0]<(l_p/2.0+r_c))
    {
      float d_r = sqrt(r[3*i_p+1]*r[3*i_p+1]+r[3*i_p+2]*r[3*i_p+2]);
      if( d_r>(r_p-r_c)){ accept = 0;}
      if( r[3*i_p+0]<0){ accept = 0;}
    }
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
      n_failures = 0;
    }
    else
    {
      n_failures++;
      if( n_failures>8)
      {
        dir_old[0] = 1.0;
        dir_old[1] = 0.0;
        dir_old[2] = 0.0;
        i_p = 0;
      }
      else
      {
        i_p--;
      }
    }
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
void calc_extern_f_LANGEVIN( int N, float eta, float *v, float *f_c, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      f[3*i_p+i_c] = f_c[3*i_p+i_c];

      f[3*i_p+i_c] += -eta*v[3*i_p+i_c];
    }
  }
}

__global__
void calc_extern_f_BROWNIAN( int N, float *f_c, float *f)
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
void calc_pore_f( int N, float *r, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    int side;
    float d_x, d_r;
    float dpp[3] = {0.0, 0.0, 0.0};
    side = ((r[3*i_p+0]>0)-(r[3*i_p+0]<0));
    d_x = side*r[3*i_p+0]-l_p/2.0;
    if( d_x<r_c)
    {
      d_r = sqrt(r[3*i_p+1]*r[3*i_p+1]+r[3*i_p+2]*r[3*i_p+2]);
      if( d_r<r_p)
      {
        dpp[1] = (-r[3*i_p+1]/d_r)*(r_p-d_r);
        dpp[2] = (-r[3*i_p+2]/d_r)*(r_p-d_r);
      }
      if( d_x<0)
      {
        f[3*i_p+0] += -q_e*V/l_p;
      }
      else
      {
        dpp[0]=side*d_x;
      }
      float k_LJ;
      float d2 = 0.0;
      for( int i_c = 0; i_c<3; i_c++)
      {
        d2 += dpp[i_c]*dpp[i_c];
      }
      if( d2<r_c*r_c)
      {
        k_LJ = 48.0/(d2*d2*d2*d2*d2*d2*d2)-24.0/(d2*d2*d2*d2);
        for( int i_c = 0; i_c<3; i_c++)
        {
          f[3*i_p+i_c] += k_LJ*dpp[i_c];
        }
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
void calc_exclvol_f_SLOW( int N, float *r, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int j_p = i_p-3; j_p>=0; j_p -= 1)
    {
      calc_LJ_f(N,r,i_p,j_p,f);
    }
    for( int j_p = i_p+3; j_p<N; j_p += 1)
    {
      calc_LJ_f(N,r,i_p,j_p,f);
    }
  }
}

__global__
void RK_stage_1_LANGEVIN( int N, float *r_1, float *r_2, float *v_1, float *v_2, float *f_2, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      r_1[3*i_p+i_c] = r_2[3*i_p+i_c]+v_2[3*i_p+i_c]*dt;

      v_1[3*i_p+i_c] = v_2[3*i_p+i_c]+f_2[3*i_p+i_c]*dt/m+nrn[3*i_p+i_c]/m;
    }
  }
}

__global__
void RK_stage_1_BROWNIAN( int N, float *r_1, float *r_2, float *f_2, float *nrn)
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
void RK_stage_2_LANGEVIN( int N, float *r_2, float *v_1, float *v_2, float *f_1, float *f_2, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      r_2[3*i_p+i_c] = r_2[3*i_p+i_c]+0.5*(v_1[3*i_p+i_c]+v_2[3*i_p+i_c])*dt;

      v_2[3*i_p+i_c] = v_2[3*i_p+i_c]+0.5*(f_1[3*i_p+i_c]+f_2[3*i_p+i_c])*dt/m+nrn[3*i_p+i_c]/m;
    }
  }
}

__global__
void RK_stage_2_BROWNIAN( int N, float *r_2, float *f_1, float *f_2, float *nrn)
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
  if( argc!=1){ printf("Too many arguments.\n"); exit(-1);}

  clock_t t_s;
  clock_t t_e;

  float t[3][sim_num];

  rand_var_props t_p[3];

  srand(time(NULL));

  //Simulation variables

  float c_rn = sqrt(2.0*xi*k_B*T*dt);

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

  printf("#%-4s","N");
  printf("%-20s","IDEAL");
  printf("%-20s","EXCLVOL");
  printf("%-20s","EXCLVOL SLOW");
  printf("\n");

  for( int N = 32; N<=2048; N += 32)
  {
    printf("%4d",N);

    size_t threads_block = 256;
    size_t n_blocks = (N+threads_block-1)/threads_block;
    size_t n_threads = n_blocks*threads_block;

    //Memory allocation

    cuda_check( cudaMallocManaged( &r_2, 3*N*sizeof(float)));
    cuda_check( cudaMallocManaged( &r_1, 3*N*sizeof(float)));

    cuda_check( cudaMallocManaged( &v_2, 3*N*sizeof(float)));
    cuda_check( cudaMallocManaged( &v_1, 3*N*sizeof(float)));

    cuda_check( cudaMallocManaged( &f_2, 3*N*sizeof(float)));
    cuda_check( cudaMallocManaged( &f_1, 3*N*sizeof(float)));

    cuda_check( cudaMallocManaged( &nrn, 3*n_threads*sizeof(float)));
    cuda_check( cudaMallocManaged( &state, n_threads*sizeof(PRNGstate)));

    cuda_check( cudaMallocManaged( &b, 3*(N+3)*sizeof(float)));
    cuda_check( cudaMallocManaged( &invlen, (N+3)*sizeof(float)));
    cuda_check( cudaMallocManaged( &cosine, (N+4)*sizeof(float)));

    cuda_check( cudaMallocManaged( &f_c, 3*N*sizeof(float)));

    //Exceptions for the polymer ends

    b[0] = b[1] = b[2] = 0.0;
    b[3] = b[4] = b[5] = 0.0;
    b[3*(N+2)+0] = b[3*(N+2)+1] = b[3*(N+2)+2] = 0.0;
    b[3*(N+1)+0] = b[3*(N+1)+1] = b[3*(N+1)+2] = 0.0;

    invlen[0] = invlen[1] = 0.0;
    invlen[N+2] = invlen[N+1] = 0.0;

    cosine[0] = cosine[1] = cosine[2] = 0.0;
    cosine[N+3] = cosine[N+2] = cosine[N+1] = 0.0;

    //Constant force

    for( int i_p = 0; i_p<N; i_p++)
    {
      for( int i_c = 0; i_c<3; i_c++)
      {
        f_c[3*i_p+i_c] = 0.0;
      }
    }

    //PRNG initialization

    setup_PRNG<<<n_blocks,threads_block>>>(time(NULL),state);
    cuda_check( cudaDeviceSynchronize());

    //Simulations

    for( int sim_idx = 0; sim_idx<sim_num; sim_idx++)
    {
      //IDEAL

      generate_initial_configuration(N,r_2,v_2);

      t_s = clock();

      for( int i_s = 0; i_s<n_s; i_s++)
      {
        call_PRNG<<<n_blocks,threads_block>>>(c_rn,nrn,state);

        calc_extern_f_LANGEVIN<<<n_blocks,threads_block>>>(N,xi,v_2,f_c,f_2);
        calc_bonds<<<n_blocks,threads_block>>>(N,r_2,b,invlen);
        calc_cosines<<<n_blocks,threads_block>>>(N,b,invlen,cosine);
        calc_intern_f<<<n_blocks,threads_block>>>(N,b,invlen,cosine,f_2);

        RK_stage_1_LANGEVIN<<<n_blocks,threads_block>>>(N,r_1,r_2,v_1,v_2,f_2,nrn);

        calc_extern_f_LANGEVIN<<<n_blocks,threads_block>>>(N,xi,v_1,f_c,f_1);
        calc_bonds<<<n_blocks,threads_block>>>(N,r_1,b,invlen);
        calc_cosines<<<n_blocks,threads_block>>>(N,b,invlen,cosine);
        calc_intern_f<<<n_blocks,threads_block>>>(N,b,invlen,cosine,f_1);

        RK_stage_2_LANGEVIN<<<n_blocks,threads_block>>>(N,r_2,v_1,v_2,f_1,f_2,nrn);
      }

      cuda_check( cudaDeviceSynchronize());

      t_e = clock();

      t[0][sim_idx] = (1000.0/n_s)*((double)(t_e-t_s))/CLOCKS_PER_SEC;

      //EXCLVOL

      generate_initial_configuration(N,r_2,v_2);

      t_s = clock();

      for( int i_s = 0; i_s<n_s; i_s++)
      {
        call_PRNG<<<n_blocks,threads_block>>>(c_rn,nrn,state);

        calc_extern_f_LANGEVIN<<<n_blocks,threads_block>>>(N,xi,v_2,f_c,f_2);
        calc_bonds<<<n_blocks,threads_block>>>(N,r_2,b,invlen);
        calc_cosines<<<n_blocks,threads_block>>>(N,b,invlen,cosine);
        calc_intern_f<<<n_blocks,threads_block>>>(N,b,invlen,cosine,f_2);
        calc_exclvol_f<<<n_blocks,threads_block>>>(N,r_2,f_2);

        RK_stage_1_LANGEVIN<<<n_blocks,threads_block>>>(N,r_1,r_2,v_1,v_2,f_2,nrn);

        calc_extern_f_LANGEVIN<<<n_blocks,threads_block>>>(N,xi,v_1,f_c,f_1);
        calc_bonds<<<n_blocks,threads_block>>>(N,r_1,b,invlen);
        calc_cosines<<<n_blocks,threads_block>>>(N,b,invlen,cosine);
        calc_intern_f<<<n_blocks,threads_block>>>(N,b,invlen,cosine,f_1);
        calc_exclvol_f<<<n_blocks,threads_block>>>(N,r_1,f_1);

        RK_stage_2_LANGEVIN<<<n_blocks,threads_block>>>(N,r_2,v_1,v_2,f_1,f_2,nrn);
      }

      cuda_check( cudaDeviceSynchronize());

      t_e = clock();

      t[1][sim_idx] = (1000.0/n_s)*((double)(t_e-t_s))/CLOCKS_PER_SEC;

      //EXCLVOL SLOW

      generate_initial_configuration(N,r_2,v_2);

      t_s = clock();

      for( int i_s = 0; i_s<n_s; i_s++)
      {
        call_PRNG<<<n_blocks,threads_block>>>(c_rn,nrn,state);

        calc_extern_f_LANGEVIN<<<n_blocks,threads_block>>>(N,xi,v_2,f_c,f_2);
        calc_bonds<<<n_blocks,threads_block>>>(N,r_2,b,invlen);
        calc_cosines<<<n_blocks,threads_block>>>(N,b,invlen,cosine);
        calc_intern_f<<<n_blocks,threads_block>>>(N,b,invlen,cosine,f_2);
        calc_exclvol_f_SLOW<<<n_blocks,threads_block>>>(N,r_2,f_2);

        RK_stage_1_LANGEVIN<<<n_blocks,threads_block>>>(N,r_1,r_2,v_1,v_2,f_2,nrn);

        calc_extern_f_LANGEVIN<<<n_blocks,threads_block>>>(N,xi,v_1,f_c,f_1);
        calc_bonds<<<n_blocks,threads_block>>>(N,r_1,b,invlen);
        calc_cosines<<<n_blocks,threads_block>>>(N,b,invlen,cosine);
        calc_intern_f<<<n_blocks,threads_block>>>(N,b,invlen,cosine,f_1);
        calc_exclvol_f_SLOW<<<n_blocks,threads_block>>>(N,r_1,f_1);

        RK_stage_2_LANGEVIN<<<n_blocks,threads_block>>>(N,r_2,v_1,v_2,f_1,f_2,nrn);
      }

      cuda_check( cudaDeviceSynchronize());

      t_e = clock();

      t[2][sim_idx] = (1000.0/n_s)*((double)(t_e-t_s))/CLOCKS_PER_SEC;
    }

    calc_statistics(sim_num,t[0],&t_p[0]);
    printf("%10.6lf%10.6lf",t_p[0].m_1,t_p[0].sem);

    calc_statistics(sim_num,t[1],&t_p[1]);
    printf("%10.6lf%10.6lf",t_p[1].m_1,t_p[1].sem);

    calc_statistics(sim_num,t[2],&t_p[2]);
    printf("%10.6lf%10.6lf",t_p[2].m_1,t_p[2].sem);

    printf("\n");

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
  }

  return 0;
}
