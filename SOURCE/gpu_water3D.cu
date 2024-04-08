
/****************************************************************/
/*								*/
/*			Bulk FS Water				*/
/*								*/
/*								*/
/****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <curand.h>
#include <curand_kernel.h>
#include <sys/time.h>
#include <vector_types.h>
#include <string.h>


#include "mersenne.h"
#include "chessboard.h"

#define delta(X, Y)  ((X) == (Y)) // ? 1 : 0)
#define theta(X, Y)  ((X) <= (Y)) // ? 1 : 0)

//#define CPU
#define GPU

/* Lattice NX*NY*NZ, adjust */
#define NX 32            
#define NY 32
#define NZ 32

/* Size of the random number arrays (GPU calculation) */
#define N_RANDOM 500 

/* How many iterations for each measurements */
int N_MEASURE;

int NLOOPS;

int NEQUILIBRIUM;

/* 1 MC step = NMETROPOLIS + NCLUSTER */
int NMETROPOLIS, NCLUSTER;

void MonteCarlo_Step(int time, int &loop, int &offset, int &loop_chess, int &offset_chess, int nblocks, int nthreads, int *order,
     FILE* logfile, int *cluster_size, FILE* logsize);


/* Model parameters */
#define J 		1.f
#define J_hb 		0.5f
#define J_sig 	0.08f
#define v_hb 		0.6f

#define q	6
#define NHBMAX 	4
#define ARMS		6

/* Number of parallel threads per block */
//#define NTHREADS 1024
#define NTHREADS 512
//#define NTHREADS 64
//#define NTHREADS 36
//#define NTHREADS 1

int flag_config;
int flag_distances;

float T, P;  //Temperature, Pressure
double V;    //Volume

//   xorshift1024*Phi pseudorandom number generator
uint64_t xor_s[16]; 
int xor_p;
void initRand(uint64_t  seed);
uint64_t nextRand(void);
double uniformDoubleRand(void);

/// This macro activates the check for errors with the GPU (may influence performance)
#define CUDA_ERROR_CHECK
#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

typedef uint8_t byte;

int deviceID;
void __cudaSafeCall( cudaError err, const char *file, const int line );
void __cudaCheckError( const char *file, const int line );

__global__ void gpu_RNG_setup (curandState * state, uint * seed);
__global__ void gpu_RNG_generate ( curandState* globalState, float * Rnd, int n_rand);
__global__ void gpu_update (byte * dev_s, byte * dev_active_bond, byte * dev_nhb,
   float * dev_Rnd, float * dev_Rnd_spin, int indices, float P, float T, float V, int offset, byte random_index, int flag_chessboard);
   
void cpu_update();
   
//#define cpu_METROPOLIS
#define gpu_METROPOLIS

__global__ static void _rand(float* vec, uint* z1, uint* z2, uint* z3, uint* z4);
__device__ static float HybridTaus(uint& z1, uint& z2, uint& z3, uint& z4);

void compute_indices(int *a,int *b);

curandState* devStates;

/* van der Waals interacion */
float * distance;
int * frequency;
float distance_PBC(int x1, int y1, int z1, int x2, int y2, int z2);
int num_distances;
float R_cutoff, LJ_InfiniteBarrier=0;

void volume_step(double * r_cell, double * energy, double * volume);
void calculate_energy_vdW (double * energy);

/* sigma variables storage */
byte * s, * dev_s, * nhb, * dev_nhb;
double beta;

/* eta variables storage */
byte * active_bond, * dev_active_bond;
int * rand_Chessboard;

#define ACTIVE 1
#define NON_ACTIVE 0

/* Neighbour index arrays, to be filled at the beginning  */
int * xup, * yup, * zup, 
    * xdn, * ydn, * zdn, 
    * dev_xup, * dev_yup, * dev_zup,
    * dev_xdn, * dev_ydn, * dev_zdn;
float * dev_Rnd, * Rnd;
float * dev_Rnd_spin, * Rnd_spin;
float * dev_Rnd_array, * Rnd_array;
float * acc, * dev_acc;

void update(void);
void measure(int iter, FILE *f);

uint * seed_array, * dev_seed_array;

double r_cell = 1.0001;
double energy_vdW = 0;  //van der Waals energy

double E = 0, E2 = 0, M = 0, M2 = 0;

// save configurations
//int flag_save_config;
//int * logaritmic_times;
//int max_power, num_of_times;
//bool is_logaritmic_time(int);

/* autocorrelation function calculation */

int flag_correlation;

#define numcorrelators 32
#define p 16
#define m_ 2

// correlator
struct correlator {
	double ** shift;
	double ** correlation;
	long int ** ncorrelation;
	double * accumulator;
	int * naccumulator;
	int * insertindex;
	double * t;
	double * f;
	int npcorr;
	int kmax;
	double accval;
};

void initialize(int N, struct correlator mol[]);
void add(struct correlator mol[], int i, double w, int k);
void evaluate(struct correlator mol[], int N, int norm);
void printCorr(struct correlator corr[], char filename[],int total_correlators);
//void HistogramTest(struct correlator mol[], int N, int norm);

struct correlator *Corr;
struct correlator *corr;

int *histo_sigma;

/* Wolff CPU cluster algorithm */
//#define WOLFF
double Jeff;
void cluster_step();
int site (int arm, int cell);
void shuffle ();
int cluster_poke ();
int add_to_cluster (int arm, int cell);
void update_cluster (int arm, int cell);
int * order, * neigh_arm, * neighbor;
bool * is_cluster;
int cluster_size;
float pJ_hb, pJ_s;

/* Swendsen-Wang GPU cluster algorithm */
#define SWENDSEN_WANG
#define SW_links_per_cell 18  
   // xup, yup, zup + 15 sigma bonds = 18
   // in each call to gpu cluster kernels, a cell works on 18 links between spins
float * dev_Rnd_cluster, * dev_Rnd_update;
byte * dev_delta, * cpu_delta, * dev_converges, *converges;
int * dev_label, * cpu_label, *dev_prev_label;

__global__ void gpu_initialize_cluster_variables (byte *dev_delta, int *dev_label, int *dev_prev_label);
__global__ void gpu_create_cluster( float* dev_Rnd_cluster, byte* dev_delta, byte* dev_s, 
                                    byte * dev_active_bond, float pJ_hb, float pJ_s, float Jeff);
__global__ void gpu_cluster_scanning_covalent (byte * dev_delta, int * dev_label);
__global__ void gpu_cluster_scanning_sigma (byte * dev_delta, int * dev_label);
__global__ void gpu_cluster_analysis(int *dev_label, int spin);
__global__ void gpu_update_cluster(int * dev_label, byte * dev_s, int offset_update, float * dev_Rnd_update);
__global__ void gpu_convergence_test(int * dev_label, int * dev_prev_label, byte * converges);

int offset_update, loop_update;
bool cluster_converges;
int max_SW_scans=-1;
int max_cluster_size, max_cluster_size2, tmp_max_cluster_size, num_clusters;
int *SW_cluster_size;
float average_cluster_size;

/* Swendsen-Wang CPU algorithm */
//#define cpu_SWENDSEN_WANG
int ** SW_bonded;
int * cpu_SW_L, * cpu_SW_N, * neighbor_spin, * SW_new_spin, * HK_label;// * cpu_SW_L2, * cpu_SW_L3;

void cpu_SW_initialize();
void cpu_SW_neighboring_spins(int i);
void cpu_SW_update(int* label);

//void dfs_rec ( int u, int this_label );
//int *visited;

/* chessboard / checkerboard algorithm: update of eta variables */
int flag_chessboard;

int * devChessboard_state, * devChessboard_edge, * devChessboard_vertex, * devSolutions;
__global__ void gpu_set_cubes(int * devChessboard_edge, int * devChessboard_vertex);

__global__ void gpu_Chessboard_start_cubes(int N, float * dev_Rnd_cubeFlip,
 int * devSolutions, byte * dev_active_bond, int * devChessboard_edge,
 int * devChessboard_vertex, int * devChessboard_state );
 
__global__ void gpu_Chessboard_set_state( int N, byte * dev_active_bond,
 int * devChessboard_vertex, int * devChessboard_state );

__global__ void gpu_Chessboard_set_active_bonds (byte * dev_s, int N, double Jeff, double beta,
     float * rand0, float * rand1, int offset_chess,
     int * devSolutions, byte * dev_active_bond, int * devChessboard_edge,
     int * devChessboard_vertex, int * devChessboard_state);

float * dev_Rnd_chess, * dev_Rnd_cubeFlip;
int offset_chess, loop_chess;

//cluster statisitcs

void print_cluster_statistics(int time, FILE * logsize, int* label);
int * cluster_size_;
byte * Xcoord, *Ycoord, *Zcoord; 
int is_percolating_cluster(int *label, int i);


//////////////////////
///////  MAIN
//////////////////////

int main ( int argc, char * argv[] ) 
{
  int i,j,n_loops=0,loop=0,offset=0,n_equilibrium=0;
  uint64_t seed;
  char input_data[200];
  char output_data[200];
  char input_config[200];
  char output_config[200];
  char log_file[200];
  FILE *idata, *odata, *iconfig, *oconfig, *logfile;

  int total_global_memory = 0;

  /* Read the input */
     sprintf(input_data,  "input_data");
     idata = fopen(input_data, "r");
     if(idata == NULL){
       fprintf(stderr,"Can't open input file\n");
       exit(1);
     }

     char parameter[50];
     
     if (fscanf(idata,"%s %f", &parameter [0], &P) != 2){
          fprintf(stderr,"Error [Pressure]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %f", &parameter [0], &T) != 2){
          fprintf(stderr,"Error [Temperature]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %ld", &parameter [0], &seed) != 2){
          fprintf(stderr,"Error[Seed]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &NMETROPOLIS) != 2){
          fprintf(stderr,"Error[Metropolis_Steps]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &NCLUSTER) != 2){
          fprintf(stderr,"Error[Cluster_Steps]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &NEQUILIBRIUM) != 2){
          fprintf(stderr,"Error[Equilibration_Steps]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &NLOOPS) != 2){
          fprintf(stderr,"Error[Sampling_Steps]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &N_MEASURE) != 2){
          fprintf(stderr,"Error[Sampling_Interval]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %f", &parameter [0], &R_cutoff) != 2){
          fprintf(stderr,"Error[R_cutoff]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &flag_config) != 2){
          fprintf(stderr,"Error[Flag_config]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &flag_distances) != 2){
          fprintf(stderr,"Error[Flag_distances]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &flag_correlation) != 2){
          fprintf(stderr,"Error[Flag_correlation]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &flag_chessboard) != 2){
          fprintf(stderr,"Error[Flag_chessboard]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     if (fscanf(idata,"%s %d", &parameter [0], &deviceID) != 2){
          fprintf(stderr,"Error[Device_ID]: Can't read input file\n");
          fclose(idata);
          exit(1);
     }

     fclose(idata);

  beta = 1./T;
  n_loops = NLOOPS;
  n_equilibrium = NEQUILIBRIUM;
  Jeff = J_hb - P*v_hb;

#ifdef CPU
  sprintf(output_data,  "test_cpu.out");
#endif
#ifdef GPU
  sprintf(output_data,  "data.out");
  sprintf(output_config, "config.out");
  sprintf(input_config, "input_config");
#endif

/* Verify that only one, cpu or gpu Metropolis, is enabled */

#ifdef gpu_METROPOLIS
  #ifdef cpu_METROPOLIS
       fprintf(stderr,"ERROR!! Both Metropolis CPU and GPU algorithms are active!!\n");
       fprintf(stderr,"Please check source file\n");
       exit(1);
  #endif
#endif

/* Verify that only one, Wolff, cpu SW, or gpu SW, is enabled */
#ifdef WOLFF
   #ifdef SWENDSEN_WANG
       fprintf(stderr,"ERROR!! Both Swendsen-Wang and Wolff cluster algorithms are active!!\n");
       fprintf(stderr,"Please check source file\n");
       exit(1);
   #endif
   #ifdef cpu_SWENDSEN_WANG
       fprintf(stderr,"ERROR!! Both Swendsen-Wang and Wolff cluster algorithms are active!!\n");
       fprintf(stderr,"Please check source file\n");
       exit(1);
   #endif
#endif

#ifdef cpu_SWENDSEN_WANG
   #ifdef SWENDSEN_WANG
       fprintf(stderr,"ERROR!! Both Swendsen-Wang CPU and GPU algorithms are active!!\n");
       fprintf(stderr,"Please check source file\n");
       exit(1);
   #endif
#endif

  odata = fopen(output_data,"w");       /* open output file */
  if(odata == NULL){
     fprintf(stderr,"Can't create output file\n");
     exit(1);
  }

  sprintf(log_file,  "water3D.log");     /* open log file */
  logfile = fopen(log_file, "w");  
  if(logfile == NULL){
     fprintf(stderr,"Can't create log file\n");
     exit(1);
  }

  fprintf(logfile, "++++ Welcome to water3D log file ++++\n");

  seed_mersenne( seed );    //initialize mersenne
  initRand( seed );         //initialize xorsifht1024*Phi
  
  if ((NX*NY*NZ) %  NTHREADS != 0){
    fprintf(stderr,"Error, NX*NY*NZ %d is not a multiple of NTHREADS %d \n",NX*NY*NZ, NTHREADS);
    return 1;
  }

#ifdef GPU

    int total_devices;  /* set GPU device */
    CudaSafeCall(cudaGetDeviceCount(&total_devices));
    if ( deviceID < 0 || deviceID >=total_devices ){
       fprintf(stderr,"Error: Tried to execute in GPU deive ID = %d\n",deviceID);
       fprintf(stderr,"\tTotal number of devices is %d\n",total_devices);
       return 1;
    }
    CudaSafeCall(cudaSetDevice(deviceID));

#endif
    int maxMB = 5000;     // Avoid excess of global memory storage (device). The limit is hardware dependent and may be modified.

  int total_RNG_memory = (5*N_RANDOM+SW_links_per_cell)*NX*NY*NZ;  
           //4 vectors of sizes N_RANDOM (Metropolis & Chessboard) + 1 vector of size 18 (Swendsen-Wang)
  
  fprintf(logfile,"RNG allocated memory: %d MB\n", (total_RNG_memory*sizeof(float))/(1024*1024));
  
  if((total_RNG_memory*sizeof(float))/(1024*1024) > maxMB) {
    fprintf(stderr,"Error, memory allocated > %d MB.\n",maxMB);
    return 1;
  }
  
  // NX*NY*NZ of the lattice and initial volume
  V = NX*NY*NZ*r_cell*r_cell*r_cell;

  fprintf(logfile," ++++++++++++++++++++++++++++++++++++++++++\n");
  fprintf(logfile,"Source file: %s\n",__FILE__);
  fprintf(logfile," 3D water model, %d x %d x %d lattice, N = %d, P = %1.2f, T = %1.2f\n",NX,NY,NZ,NX*NY*NZ,P,T);
  fprintf(logfile," J_hb = %1.2f, J_sigma = %1.2f, v_HB = %1.2f\n",J_hb,J_sig,v_hb);
  fprintf(logfile," Flag configuration = %d\n",flag_config);
  fprintf(logfile," Monte Carlo Step: %d Metropolis + %d Cluster\n",NMETROPOLIS,NCLUSTER);
#ifdef gpu_METROPOLIS
  fprintf(logfile," Metropolis: Parallel GPU algorithm\n");
#endif
#ifdef cpu_METROPOLIS
  fprintf(logfile," Metropolis: Sequential CPU algorithm\n");
#endif
#ifdef WOLFF
  fprintf(logfile," Cluster: Wolff Recursive CPU algorithm\n");
#endif
#ifdef SWENDSEN_WANG
  fprintf(logfile," Cluster: Swendsen-Wang Parallel GPU algorithm\n");
#endif
#ifdef cpu_SWENDSEN_WANG
  fprintf(logfile," Cluster: Swendsen-Wang Sequential CPU algorithm\n");
#endif
  fprintf(logfile," Device = %d\n",deviceID);
  fprintf(logfile," %d equilibration, %d updates/measurement, %d updates\n",
	 n_equilibrium, N_MEASURE, n_loops);
  fprintf(logfile," Output file %s\n",output_data);
  fprintf(logfile," Random seed %ld\n", seed );

  int dist_index=0;      //index to build distance and frequency
  num_distances= (int) 50;    //guess value for the total number of distances
  distance = (float*)calloc(num_distances,sizeof(float));
  frequency = (int*)calloc(num_distances,sizeof(int));

#ifdef GPU
  int nthreads = NTHREADS;
  int nblocks = NX*NY*NZ/nthreads; /// Synchronous blocks

  if (nthreads > 1024){   /* Hardware dependent. May be modified */
     fprintf(stderr,"Excessive number of threads. Reduce NTHREADS.\n");
     exit(-1);
  }
  if (nblocks > 65000){   //* Hardware dependent. May be modified */
     fprintf(stderr,"Excessive number of blocks. Increase the number of threads.\n");
     exit(-1);
  }
#endif  

  s = (byte*)calloc(ARMS*NX*NY*NZ,sizeof(byte));	// sigma variables
  
  nhb = (byte*)calloc(NX*NY*NZ,sizeof(byte));		// controls nhb per cell (if chessboard OFF)
 
  seed_array = (uint*)calloc((NX*NY*NZ),sizeof(uint));  // array of seed for random number generator

  active_bond = (byte*)calloc(ARMS*NX*NY*NZ,sizeof(byte));  // eta variables

  srand(seed);  
  
  for (i=0; i<(NX*NY*NZ); ++i) {
    if(i == 0)
      seed_array[i] = rand();
    else
      seed_array[i] = rand() ^ seed_array[i-1];
  }

#ifdef GPU  

  /* Allocate arrays in the GPU */

  total_global_memory += (NX*NY*NZ)*sizeof(curandState);
  CudaSafeCall(cudaMalloc( (void**)&devStates, (NX*NY*NZ)*sizeof(curandState)));

  total_global_memory += (NX*NY*NZ)*sizeof(uint);
  CudaSafeCall(cudaMalloc( (void**)&dev_seed_array, (NX*NY*NZ)*sizeof(uint)));

          cudaError_t err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Allocated memory in GPU: %d MB\n",total_global_memory/(1024*1024));
             fprintf(stderr, "Failed to launch cudaMalloc (part 1) (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }

  total_global_memory += (1+SW_links_per_cell+2*ARMS)*(NX*NY*NZ)*sizeof(byte);
  CudaSafeCall(cudaMalloc( (void**)&dev_s, ARMS*NX*NY*NZ*sizeof(byte) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_active_bond, ARMS*NX*NY*NZ*sizeof(byte) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_nhb, NX*NY*NZ*sizeof(byte) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_delta, SW_links_per_cell*NX*NY*NZ*sizeof(byte) ));
 
          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Allocated memory in GPU: %d MB\n",total_global_memory/(1024*1024));
             fprintf(stderr, "Failed to launch cudaMalloc (part 2) (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }

  total_global_memory += (NX*NY*NZ)*(2*ARMS*sizeof(int) + sizeof(byte));
  CudaSafeCall(cudaMalloc( (void**)&dev_label,ARMS*NX*NY*NZ*sizeof(int) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_prev_label,ARMS*NX*NY*NZ*sizeof(int) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_converges,NX*NY*NZ*sizeof(byte) ));

  total_global_memory += (5*N_RANDOM+SW_links_per_cell)*NX*NY*NZ*sizeof(float);
  CudaSafeCall(cudaMalloc( (void**)&dev_Rnd, (NX*NY*NZ)*N_RANDOM*sizeof(float) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_Rnd_spin, (NX*NY*NZ)*N_RANDOM*sizeof(float) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_Rnd_cubeFlip, (NX*NY*NZ)*N_RANDOM*sizeof(float) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_Rnd_chess, (NX*NY*NZ)*N_RANDOM*sizeof(float) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_Rnd_cluster, (NX*NY*NZ)*SW_links_per_cell*sizeof(float) ));
  CudaSafeCall(cudaMalloc( (void**)&dev_Rnd_update, (NX*NY*NZ)*N_RANDOM*sizeof(float) ));

          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Allocated memory in GPU: %d MB\n",total_global_memory/(1024*1024));
             fprintf(stderr, "Failed to launch cudaMalloc (part 3) (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }


  CudaSafeCall(cudaMemcpy( dev_seed_array, seed_array, (NX*NY*NZ)*sizeof(uint), cudaMemcpyHostToDevice ));

          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Failed to launch cudaMemcpy dev_seed_array (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }

  gpu_RNG_setup <<<nblocks,nthreads>>> (devStates, dev_seed_array);

          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Failed to launch gpu_RNG_setup (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }

  
  gpu_RNG_generate <<<nblocks,nthreads>>> ( devStates, dev_Rnd, N_RANDOM);
  gpu_RNG_generate <<<nblocks,nthreads>>> ( devStates, dev_Rnd_spin, N_RANDOM);
  gpu_RNG_generate <<<nblocks,nthreads>>> ( devStates, dev_Rnd_cubeFlip, N_RANDOM);
  gpu_RNG_generate <<<nblocks,nthreads>>> ( devStates, dev_Rnd_chess, N_RANDOM);
  gpu_RNG_generate <<<nblocks,nthreads>>> ( devStates, dev_Rnd_cluster, SW_links_per_cell);
  gpu_RNG_generate <<<nblocks,nthreads>>> ( devStates, dev_Rnd_update, N_RANDOM);

          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Failed to launch gpu_RNG_generate (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }
#endif

////
////  CHESSBOARD ALGORITHM
////

  if (flag_chessboard == 1){

      fprintf(logfile, "Starting Chessboard\n");

#ifdef CPU
      rand_Chessboard = (int*) calloc((NX*NY*NZ)/4,sizeof(int));

      initialize_Chessboard (NX, NY, NZ);

      for (int i =0; i<(NX*NY*NZ)/4; i++)
        rand_Chessboard[i] = (int) (NSOLUTIONS*mersenne());

      Chessboard_set_active_bonds_Typewriter(0, s, ARMS, NX*NY*NZ, Jeff , beta, rand_Chessboard, active_bond);
#endif

#ifdef GPU
      if (NX != NY || NX != NZ || NY != NZ ) { fprintf(stderr,"Error: can't build chessboard.\n\tThe system is not a cube.\n"); exit(1);} 
		//Maybe it could be solved for rectangles

      if ( NX % 4 != 0 ) { fprintf(stderr,"Error: can't build chessboard.\n\tUse NX muliply of 4.\n"); exit(1);} 

      cudaError_t err = cudaSuccess;

      total_global_memory += (1+3*EDGES)*(NX*NY*NZ/4)*sizeof(int);

      CudaSafeCall(cudaMalloc( (void**)&devChessboard_vertex, EDGES*((NX*NY*NZ)/4)*2*sizeof(int) ));
      CudaSafeCall(cudaMalloc( (void**)&devChessboard_edge, EDGES*((NX*NY*NZ)/4)*sizeof(int) ));
      CudaSafeCall(cudaMalloc( (void**)&devChessboard_state, ((NX*NY*NZ)/4)*sizeof(int) ));

          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Allocated memory in GPU: %d MB\n",total_global_memory/(1024*1024));
             fprintf(stderr, "Failed to launch cudaMalloc devChessboard_* (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }

      gpu_set_cubes <<<nblocks,nthreads>>> (devChessboard_edge, devChessboard_vertex);
          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Failed to launch gpu_set_cubes kernel (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }

      total_global_memory += EDGES*NSOLUTIONS*sizeof(int);
      CudaSafeCall(cudaMalloc( (void**)&devSolutions, EDGES*NSOLUTIONS*sizeof(int)));

          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Allocated memory in GPU: %d MB\n",total_global_memory/(1024*1024));
             fprintf(stderr, "Failed to launch cudaMalloc devSolutions (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }

      int * cpuSolutions;
      cpuSolutions = (int*) malloc(EDGES*NSOLUTIONS*sizeof(int));

      setSolutions_Typewriter(cpuSolutions); // fill vector in host, then copy to device

      CudaSafeCall(cudaMemcpy( devSolutions, cpuSolutions, EDGES*NSOLUTIONS*sizeof(int), cudaMemcpyHostToDevice ));
          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Failed to launch cudaMemcpy devSolutions HostToDevice (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }
    
      gpu_Chessboard_start_cubes <<<nblocks,nthreads>>> (NX*NY*NZ, dev_Rnd_cubeFlip, devSolutions, dev_active_bond, 
          devChessboard_edge, devChessboard_vertex, devChessboard_state);

          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Failed to launch gpu_Chessboard_start_cubes kernel (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }

      CudaSafeCall(cudaMemcpy( active_bond, dev_active_bond, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));

#endif
  } else {   //Chessboard algorithm is not active
     for (int i=0; i<ARMS*NX*NY*NZ; i++)
         active_bond[i] = 1;
  }
  
////
////  GENERATE INITIAL CONFIGURATION
////

  if (flag_config == 1){ /* Initial spins - read from file */
     iconfig = fopen(input_config,"r");
     if (iconfig == NULL){
        fprintf(stderr,"Can't open the initial configuration file\n");
        fprintf(stderr,"File named 'input_config' is expected\n");
        exit(1);
     }

     byte sp0, sp1, sp2, sp3, sp4, sp5;
     int ab0, ab1, ab2, ab3, ab4, ab5; //variables to read from file

     int i=0;
     while(fscanf(iconfig, "%d %d %d %d %d %d %d %d %d %d %d %d",&sp0,&sp1,&sp2,&sp3,&sp4,&sp5,
		&ab0,&ab1,&ab2,&ab3,&ab4,&ab5) == 12)  {

       s[0*NX*NY*NZ + i] = sp0;
       s[1*NX*NY*NZ + i] = sp1;
       s[2*NX*NY*NZ + i] = sp2;
       s[3*NX*NY*NZ + i] = sp3;
       s[4*NX*NY*NZ + i] = sp4;
       s[5*NX*NY*NZ + i] = sp5;

       active_bond[0*NX*NY*NZ + i] = ab0;
       active_bond[1*NX*NY*NZ + i] = ab1;
       active_bond[2*NX*NY*NZ + i] = ab2;
       active_bond[3*NX*NY*NZ + i] = ab3;
       active_bond[4*NX*NY*NZ + i] = ab4;
       active_bond[5*NX*NY*NZ + i] = ab5;

       for (int j=0; j<ARMS; j++){

         if (s[j*NX*NY*NZ + i] < 0 || s[j*NX*NY*NZ + i] > 5){
    	     fprintf(stderr,"Initial configuration file error\n");
             fprintf(stderr,"Unexpected value for spin (%d) appeared in line %d\n", s[j*NX*NY*NZ + i], i);
             fclose(iconfig);
             exit(1);
         }

         if (active_bond[j*NX*NY*NZ + i] < 0 || active_bond[j*NX*NY*NZ + i] > 1){
    	     fprintf(stderr,"Initial configuration file error\n");
             fprintf(stderr,"Unexpected value for active_bond (%d) appeared in line %d\n", active_bond[j*NX*NY*NZ + i], i);
             fclose(iconfig);
             exit(1);
         }
       } 

       i++;
       if (i > (NX*NY*NZ)){
         fprintf(stderr,"Initial configuration file error\n");
         fprintf(stderr,"There are more rows than cells in the file\n");
         fclose(iconfig);
         exit(1);
       }
  
     }
       
     if ( i < (NX*NY*NZ) ){
        fprintf(stderr,"Initial configuration file error\n");
        fprintf(stderr,"Unexpected ending. Can't fill the spin array\n");
        fclose(iconfig);
        exit(1);
     }
     fclose(iconfig);
     
     if ( flag_chessboard == 1 ){  //change cube states according to input active bonds
#ifdef GPU
     CudaSafeCall(cudaMemcpy( dev_active_bond, active_bond, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyHostToDevice ));

  //   cudaError_t err = cudaSuccess;

     gpu_Chessboard_set_state <<<nblocks,nthreads>>> ( NX*NY*NZ, dev_active_bond, devChessboard_vertex, devChessboard_state ); 

          err = cudaGetLastError();
          if (err != cudaSuccess)
          {
             fprintf(stderr, "Failed to launch gpu_set_state kernel (error code: %s)!\n", cudaGetErrorString(err));
             exit(EXIT_FAILURE);
          }

#endif
     }

  } else if (flag_config == 0) {   /* Initial spins - completely disordered */

    
    for (i=0; i<NX*NY*NZ; ++i)
      for (j=0; j<ARMS; ++j)
        s[j*NX*NY*NZ + i] = (byte)(q*mersenne());

  } else if (flag_config == 2) {  // Initial spins - completely ordered

    byte initial_spin = (byte)(q*mersenne());

    for (i=0; i<NX*NY*NZ; ++i)
      for (j=0; j<ARMS; ++j)
        s[j*NX*NY*NZ + i] = initial_spin;

  } else {
    fprintf(stderr,"Error: Wrong value for flag_config = %d\n",flag_config);
    fprintf(stderr,"Flag_config=0  Initial spins completely disordered\n");
    fprintf(stderr,"Flag_config=1  Initial spins read from file\n");
    fprintf(stderr,"Flag_config=2  Initial spins completely ordered (T=0)\n");
    exit(1);
  }

  int x,y,z;
  for (i=0; i < (NX*NY*NZ); ++i) {  //setting nhb

    x = i % NX;
    y = (i / NX)%NY;
    z = i / (NX*NY);

    int xup = (x+1)%NX + y*NX + z*NX*NY;
    int yup = x + ((y+1)%NY)*NX + z*NX*NY;
    int zup = x + y*NX + ((z+1)%NZ)*NX*NY;
    
    int xdn = (x-1+NX)%NX + y*NX + z*NX*NY;
    int ydn = x + ((y-1+NY)%NY)*NX + z*NX*NY;
    int zdn = x + y*NX + ((z-1+NZ)%NZ)*NX*NY;
    
    byte acc = 0;
    
    acc += active_bond[0*NX*NY*NZ + i] * delta(s[0*NX*NY*NZ + i],s[1*NX*NY*NZ + xup]);
    acc += active_bond[1*NX*NY*NZ + i] * delta(s[1*NX*NY*NZ + i],s[0*NX*NY*NZ + xdn]);
    
    acc += active_bond[2*NX*NY*NZ + i] * delta(s[2*NX*NY*NZ + i],s[3*NX*NY*NZ + yup]);
    acc += active_bond[3*NX*NY*NZ + i] * delta(s[3*NX*NY*NZ + i],s[2*NX*NY*NZ + ydn]);
    
    acc += active_bond[4*NX*NY*NZ + i] * delta(s[4*NX*NY*NZ + i],s[5*NX*NY*NZ + zup]);
    acc += active_bond[5*NX*NY*NZ + i] * delta(s[5*NX*NY*NZ + i],s[4*NX*NY*NZ + zdn]);
    
    nhb[i] = acc;
  }

  //setting distance and frequency
  if ( flag_distances == 0 ){
       x = 0 % NX;       // x, y, z for molecule 0
       y = (0 / NX)%NY;
       z = 0 / (NX*NY);

       for (int j=1; j < (NX*NY*NZ); ++j){   //calculate distances between cells and molecule 0
         int x2, y2, z2;
         float dist;
         int flag;

         x2 = j % NX;
         y2 = (j / NX)%NY;
         z2 = j / (NX*NY);

         dist = distance_PBC(x,y,z,x2,y2,z2);

         if (dist < R_cutoff) {
 
           flag = 1;
 
           for (int k=0; k<dist_index; ++k)
              if (dist == distance[k]){
                 flag = 0;
                 frequency[k] ++;
                 break;
              }

           if (flag){

              distance[dist_index] = dist;
              frequency[dist_index] = 1;
              dist_index ++;

              if (dist_index == num_distances){    //reallocation of distance and frequency
                 num_distances += int( cbrt( float(NX*NY*NZ) ) );
                 distance = (float*)realloc(distance,(num_distances)*sizeof(float));
                 frequency = (int*)realloc(frequency,(num_distances)*sizeof(int));
                 for (int l = dist_index; l<num_distances; ++l){
                     distance[l] = 0;
                     frequency[l] = 0;
                 }
              }
  
           }
         }  //end if dist>R_cutoff
       }  //end for j=i+1

       for(int l=0; l<dist_index; ++l){  // final result for full lattice
            frequency[l] /= 2;
            frequency[l] *= (NX*NY*NZ);
       }

      // writing into file
      char dist_file[200];
      FILE *odist;

      sprintf(dist_file,  "distances");
      odist = fopen(dist_file, "w");
      if(odist == NULL){
         fprintf(stderr,"Can't create distances file\n");
         exit(1);
      }

      fprintf(odist,"%d %d %d\n",NX,NY,NZ);
      for( int i=0; i<dist_index; ++i)
         fprintf(odist,"%f %d\n",distance[i],frequency[i]);

      fclose(odist);

  } else {     // reading distances from file

    char dist_file[200];
    FILE *idist;
    sprintf(dist_file,  "distances");
    idist = fopen(dist_file, "r");
    if(idist == NULL){
       fprintf(stderr,"Can't open distances file\n");
       exit(1);
    }

    int LX, LY, LZ;

    if (fscanf(idist,"%d %d %d", &LX, &LY, &LZ) != 3){
        fprintf(stderr,"Error: Can't read distances file\n");
        fclose(idist);
        exit(1);
    }


    if (LX != NX){
        fprintf(stderr,"Error: Distances file. X rank mismatch\n");
        fclose(idist);
        exit(1);
    }

    if ( LY != NY){
        fprintf(stderr,"Error: Distances file. Y rank mismatch\n");
        fclose(idist);
        exit(1);
    }

    if ( LZ!= NZ){
        fprintf(stderr,"Error: Distances file. Z rank mismatch\n");
        fclose(idist);
        exit(1);
    }

    dist_index=0;
    while(fscanf(idist,"%f %d", &distance[dist_index], &frequency[dist_index])==2){
       dist_index++;
       if (dist_index == num_distances){    //reallocation of distance and frequency
                 num_distances += int( cbrt( float(NX*NY*NZ) ) );
                 distance = (float*)realloc(distance,(num_distances)*sizeof(float));
                 frequency = (int*)realloc(frequency,(num_distances)*sizeof(int));
                 for (int l = dist_index; l<num_distances; ++l){
                     distance[l] = 0;
                     frequency[l] = 0;
                 }
       }
    }
  }

  num_distances = dist_index;
  distance = (float*)realloc(distance,num_distances*sizeof(float));  //free memory
  frequency = (int*)realloc(frequency,num_distances*sizeof(int));

// Cluster Probabilities (Wolff & Swendsen-Wang)
  pJ_hb = 1-exp(-beta*fabs(Jeff));   // Jeff = J_hb - P*v_hb where J_hb > 0 is attractive
  pJ_s = 1-exp(-beta*fabs(J_sig));   // J_sig > 0 is attractive
                                     // E_HB = -J_hb*N_HB - J_sig*N_sig

// Vectors for cluster statistics
  cluster_size_ = (int*) calloc(ARMS*NX*NY*NZ,sizeof(int));
  Xcoord = (byte*) calloc(NX,sizeof(byte));
  Ycoord = (byte*) calloc(NY,sizeof(byte));
  Zcoord = (byte*) calloc(NZ,sizeof(byte));
  
// initialize Wolff cluster algorithm variables
#ifdef WOLFF
  order = (int*)calloc(ARMS,sizeof(int));
  for ( i=0; i<ARMS; i++ )
      order[i] = i;

  neigh_arm = (int*)calloc(ARMS,sizeof(int));
  neigh_arm[0] = 1;
  neigh_arm[1] = 0;
  neigh_arm[2] = 3;
  neigh_arm[3] = 2;
  neigh_arm[4] = 5;
  neigh_arm[5] = 4;

  is_cluster = (bool*) calloc(ARMS*NX*NY*NZ,sizeof(bool));
  for(i =0; i<ARMS*NX*NY*NZ; i++)
    is_cluster[i] = false;

  neighbor = (int*) calloc(ARMS*NX*NY*NZ,sizeof(int));  
  for(i=0; i<ARMS*NX*NY*NZ; i++){
     int arm = i/(NX*NY*NZ);
     int cell = i%(NX*NY*NZ);

     int x = cell % NX;
     int y = (cell/NX)%NY;
     int z = cell/(NX*NY);

    switch(arm){
      case 0: // forward X
        neighbor[i] = (x+1)%NX + y*NX + z*NX*NY;
        break;
      case 1: // backward X
        neighbor[i] = (x-1+NX)%NX + y*NX + z*NX*NY;
        break;
      case 2: // forward Y
        neighbor[i] = x + ((y+1)%NY)*NX + z*NX*NY;
        break;
      case 3: // backward Y
        neighbor[i] = x + ((y - 1 + NY)%NY)*NX + z*NX*NY;
        break;
      case 4: // forward Z
        neighbor[i] = x + y*NX + ((z+1)%NZ)*NX*NY;
        break;
      case 5: // backward Z
        neighbor[i] = x + y*NX + ((z-1+NZ)%NZ)*NX*NY;
        break;
    }
     
  }

  cluster_size = 0;
#endif

#ifdef SWENDSEN_WANG
  cpu_label = (int*) calloc(ARMS*NX*NY*NZ,sizeof(int));  
  //prev_label = (int*) calloc(ARMS*NX*NY*NZ,sizeof(int));

  cpu_delta = (byte*) calloc(SW_links_per_cell*NX*NY*NZ,sizeof(byte));
  converges = (byte*) calloc(NX*NY*NZ,sizeof(byte));
  //device SW variables are allocated elsewhere and initialize at the beginning of each step with a GPU kernel
  //check if cuda Atomic function could be used for checking convergence instead of this method

#endif

#ifdef cpu_SWENDSEN_WANG

  SW_bonded = (int**) calloc(ARMS*NX*NY*NZ,sizeof(int*));
  for (int i=0; i<ARMS*NX*NY*NZ; i++)
    SW_bonded[i] = (int*) calloc(ARMS,sizeof(int));
   
  //Hoshen Kopelman labels
  cpu_SW_L = (int*) calloc (ARMS*NX*NY*NZ,sizeof(int));
  cpu_SW_N = (int*) calloc (ARMS*NX*NY*NZ,sizeof(int));
  HK_label = (int*) calloc (ARMS*NX*NY*NZ,sizeof(int));
  
 // cpu_SW_L2 = (int*) calloc (ARMS*NX*NY*NZ,sizeof(int));
  
 // cpu_SW_L3 = (int*) calloc (ARMS*NX*NY*NZ,sizeof(int));
 // visited = (int*) calloc (ARMS*NX*NY*NZ,sizeof(int));
  
  neighbor_spin = (int*) calloc (ARMS,sizeof(int));
  SW_new_spin = (int*) calloc (ARMS*NX*NY*NZ,sizeof(int));

#endif

////
////  ALLOCATE AND INITIALZE CORRELATORS
////

// Observable M = fraction of spins at the most populated state

  int total_correlators = 1;    // change to NX*NY*NZ to recover correlation previous to 3.1
  if ( flag_correlation == 1 ){

   Corr = (struct correlator*) malloc(sizeof(struct correlator)*total_correlators);

   for(i = 0; i < total_correlators; i ++) {
      Corr[i].shift = (double **)malloc(sizeof(double *)*numcorrelators);
      for(j = 0; j < numcorrelators; j ++)
          Corr[i].shift[j] = (double *)calloc(p,sizeof(double));
       
      Corr[i].correlation = (double **)malloc(sizeof(double *)*numcorrelators);
      for(j = 0; j < numcorrelators; j ++)
          Corr[i].correlation[j] = (double *)calloc(p,sizeof(double));
       
      Corr[i].ncorrelation = (long int **)malloc(sizeof(long int *)*numcorrelators);
      for(j = 0; j < numcorrelators; j ++)
          Corr[i].ncorrelation[j] = (long int *)calloc(p,sizeof(long int));
        
      Corr[i].accumulator = (double *)calloc(numcorrelators,sizeof(double));
      Corr[i].naccumulator = (int *)calloc(numcorrelators,sizeof(int));
      Corr[i].insertindex = (int *)calloc(numcorrelators,sizeof(int));
      Corr[i].t = (double *)calloc(numcorrelators*p,sizeof(double));
      Corr[i].f = (double *)calloc(numcorrelators*p,sizeof(double));
    }

   initialize(total_correlators,Corr);

   histo_sigma = (int*) calloc(q,sizeof(int));
   
  }

////
////  FIRST CALCULATION OF VDW ENERGY
////

  calculate_energy_vdW ( &energy_vdW );
  
////
////  COPY SYSTEM FROM HOST TO DEVICE
////

#ifdef GPU  
  CudaSafeCall(cudaMemcpy( dev_s, s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyHostToDevice ));

  CudaSafeCall(cudaMemcpy( dev_active_bond, active_bond, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyHostToDevice ));

  CudaSafeCall(cudaMemcpy( dev_nhb, nhb, NX*NY*NZ*sizeof(byte), cudaMemcpyHostToDevice ));
#endif

  offset = 0;  //Metropolis random number
  loop = 0;
  offset_chess = 0; //Chessboard random number
  loop_chess = 0;
  offset_update = 0; //Swendsen-Wang random number
  loop_update = 0;

 /* {  // write first config to compare to config.out //
        FILE * oconfig_time;
        char output_config_time[200];
        sprintf(output_config_time,"config0");
        oconfig_time = fopen(output_config_time,"w");
        for (int k=0; k<NX*NY*NZ; ++k){
           fprintf(oconfig_time,"%d %d %d %d %d %d ",s[0*NX*NY*NZ + k],s[1*NX*NY*NZ + k],s[2*NX*NY*NZ + k],s[3*NX*NY*NZ + k],s[4*NX*NY*NZ + k],s[5*NX*NY*NZ + k]);
           fprintf(oconfig_time,"%d %d %d %d %d %d\n",active_bond[0*NX*NY*NZ + k],active_bond[1*NX*NY*NZ + k],active_bond[2*NX*NY*NZ + k],active_bond[3*NX*NY*NZ + k],active_bond[4*NX*NY*NZ + k],active_bond[5*NX*NY*NZ + k]);
        }
        fclose(oconfig_time);
  } */

  max_cluster_size = 0;
  max_cluster_size2 = 0;
  SW_cluster_size = (int*) calloc(ARMS*NX*NY*NZ,sizeof(int));

  FILE *logsize = fopen("cluster.out", "w");
  if(logsize == NULL){
    fprintf(stderr,"Can't open input file\n");
    exit(1);
  }
  
////
////  EQULIIBRATION STEPS
////
  
  /* the update/equilibrate loop */
  for (i=0; i<n_equilibrium; ++i){
  
   //  printf("time=%d\n",i);

 /*    clock_t start, end;
     double cpu_time_used;
    
     if(i%2000 == 0)
        start = clock();*/

    MonteCarlo_Step(i, loop, offset, loop_chess, offset_chess, nblocks, nthreads, order, logfile, SW_cluster_size, logsize);

  /*   if((i+1)%2000 == 0){
         end = clock(); 
         cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
         printf("Equil. step %d, lasted %f\n",i,cpu_time_used);
     }*/

  }

  max_cluster_size = 0;
  max_cluster_size2 = 0;
  
////
////  PRODUCTION RUN
////
          
  /* and the update/measure loop */
  for (i=0; i<n_loops; ++i) {

 //   if(i % (n_loops/10) == 0)
 //       fprintf(logfile,"%.1f%% elapsed\n",10*i/(n_loops/10.));
    

  /*   clock_t start, end;
     double cpu_time_used;
    
     if(i%2000 == 0)
        start = clock(); */
 
 //printf("time=%d\n",i);

    MonteCarlo_Step(i+n_equilibrium, loop, offset, loop_chess, offset_chess, nblocks, nthreads, order, logfile, SW_cluster_size,logsize);

   /*  if((i+1)%2000 == 0){
         end = clock(); 
         cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
         printf("Sampl. step %d, lasted %f\n",i,cpu_time_used); 

     }*/
    
    if (i % N_MEASURE == 0){
#ifdef GPU
      CudaSafeCall(cudaMemcpy( s, dev_s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
      CudaSafeCall(cudaMemcpy( active_bond, dev_active_bond, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
      CudaSafeCall(cudaMemcpy( nhb, dev_nhb, NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
#endif

#ifdef SWENDSEN_WANG 
     CudaSafeCall(cudaMemcpy( cpu_label, dev_label, ARMS*NX*NY*NZ*sizeof(int), cudaMemcpyDeviceToHost ));
     print_cluster_statistics(i,logsize, cpu_label);
#endif

#ifdef cpu_SWENDSEN_WANG
      print_cluster_statistics(i,logsize, HK_label);
#endif

      measure(i,odata);
    }

    if ( flag_correlation == 1 ){
#ifdef GPU
      CudaSafeCall(cudaMemcpy( s, dev_s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
#endif
      for (int k=0;k<total_correlators; k++){
//        double sigma = (double) (s[0*NX*NY*NZ + k]+s[1*NX*NY*NZ + k]+s[2*NX*NY*NZ + k]+s[3*NX*NY*NZ + k]+s[4*NX*NY*NZ + k]+s[5*NX*NY*NZ + k]);
//        sigma /= 6.0;
        for (int l=0; l<q; l++)
           histo_sigma[l] = 0;

        for (int l=0; l<ARMS*NX*NY*NZ; l++){

           if ( s[l] < 0 || s[l] >= q ) { 
              fprintf(stderr,"ERROR while correlation estimate: s[%d]=%d out of range\n",
                            l,s[l]);
              exit(1);
           }
           histo_sigma[s[l]] ++;
        }

        double max_mass = -2E10;
        for (int l=0; l<q; l++)
          if (histo_sigma[l] > max_mass) max_mass = histo_sigma[l];

        max_mass /= (double) (ARMS*NX*NY*NZ);

        add(Corr,k,max_mass,0);
      }
    }
  } //end of main loop

  fclose(odata);
  fclose(logsize);
  
////
////  STORE DATA AND FINAL CALCULATIONS
////

  M /= n_loops*1./N_MEASURE;
  M2 /= n_loops*1./N_MEASURE;
  E /= n_loops*1./N_MEASURE;
  E2 /= n_loops*1./N_MEASURE;

  double dM = sqrt(M2 - M*M);
  double dE = sqrt(E2 - E*E);

  fprintf(logfile,"Magnetization:\t%1.6lf\t+/-\t%1.6lf\n", M, dM);
  fprintf(logfile,"Energy:\t\t%1.6lf\t+/-\t%1.6lf\n", E, dE);
  
  /* Final spins - write output file */
  oconfig = fopen(output_config,"w");
#ifdef GPU
  CudaSafeCall(cudaMemcpy( s, dev_s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
  CudaSafeCall(cudaMemcpy( active_bond, dev_active_bond, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
#endif
  for (i=0; i<NX*NY*NZ; ++i){
      fprintf(oconfig,"%d %d %d %d %d %d ",s[0*NX*NY*NZ + i],s[1*NX*NY*NZ + i],s[2*NX*NY*NZ + i],s[3*NX*NY*NZ + i],s[4*NX*NY*NZ + i],s[5*NX*NY*NZ + i]);
      fprintf(oconfig,"%d %d %d %d %d %d\n",active_bond[0*NX*NY*NZ + i],active_bond[1*NX*NY*NZ + i],active_bond[2*NX*NY*NZ + i],active_bond[3*NX*NY*NZ + i],active_bond[4*NX*NY*NZ + i],active_bond[5*NX*NY*NZ + i]);
  }
  fclose(oconfig);

  /* Correlation evaluation */
  if (flag_correlation == 1) {
    evaluate(Corr,total_correlators,1);

    char filename[100];
    sprintf(filename,"correlation.dat");
    printCorr(Corr,filename,total_correlators);
  }

  fprintf(logfile,"Max_SW_scans = %d\n",max_SW_scans);
  fprintf(logfile," ** simulation done\n");
  fclose(logfile);
  
  return 0;
}

float distance_PBC(int x1, int y1, int z1, int x2, int y2, int z2){
   int x_PBC, y_PBC, z_PBC;

   x_PBC = abs(x1 - x2);
   if ( x_PBC > NX/2 ) x_PBC = NX - x_PBC;

   y_PBC = abs(y1 - y2);
   if ( y_PBC > NY/2 ) y_PBC = NY - y_PBC;

   z_PBC = abs(z1 - z2);
   if ( z_PBC > NZ/2 ) z_PBC = NZ - z_PBC;

   return sqrt(x_PBC*x_PBC + y_PBC*y_PBC + z_PBC*z_PBC);

}

////
////  MONTE CARLO STEP
////	Update eta variables (if chessboard enabled)
////	Update r_cell
////	Update sigma variables (Metropolis, SW, or Wolff)
////

void MonteCarlo_Step(int time, int &loop, int &offset, int &loop_chess, int &offset_chess, int nblocks, int nthreads, int *order, FILE* logfile, int *SW_cluster_size, FILE *logsize){

    if (flag_chessboard == 1){

#ifdef CPU
    	for (int i =0; i<(NX*NY*NZ)/4; i++)
   	   rand_Chessboard[i] = (int) (9*mersenne());  //warning: I haven't checked if mersenne() returns [0,1) or (0,1]
 
	Chessboard_set_active_bonds_Typewriter(1,s, ARMS, NX*NY*NZ, Jeff , beta, rand_Chessboard, active_bond);
#endif

#ifdef GPU
        offset_chess = ((NX*NY*NZ)/4)*loop_chess;  

        cudaError_t err2 = cudaSuccess;

        gpu_Chessboard_set_active_bonds <<<nblocks,nthreads>>> (dev_s, NX*NY*NZ, Jeff, beta,
                                       dev_Rnd_cubeFlip, dev_Rnd_chess, offset_chess, devSolutions, dev_active_bond,
                                       devChessboard_edge, devChessboard_vertex, devChessboard_state);
          err2 = cudaGetLastError();
          if (err2 != cudaSuccess)
          {
             fprintf(stderr, "Failed to launch gpu_Chessboard_set_active_bonds kernel (error code: %s)!\n", cudaGetErrorString(err2));
             exit(EXIT_FAILURE);
          }

        loop_chess ++;

        if(loop_chess == 4*N_RANDOM) { 
          loop_chess = 0;
          gpu_RNG_generate <<<nblocks,nthreads>>> (devStates, dev_Rnd_cubeFlip, N_RANDOM);
          gpu_RNG_generate <<<nblocks,nthreads>>> (devStates, dev_Rnd_chess, N_RANDOM);
        } 

        CudaSafeCall(cudaMemcpy( active_bond, dev_active_bond,ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
#endif
    }

    volume_step(&r_cell, &energy_vdW, &V);

#ifdef gpu_METROPOLIS
#ifdef GPU
    cudaError_t err3 = cudaSuccess;
    for (int j=0; j<NMETROPOLIS; j++){

      for(int k=0; k < 2*ARMS; ++k) {

        byte spin = (byte) ARMS*mersenne();
        int indices = (int) 2*mersenne();

        offset = (NX*NY*NZ)*loop;

        gpu_update <<<nblocks,nthreads>>> (dev_s,dev_active_bond,dev_nhb,dev_Rnd,
                                           dev_Rnd_spin,indices,P, T, V, offset, spin,flag_chessboard);
          err3 = cudaGetLastError();
          if (err3 != cudaSuccess)
          {
             fprintf(stderr, "Failed to launch gpu_update kernel (error code: %s)!\n", cudaGetErrorString(err3));
             exit(EXIT_FAILURE);
          }

        ++loop;

        if(loop == N_RANDOM) {
          loop = 0;
          gpu_RNG_generate <<<nblocks,nthreads>>> ( devStates, dev_Rnd, N_RANDOM);
          gpu_RNG_generate <<<nblocks,nthreads>>> ( devStates, dev_Rnd_spin, N_RANDOM);
        }
      }
    }
#endif
#endif

#ifdef cpu_METROPOLIS
    for (int j=0; j<NMETROPOLIS; j++)
       cpu_update();
       
  #ifdef GPU
        CudaSafeCall(cudaMemcpy( dev_s, s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyHostToDevice ));
  #endif
#endif

#ifdef WOLFF
    if (NCLUSTER > 0){
#ifdef GPU
        CudaSafeCall(cudaMemcpy( s, dev_s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost )); 
        CudaSafeCall(cudaMemcpy( active_bond, dev_active_bond,ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
#endif
        for (int j=0; j<NCLUSTER; j++)
           cluster_step();

     //   fprintf(logsize,"%d %d %f %d\n",time,tmp_max_cluster_size,average_cluster_size,num_clusters);
#ifdef GPU
        CudaSafeCall(cudaMemcpy( dev_s, s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyHostToDevice ));
#endif
    }
#endif  //end ifdef Wolff

#ifdef SWENDSEN_WANG
#ifdef GPU
    cudaError_t err4;
    if (NCLUSTER > 0){

        for (int j=0; j<NCLUSTER; j++){

            tmp_max_cluster_size = 0;

//            printf("Inicio paso SW\n");

            gpu_initialize_cluster_variables <<<nblocks,nthreads>>> (dev_delta, dev_label, dev_prev_label);
                err4 = cudaGetLastError();
                if (err4 != cudaSuccess)
                {
                   fprintf(stderr, "Failed to launch gpu_initialize_cluster_variables kernel (error code: %s)!\n", cudaGetErrorString(err4));
                   exit(EXIT_FAILURE);
                }

       /*     CudaSafeCall(cudaMemcpy( active_bond, dev_active_bond, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
            fprintf(logfile,"Initial ab[]\n");
            for (int jj=0; jj<ARMS*NX*NY*NZ; jj++)
              fprintf(logfile,"\tI ab[%d] = %d\n",jj,active_bond[jj]);

*/
            gpu_create_cluster <<<nblocks,nthreads>>> (dev_Rnd_cluster, dev_delta, dev_s, dev_active_bond, pJ_hb, pJ_s, Jeff);
                err4 = cudaGetLastError();
                if (err4 != cudaSuccess)
                {
                   fprintf(stderr, "Failed to launch gpu_create_cluster kernel (error code: %s)!\n", cudaGetErrorString(err4));
                   exit(EXIT_FAILURE);
                }

            cluster_converges = false;
            int loop_counter=0;

            while ( ! cluster_converges ){

                loop_counter ++;

    //            fprintf(logfile,"Step %d Scanning\n",loop_counter); 

                gpu_cluster_scanning_covalent <<<nblocks,nthreads>>> (dev_delta, dev_label);
                   err4 = cudaGetLastError();
                   if (err4 != cudaSuccess)
                   {
                      fprintf(stderr, "Failed to launch gpu_cluster_scanning kernel (error code: %s)!\n",
                                      cudaGetErrorString(err4));
                      exit(EXIT_FAILURE);
                   }

      /*          CudaSafeCall(cudaMemcpy( cpu_label, dev_label, ARMS*NX*NY*NZ*sizeof(int), cudaMemcpyDeviceToHost ));
                fprintf(logfile,"Scanning_Covalent Step %d Label[]\n",loop_counter);
                for (int jj=0; jj<ARMS*NX*NY*NZ; jj++)
                   fprintf(logfile,"\tScanning_Covalent %d Label[%d] = %d\n",loop_counter,jj,cpu_label[jj]);
        */    

                gpu_cluster_scanning_sigma <<<nblocks,nthreads>>> (dev_delta, dev_label);
                   err4 = cudaGetLastError();
                   if (err4 != cudaSuccess)
                   {
                      fprintf(stderr, "Failed to launch gpu_cluster_scanning kernel (error code: %s)!\n",
                                      cudaGetErrorString(err4));
                      exit(EXIT_FAILURE);
                   }

          /*      CudaSafeCall(cudaMemcpy( cpu_label, dev_label, ARMS*NX*NY*NZ*sizeof(int), cudaMemcpyDeviceToHost ));
                fprintf(logfile,"Scanning_Sigma Step %d Label[]\n",loop_counter);
                for (int jj=0; jj<ARMS*NX*NY*NZ; jj++)
                   fprintf(logfile,"\tScanning_Sigma %d Label[%d] = %d\n",loop_counter,jj,cpu_label[jj]);
            */

                for (int ii=0; ii<ARMS; ii++){

                   gpu_cluster_analysis <<<nblocks,nthreads>>> (dev_label, ii);
                      err4 = cudaGetLastError();
                      if (err4 != cudaSuccess)
                      {
                         fprintf(stderr, "Failed to launch gpu_cluster_analysis kernel (error code: %s)!\n",
                                         cudaGetErrorString(err4));
                         exit(EXIT_FAILURE);
                      }

                }

              /*  CudaSafeCall(cudaMemcpy( cpu_label, dev_label, ARMS*NX*NY*NZ*sizeof(int), cudaMemcpyDeviceToHost ));
                fprintf(logfile,"Analysis Step %d Label[]\n",loop_counter);
                for (int jj=0; jj<ARMS*NX*NY*NZ; jj++)
                   fprintf(logfile,"\tAnalysis %d Label[%d] = %d\n",loop_counter,jj,cpu_label[jj]);
            */

                /// CONVERGENCE TEST 
                gpu_convergence_test <<<nblocks,nthreads>>> (dev_label, dev_prev_label, dev_converges);
                   err4 = cudaGetLastError();
                   if (err4 != cudaSuccess)
                   {
                      fprintf(stderr, "Failed to launch gpu_convergence_test kernel (error code: %s)!\n",
                                      cudaGetErrorString(err4));
                      exit(EXIT_FAILURE);
                   }
                CudaSafeCall(cudaMemcpy( converges, dev_converges, NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));

              //  fprintf(logfile,"Test Convergence Step %d\n",loop_counter);
                cluster_converges = true;
                for (int ii=0; ii<NX*NY*NZ; ii++)
                    if ( converges[ii] != 0 ) {
                //       fprintf(logfile,"\tCell %d has not converged. Conv=%d\n",ii,converges[ii]);
                       cluster_converges = false;
                       break;
                    }

            }   

        //    fprintf(logfile,"Cluster has converged\n");
            loop_counter --;
            if (loop_counter > max_SW_scans) max_SW_scans = loop_counter;

            offset_update = (NX*NY*NZ)*loop_update;

            gpu_update_cluster <<<nblocks,nthreads>>> (dev_label, dev_s, offset_update, dev_Rnd_update);
                err4 = cudaGetLastError();
                if (err4 != cudaSuccess)
                {
                   fprintf(stderr, "Failed to launch gpu_update_cluster kernel (error code: %s)!\n", cudaGetErrorString(err4));
                   exit(EXIT_FAILURE);
                }


          /*  CudaSafeCall(cudaMemcpy( cpu_label, dev_label, ARMS*NX*NY*NZ*sizeof(int), cudaMemcpyDeviceToHost ));
	    fprintf(logfile,"Final Label[]\n");
            for (int jj=0; jj<ARMS*NX*NY*NZ; jj++)
		fprintf(logfile,"\tF Label[%d] = %d\n",jj,cpu_label[jj]);

            CudaSafeCall(cudaMemcpy( s, dev_s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
            fprintf(logfile,"Final s[]\n");
            for (int jj=0; jj<ARMS*NX*NY*NZ; jj++)
              fprintf(logfile,"\tF s[%d] = %d\n",jj,s[jj]);
*/
            loop_update ++;

            if (loop_update == N_RANDOM){
                loop_update = 0;
                gpu_RNG_generate <<<nblocks,nthreads>>> (devStates, dev_Rnd_update, N_RANDOM);
            }

            gpu_RNG_generate <<<nblocks,nthreads>>> (devStates, dev_Rnd_cluster, SW_links_per_cell);
                err4 = cudaGetLastError();
                if (err4 != cudaSuccess)
                {
                   fprintf(stderr, "Failed to launch gpu_RNG_generate SW cluster kernel (error code: %s)!\n", cudaGetErrorString(err4));
                   exit(EXIT_FAILURE);
                }

       //     printf("End SW step\n");

        }       
    }
        
#endif
#endif //end ifdef SWENDSEN_WANG

#ifdef cpu_SWENDSEN_WANG
#ifdef GPU
        CudaSafeCall(cudaMemcpy( s, dev_s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost )); 
        CudaSafeCall(cudaMemcpy( active_bond, dev_active_bond,ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyDeviceToHost ));
#endif

   if (NCLUSTER > 0){

        for (int j=0; j<NCLUSTER; j++){

          cpu_SW_initialize();
          
          for ( int k=NX*NY*NZ; k<ARMS*NX*NY*NZ; k++){  //create matrix of bonds
          
             cpu_SW_neighboring_spins(k);  // set vector neighbor_spin[]
             
             // covalent bonds
             if ( neighbor_spin[0] < k ){  //avoid double attempt create same bond
             
               if ( Jeff > 0 ){ 
               
                  if ( active_bond[k] != active_bond[neighbor_spin[0]] ){
                     fprintf(stderr,"cpu_SW detecta un error en active_bond on el calculo del vecino.\n");
                     exit(1);
                  }
               
                  if ( active_bond[k] ){
                        if ( s[k]  == s[neighbor_spin[0]] ){
                            if ( mersenne() < pJ_hb )
                                SW_bonded[k][0] = 1; 
                        }
                      
                  }
               } else {
                  if ( s[k]  != s[neighbor_spin[0]] ){
                        if ( mersenne() < pJ_hb )
                            SW_bonded[k][0] = 1; 
                  }      
               }
               
             }
               
             // cooperative bonds
             for ( int kk=1; kk<ARMS; kk++){
                if ( neighbor_spin[kk] < k ){ //avoid double attempt create same bond
                  if ( s[k] == s[neighbor_spin[kk]] ){
                        if ( mersenne() < pJ_s )
                            SW_bonded[k][kk] = 1; 
                  }             
                }
               
             }
    
          } // end matrix of bonds 
       
          
          //Hoshen Kopelman Algorithm
          int cpu_SW_converges = 0;
          while ( cpu_SW_converges == 0 ){
          
             cpu_SW_converges = 1;
              
             for ( int k=NX*NY*NZ; k<ARMS*NX*NY*NZ; k++){
    
                 cpu_SW_neighboring_spins(k);
               
                 int nbonds = 0;
                 int min_label = 2*ARMS*NX*NY*NZ;
                 for ( int kk=0 ; kk<ARMS; kk++){    // estimate min_label
               
                    if ( neighbor_spin[kk] < k )
                      if ( SW_bonded[k][kk] == 1 ){
                         nbonds ++;
                       
                         if ( cpu_SW_L[neighbor_spin[kk]] < min_label )
                            min_label = cpu_SW_L[neighbor_spin[kk]];
                      }
                 }
               
                 if ( nbonds >= 1 )
                     cpu_SW_L[k] = min_label;
                 
              
                 if ( nbonds > 1 ){
              
                     int new_N = -2*ARMS*NX*NY*NZ;

                     for ( int kk=0 ; kk<ARMS; kk++){  //find max N among bonded spins
               
                       if ( neighbor_spin[kk] < k )
                         if ( SW_bonded[k][kk] == 1 ){  
                             if ( new_N < cpu_SW_N[cpu_SW_L[neighbor_spin[kk]]] )
                                new_N = cpu_SW_N[cpu_SW_L[neighbor_spin[kk]]];
                                
                         }
                    
                     }
                  
                     for ( int kk=0 ; kk<ARMS; kk++){
               
                       if ( neighbor_spin[kk] < k )
                         if ( SW_bonded[k][kk] == 1 ){
                         
                           if ( cpu_SW_N[ cpu_SW_L[neighbor_spin[kk]] ] != new_N){ 
                              cpu_SW_converges = 0;
                              cpu_SW_N[ cpu_SW_L[neighbor_spin[kk]] ] = new_N;   
                           }
                           
                         }
                        
                     }
                 }             
             }
         
          }  //end Hoshen Kopelman
          
          cpu_SW_update(HK_label);

        }
   }
   
   
#ifdef GPU
        CudaSafeCall(cudaMemcpy( dev_s, s, ARMS*NX*NY*NZ*sizeof(byte), cudaMemcpyHostToDevice ));
#endif

#endif //end ifdef cpu_SWENDSEN_WANG

}// end subroutine Monte Carlo Step


////
////  VDW VOLUME STEP
////

void volume_step(double * r_cell, double * energy, double * volume)
{

    double trial_r_cell;
    double volume_step_epsilon = 0.01;
    
    do {
        trial_r_cell = (*r_cell) + volume_step_epsilon* 2 *( uniformDoubleRand() - 0.5 );
    } while ( trial_r_cell < LJ_InfiniteBarrier );
    
    double trial_volume = trial_r_cell*trial_r_cell*trial_r_cell * NX*NY*NZ;
    double dV = trial_volume - (*volume);         
    
    double trial_energy = 0, tmp, attractive, repulsive, shift;
    
    // shift: correction on the potential to avoid step at the cutoff
    tmp = 1.0/(R_cutoff);
    tmp = tmp * tmp * tmp;
    attractive = tmp * tmp;
    repulsive = attractive * attractive;
    shift = (repulsive - attractive);

    for ( int i = 0 ; i < num_distances ; ++i )  
      if ( trial_r_cell * distance[i] < R_cutoff ){
          tmp = 1.0/( trial_r_cell * distance[i] );
          tmp = tmp * tmp * tmp;
          attractive = tmp * tmp;
          repulsive = attractive * attractive;
          trial_energy += frequency[i] * (repulsive - attractive - shift);
      }

    
    double dE = ( trial_energy - (*energy) ) - NX*NY*NZ * T * log( trial_volume/(*volume) );

    if ( uniformDoubleRand() < exp( - beta *( dE + P* dV )) )
    {
        * r_cell = trial_r_cell;
        * volume = trial_volume;
        * energy = trial_energy;
    }

}


void calculate_energy_vdW (double * energy){

    double tmp, attractive, repulsive, ret=0, shift;
    
    // shift: correction on the potential to avoid step at the cutoff
    tmp = 1.0/(R_cutoff);
    tmp = tmp * tmp * tmp;
    attractive = tmp * tmp;
    repulsive = attractive * attractive;
    shift = (repulsive - attractive);

    for ( int i = 0 ; i < num_distances ; ++i )  //energy of one cell
						 //energy of the lattice when multiplying by frequency[i]
       if ( r_cell * distance[i] < R_cutoff ){
           tmp = 1.0/( r_cell * distance[i] );
           tmp = tmp * tmp * tmp;
           attractive = tmp * tmp;
           repulsive = attractive * attractive;
           ret += frequency[i] * (repulsive - attractive - shift);
       }
     * energy = ret;
}

////
////  CPU METROPOLIS ALGORITHM
////

void cpu_update() 
{
   
  int this_site, cell, arm, neigh_cell, neighbor_spin,x,y,z;
  int new_spin, old_spin;
  float dE;
  for (int j=0; j < ARMS*NX*NY*NZ; ++j) {
  
    this_site = (int) (mersenne()*ARMS*NX*NY*NZ);
    
    old_spin = s[this_site];
    new_spin = (int) (q*mersenne());
    
    if ( old_spin == new_spin ) continue;
    
    cell = this_site%(NX*NY*NZ);
    arm = this_site/(NX*NY*NZ);
    
    x = cell % NX;
    y = (cell / NX) % NY;
    z = cell / (NX*NY);
    
    switch (arm){
   
      case 0:
        neigh_cell = (x+1)%NX + y*NX + z*NX*NY;
        neighbor_spin = 1*NX*NY*NZ+neigh_cell;
        break;
        
      case 1:
        neigh_cell = (x-1+NX)%NX + y*NX + z*NX*NY;
        neighbor_spin = 0*NX*NY*NZ+neigh_cell;
        break;
        
      case 2:
        neigh_cell = x + ((y+1)%NY)*NX + z*NX*NY;
        neighbor_spin = 3*NX*NY*NZ+neigh_cell;
        break;
        
      case 3:
        neigh_cell = x + ((y-1+NY)%NY)*NX + z*NX*NY;
        neighbor_spin = 2*NX*NY*NZ+neigh_cell;
        break;
        
      case 4:
        neigh_cell = x + y*NX + ((z+1)%NZ)*NX*NZ;
        neighbor_spin = 5*NX*NY*NZ+neigh_cell;
        break;
        
      case 5:
        neigh_cell = x + y*NX + ((z-1+NZ)%NZ)*NX*NY;
        neighbor_spin = 4*NX*NY*NZ+neigh_cell;
        break; 
        
      default:
        break;  
    }
    
    dE = - Jeff * active_bond[this_site] * (delta(new_spin,s[neighbor_spin]) -  delta(old_spin,s[neighbor_spin]));
    
    for ( int k=0; k<ARMS; k++ ){
       if ( k == arm ) continue; 
       neighbor_spin = k*NX*NY*NZ + cell;
       dE -= J_sig * ( delta(new_spin,s[neighbor_spin]) - delta(old_spin,s[neighbor_spin]) );
    }
    
    if ( dE <= 0 )
       s[this_site] = new_spin;
    else if ( mersenne() < exp(-beta*dE) )
       s[this_site] = new_spin;
      
  }   

}

////
////  GPU METROPOLIS ALGORITHM
////

__global__ void gpu_update (
  byte * dev_s,
  byte * dev_active_bond,
  byte * dev_nhb,
  float * dev_Rnd,
  float * dev_Rnd_spin,
  int indices,
  float P,
  float T,
  float V,
  int offset,
  byte spin_index,
  int flag_chessboard)
{
  /* 1D grid of 1D blocks */
    uint tid = blockIdx.x*blockDim.x + threadIdx.x;
	
    uint x = tid % NX;
    uint y = (tid/NX) % NY;
    uint z = tid/(NX*NY);

    // water model algorithm
      
    if(spin_index == 0 || spin_index == 1) { // update Y-Z planes
      if(indices && x%2)
        return;
      else if(!indices && !(x%2))
        return;
    } else if(spin_index == 2 || spin_index == 3) { // update X-Z planes
      if(indices && y%2)
        return;
      else if(!indices && !(y%2))
        return;
    } else if(spin_index == 4 || spin_index == 5) { // update X-Y planes
      if(indices && z%2)
        return;
      else if(!indices && !(z%2))
        return;
    }
        
    /// attempt to coalesced reads and writes
    // coalesced read
    uint neighbor;
    byte neig_spin;
    
    //  NO warp divergence in this block, spin_index is the same for all threads
    switch(spin_index){
      case 0: // forward X
        neighbor = (x+1)%NX + y*NX + z*NX*NY;
        neig_spin = dev_s[1*NX*NY*NZ + neighbor];
        break;
      case 1: // backward X
        neighbor = (x-1+NX)%NX + y*NX + z*NX*NY;
        neig_spin = dev_s[0*NX*NY*NZ + neighbor];
        break;
      case 2: // forward Y
        neighbor = x + ((y+1)%NY)*NX + z*NX*NY;
        // nearly coalesced
        neig_spin = dev_s[3*NX*NY*NZ + neighbor];
        break;
      case 3: // backward Y
        neighbor = x + ((y - 1 + NY)%NY)*NX + z*NX*NY;
        neig_spin = dev_s[2*NX*NY*NZ + neighbor];
        break;
      case 4: // forward Z
        neighbor = x + y*NX + ((z+1)%NZ)*NX*NY;
        neig_spin = dev_s[5*NX*NY*NZ + neighbor];
        break;
      case 5: // backward Z
        neighbor = x + y*NX + ((z-1+NZ)%NZ)*NX*NY;
        neig_spin = dev_s[4*NX*NY*NZ + neighbor];
        break;
    }
    
    byte old_spin = dev_s[spin_index*NX*NY*NZ + tid];
    
    // attempt spin change
    //  byte new_spin = (old_spin + 1 + (byte)((q-1)*dev_Rnd_spin[tid + offset])) % q;
      byte new_spin = (byte) (q*dev_Rnd_spin[tid + offset]);
      new_spin %= q;    // avoid event new_spin = 6

    // compute new energy
    int8_t deltaNHB = (int8_t)(delta(new_spin,neig_spin) - delta(old_spin,neig_spin));
    
    byte nhbi = dev_nhb[tid];
    byte nhbj = dev_nhb[neighbor];
    
    float deltaE = 0;
    if ( flag_chessboard == 1 ){
       deltaE += delta(new_spin,neig_spin) * dev_active_bond[spin_index*NX*NY*NZ + tid];
       deltaE -= delta(old_spin,neig_spin) * dev_active_bond[spin_index*NX*NY*NZ + tid];
    }else {
       deltaE += (theta(nhbi+deltaNHB,NHBMAX) * theta(nhbj+deltaNHB,NHBMAX)) * delta(new_spin,neig_spin);
       deltaE -= (theta(nhbi,NHBMAX)          * theta(nhbj,NHBMAX))          * delta(old_spin,neig_spin);
    }
        
    deltaE *= (-J_hb + P*v_hb);
    
    for (uint i=0; i < ARMS; ++i)
    {
	byte spin = dev_s[i*NX*NY*NZ + tid];
        // coalesced reads
        deltaE += -J_sig * (delta(new_spin,spin) - delta(old_spin,spin));
    }
    
    // avoid double counting of same spin
    deltaE -= J_sig;    //warning: if (new_spin == old_spin) deltaE = -J_sig, but it should be deltaE = 0
			//This mistake has no effect in simulation because it doesn't change resulting dev_s[] and dev_nhb[]    

    // metropolis acceptance rule
    if (dev_Rnd[tid + offset] > expf(-deltaE/T)) {
      new_spin = old_spin;
      deltaNHB = 0;
    }
        
    dev_s[spin_index*NX*NY*NZ + tid] = new_spin;
    dev_nhb[tid] += deltaNHB;
    dev_nhb[neighbor] += deltaNHB;
}

////
////  MEASUREMENT SUBROUTINE
////

void measure(int iter,FILE *f)
{
  /* measure energy and magnetization, and print out
   */
  int i,j,k,x,y,z,xup,yup,zup,neighbor,neig_spin;
  double e,h,m,v;
  double n_hb = 0;
  
  /* sum over the neighbour sites - typewriter fashion */
  e = m = 0;

  int count_nhb[NX*NY*NZ];
  for (i=0; i < (NX*NY*NZ); ++i)
    count_nhb[i] = 0;
  
  for (i=0; i < (NX*NY*NZ); ++i) {
    x = i % NX;
    y = (i / NX)%NY;
    z = i / (NX*NY);
    xup = (x+1)%NX + y*NX + z*NX*NY;
    yup = x + ((y+1)%NY)*NX + z*NX*NY;
    zup = x + y*NX + ((z+1)%NZ)*NX*NY;

    if (flag_chessboard == 1){

      // forward X
      neighbor = xup;
      neig_spin = s[1*NX*NY*NZ + neighbor];
      if ( delta(s[0*NX*NY*NZ + i], neig_spin) ){
         count_nhb[i] += active_bond[0*NX*NY*NZ + i];
         count_nhb[neighbor] += active_bond[1*NX*NY*NZ + neighbor];
      }

      // forward Y
      neighbor = yup;
      neig_spin = s[3*NX*NY*NZ + neighbor];
      if ( delta(s[2*NX*NY*NZ + i], neig_spin) ){
         count_nhb[i] += active_bond[2*NX*NY*NZ + i];
         count_nhb[neighbor] += active_bond[3*NX*NY*NZ + neighbor];
      }

      // forward Z
      neighbor = zup;
      neig_spin = s[5*NX*NY*NZ + neighbor];
      if ( delta(s[4*NX*NY*NZ + i], neig_spin) ){
         count_nhb[i] += active_bond[4*NX*NY*NZ + i];
         count_nhb[neighbor] += active_bond[5*NX*NY*NZ + neighbor];
      }

    }  // end count nhb with chessboard
    else {

      // forward X
      neighbor = xup;
      neig_spin = s[1*NX*NY*NZ + neighbor];
      if(delta(s[0*NX*NY*NZ + i], neig_spin)) {
        n_hb ++;
        if(count_nhb[i] < NHBMAX && count_nhb[neighbor] < NHBMAX) {
          count_nhb[i] ++;
          count_nhb[neighbor] ++;
        }
      }
    
      // forward Y
      neighbor = yup;  
      neig_spin = s[3*NX*NY*NZ + neighbor];
      if(delta(s[2*NX*NY*NZ + i], neig_spin)) {
        n_hb ++;
        if(count_nhb[i] < NHBMAX && count_nhb[neighbor] < NHBMAX) {
          count_nhb[i] ++;
          count_nhb[neighbor] ++;
        }
      }
      
      // forward Z
      neighbor = zup;
      neig_spin = s[5*NX*NY*NZ + neighbor];
      if(delta(s[4*NX*NY*NZ + i], neig_spin)) {
        n_hb ++;
        if(count_nhb[i] < NHBMAX && count_nhb[neighbor] < NHBMAX) {
          count_nhb[i] ++;
          count_nhb[neighbor] ++;
        }
   
      }
    }  // end count nhb without chessboard

    for(j=0; j < ARMS-1; ++j)
      for(k=j+1; k < ARMS; ++k) 
        m += delta(s[j*NX*NY*NZ + i], s[k*NX*NY*NZ + i]);
	
  } // end loop over cells

 // check count_nhb
 /* 
 for (i=0; i < (NX*NY*NZ); ++i) {
        if (count_nhb[i] < 0 || count_nhb[i] > NHBMAX ){
              printf("ERROR! nhb(%d)=%d\n",i, count_nhb[i]);
              exit(1);
        }
  }*/
  
  // check the distribution of sigma states across the systems
  /*
  if(iter%(NLOOPS/10)==0) {
    printf("iter %d\n",iter);
    printf("sigmas\n");
    for(j=0; j < q; ++j)
      printf("%d: %1.2f%%\n",j,sigma_dist[j]*100./(NX*NY*NZ*q));
    printf("hbs\n");
    for(j=0; j < 1+ARMS; ++j)
      printf("%d: %1.2f%%\n",j,nhb_dist[j]*100./(NX*NY*NZ));
    
    
    int check_nhb=0;
    for(i=0; i < NX*NY*NZ; ++i) {
      if(nhb[i] < 0)
        printf("warning! nhb(%d)=%d\n",i, nhb[i]);
      check_nhb += nhb[i];
    }  
    if(check_nhb != (int)(2*n_hb))
      printf("warning! check nhb: %.0f %d\n",2*n_hb,check_nhb);
  }*/

  calculate_energy_vdW ( &e );
  
  v = NX*NY*NZ*r_cell*r_cell*r_cell;

  n_hb = 0;  //now it is the total hb taking into account the restriction NHBMAX = 4
  for (i=0; i < (NX*NY*NZ); ++i){
    v += 0.5*count_nhb[i] * v_hb;
    e += 0.5*count_nhb[i] * (-J_hb);
    n_hb += 0.5*count_nhb[i];
  }

  e += (-J_sig)*m;  
  h = e + P*v;  // enthalpy

  m /= 1.*NX*NY*NZ*(ARMS*(ARMS-1)/2.);
  
  E += e;
  E2 += e*e;
  M += m;
  M2 += m*m;
  
  fprintf(f,"%d %1.7e %1.7e %1.7e %1.7e %1.7e %1.7e\n", iter, e, h, v, n_hb, m, r_cell);
}

/////////////////  SAVE CONFIGS TO CALCULATE CORRELATIONS   /////////////////////

/*
void set_logaritmic_tiimes()
{
        int max_power = 10;
        int num_of_times = 9*max_power + 1;
        int logaritmic_times[num_of_times];
        int power = 1;
        int k = 0;
        for(int i=0; i<max_power; i++){
                for(int j=1; j<=9; j++){
                        logaritmic_times[k] = j*power;
                        k++;
                }
                power *= 10;
        }
        logaritmic_times[k] = power;
}


bool is_logaritmic_time (int time)
{
        bool ret = false;
        int i = 0;
        while (!ret && i<num_of_times){
                if(time == logaritmic_times[i]) ret = true;
                i++;
        }
        return ret;
}
*/

////
////  WOLF CLUSTER ALGORITHM
////

void cluster_step()
{
   int cumm_cluster_size = 0;

   tmp_max_cluster_size = 0;
   num_clusters = 0;
   int cs;
   while (cumm_cluster_size < ARMS*NX*NY*NZ){
    //     cumm_cluster_size += cluster_poke();
         cs = cluster_poke();
         cumm_cluster_size += cs;
         num_clusters ++;
         if ( cs > tmp_max_cluster_size ) tmp_max_cluster_size = cs;
   }

   average_cluster_size = (float) cumm_cluster_size / (float) num_clusters;
  
}

int site (int arm, int cell)
{ 
   if (arm*NX*NY*NZ + cell >= NX*NY*NZ*6){
	printf("ERROR (1) SITE \n");
        printf("arm=%d/6 ; cell=%d/%d; res=%d/%d\n",arm, cell,NX*NY*NZ,arm*NX*NY*NZ + cell,NX*NY*NZ*6);
        exit(-1);
   } else if (arm*NX*NY*NZ + cell < 0){
	printf("ERROR (2) SITE \n");
        exit(-1);
   }

   return arm*NX*NY*NZ + cell; 
}

int cluster_poke ()
{

  int cell = (int) NX*NY*NZ*mersenne();
  int arm = (int) q*mersenne();
  int this_spin = site(arm,cell);

  byte old_spin = s[this_spin];

  shuffle();

  byte new_spin = (byte) ( (q-1)*mersenne() + 1 );  // integer between 1 and q-1
        // I will rotate the cluster in update_cluster by adding new_spin to old_spin

  if ( new_spin <=0 || new_spin >= q){
	fprintf(stderr,"ERROR Cluster: new_spin out of range\n ns=%d\n",new_spin);
        exit(-1);	
  }

  cluster_size = add_to_cluster( arm, cell );

  int neigh_spin = site(neigh_arm[arm], neighbor[this_spin]);
  if ( ! is_cluster[ neigh_spin ] ){

      if ( Jeff > 0 )
      {
         if (  ( active_bond[ neigh_spin ] == ACTIVE ) &&
               ( s[ this_spin ] == s[ neigh_spin ] ) &&
               ( mersenne() < pJ_hb )   ){

                cluster_size += add_to_cluster( neigh_arm[arm] , neighbor[this_spin] );
                update_cluster (neighbor[this_spin],new_spin);
         }
      }
      else if ( Jeff < 0 )
      {
         if (s[ this_spin ] != s[ neigh_spin ] &&
             mersenne() < pJ_hb ){
                cluster_size += add_to_cluster( neigh_arm[arm] , neighbor[this_spin] );
                update_cluster (neighbor[this_spin],new_spin);
         }          
      }

  }

  update_cluster(cell, new_spin);

  return cluster_size;
}

void shuffle()
{

    int temp, pos;
    for ( int i = 0 ; i < q ; i++ )
    {
        temp = order[i];
        pos = (int) q*mersenne();
        order[i] = order[pos];
        order[pos] = temp;
    }
}

int add_to_cluster( int arm, int cell ){

   int this_spin = site(arm,cell);  //seed spin of this (sub)cluster
 
   is_cluster[this_spin] = true;
   int local_size = 1;

   for (int i=0; i<6; i++){
      int next_arm = order[i];
      int next_spin = site(next_arm,cell);

      if ( ! is_cluster[next_spin] && s[next_spin] == s[this_spin] && mersenne() < pJ_s){

         local_size += add_to_cluster(next_arm, cell);

         if ( active_bond[ next_spin ] == NON_ACTIVE ) continue;

         int neigh_spin = site(neigh_arm[next_arm] , neighbor[next_spin] );
         if ( ! is_cluster[ neigh_spin ] ){

               if ( Jeff > 0 )
               {
                  if (s[ next_spin ] == s[ neigh_spin ] &&
                      mersenne() < pJ_hb )

                         local_size += add_to_cluster( neigh_arm[next_arm] , neighbor[next_spin] );
               }
               else if ( Jeff < 0 )
               {
                  if (s [next_spin ] != s[ neigh_spin ] &&
                      mersenne() < pJ_hb )
                         local_size += add_to_cluster( neigh_arm[next_arm] , neighbor[next_spin] );         
               }

         }
         
      }

   }

   return local_size;
}

void update_cluster( int cell, int new_spin ){

  for (int i=0; i<6; i++){
    int arm = order[i];
    int this_spin = site(arm,cell);

    if ( is_cluster[this_spin] )
    {
       s[this_spin] = ( s[this_spin] + new_spin ) % q;
       is_cluster[this_spin] = false;

       int neigh_spin = site(neigh_arm[arm] , neighbor[this_spin]);
       if ( is_cluster[ neigh_spin ] ){
          update_cluster(neighbor[this_spin],new_spin);
       }     
    }

  }

}

////
////  GPU SWENDSEN WANG CLUSTER ALGORITHM
////

// Follows algorithm described by Y. Komura and Y. Okabe in https://doi.org/10.1016/j.cpc.2012.01.017

__global__ void gpu_initialize_cluster_variables (byte *dev_delta, int *dev_label, int *dev_prev_label)
{

    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    int sp;

    for ( int i=0; i<SW_links_per_cell; i++ )
      dev_delta[i*NX*NY*NZ+tid] = 0;

    for ( int i=0; i<ARMS; i++ ){
      sp = i*NX*NY*NZ+tid;
      dev_label[sp] = sp;
      dev_prev_label[sp] = sp;
    }

}

__global__ void gpu_create_cluster ( float * dev_Rnd_cluster, byte * dev_delta, byte * dev_s,
   byte * dev_active_bond, float pJ_hb, float pJ_s, float Jeff)
{

    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    uint x = tid % NX;
    uint y = (tid/NX) % NY;
    uint z = tid/(NX*NY);

    uint xup = (x+1)%NX + y*NX + z*NX*NY;
    uint yup = x + ((y+1)%NY)*NX + z*NX*NY;
    uint zup =  x + y*NX + ((z+1)%NZ)*NX*NZ;


    if ( Jeff > 0 ) //no thread divergence, Jeff is the same for all threads
    {

      if ( dev_active_bond[0*NX*NY*NZ+tid] && delta(dev_s[0*NX*NY*NZ+tid],dev_s[1*NX*NY*NZ+xup]) &&
           dev_Rnd_cluster[0*NX*NY*NZ+tid] < pJ_hb )
        dev_delta[0*NX*NY*NZ+tid] = 1;

      if ( dev_active_bond[2*NX*NY*NZ+tid] && delta(dev_s[2*NX*NY*NZ+tid],dev_s[3*NX*NY*NZ+yup]) && 
           dev_Rnd_cluster[1*NX*NY*NZ+tid] < pJ_hb )
        dev_delta[1*NX*NY*NZ+tid] = 1;

      if ( dev_active_bond[4*NX*NY*NZ+tid] && delta(dev_s[4*NX*NY*NZ+tid],dev_s[5*NX*NY*NZ+zup]) &&
           dev_Rnd_cluster[2*NX*NY*NZ+tid] < pJ_hb )
        dev_delta[2*NX*NY*NZ+tid] = 1;

    } else if ( Jeff < 0) 
    {

      if ( dev_s[0*NX*NY*NZ+tid] != dev_s[1*NX*NY*NZ+xup] && dev_Rnd_cluster[0*NX*NY*NZ+tid] < pJ_hb )
        dev_delta[0*NX*NY*NZ+tid] = 1;

      if ( dev_s[2*NX*NY*NZ+tid] != dev_s[3*NX*NY*NZ+yup] && dev_Rnd_cluster[1*NX*NY*NZ+tid] < pJ_hb )
        dev_delta[1*NX*NY*NZ+tid] = 1;

      if ( dev_s[4*NX*NY*NZ+tid] != dev_s[5*NX*NY*NZ+zup] && dev_Rnd_cluster[2*NX*NY*NZ+tid] < pJ_hb )
        dev_delta[2*NX*NY*NZ+tid] = 1;

    }

    int link_counter = 3;
    for ( int ii=0; ii<ARMS-1; ii ++){
       for ( int jj=ii+1; jj<ARMS; jj++){

          if ( delta(dev_s[ii*NX*NY*NZ+tid],dev_s[jj*NX*NY*NZ+tid]) && dev_Rnd_cluster[link_counter*NX*NY*NZ+tid] < pJ_s )
             dev_delta[link_counter*NX*NY*NZ+tid] = 1;

          link_counter ++;
       }
    }

}
 
__global__ void gpu_cluster_scanning_covalent (byte * dev_delta, int * dev_label)
{

    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    uint x = tid % NX;
    uint y = (tid/NX) % NY;
    uint z = tid/(NX*NY);

    uint xup = (x+1)%NX + y*NX + z*NX*NY;
    uint yup = x + ((y+1)%NY)*NX + z*NX*NY;
    uint zup =  x + y*NX + ((z+1)%NZ)*NX*NY;

    //scan xup
    if( dev_delta[0*NX*NY*NZ+tid] ){
      if ( dev_label[0*NX*NY*NZ+tid] < dev_label[1*NX*NY*NZ+xup] )
         dev_label[1*NX*NY*NZ+xup] = dev_label[0*NX*NY*NZ+tid];
      else
         dev_label[0*NX*NY*NZ+tid] = dev_label[1*NX*NY*NZ+xup];
    }

    //scan yup
    if( dev_delta[1*NX*NY*NZ+tid] ){
      if ( dev_label[2*NX*NY*NZ+tid] < dev_label[3*NX*NY*NZ+yup] )
         dev_label[3*NX*NY*NZ+yup] = dev_label[2*NX*NY*NZ+tid];
      else
         dev_label[2*NX*NY*NZ+tid] = dev_label[3*NX*NY*NZ+yup];
    }

    //scan zup
    if( dev_delta[2*NX*NY*NZ+tid] ){
      if ( dev_label[4*NX*NY*NZ+tid] < dev_label[5*NX*NY*NZ+zup] )
         dev_label[5*NX*NY*NZ+zup] = dev_label[4*NX*NY*NZ+tid];
      else
         dev_label[4*NX*NY*NZ+tid] = dev_label[5*NX*NY*NZ+zup];
    }

}

__global__ void gpu_cluster_scanning_sigma (byte * dev_delta, int * dev_label)
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    //scan sigma bonds
    int link_counter = 3;
    for ( int ii=0; ii<ARMS-1; ii ++){
       for ( int jj=ii+1; jj<ARMS; jj++){

          if ( dev_delta[link_counter*NX*NY*NZ+tid] ){

             if ( dev_label[ii*NX*NY*NZ+tid] < dev_label[jj*NX*NY*NZ+tid])
                  dev_label[jj*NX*NY*NZ+tid] = dev_label[ii*NX*NY*NZ+tid];
             else
                  dev_label[ii*NX*NY*NZ+tid] = dev_label[jj*NX*NY*NZ+tid];

          }

          link_counter ++;
       }
    }

}

// This implementation suffers from race conditions. It should not be a problem as they
// will eventually be solved correctly in subsequently calls of the scanning and analysis kernels.
// Splitting analysis into six calls from sp=0 to sp=5 minimize the problem.
__global__ void gpu_cluster_analysis(int *dev_label, int sp)
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    int idx = sp*NX*NY*NZ+tid;

     if ( dev_label[dev_label[idx]] != dev_label[idx] )
        dev_label[idx] = dev_label[dev_label[idx]];

}

__global__ void gpu_convergence_test(int *dev_label, int * dev_prev_label, byte * dev_converges)
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    int sp;
    dev_converges[tid]=0;
    for (int k=0; k<ARMS; k++){
      sp = k*NX*NY*NZ+tid;
      if ( dev_prev_label[sp] != dev_label[sp] ) dev_converges[tid] = 1;
      dev_prev_label[sp] = dev_label[sp];  //update for the next step in the loop
    }

}

__global__ void gpu_update_cluster(int * dev_label, byte * dev_s, int offset_update, float * dev_Rnd_update)
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    int sp;
    for ( int arm=0; arm<ARMS; arm++ ){
       sp = arm*NX*NY*NZ+tid;
       dev_s[sp] += (byte) (q*dev_Rnd_update[dev_label[sp]+offset_update]);
       dev_s[sp] = dev_s[sp]%q;  
    }

}

////
////  SEQUENTIAL CPU SWENDSEN WANG CLUSTER ALGORITHM
////

void cpu_SW_initialize(){

   // matrix of bonds
   for (int i=0; i<ARMS*NX*NY*NZ; i++){
      for (int j=0; j<ARMS; j++)
        SW_bonded[i][j] = 0; 
   }
   
   //Hoshen Kopelman labels
   for (int i=0; i<ARMS*NX*NY*NZ; i++){
      cpu_SW_L[i] = i;
      cpu_SW_N[i] = -i;
 /*     cpu_SW_L2[i] = i;
      cpu_SW_L3[i] = -1;
      visited[i] = 0;*/
   }

}

void cpu_SW_neighboring_spins(int i){

   int cell = i % (NX*NY*NZ);
   int arm = i / (NX*NY*NZ);
   
   int x = cell % NX;
   int y = (cell/NX) % NY;
   int z = cell/(NX*NY);
   
   int neigh_cell;
   
   switch (arm){
   
      case 0:
        neigh_cell = (x+1)%NX + y*NX + z*NX*NY;
        neighbor_spin[0] = 1*NX*NY*NZ+neigh_cell;
        break;
        
      case 1:
        neigh_cell = (x-1+NX)%NX + y*NX + z*NX*NY;
        neighbor_spin[0] = 0*NX*NY*NZ+neigh_cell;
        break;
        
      case 2:
        neigh_cell = x + ((y+1)%NY)*NX + z*NX*NY;
        neighbor_spin[0] = 3*NX*NY*NZ+neigh_cell;
        break;
        
      case 3:
        neigh_cell = x + ((y-1+NY)%NY)*NX + z*NX*NY;
        neighbor_spin[0] = 2*NX*NY*NZ+neigh_cell;
        break;
        
      case 4:
        neigh_cell = x + y*NX + ((z+1)%NZ)*NX*NZ;
        neighbor_spin[0] = 5*NX*NY*NZ+neigh_cell;
        break;
        
      case 5:
        neigh_cell = x + y*NX + ((z-1+NZ)%NZ)*NX*NY;
        neighbor_spin[0] = 4*NX*NY*NZ+neigh_cell;
        break; 
        
      default:
        break;  
   }
   
   int idx = 1;
   for ( int j=0; j<ARMS; j++ ){
   
       if ( j == arm ) continue;
       
       neighbor_spin[idx] = j*NX*NY*NZ+cell;
  
       idx ++;
   }
   
}


void cpu_SW_update(int * HK_label){

   for ( int i=0; i<ARMS*NX*NY*NZ ; i++)
      SW_new_spin[i] = -1;

   int cluster_label;
   for ( int i=0; i<ARMS*NX*NY*NZ ; i++){
   
      cluster_label = i;
      
      
      do{   //N is not enough to label the clusters. This loop looks for N[i] = -i
         cluster_label = - cpu_SW_N[cpu_SW_L[cluster_label]];
      }while( cluster_label !=  - cpu_SW_N[cpu_SW_L[cluster_label]] );  
      
      cpu_SW_N[cpu_SW_L[i]] = - cluster_label;  // update N to its correct value. 
                                       // Shortcut for future appearences of the same label
      
      if ( - cpu_SW_N[cluster_label] != cpu_SW_L[cluster_label] ){
   //      fprintf(stderr,"Time %d. ERROR Seq. Swendsen Wang: An error in Hoshen Kopelman algorithm occured\n",time);
         fprintf(stderr,"L[%d] = %d , L[%d] = %d\n",i,cpu_SW_L[i],cluster_label,cpu_SW_L[cluster_label]);
         fprintf(stderr,"N[%d] = %d , N[%d] = %d\n",cpu_SW_L[i],cpu_SW_N[cpu_SW_L[i]],cluster_label,cpu_SW_N[cluster_label]);
         exit(1);
      }
      
      if ( SW_new_spin[cluster_label] == -1 )
         SW_new_spin[cluster_label] = (int) (q*mersenne());
         
      if ( SW_new_spin[cluster_label] < 0 || SW_new_spin[cluster_label] > 5 ){
         fprintf(stderr,"Error: SW new spin mal calculado\n");
         exit(1);
      }
      
      s[i] = (byte) ( (s[i] + SW_new_spin[cluster_label]) % q );
      
      if ( s[i] < 0 || s[i] > 5 ){
         fprintf(stderr,"Error: SW new spin mal calculado\n");
         exit(1);
      }
      
      HK_label[i] = cluster_label;
      
   }

}

////
//// CLUSTER ANALYSIS SUBROUTINES
////

void print_cluster_statistics (int time, FILE * logsize, int * label){

   for ( int i=0; i<ARMS*NX*NY*NZ ; i++)  //loop over labels
      cluster_size_[i] = 0;
   
   for ( int i=0; i<ARMS*NX*NY*NZ ; i++){ //loop over sigma_ij variables
   
      cluster_size_[label[i]] ++;
      
   }
   
   int largest_size=0, second_size=0,num_clusters=0,num_percolating_clusters = 0;
   float average=0, average2=0;
   
   for ( int i=0; i<ARMS*NX*NY*NZ ; i++){ //loop over labels
   
       if ( cluster_size_[i] == 0 ) continue;
       
       num_clusters ++;
       average += cluster_size_[i];
       average2 += cluster_size_[i]*cluster_size_[i];
       
       if ( cluster_size_[i] > largest_size ){
          second_size = largest_size;
          largest_size = cluster_size_[i];
       } else if ( cluster_size_[i] > second_size ){
          second_size = cluster_size_[i];
       }
             
       // minimum size of a percolating cluster is 2L (chain from 0 to L, including up and dn spins)
       if ( cluster_size_[i] >= 2*NX || cluster_size_[i] >= 2*NY || cluster_size_[i] >= 2*NZ )
           num_percolating_clusters += is_percolating_cluster(label,i);  // returns 1 if percolates, 0 otherwise
       
   }
      
   average /= (float) num_clusters;
   average2 /= (float) num_clusters;
   
   average2 -= average*average;
    
   fprintf(logsize,"%d %d %d %d %f %f %d\n",time,largest_size,second_size,num_percolating_clusters,
                                         average,average2,num_clusters);  
   

}

int is_percolating_cluster(int * label, int id){
 
   for ( int i=0; i<NX; i++ )
       Xcoord[i] = 0;
   for ( int i=0; i<NY; i++ )
       Ycoord[i] = 0;
   for ( int i=0; i<NZ; i++ )
       Zcoord[i] = 0;
     
   int counter=0;
   for ( int i=0; i<ARMS*NX*NY*NZ; i++ ){  // loop over sigma_ij variables

         if ( label[i] == id ){
               
            int cell = i%(NX*NY*NZ);   
            int x = cell % NX;
            int y = (cell/NX) % NY;
            int z = cell/(NX*NY);
            
            Xcoord[x] = 1;
            Ycoord[y] = 1;
            Zcoord[z] = 1;
            
            counter ++;
         
         }  
         
         if ( counter == cluster_size_[id] ) //the whole cluster has been ckecked
            break;
            
   }
   
   int flag_x=1,flag_y=1,flag_z=1;
   
   for ( int i=0; i<NX; i++ )
      if  ( Xcoord[i] == 0 ){
         flag_x = 0;
         break;
      }
      
   for ( int i=0; i<NY; i++ )
      if  ( Ycoord[i] == 0 ){
         flag_y = 0;
         break;
      }
      
   for ( int i=0; i<NZ; i++ )
      if  ( Zcoord[i] == 0 ){
         flag_z = 0;
         break;
      }
   
   if ( flag_x == 1 || flag_y == 1 || flag_z==1 )
      return 1;
   else
      return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////////

////
//// UTILITES
////

void compute_indices (int * indicesA,int * indicesB)
{
  int i,x,y,z,iA=0,iB=0;
  
  for (i=0; i < NX*NY*NZ; ++i) {
    
    x = i % NX;
    y = (i / NX) % NY;
    z = i / (NX*NY);

    if (z % 2) {
      if (y % 2) { 
        if (x % 2)
          indicesA[iA++] = i;
        else
          indicesB[iB++] = i;
      }
      else {
        if (x % 2)
          indicesB[iB++] = i;
        else
          indicesA[iA++] = i;
      }
    }
    else {
      if (y % 2) { 
        if (x % 2)
          indicesB[iB++] = i;
        else
          indicesA[iA++] = i;
      }
      else {
        if (x % 2)
          indicesA[iA++] = i;
        else
          indicesB[iB++] = i;
      }
    }
  }
}

////
////  GPU RANDOM NUMBER GENERATOR
////

__global__ void gpu_RNG_setup ( 
	curandState * state, 
	uint * seed)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    
    //curand_init( (seed << 20) + id, 0, 0, &state[id]);
    curand_init( seed[id], 0, 0, &state[id]);
}

__global__ void gpu_RNG_generate ( 
	curandState* globalState, 
	float * Rnd,
        int n_rand) 
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;
    
    curandState localState = globalState[ind];
    
    while(ind < (NX*NY*NZ) *n_rand) {
	    
	   // __syncthreads();
	    
	    Rnd[ind] = curand_uniform( &localState );   //returns a random number within the range (0, 1]
	    
	    ind += blockDim.x*gridDim.x;
    }
    
    ind = blockIdx.x * blockDim.x + threadIdx.x;
    
    globalState[ind] = localState; 
}


inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
	if ( cudaSuccess != err )
	{
		fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
		file, line, cudaGetErrorString( err ) );
		exit( -1 );
	}
#endif
	
	return;
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
	cudaError err = cudaGetLastError();
	if ( cudaSuccess != err )
	{
		fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
		file, line, cudaGetErrorString( err ) );
		exit( -1 );
	}
	 
	// More careful checking. However, this will affect performance.
	// Comment away if needed.
//	err = cudaDeviceSynchronize();
//	if( cudaSuccess != err )
//	{
//		fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
//		file, line, cudaGetErrorString( err ) );
//		exit( -1 );
//	}
#endif
	 
	return;
}




//////////////////////// FAST CUDA RANDOM GENERATOR ////////////////////////////



// S1, S2, S3, and M are all constants, z is the inner state  
__device__
static uint TausStep(uint &z, int S1, int S2, int S3, uint M) {  
  uint b=(((z << S1) ^ z) >> S2); 
  return z = (((z & M) << S3) ^ b);  
} 

__device__
// A and C are constants 
static uint LCGStep(uint &z, uint A, uint C) {  
  return z=(A*z+C);  
} 


__device__
static float HybridTaus(uint& z1, uint& z2, uint& z3, uint& z4) {  
  // Combined period is lcm(p1,p2,p3,p4)~ 2^121
  float randval;

  //return 2.3283064365387e-10f*LCGStep(z4, 1664525, 1013904223UL);

  //do { 
   randval = 2.3283064365387e-10f * (          // Periods  
    TausStep(z1, 13, 19, 12, 4294967294UL) ^  // p1=2^31-1  
    TausStep(z2, 2, 25, 4, 4294967288UL) ^    // p2=2^30-1  
    TausStep(z3, 3, 11, 17, 4294967280UL) ^   // p3=2^28-1  
    LCGStep(z4, 1664525, 1013904223UL)        // p4=2^32  
   );
  //} while (!(randval > 0.0f && randval < 1.0f));
  return randval;
}  

__global__
static void _rand(float* vec, uint* z1, uint* z2, uint* z3, uint* z4) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  vec[i] = HybridTaus(z1[i],z2[i],z3[i],z4[i]);
}


//////////////////// xorshift1024*Phi Random Generator   /////////////

void initRand(uint64_t  seed)
    {
        xor_s[0] = seed;
        xor_s[1] = seed + seed;
        
        for ( int i = 0 ; i < 16 ; ++i )
        {
            xor_s[ i ] = nextRand();
        }
        xor_p = nextRand() % 15;
    }

uint64_t nextRand(void) {
	const uint64_t s0 = xor_s[xor_p];
	uint64_t s1 = xor_s[xor_p = (xor_p + 1) & 15];
	s1 ^= s1 << 31; // a
	xor_s[xor_p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
	return xor_s[xor_p] * 0x9e3779b97f4a7c13;
}

double uniformDoubleRand(void){
     return (double) nextRand() / (double) UINT64_MAX;
}


////
////  CORRELATOR SUBROUTINES
////

void initialize(int N, struct correlator mol[]) {
	
	int i,j,k,length=numcorrelators*p;
	
	for (k = 0; k < N; k ++) {
		for (j=0;j<numcorrelators;++j) {
			for (i=0;i<p;++i) {
				mol[k].shift[j][i] = -2E10;
				mol[k].correlation[j][i] = 0;
				mol[k].ncorrelation[j][i] = 0;
			}
			mol[k].accumulator[j] = 0.0;
			mol[k].naccumulator[j] = 0;
			mol[k].insertindex[j] = 0;
		}
	
		for (i=0;i<length;++i) {
			mol[k].t[i] = 0;
			mol[k].f[i] = 0;
		}
		
		mol[k].npcorr = 0;
		mol[k].kmax=0;
		mol[k].accval=0;
	}
}

void add(struct correlator mol[], int i, double w, int k) {
	int j, dmin = p/m_;
	
	/// If we exceed the correlator side, the value is discarded
	if (k == numcorrelators)
		return;
	if (k > mol[i].kmax) 
		mol[i].kmax = k;

	/// Insert new value in shift array
	mol[i].shift[k][mol[i].insertindex[k]] = w;

	/// Add to average value
	if (k==0)
		mol[i].accval += w;

	/// Add to accumulator and, if needed, add to next correlator
	mol[i].accumulator[k] += w;
	++mol[i].naccumulator[k];
	if (mol[i].naccumulator[k]==m_) {
		add(mol,i,mol[i].accumulator[k]/m_, k+1);
		mol[i].accumulator[k]=0;
		mol[i].naccumulator[k]=0;
	}

	/// Calculate correlation function
	int ind1 = mol[i].insertindex[k];
	
	if (k==0) { /// First correlator is different
		int ind2 = ind1;
		
		for (j=0;j<p;++j) {
			if (mol[i].shift[k][ind2] > -1e10) {
				mol[i].correlation[k][j] += mol[i].shift[k][ind1] * mol[i].shift[k][ind2];
				++mol[i].ncorrelation[k][j];
			}
			--ind2;
			if (ind2<0)
				ind2+=p;
		}
	}
	else {
		int ind2=ind1-dmin;
		for (j=dmin;j<p;++j) {
			if (ind2<0) 
				ind2+=p;
			if (mol[i].shift[k][ind2] > -1e10) {
				mol[i].correlation[k][j] += mol[i].shift[k][ind1] * mol[i].shift[k][ind2];
				++mol[i].ncorrelation[k][j];
			}
			--ind2;
		}
	}

	++ mol[i].insertindex[k];
	if (mol[i].insertindex[k]==p) 
		mol[i].insertindex[k]=0;
}

void evaluate(struct correlator mol[], int N, int norm) {
	int i,j,k,im,dmin=p/m_;

	double aux;
	
	for (j = 0; j < N; j ++) {
		
		aux = 0;
		im = 0;
		
		if (norm)
			aux = (mol[j].accval/mol[j].ncorrelation[0][0])*(mol[j].accval/mol[j].ncorrelation[0][0]);
	
		// First correlator
		for (i=0;i<p;++i) {
			if (mol[j].ncorrelation[0][i] > 0) {
				mol[j].t[im] = i;
				mol[j].f[im] = mol[j].correlation[0][i]/mol[j].ncorrelation[0][i] - aux;
				++im;
			}
		}
	
		// Subsequent correlators
		for (k=1;k<mol[j].kmax;++k) {
			for (i=dmin;i<p;++i) {
				if (mol[j].ncorrelation[k][i]>0) {
					mol[j].t[im] = i * pow((double)m_, k);
					mol[j].f[im] = mol[j].correlation[k][i] / mol[j].ncorrelation[k][i] - aux;
					++im;
				}
			}
		}
	
		mol[j].npcorr = im;
	
	}
}


void printCorr(struct correlator corr[], char filename[], int N) {
	
	FILE * fp;
//	char filename[100];
	int i,j;
	
//	sprintf(filename,"autocorr/autocorrelation-N%d-P%1.3f-%1.3f.dat",N,P,T);
	fp = fopen(filename,"w");
	
	double meanCorr[corr[0].npcorr];
	
	for (i = 0; i < corr[0].npcorr; i ++) {
		meanCorr[i] = 0;
		for(j = 0; j < N; j ++)
			meanCorr[i] += corr[j].f[i];
		meanCorr[i] /= (double) N;
		fprintf(fp,"%.0f\t%f\n",corr[0].t[i],meanCorr[i]/meanCorr[0]);
	}
	
	fclose(fp);
}


////
////  CHESSBOARD GPU ALGORITHM 
////


/****   EDGES AND VERTEX DEFINITION IN A CUBE

VERTEX DEFINITION

	Bottom layer (Z=0)       Y direction
                                  ^
    v3 [----------------] v8      [
       [                ]         [
       [                ]         [
       [                ]         [
       [                ]         [
       [                ]         [
    v1 [----------------] v2      [---------------> X direction

                                 v1 is placed at the origin of coordintates (X=0,Y=0,Z=0)
	Top layer (Z=1)

     v6	[----------------] v5
        [                ]
        [                ]
        [                ]
        [                ]
        [                ]
     v4 [----------------] v7

EDGES DEFINITION
	0	v1 -> v2
	1	v1 -> v3
	2	v1 -> v4
	3	v2 -> v7
	4	v2 -> v8
	5	v3 -> v6
	6	v3 -> v8
	7	v4 -> v6
	8	v4 -> v7
	9	v6 -> v5
	10	v7 -> v5
	11	v8 -> v5



****/


__global__ void gpu_set_cubes (int * devChessboard_edge, int * devChessboard_vertex)
{
    uint tid = blockIdx.x*blockDim.x + threadIdx.x;
	
    uint x = tid % NX;
    uint y = (tid/NX) % NY;
    uint z = tid/(NX*NY);


    if ( z%2 == 0 ){          // only survive threads with x, y, z all even or all odd
        if ( (x%2 != 0) || (y%2 != 0) ){
           return;
        }
    } else if ( z%2 == 1 ){
        if ( (x%2 != 1) || (y%2 != 1) ){
           return;
        }
    }

    int cube = (x/2) + (NX/2)*(y/2) + ((NX*NY)/4)*z;
    int NCUBES = (NX*NY*NZ/4);

//    if ( cube > NX*NY*NZ/4 )
//       return;

    for (int i=0; i < EDGES; i++)
             devChessboard_edge[cube+i*NCUBES] = 0;
    
    ////// set Vertex
     int v1, v2, v3, v4, v5, v6, v7, v8;

     v1 = x + y*NX + z*NX*NY;				// v1 = tid
     v2 = (x+1)%NX + y*NX + z*NX*NY;      		// v2 = xup ( v1 )
     v3 = x + ((y+1)%NY)*NX + z*NX*NY;    		// v3 = yup ( v1 )
     v4 = x + y*NX + ((z+1)%NZ)*NX*NY;    		// v4 = zup ( v1 )
     v6 = x + ((y+1)%NY)*NX + ((z+1)%NZ)*NX*NY;    	// v6 = zup ( v3 ) = zup ( yup ( v1 ) )
     v7 = (x+1)%NX + y*NX + ((z+1)%NZ)*NX*NY;		// v7 = zup ( v2 ) = zup ( xup ( v1 ) )
     v8 = (x+1)%NX + ((y+1)%NY)*NX + z*NX*NY;		// v8 = yup ( v2 ) = yup ( xup ( v1 ) )
     v5 = (x+1)%NX + ((y+1)%NY)*NX + ((z+1)%NZ)*NX*NY;  // v5 = xup ( v6 ) = xup ( zup ( yup ( v1 ) ) )
 
     devChessboard_vertex[cube+0*NCUBES+0*NCUBES*EDGES] = v1;
     devChessboard_vertex[cube+0*NCUBES+1*NCUBES*EDGES] = v2;

     devChessboard_vertex[cube+1*NCUBES+0*NCUBES*EDGES] = v1;
     devChessboard_vertex[cube+1*NCUBES+1*NCUBES*EDGES] = v3;

     devChessboard_vertex[cube+2*NCUBES+0*NCUBES*EDGES] = v1;
     devChessboard_vertex[cube+2*NCUBES+1*NCUBES*EDGES] = v4;

     devChessboard_vertex[cube+3*NCUBES+0*NCUBES*EDGES] = v2;
     devChessboard_vertex[cube+3*NCUBES+1*NCUBES*EDGES] = v7;

     devChessboard_vertex[cube+4*NCUBES+0*NCUBES*EDGES] = v2;
     devChessboard_vertex[cube+4*NCUBES+1*NCUBES*EDGES] = v8;

     devChessboard_vertex[cube+5*NCUBES+0*NCUBES*EDGES] = v3;
     devChessboard_vertex[cube+5*NCUBES+1*NCUBES*EDGES] = v6;

     devChessboard_vertex[cube+6*NCUBES+0*NCUBES*EDGES] = v3;
     devChessboard_vertex[cube+6*NCUBES+1*NCUBES*EDGES] = v8;

     devChessboard_vertex[cube+7*NCUBES+0*NCUBES*EDGES] = v4;
     devChessboard_vertex[cube+7*NCUBES+1*NCUBES*EDGES] = v6;

     devChessboard_vertex[cube+8*NCUBES+0*NCUBES*EDGES] = v4;
     devChessboard_vertex[cube+8*NCUBES+1*NCUBES*EDGES] = v7;

     devChessboard_vertex[cube+9*NCUBES+0*NCUBES*EDGES] = v6;
     devChessboard_vertex[cube+9*NCUBES+1*NCUBES*EDGES] = v5;

     devChessboard_vertex[cube+10*NCUBES+0*NCUBES*EDGES] = v7;
     devChessboard_vertex[cube+10*NCUBES+1*NCUBES*EDGES] = v5;

     devChessboard_vertex[cube+11*NCUBES+0*NCUBES*EDGES] = v8;
     devChessboard_vertex[cube+11*NCUBES+1*NCUBES*EDGES] = v5;
    //////

}

__global__ void gpu_Chessboard_start_cubes( int N, float * dev_Rnd_cubeFlip, int * devSolutions, byte * dev_active_bond,
  int * devChessboard_edge, int * devChessboard_vertex, int * devChessboard_state )   
// output ret in typewriter fashion
{

    uint tid = blockIdx.x*blockDim.x + threadIdx.x;
	
    uint x = tid % NX;
    uint y = (tid/NX) % NY;
    uint z = tid/(NX*NY);
 

    if ( z%2 == 0 ){          // only survive threads with x, y, z all even or all odd
        if ( (x%2 != 0) || (y%2 != 0) ){
           return;
        }
    }else if ( z%2 == 1 ){
        if ( (x%2 != 1) || (y%2 != 1) ){
           return;
        }
    }

   int cube = (x/2) + (NX/2)*(y/2) + ((NX*NY)/4)*z;
   int NCUBES = (NX*NY*NZ/4);

   int r = (int) (dev_Rnd_cubeFlip[cube]*NSOLUTIONS);
   r %= NSOLUTIONS;

   devChessboard_state[cube] = r;

    for(int i=0; i<EDGES; i++)
        devChessboard_edge[cube+i*NCUBES] = devSolutions[r*EDGES+i];
    
// Translate from chessboard to water bonding vector

    for (int j=0; j<EDGES; j ++){

        int v0 = devChessboard_vertex[cube+j*NCUBES+0*NCUBES*EDGES];
        int v1 = devChessboard_vertex[cube+j*NCUBES+1*NCUBES*EDGES];

        if ( j == 0 || j == 6 || j == 8 || j == 9 ){   // xup edges
              dev_active_bond[0*N+v0] = devChessboard_edge[cube+j*NCUBES];
              dev_active_bond[1*N+v1] = devChessboard_edge[cube+j*NCUBES];
        } else if ( j == 1 || j == 4 || j == 7 || j == 10 ){    // yup edges
              dev_active_bond[2*N+v0] = devChessboard_edge[cube+j*NCUBES];
              dev_active_bond[3*N+v1] = devChessboard_edge[cube+j*NCUBES];
        } else if ( j == 2 || j == 3 || j == 5 || j == 11 ){   // zup edges
              dev_active_bond[4*N+v0] = devChessboard_edge[cube+j*NCUBES];
              dev_active_bond[5*N+v1] = devChessboard_edge[cube+j*NCUBES];
        } 

    }

}

__global__ void gpu_Chessboard_set_state( int N, byte * dev_active_bond,
                int * devChessboard_vertex, int * devChessboard_state )   
// output ret in typewriter fashion
// sets cube state according to active bonds
{

    uint tid = blockIdx.x*blockDim.x + threadIdx.x;
	
    uint x = tid % NX;
    uint y = (tid/NX) % NY;
    uint z = tid/(NX*NY);
 

    if ( z%2 == 0 ){          // only survive threads with x, y, z all even or all odd
        if ( (x%2 != 0) || (y%2 != 0) ){
           return;
        }
    }else if ( z%2 == 1 ){
        if ( (x%2 != 1) || (y%2 != 1) ){
           return;
        }
    }

   int cube = (x/2) + (NX/2)*(y/2) + ((NX*NY)/4)*z;
   int NCUBES = (NX*NY*NZ/4);

   //cube vertex (cell id)
   int v1 = devChessboard_vertex[cube+0*NCUBES+0*NCUBES*EDGES];
   int v2 = devChessboard_vertex[cube+3*NCUBES+0*NCUBES*EDGES];
   int v3 = devChessboard_vertex[cube+5*NCUBES+0*NCUBES*EDGES];
   int v4 = devChessboard_vertex[cube+7*NCUBES+0*NCUBES*EDGES];
   int v5 = devChessboard_vertex[cube+7*NCUBES+1*NCUBES*EDGES];
   int v6 = devChessboard_vertex[cube+5*NCUBES+1*NCUBES*EDGES];
   int v7 = devChessboard_vertex[cube+3*NCUBES+1*NCUBES*EDGES];
   int v8 = devChessboard_vertex[cube+4*NCUBES+1*NCUBES*EDGES];

   if ( dev_active_bond[0*N+v1] == 0 ){

      if ( dev_active_bond[0*N+v6] == 0 ) devChessboard_state[cube] = 2;
      else if ( dev_active_bond[4*N+v3] == 0 ) devChessboard_state[cube] = 3;
      else devChessboard_state[cube] = 8;

   } else if ( dev_active_bond[2*N+v1] == 0 ){

      if ( dev_active_bond[4*N+v2] == 0 ) devChessboard_state[cube] = 4;
      else if ( dev_active_bond[2*N+v7] == 0 ) devChessboard_state[cube] = 5;
      else devChessboard_state[cube] = 7;
 
  } else {

      if ( dev_active_bond[0*N+v3] == 0 ) devChessboard_state[cube] = 0;
      else if ( dev_active_bond[4*N+v8] == 0 ) devChessboard_state[cube] = 1;
      else devChessboard_state[cube] = 6;

  }


}


__global__ void gpu_Chessboard_set_active_bonds(
   byte * dev_s, 
   int N, 
   double Jeff, 
   double beta, 
   float * dev_Rnd_cubeFlip, 
   float * dev_Rnd_chess,
   int offset_chess, 
   int * devSolutions,
   byte * dev_active_bond,
   int * devChessboard_edge,
   int * devChessboard_vertex,
   int * devChessboard_state )
{

    uint tid = blockIdx.x*blockDim.x + threadIdx.x;
	
    uint x = tid % NX;
    uint y = (tid/NX) % NY;
    uint z = tid/(NX*NY);

    if ( z%2 == 0 ){          // only survive threads with x, y, z all even or all odd
        if ( (x%2 != 0) || (y%2 != 0) ){
           return;
        }
    }else if ( z%2 == 1 ){
        if ( (x%2 != 1) || (y%2 != 1) ){
           return;
        }
    }

    int cube = (x/2) + (NX/2)*(y/2) + ((NX*NY)/4)*z;
    int NCUBES = NX*NY*NZ/4;

    int old_state = devChessboard_state[cube];

//    int new_state =  ( old_state + 1 + (int) ((NSOLUTIONS-1)*dev_Rnd_cubeFlip[cube+offset_chess]) ) % NSOLUTIONS; 
    int new_state = (int) (dev_Rnd_cubeFlip[cube+offset_chess]*NSOLUTIONS);
    new_state %= NSOLUTIONS;

    double delta_H=0;
			// Chessboard Metropolis
    for (int j=0; j<EDGES; j++){

        int node0 = devChessboard_vertex[cube+j*NCUBES+0*EDGES*NCUBES];
        int node1 = devChessboard_vertex[cube+j*NCUBES+1*EDGES*NCUBES];

        int arm0, arm1;
        if ( j == 0 || j == 6 || j == 8 || j == 9 ){   // xup edges
              arm0 = 0;
              arm1 = 1;
        } else if ( j == 1 || j == 4 || j == 7 || j == 10 ){    // yup edges
              arm0 = 2;
              arm1 = 3;
        } else if ( j == 2 || j == 3 || j == 5 || j == 11 ){   // zup edges
	      arm0 = 4;
              arm1 = 5;
        } 

        delta_H += ( devSolutions[new_state*EDGES+j] - devSolutions[old_state*EDGES+j] ) * delta(dev_s[arm0*N+node0],dev_s[arm1*N+node1]);
    }
    
    delta_H *= Jeff;

    if ( dev_Rnd_chess[cube+offset_chess] > expf(beta*delta_H) ){
              new_state = old_state;
    }
    
    devChessboard_state[cube] = new_state;
    for(int j=0; j<EDGES; j++)
        devChessboard_edge[cube+j*NCUBES] = devSolutions[new_state*EDGES+j];
    
// Translate from chessboard to water bonding vector

    for (int j=0; j<EDGES; j ++){

        int v0 = devChessboard_vertex[cube+j*NCUBES+0*EDGES*NCUBES];
        int v1 = devChessboard_vertex[cube+j*NCUBES+1*EDGES*NCUBES];

        if ( j == 0 || j == 6 || j == 8 || j == 9 ){   // xup edges
              dev_active_bond[0*N+v0] = devChessboard_edge[cube+j*NCUBES];
              dev_active_bond[1*N+v1] = devChessboard_edge[cube+j*NCUBES];
        } else if ( j == 1 || j == 4 || j == 7 || j == 10 ){    // yup edges
              dev_active_bond[2*N+v0] = devChessboard_edge[cube+j*NCUBES];
              dev_active_bond[3*N+v1] = devChessboard_edge[cube+j*NCUBES];
        } else if ( j == 2 || j == 3 || j == 5 || j == 11 ){   // zup edges
              dev_active_bond[4*N+v0] = devChessboard_edge[cube+j*NCUBES];
              dev_active_bond[5*N+v1] = devChessboard_edge[cube+j*NCUBES];
        } 

    }

}
