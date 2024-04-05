#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#define EDGES 12
#define NSOLUTIONS 9

typedef uint8_t byte;

//#define CHECK

/****   EDGES AND VERTEX DEFINITION IN A CUBE

VERTEX DEFINITION

	Bottom layer           Y direction
                                  ^
    v3 [----------------] v8      [
       [                ]         [
       [                ]         [
       [                ]         [
       [                ]         [
       [                ]         [
    v1 [----------------] v2      [---------------> X direction

                                 v1 is placed at the origin of coordintates (X=0,Y=0,Z=0)
	Top layer 

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


struct cube {

  int * edge;
  int ** vertex;

};

struct gpu_cube {

  int * edge;
  int * vertex; // GPU does not allow 2D arrays
  int state;    // state identifies which of the nine possible solutions holds for the cube
		// it is introduced to help in the GPU implementation of the algorithm

};

struct cube *Chessboard;
int ** Solutions;
int * next_cube;
int NUM_CUBES;

int function_xup ( int i, int NX, int NY ){ 
   
   int x = i % NX;
   int y = (i / NX)%NY;
   int z = i / (NX*NY);

   return (x+1)%NX + y*NX + z*NX*NY; 
}

int function_yup ( int i, int NX, int NY ){ 
   
   int x = i % NX;
   int y = (i / NX)%NY;
   int z = i / (NX*NY);

   return x + ((y+1)%NY)*NX + z*NX*NY;
}

int function_zup ( int i, int NX, int NY , int NZ ){ 
   
   int x = i % NX;
   int y = (i / NX)%NY;
   int z = i / (NX*NY);

   return x + y*NX + ((z+1)%NZ)*NX*NY;
}


void setVertex ( int i, int v1, int NX, int NY, int NZ )
{
     int v2, v3, v4, v5, v6, v7, v8;

     v2 = function_xup( v1 , NX, NY );
     v3 = function_yup( v1 , NX, NY );
     v4 = function_zup( v1 , NX, NY , NZ );
     v6 = function_zup( v3 , NX, NY , NZ );
     v7 = function_zup( v2 , NX, NY , NZ );
     v8 = function_yup( v2 , NX, NY );
     v5 = function_xup( v6 , NX, NY );

     Chessboard[i].vertex[0][0] = v1;
     Chessboard[i].vertex[0][1] = v2;

     Chessboard[i].vertex[1][0] = v1;
     Chessboard[i].vertex[1][1] = v3;

     Chessboard[i].vertex[2][0] = v1;
     Chessboard[i].vertex[2][1] = v4;

     Chessboard[i].vertex[3][0] = v2;
     Chessboard[i].vertex[3][1] = v7;

     Chessboard[i].vertex[4][0] = v2;
     Chessboard[i].vertex[4][1] = v8;

     Chessboard[i].vertex[5][0] = v3;
     Chessboard[i].vertex[5][1] = v6;

     Chessboard[i].vertex[6][0] = v3;
     Chessboard[i].vertex[6][1] = v8;

     Chessboard[i].vertex[7][0] = v4;
     Chessboard[i].vertex[7][1] = v6;

     Chessboard[i].vertex[8][0] = v4;
     Chessboard[i].vertex[8][1] = v7;

     Chessboard[i].vertex[9][0] = v6;
     Chessboard[i].vertex[9][1] = v5;

     Chessboard[i].vertex[10][0] = v7;
     Chessboard[i].vertex[10][1] = v5;

     Chessboard[i].vertex[11][0] = v8;
     Chessboard[i].vertex[11][1] = v5;

}

void setSolutions()
{

   for (int i=0; i<NSOLUTIONS; i++)
      for (int j=0; j<EDGES; j++)
         Solutions[i][j] = 1;			//default value: all edges are active

                                    // Now, I write the nine solutions (deactivate 4 edges in each of them)
  Solutions[0][2] = 0;
  Solutions[0][3] = 0;
  Solutions[0][6] = 0;
  Solutions[0][9] = 0;

  Solutions[1][2] = 0;
  Solutions[1][3] = 0;
  Solutions[1][5] = 0;
  Solutions[1][11] = 0;

  Solutions[2][0] = 0;
  Solutions[2][6] = 0;
  Solutions[2][8] = 0;
  Solutions[2][9] = 0;

  Solutions[3][0] = 0;
  Solutions[3][5] = 0;
  Solutions[3][8] = 0;
  Solutions[3][11] = 0;

  Solutions[4][1] = 0;
  Solutions[4][3] = 0;
  Solutions[4][7] = 0;
  Solutions[4][11] = 0;

  Solutions[5][1] = 0;
  Solutions[5][4] = 0;
  Solutions[5][7] = 0;
  Solutions[5][10] = 0;

  Solutions[6][2] = 0;
  Solutions[6][4] = 0;
  Solutions[6][5] = 0;
  Solutions[6][10] = 0;

  Solutions[7][1] = 0;
  Solutions[7][4] = 0;
  Solutions[7][8] = 0;
  Solutions[7][9] = 0;

  Solutions[8][0] = 0;
  Solutions[8][6] = 0;
  Solutions[8][7] = 0;
  Solutions[8][10] = 0;
}


void setSolutions_Typewriter(int * cpuSolutions)
{

   for (int i=0; i<NSOLUTIONS*EDGES; i++)
         cpuSolutions[i] = 1;			//default value: all edges are active

                                    // Now, I write the nine solutions (deactivate 4 edges in each of them)
  cpuSolutions[0*EDGES+2] = 0;
  cpuSolutions[0*EDGES+3] = 0;
  cpuSolutions[0*EDGES+6] = 0;
  cpuSolutions[0*EDGES+9] = 0;

  cpuSolutions[1*EDGES+2] = 0;
  cpuSolutions[1*EDGES+3] = 0;
  cpuSolutions[1*EDGES+5] = 0;
  cpuSolutions[1*EDGES+11] = 0;

  cpuSolutions[2*EDGES+0] = 0;
  cpuSolutions[2*EDGES+6] = 0;
  cpuSolutions[2*EDGES+8] = 0;
  cpuSolutions[2*EDGES+9] = 0;

  cpuSolutions[3*EDGES+0] = 0;
  cpuSolutions[3*EDGES+5] = 0;
  cpuSolutions[3*EDGES+8] = 0;
  cpuSolutions[3*EDGES+11] = 0;

  cpuSolutions[4*EDGES+1] = 0;
  cpuSolutions[4*EDGES+3] = 0;
  cpuSolutions[4*EDGES+7] = 0;
  cpuSolutions[4*EDGES+11] = 0;

  cpuSolutions[5*EDGES+1] = 0;
  cpuSolutions[5*EDGES+4] = 0;
  cpuSolutions[5*EDGES+7] = 0;
  cpuSolutions[5*EDGES+10] = 0;

  cpuSolutions[6*EDGES+2] = 0;
  cpuSolutions[6*EDGES+4] = 0;
  cpuSolutions[6*EDGES+5] = 0;
  cpuSolutions[6*EDGES+10] = 0;

  cpuSolutions[7*EDGES+1] = 0;
  cpuSolutions[7*EDGES+4] = 0;
  cpuSolutions[7*EDGES+8] = 0;
  cpuSolutions[7*EDGES+9] = 0;

  cpuSolutions[8*EDGES+0] = 0;
  cpuSolutions[8*EDGES+6] = 0;
  cpuSolutions[8*EDGES+7] = 0;
  cpuSolutions[8*EDGES+10] = 0;
}


void initialize_Chessboard(int NX, int NY, int NZ){

    if (NX != NY || NX != NZ || NY != NZ ) {printf("Error: can't build chessboard.\n\tThe system is not a cube.\n"); exit(1);} 
	//Maybe it could be solved form rectangles

    int N = NX*NY*NZ;
    int L = NX;

    if ( L % 4 != 0 ) {printf("Error: can't build chessboard.\n\tUse L muliply of 4.\n"); exit(1);} 

    NUM_CUBES = N/4;

    Chessboard = (struct cube*) malloc(sizeof(struct cube)*NUM_CUBES);

    for ( int i=0; i < NUM_CUBES; i++ ){
          Chessboard[i].edge = (int *) calloc(EDGES,sizeof(int));
          Chessboard[i].vertex = (int **) malloc(sizeof(int*)*EDGES);

          for ( int j=0; j < EDGES; j++){
             Chessboard[i].edge[j] = 0;
             Chessboard[i].vertex[j] = (int *) calloc(2,sizeof(int));   //The edge j joins vertices vertex[j][0] and vertex[j][1]
          }
    } 

    int cube_index = 0;
    for ( int i=0; i<N ; i++ ){

      int x = i % NX;
      int y = (i / NX)%NY;
      int z = i / (NX*NY);

      if ( z % 2 == 0 )
        if ( y % 2 == 0 )
           if ( x % 2 == 0 ){
               if (cube_index >= N/4) {printf("Fatal error. Attempted to bild extra cube in chessboard\n"); exit(1);}
               setVertex(cube_index,i,NX,NY,NZ);
               cube_index ++;
           }
      if ( z % 2 != 0 )
         if ( y % 2 != 0 )
           if ( x % 2 != 0 ){
               if (cube_index >= N/4) {printf("Fatal error. Attempted to bild extra cube in chessboard\n"); exit(1);}
                setVertex(cube_index,i,NX,NY,NZ);
               cube_index ++;
           } 

    }

    if (cube_index != NUM_CUBES){printf("Fatal error. Chessboard has not been built correctly\n"); exit(1);}

    Solutions = (int**) malloc(sizeof(int*)*NSOLUTIONS);

    for(int i=0; i<NSOLUTIONS; i++)
       Solutions[i] = (int*) calloc(EDGES,sizeof(int));

    setSolutions();

    next_cube = (int*) calloc(EDGES,sizeof(int));
}  

//intput s in typewriter fashion
void Chessboard_Metropolis_Typewriter(int index_cube, int r, double J, double beta, int N, byte * s, int * ret)   // output -> ret
{    // J = J_HB-Press*v_HB 	where J_HB > 0    (Following Valentino's criteria)

     int *previous, *next ;
     previous = (int *) calloc(EDGES,sizeof(int));
     next = (int *) calloc(EDGES,sizeof(int));

     for (int j=0; j<EDGES; j++){
        previous[j] = Chessboard[index_cube].edge[j];
        next[j] = Solutions[r][j];
     }

     double delta_H=0;
     
     for (int j=0; j<EDGES; j++){

        if (previous[j] == next[j]) continue;   //no change form active to not active

        int node0 = Chessboard[index_cube].vertex[j][0];
        int node1 = Chessboard[index_cube].vertex[j][1];

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

        int sp0 = s[arm0*N+node0];
        int sp1 = s[arm1*N+node1];

        if (sp0 != sp1) continue;  // no HB is created nor destroyed

        if ( (previous[j] == 1) ){   // destroy HB
              delta_H -= J;
        } else if ( (previous[j] == 0) ){   // create HB
              delta_H += J;
        } else {printf("ERROR: Chessboard_Metropolis\n"); exit(1);}

     }

     if ( (delta_H >= 0.) || (drand48() < exp(beta*delta_H)) )   // accept move
	for (int j=0; j<EDGES; j++)
           ret[j] = next[j];
     else							// reject move
	for (int j=0; j<EDGES; j++)
           ret[j] = previous[j];

     free(next);
     free(previous);
}


void Chessboard_set_active_bonds_Typewriter(int flag, byte * s, int ARMS, int N, double J, double beta, int * rand_vector, int * ret)   
// output ret in typewriter fashion
{

    for(int i=0; i<NUM_CUBES; i++){

       int r = rand_vector[i];

	#ifdef CHECK
       if ( r < 0 || r >= NSOLUTIONS ){
	 fprintf(stderr,"Set_active_bonds Error: Wrong value at input random vector\n");
	 exit(1);
       }
	#endif

       if ( flag == 1 ){
           Chessboard_Metropolis_Typewriter(i, r, J, beta, N, s, next_cube);  //Accept or reject change
       }
       else{
          for (int j=0; j<EDGES; j++)
             next_cube[j] = Solutions[r][j];	// First MC step, before initializing vector s
       }

       for(int j=0; j<EDGES; j++)
           Chessboard[i].edge[j] = next_cube[j];

    }


// Translate from chessboard to water bonding vector

#ifdef CHECK
    int * nActive;    // used as control
    nActive = (int*) calloc(N,sizeof(int));
    for (int i=0; i<N; i++)
       nActive[i]=0;    
#endif

    for (int i=0; i<NUM_CUBES; i++){
     for (int j=0; j<EDGES; j ++){
        int v0 = Chessboard[i].vertex[j][0];
        int v1 = Chessboard[i].vertex[j][1];

        if ( j == 0 || j == 6 || j == 8 || j == 9 ){   // xup edges
              if ( Chessboard[i].edge[j] == 0 ){
                 ret[0*N+v0] = 0;
                 ret[1*N+v1] = 0;
              } else if ( Chessboard[i].edge[j] == 1 ){
                 ret[0*N+v0] = 1;
                 ret[1*N+v1] = 1;
		#ifdef CHECK
                 nActive[v0] ++;
                 nActive[v1] ++;
		#endif
              }
        } else if ( j == 1 || j == 4 || j == 7 || j == 10 ){    // yup edges
              if ( Chessboard[i].edge[j] == 0 ){
                 ret[2*N+v0] = 0;
                 ret[3*N+v1] = 0;
              } else if ( Chessboard[i].edge[j] == 1 ){
                 ret[2*N+v0] = 1;
                 ret[3*N+v1] = 1;
		#ifdef CHECK
                 nActive[v0] ++;
                 nActive[v1] ++;
		#endif
              }
        } else if ( j == 2 || j == 3 || j == 5 || j == 11 ){   // zup edges
              if ( Chessboard[i].edge[j] == 0 ){
                 ret[4*N+v0] = 0;
                 ret[5*N+v1] = 0;
              } else if ( Chessboard[i].edge[j] == 1 ){
                 ret[4*N+v0] = 1;
                 ret[5*N+v1] = 1;
		#ifdef CHECK
                 nActive[v0] ++;
                 nActive[v1] ++;
		#endif
              }
        } 

      }
    }

#ifdef CHECK
    for(int i=0; i<N; i++)
      if(nActive[i] != 4){
	 fprintf(stderr,"Set_active_bonds Error: Number of active bonds != 4 for a cell\n");
	 exit(1);
      }

    free(nActive);
#endif

}

