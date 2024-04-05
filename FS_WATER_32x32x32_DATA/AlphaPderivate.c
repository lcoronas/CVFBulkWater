#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

int main (){
  int num_temperatures = 100;
  int i=0;
  double *T_, *T, *V_, *V, *V2, *V2_;
  int *n, *n_;

// T vector of sorted temperatures
// T_ vector of temperatures read from file

// V2 is the variance -> V2 = <V^2> - <V>*<V>

  T = (double*) calloc(num_temperatures,sizeof(double));
  T_ = (double*) calloc(num_temperatures,sizeof(double));
  V = (double*) calloc(num_temperatures,sizeof(double));
  V_ = (double*) calloc(num_temperatures,sizeof(double));
  V2 = (double*) calloc(num_temperatures,sizeof(double));
  V2_ = (double*) calloc(num_temperatures,sizeof(double));
  n = (int*) calloc(num_temperatures,sizeof(int));
  n_ = (int*) calloc(num_temperatures,sizeof(int));

  FILE* ifile = fopen("AlphaP_input","r");

  if(ifile == NULL){
       fprintf(stderr,"Can't open input file\n");
       exit(1);
  }

  while(fscanf(ifile,"%lf %lf %lf %d", &T_[i], &V_[i], &V2_[i], &n_[i])==4){
    i++;
    if (i == num_temperatures){
       num_temperatures += 100;
       T_ = (double*)realloc(T_,num_temperatures*sizeof(double));
       T = (double*)realloc(T,num_temperatures*sizeof(double));
       V_ = (double*)realloc(V_,num_temperatures*sizeof(double));
       V = (double*)realloc(V,num_temperatures*sizeof(double));
       V2_ = (double*)realloc(V2_,num_temperatures*sizeof(double));
       V2 = (double*)realloc(V2,num_temperatures*sizeof(double));
       n_ = (int*)realloc(n_,num_temperatures*sizeof(int));
       n = (int*)realloc(n,num_temperatures*sizeof(int));
    }
  }

  // i is the resulting time
  
/// sorting algorithm
{
  int index = 0;
  int check_index;
  double trial;
  bool *check;
  check = (bool*) calloc(i,sizeof(bool));

  for (int j=0;j<i;j++){    //i is the number of data stored in T_, etc.
     T[j]=0;
     V[j]=0;
     V2[j]=0;
     n[j]=0;
     check[j] = false;
  }

  bool next = true;
  while (next){
    for (int j=0;j<i;j++){
      if (check[j] == false){
        trial = T_[j];
        check[j] = true;
        check_index = j;
        for(int k=0;k<i;k++){
          if ( k != j && check[k] == false)
            if (trial > T_[k]) {
              trial = T_[k];
              check[check_index] = false;
              check[k] = true;
              check_index = k;
            }
        }
 
        T[index] = trial;
        V[index] = V_[check_index];
	V2[index] = V2_[check_index];
	n[index] = n_[check_index];
        index ++;
      }
    }
    next = false;
    for (int j=0;j<i;j++){
       if(check[j] == false) next=true;
    }
  }

}
///// end of sorting

  FILE* ofile = fopen("AlphaP_derivate.out","w");

  if(ofile == NULL){
     fprintf(stderr,"Can't open output file\n");
     exit(1);
  }
  
  double * ap, * errorV;
  ap = (double*) calloc(i,sizeof(double));
  errorV = (double*) calloc(i,sizeof(double));
  
  for (int j=0; j<i; j++){
    
    if ( j == 0 )
       ap[j] = (V[j+1]-V[j])/(T[j+1]-T[j]);
    else if ( j == i-1 )
       ap[j] = (V[j]-V[j-1])/(T[j]-T[j-1]);
    else
       ap[j] = (V[j+1]-V[j-1])/(T[j+1]-T[j-1]);
       
    ap[j] *= 1./V[j];
    
    errorV[j] = sqrt( (V2[j]/n[j]) * 3 );   // see below

  }
  
  double error, errorDer, errorEst, dT;
  for (int j=2; j<i-2; j++){
    dT = (T[j+1] - T[j-1])/2.0;
    errorDer = fabs(V[j+2]/2.0 - V[j+1] + V[j-1] - V[j-2]/2.0)/(dT*dT*dT);  //third derivative
    errorDer *= dT*dT/6.0;
    errorEst = ( errorV[j+1] + errorV[j-1] )/(2.0*dT);
    error = sqrt( ap[j]*ap[j]*errorV[j]*errorV[j]/(V[j]*V[j]) + errorDer*errorDer/(V[j]*V[j]) + errorEst*errorEst/(V[j]*V[j])  );
    fprintf(ofile,"%lf %lf %lf\n",T[j],ap[j],error);    
  }
  
  fclose(ofile);

}

/*

err(V) = sqrt ( sigma^2/n * ( 1 + 2 * tau_A / delta_t ) )
Ref: Binder, Landau "A Guide to Monte Carlo Simulations in Statistical Mechanics"
3rd Edition, p. 32 

Assuming that in all simulations delta_t >~ tau_A, the formula is 
err(V) = sqrt ( sigma^2/n * 3 )

See Luis Coronas notebook date 06/04/2020

Modificado el programa 05/05/2022
Derivada central y error con la dervidada tercera

*/
