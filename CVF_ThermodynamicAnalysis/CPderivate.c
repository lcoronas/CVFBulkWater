#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

int main (){
  int num_temperatures = 100;
  int i;
  double *T, *H, *H2, *T_, *H_, *H2_;
  int *n, *n_;
  
  T = (double*) calloc(num_temperatures,sizeof(double));
  H = (double*) calloc(num_temperatures,sizeof(double));
  H2 = (double*) calloc(num_temperatures,sizeof(double));
  n = (int*) calloc(num_temperatures,sizeof(int));
  T_ = (double*) calloc(num_temperatures,sizeof(double));
  H_ = (double*) calloc(num_temperatures,sizeof(double));
  H2_ = (double*) calloc(num_temperatures,sizeof(double));
  n_ = (int*) calloc(num_temperatures,sizeof(int));

  FILE* ifile = fopen("CP_input","r");

  if(ifile == NULL){
       fprintf(stderr,"Can't open input file\n");
       exit(1);
  }

  i=0;
  while(fscanf(ifile,"%lf %lf %lf %d", &T_[i], &H_[i],&H2_[i],&n_[i])==4){
    i++;
    if (i == num_temperatures){
       num_temperatures += 100;
       T = (double*)realloc(T,num_temperatures*sizeof(double));
       H = (double*)realloc(H,num_temperatures*sizeof(double));
       H2 = (double*)realloc(H2,num_temperatures*sizeof(double));
       n = (int*)realloc(n,num_temperatures*sizeof(int));
       T_ = (double*)realloc(T_,num_temperatures*sizeof(double));
       H_ = (double*)realloc(H_,num_temperatures*sizeof(double));
       H2_ = (double*)realloc(H2_,num_temperatures*sizeof(double));
       n_ = (int*)realloc(n_,num_temperatures*sizeof(int));
    }
  }

/// sorting algorithm
{
  int index = 0;
  int check_index;
  double trial;
  bool *check;
  check = (bool*) calloc(i,sizeof(bool));

  for (int j=0;j<i;j++){    //i is the number of data stored in T_, etc.
     T[j]=0;
     H[j]=0;
     H2[j]=0;
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
        H[index] = H_[check_index];
        H2[index] = H2_[check_index];
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

  FILE* ofile = fopen("CP_derivate.out","w");

  if(ofile == NULL){
     fprintf(stderr,"Can't open output file\n");
     exit(1);
  }

  double *cp, *errorH;
  cp = (double *) calloc(i,sizeof(double));
  errorH = (double*) calloc(i,sizeof(double));
  
  for (int j=0; j<i; j++){
    if ( j == 0 )
       cp[j] = (H[j+1] - H[j])/(T[j+1]-T[j]);
    else if ( j== i-1 )
       cp[j] = (H[j] - H[j-1])/(T[j]-T[j-1]);
    else 
       cp[j] = (H[j+1] - H[j-1])/(T[j+1]-T[j-1]);
       
    errorH[j] = sqrt( (H2[j]/n[j]) * 3 );
    
  }
  
  double dT, errorDer,errorEst,error;
  for (int j=2; j<i-2; j++){
    dT = (T[j+1] - T[j-1])/2.0;
    errorDer = fabs(H[j+2]/2.0 - H[j+1] + H[j-1] - H[j-2]/2.0)/(dT*dT*dT);  //third derivative
    errorDer *= dT*dT/6.0;
    errorEst = ( errorH[j+1] + errorH[j-1] )/(2.0*dT);
    error = sqrt( errorDer*errorDer + errorEst*errorEst );
    fprintf(ofile,"%lf %lf %lf\n",T[j],cp[j],error);
  }
  fclose(ofile);

}
