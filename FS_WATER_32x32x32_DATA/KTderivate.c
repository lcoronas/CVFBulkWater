#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

#define KT_NOT_VALUE -10.0
#define GAS_DENS 500

void compresibilidad(double * T, double * P, double * rho, double *Erho, int j, int num_points,
                     double * KT, double * EKT);

void error_compresibilidad(double * T, double * P, double * rho, double *Erho, int j, int num_points,
                     double * KT, double * EKT);
                     
void compresibilidad_log(double * T, double * P, double * rho, double *Erho, int j, int num_points,
                     double * KT, double * EKT);

void error_compresibilidad_log(double * T, double * P, double * rho, double *Erho, int j, int num_points,
                     double * KT, double * EKT);
                     
                     
int main (){

   int num_points = 2500;  // guess value. 1902 puntos para L=32 a dia 04/05/2022
   double * T, *P, *rho, *Erho, * KT, * EKT;
   
   T = (double*) calloc(num_points,sizeof(double));
   P = (double*) calloc(num_points,sizeof(double));
   rho = (double*) calloc(num_points,sizeof(double));
   Erho = (double*) calloc(num_points,sizeof(double));  // density error
   
   FILE* ifile = fopen("KT_input","r");

   if(ifile == NULL){
       fprintf(stderr,"KTderivate.c : Can't open KT_input file\n");
       exit(1);
   }
   
   int i=0;
   while(fscanf(ifile,"%lf %lf %lf %lf", &T[i], &P[i], &rho[i], &Erho[i])==4){
    i++;
    if (i == num_points){
       num_points += 100;
       T = (double*)realloc(T,num_points*sizeof(double));
       P = (double*)realloc(P,num_points*sizeof(double));
       rho = (double*)realloc(rho,num_points*sizeof(double));
       Erho = (double*)realloc(Erho,num_points*sizeof(double));
    }
  }
  
  fclose(ifile);
  
  num_points = i;
  
//  printf("num points=%d\n",i);
  
  KT = (double*) calloc(num_points,sizeof(double));
  EKT = (double*) calloc(num_points,sizeof(double));
  
  for (int j=0; j<num_points; j++)
//     compresibilidad(T,P,rho,Erho,j,num_points,KT,EKT);
       compresibilidad_log(T,P,rho,Erho,j,num_points,KT,EKT);
       
/*  for (int j=0; j<num_points; j++){
     if ( fabs (P[j] - 0.0 ) < 1E-5 ) printf("EKT[T=%lf] = %e; erho = %e\n",T[j],EKT[j],Erho[j]);
  }*/
     
  for (int j=0; j<num_points; j++)
//     error_compresibilidad(T,P,rho,Erho,j,num_points,KT,EKT);
     error_compresibilidad_log(T,P,rho,Erho,j,num_points,KT,EKT);
  
  /*for (int j=0; j<num_points; j++){
     if ( fabs (P[j] - 0.0 ) < 1E-5 ) printf("E_total [T=%e] = %e\n",T[j],EKT[j]);
  }*/
  
  // escribir output como isobaras
  double * list_isobars;
  int n_isobars=60; //guess value   L=32 hay unas 40 isobaras aprox (04/05/22)
  
  FILE * pfile = fopen("KT_list_isobars","r");
  
  if(pfile == NULL){
       fprintf(stderr,"KTderivate.c : Can't open KT_list_isobars file\n");
       exit(1);
  }
  
  list_isobars = (double*) calloc(n_isobars,sizeof(double));
  
  i=0;
  while(fscanf(ifile,"%lf", &list_isobars[i])==1){
    i++;
    if (i == n_isobars){
       n_isobars += 10;
       list_isobars = (double*)realloc(list_isobars,n_isobars*sizeof(double));
    }
  }
  fclose(pfile);
    
  n_isobars = i;
  
  FILE* ofile;
  char filename[100];
  for ( int j=0; j<n_isobars; j++){
     
     sprintf(filename,"P%lf_KTderivate.dat",list_isobars[j]);
     ofile = fopen(filename,"w"); 
     
     if(ofile == NULL){
       fprintf(stderr,"KTderivate.c : Can't create %s file\n",filename);
       exit(1);
     }
     
     for ( int k=0; k<num_points; k++ ){
     
         if ( KT[k] == KT_NOT_VALUE ) continue;
         
         if (  fabs(P[k] - list_isobars[j] ) < 1E-5 )
             fprintf(ofile,"%lf %lf %e %lf\n",T[k],KT[k],EKT[k],P[k]);
    
     }
  
     fclose(ofile);
  }
  
  // escribir output como isotermas
  double * list_isotherms;
  int n_isotherms=60; //guess value  
  
  pfile = fopen("KT_list_isotherms","r");
  
  if(pfile == NULL){
       fprintf(stderr,"KTderivate.c : Can't open KT_list_isotherms file\n");
       exit(1);
  }
  
  list_isotherms = (double*) calloc(n_isotherms,sizeof(double));
  
  i=0;
  while(fscanf(ifile,"%lf", &list_isotherms[i])==1){
    i++;
    if (i == n_isotherms){
       n_isotherms += 10;
       list_isotherms = (double*)realloc(list_isotherms,n_isotherms*sizeof(double));
    }
  }
  fclose(pfile);
    
  n_isotherms = i;
  
  for ( int j=0; j<n_isotherms; j++){
     
     sprintf(filename,"T%lf_density_isoth.dat",list_isotherms[j]);
     ofile = fopen(filename,"w"); 
     
     if(ofile == NULL){
       fprintf(stderr,"KTderivate.c : Can't create %s file\n",filename);
       exit(1);
     }
     
     for ( int k=0; k<num_points; k++ ){
     
         if (  fabs(T[k] - list_isotherms[j] ) < 1E-5 )
             fprintf(ofile,"%lf %lf %lf %lf\n",P[k],rho[k],Erho[k],T[k]);
    
     }
  
     fclose(ofile);
  }

  
}

void compresibilidad(double * T, double *P, double *rho, double *Erho, int j, int num_points,
                       double * KT, double * EKT){

    double temp=T[j];
    double pres=P[j];
    double dens=rho[j];
    
    double Pup, Pdn;   // closest pressures up and down P[j]
    int up_idx, dn_idx;
    
    Pup = 1E12;
    Pdn = -1E12;

    up_idx=-1;
    dn_idx=-1;
    
    for ( int i=0; i<num_points; i++ ){
    
       if ( j == i ) continue; 
       
       if ( fabs(T[i]-temp) < 1E-8 ){   // busqueda a lo largo de isoterma
       
           if ( P[i] > pres )
              if ( P[i] < Pup ){ 
                 Pup = P[i];
                 up_idx = i;
              }
              
           if ( P[i] < pres )
              if ( P[i] > Pdn ){ 
                 Pdn = P[i];
                 dn_idx = i;
              }
       
       }
    
    }
    
  /*  if (up_idx == -1 )
        printf("j=%d: %lf, %lf, Limit\n",j,Pdn,pres);
    else if (dn_idx == -1 )
        printf("j=%d: Limit, %lf, %lf\n",j,pres,Pup);
    else
        printf("j=%d: %lf, %lf, %lf\n",j,Pdn,pres,Pup);*/
        
    if ( up_idx == -1 && dn_idx == -1 ){  // puede pasar si la isoterma tiene solo un punto
         
         KT[j] = KT_NOT_VALUE;
         EKT[j] = KT_NOT_VALUE;
         
    }else if ( up_idx == -1 && dn_idx != -1 ){  // Limite isoterma (P maxima)
         
         KT[j] =  -dens*(1.0/dens - 1.0/rho[dn_idx])/(pres-Pdn); 

         EKT[j] = fabs(KT[j]*Erho[j]/dens) + fabs(dens)*( Erho[j]/(dens*dens) + Erho[dn_idx]/(rho[dn_idx]*rho[dn_idx]) )/(pres-Pdn);        

    
    }else if ( up_idx != -1 && dn_idx == -1 ){  // Limite isoterma (P minima)
         
         KT[j] =  -dens*(1.0/rho[up_idx] - 1.0/dens)/(Pup-pres);
         
         EKT[j] = fabs(KT[j]*Erho[j]/dens) + fabs(dens)*( Erho[up_idx]/(rho[up_idx]*rho[up_idx]) + Erho[j]/(dens*dens) )/(Pup-pres);
         
    }else{  
    
         KT[j] =  -dens*(1.0/rho[up_idx] - 1.0/rho[dn_idx])/(Pup-Pdn);
         
         EKT[j] = fabs(KT[j]*Erho[j]/dens) + fabs(dens)*( Erho[up_idx]/(rho[up_idx]*rho[up_idx]) + Erho[dn_idx]/(rho[dn_idx]*rho[dn_idx]) )/(Pup-Pdn);
         
    }
}

void error_compresibilidad (double * T, double *P, double *rho, double *Erho, int j, int num_points,
                       double * KT, double * EKT){

    double temp=T[j];
    double pres=P[j];
    double dens=rho[j];
    
    double Pup, Pdn, Pup2, Pdn2;   // closest pressures up and down P[j]
    int up_idx, dn_idx,up_idx2,dn_idx2;
    
    Pup = 1E12;
    Pdn = -1E12;
    Pup2 = 1E12;
    Pdn2 = -1E12;

    up_idx2=-1;
    up_idx=-1;
    dn_idx=-1;
    dn_idx2=-1;
    
    for ( int i=0; i<num_points; i++ ){
    
       if ( j == i ) continue; 
       
       if ( fabs(T[i]-temp) < 1E-8 ){   // busqueda a lo largo de isoterma
       
           if ( P[i] > pres )
              if ( P[i] < Pup ){ 
                 Pup = P[i];
                 up_idx = i;
              }
              
           if ( P[i] < pres )
              if ( P[i] > Pdn ){ 
                 Pdn = P[i];
                 dn_idx = i;
              }
       
       }
    
    }
    
    for ( int i=0; i<num_points; i++ ){
    
       if ( j == i ) continue; 
       
       if ( fabs(T[i]-temp) < 1E-8 ){   // busqueda a lo largo de isoterma
       
           if ( P[i] > Pup )
              if ( P[i] < Pup2 ){ 
                 Pup2 = P[i];
                 up_idx2 = i;
              }
              
           if ( P[i] < Pdn )
              if ( P[i] > Pdn2 ){ 
                 Pdn2 = P[i];
                 dn_idx2 = i;
              }
       
       }
    
    }
    
    double errDer, dP;
    dP = (Pup-Pdn)/2.0;
    if ( dP < 1E-5 ) 
      fprintf(stderr,"Warning: Pup = %lf ~ %lf = Pdn at isotherm %lf\n",Pup,Pdn,temp);
    
       
    if ( up_idx2 != -1 && dn_idx2 != -1 ){

       errDer = fabs(1.0/(2*rho[up_idx2]) - 1.0/rho[up_idx] + 1.0/rho[dn_idx] - 1.0/(2*rho[dn_idx2]))/(6.0*dP);
       
       if ( fabs (pres - 0.0 ) < 1E-5 ) 
          if ( fabs (temp - 200.227 ) < 1E-5 ){
               printf("\tT= %lf rho_[+2] = %lf, rho_[+1]=%lf, rho_[-1]=%lf, rho_[-2]=%lf, dP=%lf \n",temp,rho[up_idx2],
                        rho[up_idx], rho[dn_idx], rho[dn_idx2], dP);
               printf("\t\t%d %lf %lf %lf\n",up_idx2,T[up_idx2],P[up_idx2],rho[up_idx2]);            
          }
   //    printf("err=sqrt(%lf**2 + %lf**2 ), P=(%lf,%lf)\n",EKT[j],errDer,Pdn,Pup);
       EKT[j] = sqrt(EKT[j]*EKT[j] + errDer*errDer );
      
    } else if ( up_idx != -1 && dn_idx != -1 ){
    
       errDer = dP*fabs(KT[dn_idx]/rho[dn_idx] - 2*KT[j]/dens + KT[up_idx]/rho[up_idx])/6.0;
       
      if ( fabs (pres - 0.0 ) < 1E-5 ) 
          if ( fabs (temp - 200.227 ) < 1E-5 )
               printf("\tT= %lf KT[-1] = %lf, rho_[-1]=%lf, KT[0]=%lf, rho[0] = %lf, KT[+1]=%lf, rho_[+1]=%lf, dP=%lf \n",
               temp,KT[dn_idx],rho[dn_idx], KT[j], dens, KT[up_idx], rho[up_idx], dP);
       
   //    printf("err=sqrt(%lf**2 + %lf**2 ), P=(%lf,%lf)\n",EKT[j],errDer, Pdn, Pup);
       EKT[j] = sqrt(EKT[j]*EKT[j] + errDer*errDer );
       
    }
    
    if ( fabs (pres - 0.0 ) < 1E-5 ) printf("E_Der [T=%f] = %e\n",temp,errDer);
}

void compresibilidad_log (double * T, double *P, double *rho, double *Erho, int j, int num_points,
                       double * KT, double * EKT){

    double temp=T[j];
    double pres=P[j];
    double dens=rho[j];
    
    double Pup, Pdn;   // closest pressures up and down P[j]
    int up_idx, dn_idx;
    
    Pup = 1E12;
    Pdn = -1E12;

    up_idx=-1;
    dn_idx=-1;
    
    for ( int i=0; i<num_points; i++ ){
    
       if ( j == i ) continue; 
       
       if ( fabs(T[i]-temp) < 1E-8 ){   // busqueda a lo largo de isoterma
       
           if ( P[i] > pres )
              if ( P[i] < Pup ){ 
                 Pup = P[i];
                 up_idx = i;
              }
              
           if ( P[i] < pres )
              if ( P[i] > Pdn ){ 
                 Pdn = P[i];
                 dn_idx = i;
              }
       
       }
    
    }
    
  /*  if (up_idx == -1 )
        printf("j=%d: %lf, %lf, Limit\n",j,Pdn,pres);
    else if (dn_idx == -1 )
        printf("j=%d: Limit, %lf, %lf\n",j,pres,Pup);
    else
        printf("j=%d: %lf, %lf, %lf\n",j,Pdn,pres,Pup);*/
        
    if ( dens < GAS_DENS ){  // gas phase
    
         KT[j] = KT_NOT_VALUE;
         EKT[j] = KT_NOT_VALUE;      
           
    }else if ( up_idx == -1 && dn_idx == -1 ){  // puede pasar si la isoterma tiene solo un punto
         
         KT[j] = KT_NOT_VALUE;
         EKT[j] = KT_NOT_VALUE;
         
    }else if ( up_idx == -1 && dn_idx != -1 ){  // Limite isoterma (P maxima)
         
         if ( rho[dn_idx] < GAS_DENS ){
         
            KT[j] = KT_NOT_VALUE;
            EKT[j] = KT_NOT_VALUE;
            
         } else{
         
            KT[j] =  (log(dens) - log(rho[dn_idx]))/(pres-Pdn); 

            EKT[j] = fabs(Erho[j]/dens + Erho[dn_idx]/rho[dn_idx] )/(pres-Pdn);   

         }
    
    }else if ( up_idx != -1 && dn_idx == -1 ){  // Limite isoterma (P minima)

         if ( rho[up_idx] < GAS_DENS ){
         
            KT[j] = KT_NOT_VALUE;
            EKT[j] = KT_NOT_VALUE;
            
         } else{
         
         
            KT[j] =  (log(rho[up_idx]) - log(dens))/(Pup-pres); 

            EKT[j] = fabs(Erho[j]/dens + Erho[up_idx]/rho[up_idx] )/(Pup-pres); 
         }
           
    }else{  
    
         if ( rho[up_idx] < GAS_DENS || rho[dn_idx] < GAS_DENS ){
         
            KT[j] = KT_NOT_VALUE;
            EKT[j] = KT_NOT_VALUE;
            
         } else{
    
            KT[j] =  (log(rho[up_idx]) - log(rho[dn_idx]))/(Pup-Pdn); 

            EKT[j] = fabs(Erho[dn_idx]/rho[dn_idx] + Erho[up_idx]/rho[up_idx] )/(Pup-Pdn);
         }
    }
}

void error_compresibilidad_log (double * T, double *P, double *rho, double *Erho, int j, int num_points,
                       double * KT, double * EKT){

    double temp=T[j];
    double pres=P[j];
    double dens=rho[j];
    
    double Pup, Pdn, Pup2, Pdn2;   // closest pressures up and down P[j]
    int up_idx, dn_idx,up_idx2,dn_idx2;
    
    Pup = 1E12;
    Pdn = -1E12;
    Pup2 = 1E12;
    Pdn2 = -1E12;

    up_idx2=-1;
    up_idx=-1;
    dn_idx=-1;
    dn_idx2=-1;
    
    for ( int i=0; i<num_points; i++ ){
    
       if ( j == i ) continue; 
       
       if ( fabs(T[i]-temp) < 1E-8 ){   // busqueda a lo largo de isoterma
       
           if ( P[i] > pres )
              if ( P[i] < Pup ){ 
                 Pup = P[i];
                 up_idx = i;
              }
              
           if ( P[i] < pres )
              if ( P[i] > Pdn ){ 
                 Pdn = P[i];
                 dn_idx = i;
              }
       
       }
    
    }
    
    for ( int i=0; i<num_points; i++ ){
    
       if ( j == i ) continue; 
       
       if ( fabs(T[i]-temp) < 1E-8 ){   // busqueda a lo largo de isoterma
       
           if ( P[i] > Pup )
              if ( P[i] < Pup2 ){ 
                 Pup2 = P[i];
                 up_idx2 = i;
              }
              
           if ( P[i] < Pdn )
              if ( P[i] > Pdn2 ){ 
                 Pdn2 = P[i];
                 dn_idx2 = i;
              }
       
       }
    
    }
    
    double errDer, dP, dP2;
    dP = (Pup-Pdn)/2.0;
    dP2 = (Pup2-Pdn2)/4.0;
    
  //  if ( dP < 1E-5 ) 
    //  fprintf(stderr,"Warning: Pup = %lf ~ %lf = Pdn at isotherm %lf\n",Pup,Pdn,temp);
    
    if ( KT[j] == KT_NOT_VALUE || EKT[j] == KT_NOT_VALUE ){
          // do nothing
    }else if ( up_idx2 != -1 && dn_idx2 != -1 ){
    
       if ( rho[up_idx2] < GAS_DENS || rho[up_idx] < GAS_DENS || rho[dn_idx2] < GAS_DENS || rho[dn_idx] < GAS_DENS ){
       
          // do nothing

       } else {
       
          errDer = fabs(0.5*log(rho[up_idx2]) - log(rho[up_idx]) + log(rho[dn_idx]) - 0.5*log(rho[dn_idx2]))/(6.0*dP2);
       
   /*    if ( fabs (pres - 0.0 ) < 1E-5 ) 
          if ( fabs (temp - 200.227 ) < 1E-5 ){
               printf("\tT= %lf rho_[+2] = %lf, rho_[+1]=%lf, rho_[-1]=%lf, rho_[-2]=%lf, dP=%lf \n",temp,rho[up_idx2],
                        rho[up_idx], rho[dn_idx], rho[dn_idx2], dP);
               printf("\t\t%d %lf %lf %lf\n",up_idx2,T[up_idx2],P[up_idx2],rho[up_idx2]);            
          }*/
   //    printf("err=sqrt(%lf**2 + %lf**2 ), P=(%lf,%lf)\n",EKT[j],errDer,Pdn,Pup);
          EKT[j] = sqrt(EKT[j]*EKT[j] + errDer*errDer );
       }
      
    } else if ( up_idx != -1 && dn_idx != -1 ){
    
        if ( KT[up_idx] == KT_NOT_VALUE || KT[dn_idx] == KT_NOT_VALUE ){
       
          // do nothing

       } else {
    
          errDer = fabs(KT[dn_idx] - 2*KT[j] + KT[up_idx])/(6.0);
       
         /* if ( fabs (pres - 0.0 ) < 1E-5 ) 
             if ( fabs (temp - 200.227 ) < 1E-5 )
                  printf("\tT= %lf KT[-1] = %lf, rho_[-1]=%lf, KT[0]=%lf, rho[0] = %lf, KT[+1]=%lf, rho_[+1]=%lf, dP=%lf \n",
                  temp,KT[dn_idx],rho[dn_idx], KT[j], dens, KT[up_idx], rho[up_idx], dP);*/
       
   //    printf("err=sqrt(%lf**2 + %lf**2 ), P=(%lf,%lf)\n",EKT[j],errDer, Pdn, Pup);
          EKT[j] = sqrt(EKT[j]*EKT[j] + errDer*errDer );
          
       }
       
    }
    
 //   if ( fabs (pres - 0.0 ) < 1E-5 ) printf("E_Der [T=%f] = %e\n",temp,errDer);
}
  
