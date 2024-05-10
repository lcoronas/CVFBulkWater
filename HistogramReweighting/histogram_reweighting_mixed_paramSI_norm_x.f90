! Programa que realiza, a tamaño fijo, un barrido en presiones, temperaturas y ángulos
! Tambien permite ajustar el parámetro B (no debería exceder 1.1 ni ser menor de 0.9, significa un error
! del +/- 10% en B (Sci Rep Bianco Franzese 2014)
! Calcula el error 
! Imprime histograma

!Version 5: 10 Marzo 2022

!Version Binder: 14 Julio 2022, a partir de la version 5 (11 Abril 2022)
!Añade un archivo de output, cumulante de Binder del parametro de orden

!Version Binder: 01 Octubre 2022, a partir de la version 14 Julio 2022
!Añade archivos de output, cumulante de Binder de energia y densidad
!Bug corregido U = 1 - <M^4>/(3*<M²>^2), antes U = 1 - <M^4>/(3*<M²>)

!Correccion 09 Octubre 2022
!Cambio el calculo de M=rho + s e por M=rho cos(alpha)+ e sen(alpha)

!10 Octubre 2022. Comprobado que  M=rho + s e por M=rho y cos(alpha)+ e sen(alpha)
!son equivalentes y que la unica solucion es s negativo,
!recupero M=rho + s e 

!11 Octubre 2022. Paso a calcular  M=rho + s e con rho en unds reales.
! e -> kJ/mol; rho -> kg/m^3; s -> (kg*mol)/(m^3*kJ)
! Este codigo toma como input T, P, rho y h en unds internas y las pasa a unds reales (e=h-P*v)
! Así evito perder precision en T y P por redondeo, lo cual puede afectar a los resultados
! Los archivos de output saldran con los valores en los dos sistemas

!28 Octubre 2022. Añado un factor de escala para que T, P, e , rho sean proximos a 1.
!Para no perder "lo ganado", e y rho se deben dividir por el mismo factor, de forma que 
!el angulo en el plano e,rho NO cambie.

!10 Noviembre 2022. Busco un factor que normalice las unidades de rho
!do forma que la distribución ajuste mejor a Ising

program main
implicit none
  integer :: i,j,k,istep,w,step,angle_step,flag, io
  integer :: indiceX,indiceY,indiceE,indiceDens     !contadores
  integer :: n_simulaciones, n_bin, n_bin_mixed, n_bin2D
  character(100) :: filename, param
  character(1) :: blank_ !blank space at input file
  real*8,dimension(:),allocatable :: That, Phat    !internal units
  real*8,dimension(:),allocatable :: T, P          !international units
  integer :: L,N
  integer,dimension(:),allocatable :: medidas
  integer,dimension(:,:), allocatable :: memo_zero
  integer :: medidas_max
  real*8 :: E, Emin, Emax, Dens, Densmin, Densmax, V, tam_binE, tam_binDens, rangeE, rangeDens
  real*8 :: tam_binE2D, tam_binDens2D
  real*8 :: X, Xmin, Xmax, Y, Ymin, Ymax, rangeX, tam_binX, rangeY, tam_binY, delta_X
  real*8 :: xx, maxlog_P, log_Z, log_hist, hist, max_arg, norma, media, var, a1, a2
  real*8,dimension(:,:,:),allocatable :: datos_originales
  real*8,dimension(:,:),allocatable :: energy_original
  real*8,dimension(:,:),allocatable :: density_original, X_original, X_integral
  real*8,dimension(:,:),allocatable :: log_P, histogram_reweighted_grid, histogram_reweighted_rotated
  real*8,dimension(:),allocatable :: ctes, ctes2, arg_exp, log_medidas, histogram_reweighted_mixed
  real*8,dimension(:),allocatable :: mixed_parameter, mixed_parameter_
  real*8 :: beta, beta1, pres1, T1,P1, S, alpha
  real*8 :: Tmin, Tmax, deltaT, Pmin, Pmax, deltaP, Smin, Smax, deltaS, bin_Factor, T_ctes, P_ctes, tmp_ctes
  real*8 :: Tmin_SI, Tmax_SI
  real*8 :: peak2D, peak3D, xs
  real*8 :: kl_2D, liu_2D, kl_3D, liu_3D
  real*8 :: kl_2Dmin, liu_2Dmin, kl_3Dmin, liu_3Dmin
  real*8 :: alpha_kl_2Dmin, alpha_liu_2Dmin, alpha_kl_3Dmin, alpha_liu_3Dmin, alpha_optim
  real*8 :: T_liu_2Dmin, T_liu_3Dmin, P_liu_2Dmin, P_liu_3Dmin, B_liu_2Dmin, B_liu_3Dmin
  real*8 :: zero, zeromin, alpha_zeromin, T_zeromin, P_zeromin, B_zeromin
  real*8 :: height, heightmin, alpha_heightmin, T_heightmin, P_heightmin, B_heightmin
  real*8 :: zero_height, zero_heightmin, alpha_zero_heightmin, T_zero_heightmin, P_zero_heightmin, B_zero_heightmin
  real*8 :: B
  real*8, dimension(:), allocatable :: err_hist_original, err_hist_reweighted
  real*8, dimension(:), allocatable :: error_hist, error_X
  real*8 :: err_distance, err_min_distance
  real*8 :: binder_result, binder_ene, binder_rho
  integer :: err_k, err_counter
  logical :: plot_Confront, file_exists, read_flag
  integer :: contador
  
  !units conversion
  real*8 :: Rhos, Rhoi, Hs, Hi, Ts, Ti, Ps, Pi 
  Ts=140.5678    ! convert to K
  Ti=185.4676
  Ps=4694.59     ! convert to bar
  Pi=-2111.5552
  Rhos=1527.3    ! convert to kg/m^3
  Rhoi=-23.1017
  Hs=11.272132   ! convert enthalpy to kJ/mol  (0.018*Ts*CPs)  12/05/2022
  Hi=9.028914    !  (0.018*Ts*CPi)

  !read input in internal units
  open(unit=10,file="input",status="old")
  read(10,'(A15,I2)') param, n_simulaciones
  read(10,'(A5,I5)') param, step
  read(10,'(A2,I2)') param, L
  read(10,'(A6,F10.8)') param, Tmin
  read(10,'(A6,F10.8)') param, Tmax
  read(10,'(A8,F12.10)') param, deltaT
  read(10,'(A6,F7.5)') param, Pmin
  read(10,'(A6,F7.5)') param, Pmax
  read(10,'(A8,F10.8)') param, deltaP
  read(10,'(A6,F8.5)') param, Smin
  read(10,'(A6,F8.5)') param, Smax
  read(10,'(A8,F9.7)') param, deltaS
  read(10,'(A8,F7.4)') param, bin_Factor
  allocate(T(n_simulaciones),P(n_simulaciones),medidas(n_simulaciones) , log_medidas(n_simulaciones))
  allocate(That(n_simulaciones),Phat(n_simulaciones))
  do i=1,n_simulaciones
    read(10,'(F9.7,A1,F6.4,A1,I9)') That(i), blank_, Phat(i), blank_, medidas(i)
  end do
  close(10)
  
  T(:) = That(:)*Ts+Ti   !T in K
  P(:) = Phat(:)*Ps+Pi  !P in bar
  
  if ( Tmax <= Tmin ) then
     print*,"Error: Tmax <= Tmin"
     call exit(1)
  end if
  
  if ( Pmax <= Pmin ) then
     print*,"Error: Pmax <= Pmin"
     call exit(1)
  end if
  
  if ( Smax <= Smin ) then
     print*,"Error: Smax <= Smin"
     call exit(1)
  end if
  
  Tmin_SI = Ts*Tmin+Ti
  Tmax_SI = Ts*Tmax+Ti

 ! P1=Phat(1)
  N = L*L*L

  log_medidas(:) = log(real(medidas(:)))

  !Abro archivos de datos y miro el rango rho-E

  flag=1
  do i=1,n_simulaciones

    if ( L<10 ) then
       if ( Phat(i) >= 0 ) then
          write(filename,'("data_L",I1,"_T",F9.7,"_P",F6.4)') L,That(i),Phat(i)
       else
          write(filename,'("data_L",I1,"_T",F9.7,"_P",F6.3)') L,That(i),Phat(i)
       end if
    else
       if ( Phat(i) >= 0 ) then 
          write(filename,'("data_L",I2,"_T",F9.7,"_P",F6.4)') L,That(i),Phat(i)
       else
          write(filename,'("data_L",I2,"_T",F9.7,"_P",F6.3)') L,That(i),Phat(i)
       end if
    end if

    open(unit=10,file=filename,status="old")
      do j=1,medidas(i)
         read(10,*) E, Dens  ! enthalpy h=H/N and density rho=N/V in internal units
         Dens = rhos*Dens+rhoi  ! Dens in kg/m^3
         E = Hs*E+Hi*That(i)-1.8*P(i)/Dens ! Energy in kJ/mol
         if (i==1 .and. j==1) then
            Emax = E
            Emin = E
            Densmax = Dens
            Densmin = Dens
         else
            if (Emax < E) then 
                Emax = E
            end if
            if (Emin > E) then 
                Emin = E
            end if
            if (Densmax < Dens) then 
                Densmax = Dens 
            end if
            if (Densmin > Dens) then 
                Densmin = Dens 
            end if
         end if
      end do
    close(10)
    if ( i == 1) then 
         medidas_max = medidas(i)
    else if (medidas_max < medidas(i)) then 
         medidas_max = medidas(i)
    end if
    if ( flag == 1 .and. (.not. medidas(i) == medidas(1)) ) then
	print*, "Warning: Histogramas de entrada con distinto numero de datos"
        flag = 0
    end if    
  end do

  rangeE = Emax - Emin
  rangeDens = Densmax - Densmin

  n_bin = int( sqrt( real(medidas_max) )/2.0 )

  tam_binE = rangeE/n_bin
  tam_binDens = rangeDens/n_bin
  
! Ampliación de histogrmas (para dejar bins vacios)	
  n_bin = n_bin+20
  rangeE = rangeE + 20*tam_binE
  Emin = Emin - 10*tam_binE
  rangeDens = rangeDens + 20*tam_binDens
  Densmin = Densmin - 10*tam_binDens
  
  n_bin2D = int(n_bin*10*bin_Factor)
 ! n_bin2D = 200
  tam_binE2D = rangeE/n_bin2D
  tam_binDens2D = rangeDens/n_bin2D
  
  !Calculo los histogramas a partir de los datos originales
 
  allocate(datos_originales(n_simulaciones,n_bin2D,n_bin2D))
!  allocate(energy_original(n_simulaciones,n_bin))
!  allocate(density_original(n_simulaciones,n_bin))
!  allocate(X_original(n_simulaciones,n_bin))
    
  datos_originales(:,:,:) = 0
!  energy_original(:,:) = 0
!  density_original(:,:) = 0
!  X_original(:,:) = 0
  
  do i=1,n_simulaciones
    if ( L<10 ) then
       if ( Phat(i) >= 0 ) then
          write(filename,'("data_L",I1,"_T",F9.7,"_P",F6.4)') L,That(i),Phat(i)
       else
          write(filename,'("data_L",I1,"_T",F9.7,"_P",F6.3)') L,That(i),Phat(i)
       end if
    else
       if ( Phat(i) >= 0 ) then 
          write(filename,'("data_L",I2,"_T",F9.7,"_P",F6.4)') L,That(i),Phat(i)
       else
          write(filename,'("data_L",I2,"_T",F9.7,"_P",F6.3)') L,That(i),Phat(i)
       end if
    end if
    open(unit=10,file=filename,status="old")
    norma = medidas(i)
    do j=1,medidas(i)
        read(10,*) E, Dens  ! enthalpy h=H/N and density rho=N/V in internal units
        Dens = rhos*Dens+rhoi  ! Dens in kg/m^3
        E = Hs*E+Hi*That(i)- 1.8*P(i)/Dens ! Energy in kJ/mol
    !   indiceX = int((Dens + Smin * E - Xmin)/tam_binX)+1
    !   indiceE = int((E-Emin)/tam_binE)+1
    !   indiceDens = int((Dens-Densmin)/tam_binDens)+1
    !   if(indiceX >= 1 .and. indiceX <= n_bin) then
    !     X_original(i,indiceX) = X_original(i,indiceX) + 1
    !   end if
    !   if(indiceE >= 1 .and. indiceE <= n_bin) then    !E=Emax queda fuera ... (igual con Densmax)
    !     energy_original(i,indiceE) = energy_original(i,indiceE) + 1
    !   end if
    !   if(indiceDens >=1 .and. indiceDens <= n_bin) then
    !     density_original(i,indiceDens) = density_original(i,indiceDens) + 1
    !   end if
       indiceE = int((E-Emin)/tam_binE2D)+1
       indiceDens = int((Dens-Densmin)/tam_binDens2D)+1
       if ( indiceE >= 1 .and. indiceE <= n_bin2D .and. indiceDens >= 1 .and. indiceDens <= n_bin2D ) then
            datos_originales(i,indiceE,indiceDens) = datos_originales(i,indiceE,indiceDens) + 1
       else
            norma = norma-1  !data out of histogram
       end if
    end do
    
    ! Normalizo el histograma de entrada. Posiblemente es innecesario...
    norma = norma*tam_binE2D*tam_binDens2D
    do j=1,n_bin2D
      do k=1,n_bin2D
        datos_originales(i,j,k) = datos_originales(i,j,k)/norma
      end do
    end do
    
 !   norma = 0
 !   do j=1,n_bin
 !       norma = norma + energy_original(i,j)
 !   end do
 !   do j=1,n_bin
 !     energy_original(i,j) = energy_original(i,j)/(norma*tam_binE)
 !   end do
    
 !   norma = 0
 !   do j=1,n_bin
 !       norma = norma + density_original(i,j)
 !   end do
 !   do j=1,n_bin
 !     density_original(i,j) = density_original(i,j)/(norma*tam_binDens)
 !   end do
    
 !   norma = 0
 !   do j=1,n_bin
 !       norma = norma + X_original(i,j)
 !   end do
 !   do j=1,n_bin
 !     X_original(i,j) = X_original(i,j)/(norma*tam_binX)
 !   end do

    close(10)
  end do

   ! do j=1,n_bin2D
   !   do k=1,n_bin2D
   !     if (datos_originales(3,j,k)>0) then
   !      print*, Emin+(j-0.5)*tam_binE2D, Densmin+(k-0.5)*tam_binDens2D, datos_originales(3,j,k)
   !     end if
   !   end do
   ! end do

!  do i=1,n_simulaciones
!    write(filename,'("Histogram_Dato_Ene_L",I2,"_T",F8.6,"_P",F4.2)') L,T(i),P(i)
!    open(unit=11,file=filename)
!    write(filename,'("Histogram_Dato_Rho_L",I2,"_T",F8.6,"_P",F4.2)') L,T(i),P(i)
!    open(unit=12,file=filename)
!    write(filename,'("Histogram_Dato_OP_L",I2,"_T",F8.6,"_P",F4.2)') L,T(i),P(i)
!    open(unit=13,file=filename)
!    do j=1,n_bin
!      write(11,*) Emin + (j-0.5)*tam_binE, energy_original(i,j)
!      write(12,*) Densmin + (j-0.5)*tam_binDens, density_original(i,j)
!      write(13,*) Xmin + (j-0.5)*tam_binX, X_original(i,j)
!    end do
!  close(11)
!  close(12)
!  close(13)
!  end do
  
 ! do i=1,n_simulaciones
 !   write(filename,'("Histogram_Dato_Ene_Rho_L",I2,"_T",F8.6,"_P",F4.2)') L,T(i),P(i)
 !   open(unit=10,file=filename)
 !   do j=1,n_bin2D
 !     do k=1,n_bin2D
 !        write(10,*) Emin + (j-0.5)*tam_binE, Densmin + (k-0.5)*tam_binDens, datos_originales(i,j,k)
 !     end do
 !     write(10,*) " "
 !   end do
 !   close(10)
 ! end do
  
 ! n_bin_mixed = int(sqrt(n_bin2D*1.0))
  
 ! allocate (X_integral(n_simulaciones,n_bin_mixed))
 ! X_integral(:,:) = 0
  
 ! tam_binX = rangeX/n_bin_mixed

 ! do i=1,n_simulaciones
 ! print*, "Integro simul",i
 ! write(filename,'("Histogram_Dato_Integral_OP_L",I2,"_T",F8.6,"_P",F4.2)') L,T(i),P(i)
 ! open(unit=14,file=filename)
 !   do j=1,n_bin2D
 !     E = Emin + (j-0.5)*tam_binE2D
 !     do k=1,n_bin2D
 !       Dens = Densmin + (k-0.5)*tam_binDens2D
 !       indiceX = int((Dens + Smin * E - Xmin)/tam_binX)+1      
 !       if ( datos_originales(i,j,k) > 0 ) then
 !         if ( indiceX >= 1 .and. indiceX <= n_bin_mixed ) then
 !            X_integral(i,indiceX) = X_integral(i,indiceX) + datos_originales(i,j,k)
 !         else
 !            print*, "Warning: Out of Histogram", indiceX, n_bin_mixed
 !         end if
 !       end if
 !     end do
 !   end do
    
 !   norma = 0
 !   do j=1,n_bin_mixed
 !       norma = norma + X_integral(i,j)
 !   end do

 !   do j=1,n_bin_mixed
 !     X_integral(i,j) = X_integral(i,j)/(norma*tam_binX)
 !     write(14,*) Xmin + (j-0.5)*tam_binX, X_integral(i,j)
 !   end do 
 ! close(14)
 ! end do


  !Histogram reweighting. Obtención de las ctes

  allocate(ctes(n_simulaciones), ctes2(n_simulaciones))
  allocate(memo_zero(n_bin2D,n_bin2D), log_P(n_bin2D,n_bin2D) )
  allocate(arg_exp(n_simulaciones)) 

  if ( L < 10 ) then
    write(filename,'("Constants_L",I1)') L
  else
    write(filename,'("Constants_L",I2)') L
  end if
  INQUIRE(FILE=filename, EXIST=file_exists)
  
  read_flag = .false.
  
  if ( file_exists ) then
    
    read_flag = .true.
    
    write(0,*) "Leyendo el archivo de constantes..."
  
    open(unit=10,file=filename,status="old")
    
    i=1
    do  ! quiero asegurarme de que el archivo tiene el mismo numero de temperaturas
  
         read(10,*,iostat=io) T_ctes, P_ctes, tmp_ctes
         
         if (io/=0) then 
             EXIT
         end if
         
         if ( i > n_simulaciones ) then
            read_flag = .false.
            EXIT
         end if
         
         ctes(i) = tmp_ctes
         ctes2(i) = tmp_ctes
                
         if ( dabs(T_ctes-That(i)) > 1E-10 ) then
             read_flag = .false.
         end if
         
         if ( dabs(P_ctes-Phat(i)) > 1E-10 ) then
             read_flag = .false.
         end if
         
         i=i+1
         
    end do
  
    close(10)
    
    if ( .not. (i-1) == n_simulaciones ) then
       read_flag = .false.
    end if
    
  end if
  
  if ( .not. read_flag  ) then
  
    if ( file_exists ) then
       print*, "WARNING: Some error occured while reading constants file"
       print*, "I recalculate the constants, but existing file will not be updated"
    end if
  
    ctes(:) = 0.0
    ctes2(:) = 0.0
  
    write(0,*) "Calculando las constantes..."

    do istep=1,step
  
      if ( mod(istep,50) == 0 ) then
         print*,"Step ", istep
      end if
  
      ctes(:) = ctes2(:)

      do k=1,n_simulaciones
           call Calcula_Histograma(T(k),P(k),log_P,log_Z)
           ctes2(k) = log_Z
      end do

      ctes2(:) = ctes2(:) - ctes2(1)

    end do
    
    do k=1,n_simulaciones
      if (dabs(ctes(k) - ctes2(k)) > 1E-07 ) then
        write(0,*) "Las constantes no han convergido"   !salida por stderr
        write(0,*) "Usar un valor de step >", step
        call exit(1)
      end if
    end do  
 
  end if

  write(0,*) "Las constantes son:"
  write(0,*) ctes2(:)

  ctes(:) = ctes2(:)
  
 
  if ( .not. file_exists ) then
  
     open(unit=10,file=filename,status="new")
     
     do i=1,n_simulaciones
     
        write(10,*) That(i), Phat(i), ctes(i)
     
     end do
     
     close(10)
  
  end if  

!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!
!  do k=1,5
!   call Calcula_Histograma(T(k),P(k),log_P,log_Z)

!   norma=0
!   do i=1,n_bin
!     do j=1,n_bin
!        if(memo_zero(i,j)==0) then
!          norma = norma + exp(log_P(i,j)-log_Z)
!       end if
!     end do
!   end do

!   write(filename,'("Histogram_Reweighted_ene_rho_L",I2,"_T",F9.5,"_P",F5.1)') L,T(k),P(k)
!   open(unit=14,file=filename)
!   do i=1,n_bin
!   do j=1,n_bin
 !     if ( memo_zero(i,j) == 0) then
  !      write(14,*) Emin+(i-1)*tam_binE, Densmin+(j-1)*tam_binDens, exp(log_P(i,j)-log_Z) / (norma*tam_binE*tam_binDens)
   !   end if
    !end do
 !  end do
  !  close(14)

!  alpha=0
 !  do while(alpha < 1.5)
!     call Calcula_Histograma_Rotado(alpha)
!     alpha=alpha+0.1
!   end do

! end do
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!

  !Cálculo de los picos de las distribuciones

  call Ising_Peak(peak2D, peak3D)

 ! pi_griego = 3.1415926536
  
  if (P(1) >= 1000.0 ) then
  
   if ( L < 10 ) then
     write(filename,'("Binder_L",I1,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,".dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
     open(unit=20,file=filename)
   else
     write(filename,'("Binder_L",I2,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,".dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
     open(unit=20,file=filename)
   end if
  
   if ( L < 10 ) then
     write(filename,'("Binder_L",I1,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_energy_density.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI
     open(unit=21,file=filename)
   else
     write(filename,'("Binder_L",I2,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_energy_density.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI
     open(unit=21,file=filename)
   end if
   
  else if (P(1) >= 0.0) then
  
   if ( L < 10 ) then
     write(filename,'("Binder_L",I1,"_P",F5.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,".dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
     open(unit=20,file=filename)
   else
     write(filename,'("Binder_L",I2,"_P",F5.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,".dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
     open(unit=20,file=filename)
   end if
  
   if ( L < 10 ) then
     write(filename,'("Binder_L",I1,"_P",F5.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_energy_density.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI
     open(unit=21,file=filename)
   else
     write(filename,'("Binder_L",I2,"_P",F5.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_energy_density.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI
     open(unit=21,file=filename)
   end if
   
  else if ( P(1) > -1000.0 ) then !P negative
  
   if ( L < 10 ) then
     write(filename,'("Binder_L",I1,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,".dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
     open(unit=20,file=filename)
   else
     write(filename,'("Binder_L",I2,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,".dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
     open(unit=20,file=filename)
   end if
  
   if ( L < 10 ) then
     write(filename,'("Binder_L",I1,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_energy_density.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI
     open(unit=21,file=filename)
   else
     write(filename,'("Binder_L",I2,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_energy_density.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI
     open(unit=21,file=filename)
   end if  
   
  else

   if ( L < 10 ) then
     write(filename,'("Binder_L",I1,"_P",F7.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,".dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
     open(unit=20,file=filename)
   else
     write(filename,'("Binder_L",I2,"_P",F7.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,".dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
     open(unit=20,file=filename)
   end if
  
   if ( L < 10 ) then
     write(filename,'("Binder_L",I1,"_P",F7.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_energy_density.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI
     open(unit=21,file=filename)
   else
     write(filename,'("Binder_L",I2,"_P",F7.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_energy_density.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI
     open(unit=21,file=filename)
   end if
  
  end if
  
  plot_Confront = .true.
  if (plot_Confront) then
  !Archivos output
   if ( L < 10 )  then

 !   write(filename,'("Confront_L",I1,"_P",F5.3,"_Tm",F9.5,"_TM",F9.5,"_Ms",F9.6,"_MS",F9.6,"Ising2D_Liu.dat")') &
 !                 L,P1,Tmin,Tmax,Smin,Smax
 !   open(unit=10,file=filename)

 !   write(filename,'("Confront_L",I1,"_P",F5.3,"_Tm",F9.5,"_TM",F9.5,"_Ms",F9.6,"_MS",F9.6,"Ising2D_KL.dat")') &
 !	  L,P1,Tmin,Tmax,Smin,Smax
 !   open(unit=11,file=filename)
     if (P(1) <= -1000.0) then
    
write(filename,'("Confront_L",I1,"_P",F7.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_Liu.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=12,file=filename)

        write(filename,'("Confront_L",I1,"_P",F7.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_KL.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=13,file=filename) 
        
     else if (P(1) < 0.0 .or. P(1) >= 1000.0) then
        write(filename,'("Confront_L",I1,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_Liu.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=12,file=filename)

        write(filename,'("Confront_L",I1,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_KL.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=13,file=filename) 

     else 
        write(filename,'("Confront_L",I1,"_P",F5.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_Liu.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=12,file=filename)

        write(filename,'("Confront_L",I1,"_P",F5.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_KL.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=13,file=filename) 
        
     end if

    else !L>10
    
 !   write(filename,'("Confront_L",I1,"_P",F5.3,"_Tm",F9.5,"_TM",F9.5,"_Ms",F9.6,"_MS",F9.6,"Ising2D_Liu.dat")') &
 !                 L,P1,Tmin,Tmax,Smin,Smax
 !   open(unit=10,file=filename)

 !   write(filename,'("Confront_L",I1,"_P",F5.3,"_Tm",F9.5,"_TM",F9.5,"_Ms",F9.6,"_MS",F9.6,"Ising2D_KL.dat")') &
 !	  L,P1,Tmin,Tmax,Smin,Smax
 !   open(unit=11,file=filename)
     if (P(1) <= -1000.0) then
    
write(filename,'("Confront_L",I2,"_P",F7.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_Liu.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=12,file=filename)

        write(filename,'("Confront_L",I2,"_P",F7.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_KL.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=13,file=filename) 
        
     else if (P(1) < 0.0 .or. P(1) >= 1000.0) then
        write(filename,'("Confront_L",I2,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_Liu.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=12,file=filename)

        write(filename,'("Confront_L",I2,"_P",F6.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_KL.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=13,file=filename) 

     else 
        write(filename,'("Confront_L",I2,"_P",F5.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_Liu.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=12,file=filename)

        write(filename,'("Confront_L",I2,"_P",F5.1,"bar_Tm",F9.5,"K_TM",F9.5,"K_mA",F8.6,"_MA",F8.6,"Ising3D_KL.dat")') &
                  L,P(1),Tmin_SI,Tmax_SI,Smin,Smax
        open(unit=13,file=filename) 
        
     end if  
     
    end if
  end if !plot Confront

  !Cálculo del histograma Repesado

  allocate(histogram_reweighted_grid(n_bin2D,n_bin2D))

  n_bin_mixed = int(sqrt(n_bin2D*1.0))
!   n_bin_mixed = n_bin2D
   
  allocate (histogram_reweighted_mixed(n_bin_mixed))
  allocate (mixed_parameter(n_bin_mixed),mixed_parameter_(n_bin_mixed))   !eje X del histograma
!  allocate (histogram_reweighted_rotated(n_bin_mixed,n_bin_mixed))

  P1 = Pmin

  liu_2Dmin = 1E12
  liu_3Dmin = 1E12
  zeromin = 1E12
  heightmin = 1E12
  zero_heightmin = 1E12

do while (P1 <= Pmax)

  T1 = Tmin

do while (T1 <= Tmax)

  call Calcula_Histograma(T1*Ts+Ti,P1*Ps+Pi,log_P,log_Z)
  
  norma = 0
  do i=1,n_bin2D
     do w=1,n_bin2D
        if (memo_zero(i,w) == 0) then
           hist = exp(log_P(i,w)-log_Z)
           histogram_reweighted_grid(i,w) = hist
           norma = norma + hist
        end if
     end do
  end do
  norma = norma*tam_binE2D*tam_binDens2D

  do i=1,n_bin2D
    do w=1,n_bin2D
      if (memo_zero(i,w) == 0) then
        histogram_reweighted_grid(i,w) = histogram_reweighted_grid(i,w)/norma
      end if
    end do
  end do
  
  call binder_energy_density(binder_ene,binder_rho)
  write(21,*) T1*Ts+Ti, P1*Ps+Pi, binder_ene, binder_rho
  
  S = Smin

  do while (S < Smax)

   !   alpha = atan(S)
   !   alpha = -alpha		! S es la tangente de la recta que une los centros de los pozos de energia libre
				! Giro un angulo -alpha para alinear la distribución al eje X

!Roto los ejes y calculo el histograma
      alpha = S
      !call Calcula_Histograma_Rotado(alpha)
      call Calcula_Histograma_Rotado_normX(alpha)
       
  ! if (alpha < 0) then
  !    write(filename,'("Histogram_Reweighted_2D_L",I2,"_T",F9.5,"_P",F5.1,"_S",F9.6)') L,T1,P1,alpha
  !  else if (alpha < 10) then
  !    write(filename,'("Histogram_Reweighted_2D_L",I2,"_T",F9.5,"_P",F5.1,"_S",F8.6)') L,T1,P1,alpha
  !  else 
  !    write(filename,'("Histogram_Reweighted_2D_L",I2,"_T",F9.5,"_P",F5.1,"_S",F9.6)') L,T1,P1,alpha
  !  end if
  !  open(unit=14,file=filename)
  !  write(14,*) "#Data simulations (Temperature): ", T(:)
  !  do i=1,n_bin_mixed
  !    do w=1,n_bin_mixed
  !    if (.not. histogram_reweighted_rotated(i,w) == 0) then
  !      write(14,*) Xmin+(i-0.5)*tam_binX, Ymin+(w-0.5)*tam_binY, histogram_reweighted_rotated(i,w)
  !    end if
  !    end do
  !  end do
  !  close(14)

 ! stop
 
      call binder ( binder_result )  ! Binder calculated over M distribution, NOT the M distribution
      	                             ! with average=0 and variance=1
      write(20,*) T1*Ts+Ti, P1*Ps+Pi, S, binder_result

      mixed_parameter_(:) = mixed_parameter(:)

     ! B = 1-0.1

 !     write(filename,'("BvsLiu3D_L",I2,"_T",F9.5,"_P",F5.1,"_S",F9.6)') L,T1,P1,alpha
 !     open(unit=20,file=filename)
      B=1
 !     do while (B < 1+0.1)
      do while ( B < 1.001)

        mixed_parameter(:) = B*mixed_parameter_(:)  !rescalado de la coordenada x=B(M-Mc)
        
        call ConfrontIsing (kl_2D, kl_3D, liu_2D, liu_3D)

 !       write(20,*) B, liu_3D
        call ZeroTest (zero)
        call HeightTest (height)
        zero_height = zero + height

        if ( plot_Confront ) then
  	!  write(10,*) T1, P1, S, liu_2D
   	!  write(11,*) T1, P1, S, kl_2D
   	  write(12,*) T1*Ts+Ti, P1*Ps+Pi, S, liu_3D
          write(13,*) T1*Ts+Ti, P1*Ps+Pi, S, kl_3D
        end if


        if ( liu_2Dmin > liu_2D) then
           liu_2Dmin = liu_2D
           alpha_liu_2Dmin = alpha
           T_liu_2Dmin = T1*Ts+Ti
           P_liu_2Dmin = P1*Ps+Pi
           B_liu_2Dmin = B
        end if

        if ( liu_3Dmin > liu_3D) then
           liu_3Dmin = liu_3D
           alpha_liu_3Dmin = alpha
           T_liu_3Dmin = T1*Ts+Ti
           P_liu_3Dmin = P1*Ps+Pi
           B_liu_3Dmin = B
        end if

        if ( zeromin > zero) then
           zeromin = zero
           alpha_zeromin = alpha
           T_zeromin = T1*Ts+Ti
           P_zeromin = P1*Ps+Pi
           B_zeromin = B
        end if
        
        if ( heightmin > height) then
           heightmin = height
           alpha_heightmin = alpha
           T_heightmin = T1*Ts+Ti
           P_heightmin = P1*Ps+Pi
           B_heightmin = B
        end if

        if ( zero_heightmin > zero_height) then
           zero_heightmin = zero_height
           alpha_zero_heightmin = alpha
           T_zero_heightmin = T1*Ts+Ti
           P_zero_heightmin = P1*Ps+Pi
           B_zero_heightmin = B
        end if

        B=B+0.01
 
      end do  !fin bucle sobre scaling factor B   

  !    close(20)    

      S = S + deltaS
       
  end do !fin bucle sobre angulos

  if (plot_Confront) then
  ! linea en blanco para contour plot
!  write(10,*) 
!  write(11,*) 
  write(12,*) 
  write(13,*) 
  end if

  T1 = T1 + deltaT

end do !fin bucle de temperaturas

  P1 = P1 + deltaP

end do !fin bucle de presiones

   if (plot_Confront) then  
   !  close(10)
   !  close(11)
     close(12)
     close(13)
   end if
   
   close(20)
   close(21)

  write(0,*) "Liu 2D=", liu_2Dmin
 ! write(0,'(A8,F9.6,A3,F9.5,A3,F5.1)')  " Alpha=", alpha_liu_2Dmin, " T=", T_liu_2Dmin, " P=", P_liu_2Dmin, " B=", B_liu_2Dmin
  write(0,'(A8,F9.6,A3,F9.5,A4,F6.1)')  " S=", alpha_liu_2Dmin, " T=", T_liu_2Dmin, "K P=", P_liu_2Dmin, "bar B=", B_liu_2Dmin

  write(0,*) "***********"

  write(0,*) "Liu 3D=", liu_3Dmin
!  write(0,'(A8,F9.6,A3,F9.5,A3,F5.1)') "  Alpha=", alpha_liu_3Dmin, " T=", T_liu_3Dmin, " P=", P_liu_3Dmin, " B=", B_liu_3Dmin
  write(0,'(A8,F9.6,A3,F9.5,A4,F6.1)') "  S=", alpha_liu_3Dmin, " T=", T_liu_3Dmin, "K P=", P_liu_3Dmin, "bar B=", B_liu_3Dmin

  write(0,*) "***********"

  write(0,*) "Zero Test=", zeromin
!  write(0,'(A8,F9.6,A3,F9.5,A3,F5.1)') "  Alpha=", alpha_zeromin, " T=", T_zeromin, " P=", P_zeromin, " B=", B_zeromin
  write(0,'(A8,F9.6,A3,F9.5,A4,F6.1)') "  S=", alpha_zeromin, " T=", T_zeromin, "K P=", P_zeromin, "bar B=", B_zeromin

  write(0,*) "***********"

  write(0,*) "Height Test=", heightmin
!  write(0,'(A8,F9.6,A3,F9.5,A3,F5.1)') "  Alpha=", alpha_heightmin, " T=", T_heightmin, " P=", P_heightmin, " B=", B_heightmin
  write(0,'(A8,F9.6,A3,F9.5,A4,F6.1)') "  S=", alpha_heightmin, " T=", T_heightmin, "K P=", P_heightmin, "bar B=", B_heightmin

  write(0,*) "***********"

  write(0,*) "Zero + Height Test=", zero_heightmin
!  write(0,'(A8,F9.6,A3,F9.5,A3,F5.1)') "  Alpha=", alpha_zero_heightmin, " T=", T_zero_heightmin, &
!         " P=", P_zero_heightmin, " B=", B_zero_heightmin
  write(0,'(A8,F9.6,A3,F9.5,A4,F6.1)') "  S=", alpha_zero_heightmin, " T=", T_zero_heightmin, &
         "K P=", P_zero_heightmin, "bar B=", B_zero_heightmin

!!!! CALCULO DE ERRORES

! punto mas cercano en P-T al punto interpolado
  do k=1,n_simulaciones
     err_distance = (T(k)-T_liu_3Dmin)*(T(k)-T_liu_3Dmin) + (P(k)-P_liu_3Dmin)*(P(k)-P_liu_3Dmin)
     if (k==1) then
        err_min_distance = err_distance
        err_k = 1
     else if ( err_distance < err_min_distance ) then
        err_min_distance = err_distance
        err_k = k
     end if
  end do
  
  call Calcula_Histograma(T(err_k),P(err_k),log_P,log_Z)	!histograma del punto de simulacion calculado via hist rew

  norma = 0
  do i=1,n_bin2D
        do w=1,n_bin2D
           if (memo_zero(i,w) == 0) then
              hist = exp(log_P(i,w)-log_Z)
              histogram_reweighted_grid(i,w) = hist
              norma = norma + hist
           end if
        end do
  end do
  norma = norma*tam_binE2D*tam_binDens2D
     
  do i=1,n_bin2D
       do w=1,n_bin2D
         if (memo_zero(i,w) == 0) then
           histogram_reweighted_grid(i,w) = histogram_reweighted_grid(i,w)/norma
         end if
       end do
  end do
       
  call Err_Calcula_Histograma_Rotado(alpha_liu_3Dmin,error_hist,error_X) 

!  print*, "ERROR_Hist= ", error_hist, "ERROR_X= ", error_X

!!!! PRINT HISTOGRAMS

!  call Print_Histogram(T_liu_2Dmin,P_liu_2Dmin,alpha_liu_2Dmin,B_liu_2Dmin)

  call Print_Histogram(T_liu_3Dmin,P_liu_3Dmin,alpha_liu_3Dmin,B_liu_3Dmin) 

!  call Print_Histogram(T_zeromin,P_zeromin,alpha_zeromin,B_zeromin)

!  call Print_Histogram(T_heightmin,P_heightmin,alpha_heightmin,B_heightmin)

!  call Print_Histogram(T_zero_heightmin,P_zero_heightmin,alpha_zero_heightmin,B_zero_heightmin)

!  call Print_Histogram_2D (T_liu_3Dmin, P_liu_3Dmin)
  
 ! call Print_Histogram_2D_Rotado (T_liu_3Dmin, P_liu_3Dmin,alpha_liu_3Dmin)

!  call Check(histogram_reweighted_mixed)

 ! open(unit=11,file="Ising2D.dat")
 ! open(unit=12,file="Ising3D.dat")
 ! xs = -3
 ! do while ( xs < 3 )
 !    write(11,*) xs, Ising2D(xs)
 !    write(12,*) xs, Ising3D(xs)
 !    xs = xs+0.001
 ! end do
 ! close(11)
 ! close(12)

  contains

  subroutine Calcula_Histograma(temp,pres,log_P,log_Z)  !Histogrma es exp(log_P(i,w)-log_Z)
                                                        !La cte es log_Z
                                                        !temp y pres de entrada: unids reales
     real*8,intent(in) :: temp, pres
     real*8,dimension(n_bin2D,n_bin2D),intent(out) :: log_P
     real*8,intent(out) :: log_Z
     integer :: i,j,w,flag, idx_i, idx_w
     real*8 :: pres1, beta_in, pres_in
    

     beta1 = 1000./(0.2293*temp)  ! temp in K , beta in mol/kJ
     pres1 = pres*100.0  ! convert P in bar to kPa
                         ! P*V in kPA * m^3/mol = kJ/mol
                         
 !    beta1 = beta1/20.0  !units close to 1
 !    pres1 = pres1/20.0
                          
     memo_zero(:,:) = 0
     do i=1,n_bin2D
       do w=1,n_bin2D
         hist=0
         do j=1,n_simulaciones
           hist = hist + datos_originales(j,i,w)
         end do
         if (.not. (hist==0)) then
            E = Emin + (i-0.5)*tam_binE2D    ! E in kJ/mol
            Dens = Densmin + (w-0.5)*tam_binDens2D   ! Dens in kg/m^3

            V = 0.018/Dens ! V in m^3/mol
            
  !          E = E/20.0
            ! DO NOT DIVIDE V. I divide by 20 beta, E y P*V. P has already been divided by 20, and so P*V 
            
            E = E*N
            V = V*N

            log_hist = log(hist)
            flag=0
            do j=1,n_simulaciones
                beta_in = 1000./(0.2293*T(j))
                pres_in = P(j)*100.0
   !             beta_in = beta_in/20.0
    !            pres_in = pres_in/20.0
                arg_exp(j) = beta1*(E+pres1*V)-beta_in*(E+pres_in*V)-ctes(j)+log_medidas(j)
                if ( ( dabs(arg_exp(j)) ) > 10) then 
                  flag=1 
                end if
            end do
            xx=0
            if (flag == 0) then
               do j=1,n_simulaciones
                  xx = xx + dexp(arg_exp(j))
               end do
               xx = log(xx)
            else
               max_arg = arg_exp(1)
               do j=2,n_simulaciones
                  if (max_arg < arg_exp(j)) then
                     max_arg = arg_exp(j)
                  end if
               end do

               do j=1,n_simulaciones
                  if ( max_arg - arg_exp(j) < 10 ) then
                   xx = xx + dexp(arg_exp(j) - max_arg)
                  end if
               end do
  
               if ( .not. xx == 0 ) then 
                  xx = log(xx)
               end if
               xx = xx + max_arg  
            end if
            log_P(i,w)=log_hist-xx
         else   ! En este caso el hist está vacío en este punto E-rho
            memo_zero(i,w) = 1
         end if
       end do
     end do

     i=1
     flag=0
     do while ( flag == 0 .and. i <= n_bin2D )
     
       do w=1,n_bin2D
         if (memo_zero(i,w)==0) then
          flag = 1
          idx_i = i
          idx_w = w
         end if
       end do
       
       i = i + 1 
     end do
     
     if(flag==0) then
         write(0,*) "Error. No hay histograma en T=",temp,"; P=",pres
         call exit(1)
     end if

     maxlog_P = log_P(idx_i,idx_w)
     do i=1,n_bin2D
       do w=1,n_bin2D
         if ( (memo_zero(i,w) == 0) .and. ( maxlog_P < log_P(i,w) ) ) then
            maxlog_P = log_P(i,w)
         end if
       end do
     end do 

     log_Z = 0
     do i=1,n_bin2D
       do w=1,n_bin2D
          if ( (memo_zero(i,w) == 0) .and. ( (maxlog_P-log_P(i,w)) < 10 ) ) then
            log_Z = log_Z + dexp(log_P(i,w)-maxlog_P)
          end if
       end do
     end do 
     if (.not. log_Z == 0) then 
        log_Z = log(log_Z) 
     end if
     log_Z = log_Z + maxlog_P
  end subroutine

  subroutine Rotate(a,b,x,y,alpha)
     real*8,intent(in) :: a,b
     real*8,intent(in) :: alpha
     real*8,intent(out) :: x,y

  !   x= a*cos(alpha)-b*sin(alpha)  !! Rotation matrix
  !   y= a*sin(alpha)+b*cos(alpha)
 
     y= a*cos(alpha)-b*sin(alpha)  !! 09/10/22. Integrate over Y instead over X. Rest of the code reads X,
                                   ! so I invert the definition.
     x= a*sin(alpha)+b*cos(alpha)  !! thus, M=sen(a)e+cos(a)rho ; s=tan(a)
   
  end subroutine

!  subroutine Mix(a,b,x,s)
!     real*8,intent(in) :: a,b
!     real*8,intent(in) :: s
!     real*8,intent(out) :: x

!     x = b+s*a
!  end subroutine

  subroutine Ising_Peak(Ising_Peak_2D, Ising_Peak_3D)
     real*8,intent(out) :: Ising_Peak_2D, Ising_Peak_3D
     logical :: flag
     real*8 :: xmin, xmax, x1, x2

     !!!!! Ising 2D
     xmin = 0.1
     xmax = 2
     x1 = xmin + (xmax - xmin)/3
     x2 = xmin + 2*(xmax - xmin)/3
     
     flag = .true.

     do while ( flag )

       if( Ising2D(x1) > Ising2D(x2) ) then
         xmax = x2
         x1 = xmin + (xmax - xmin)/3
         x2 = xmin + 2*(xmax - xmin)/3
       else if ( Ising2D(x2) > Ising2D(x1) ) then
         xmin = x1
         x1 = xmin + (xmax - xmin)/3
         x2 = xmin + 2*(xmax - xmin)/3
       else
          write(0,*) "Error en Ising_Peak"
          call exit(1)
       end if

      ! if (dabs(Ising2D(x1)-Ising2D(x2)) < 1E-06) then
       if (dabs(x1-x2) < 1E-06) then
	   flag = .false.
           if (Ising2D(x1) > Ising2D(x2) ) then
               Ising_Peak_2D = Ising2D(x1)
           else
               Ising_Peak_2D = Ising2D(x2)
           end if
       end if

     end do

     !!!!!!!!! Ising 3D
     xmin = 0.1
     xmax = 2
     x1 = xmin + (xmax - xmin)/3
     x2 = xmin + 2*(xmax - xmin)/3
     
     flag = .true.

     do while ( flag )

       if( Ising3D(x1) > Ising3D(x2) ) then
         xmax = x2
         x1 = xmin + (xmax - xmin)/3
         x2 = xmin + 2*(xmax - xmin)/3
       else if ( Ising3D(x2) > Ising3D(x1) ) then
         xmin = x1
         x1 = xmin + (xmax - xmin)/3
         x2 = xmin + 2*(xmax - xmin)/3
       else
          write(0,*) "Error en Ising_Peak"
          call exit(1)
       end if

       if (dabs(Ising3D(x1)-Ising3D(x2)) < 1E-06) then
	   flag = .false.
           if (Ising3D(x1) > Ising3D(x2) ) then
               Ising_Peak_3D = Ising3D(x1)
           else
               Ising_Peak_3D = Ising3D(x2)
           end if
       end if

     end do

 !   write(0,*) Ising_Peak_3D, Ising_Peak_2D

  end subroutine

  subroutine Calcula_Histograma_Rotado(alpha)    !resultados: histogram_reweighted_mixed()
 						 !            mixed_parameter()
						! DEBE ESTAR DEFINIDO histrogram_reweighted_grid(i,w)
    real*8,intent(in) :: alpha
    real*8 :: media_test
    
    Xmax = -1E12
!    Ymax = -1E12
    Xmin = 1E12
!    Ymin = 1E12
    
    do i=1,n_bin2D  !buscar maximo y minimo
       do w=1,n_bin2D
         if (histogram_reweighted_grid(i,w) > 0) then
            E = Emin + (i-0.5)*tam_binE2D
            Dens = Densmin + (w-0.5)*tam_binDens2D
            call Rotate(E,Dens,X,Y,alpha)
!            call Mix(E,Dens,X,alpha)
            if (Xmax < X) then 
                Xmax = X
            end if
            if (Xmin > X) then 
                Xmin = X
            end if
 !           if (Ymax < Y) then
!		Ymax = Y
!	    end if
!	    if (Ymin > Y) then
!	        Ymin = Y
!	    end if
         end if
       end do
    end do
    
    rangeX = Xmax - Xmin
    tam_binX = rangeX/n_bin_mixed
   ! rangeY = Ymax - Ymin
   ! tam_binY = rangeY/n_bin_mixed
   
    ! Bins vacios a cada lado. No se puede cambiar n_bin_mixed (alojamiento de vectores)
    rangeX = rangeX + 4*tam_binX
    Xmin = Xmin - 2*tam_binX
    Xmax = Xmax + 2*tam_binX
    tam_binX = rangeX/n_bin_mixed
    
    histogram_reweighted_mixed(:) = 0  !valores que no se tendran en cuenta
  !  histogram_reweighted_rotated(:,:) = 0

    do i=1,n_bin2D
      do w=1,n_bin2D
      !  if (memo_zero(i,w) == 0) then
        if (histogram_reweighted_grid(i,w) > 0) then
          E = Emin + (i-0.5)*tam_binE2D
          Dens = Densmin + (w-0.5)*tam_binDens2D
          call Rotate(E,Dens,X,Y,alpha)
 !         call Mix(E,Dens,X,alpha)
          indiceX = int((X-Xmin)/tam_binX)+1
!          indiceY = int((Y-Ymin)/tam_binY)+1

!Histograma en la direccion X (parametro de orden)
          if(indiceX >= 1 .and. indiceX <= n_bin_mixed) then
                 histogram_reweighted_mixed(indiceX) = histogram_reweighted_mixed(indiceX) + histogram_reweighted_grid(i,w)
          else
            if ( dabs(X - Xmin) < 1E-10 ) then
                 histogram_reweighted_mixed(1) = histogram_reweighted_mixed(1) + histogram_reweighted_grid(i,w)
            else
              print*,"Error: Puntos fuera del histograma", T1, P1, S, indiceX, n_bin_mixed
              print*,"(E, Emin, Dens, Densmin, X, Xmin):",E, Emin, Dens, Densmin,X,Xmin
              call exit(1)
            end if
          end if
!Histograma rotado en 2D
     !     if(indiceX >=1 .and. indiceX <= n_bin_mixed .and. indiceY >= 1 .and. indiceY <= n_bin_mixed) then
      !      	histogram_reweighted_rotated(indiceX,indiceY) = &
       !               histogram_reweighted_rotated(indiceX,indiceY) + histogram_reweighted_grid(i,w)
        !  else
         !     print*,"Error: Puntos fuera del histograma", T1, P1, S, indiceX, n_bin_mixed
          !    call exit(1)
       !   end if
        end if
      end do
    end do

!Normalizo el histograma en 2D
 !   norma = 0
 !   do i=1,n_bin_mixed
  !     do w=1,n_bin_mixed
   !      norma = norma + histogram_reweighted_rotated(i,w)
    !   end do
  !  end do
   ! histogram_reweighted_rotated(:,:) = histogram_reweighted_rotated(:,:)/(norma*tam_binX*tam_binY)

! Histograma sobre el parameto de orden
    do i=1,n_bin_mixed
       mixed_parameter(i) = Xmin+(i-0.5)*tam_binX 
    end do
    
!    delta_x = (mixed_parameter(1)+mixed_parameter(n_bin_mixed))/2.0   
!    mixed_parameter(:) = mixed_parameter(:) - delta_x     !traslado la distribución para centrarla en 0

    norma = 0   
    do i=1,n_bin_mixed   !hacer norma 1
      norma = norma + histogram_reweighted_mixed(i)
    end do
    histogram_reweighted_mixed(:) = histogram_reweighted_mixed(:)/(norma*tam_binX)

    media = 0
    do i=1,n_bin_mixed   !hacer media 0
       media = media + (mixed_parameter(i)*histogram_reweighted_mixed(i) * tam_binX)
    end do
    mixed_parameter(:) = mixed_parameter(:) - media   !traslado la distribución para que la media sea 0

    norma = 0
    var = 0
    do i=1,n_bin_mixed
      if (  histogram_reweighted_mixed(i) > 0 ) then
        norma = norma + histogram_reweighted_mixed(i)*tam_binX
        var = var + ( mixed_parameter(i) * mixed_parameter(i) * histogram_reweighted_mixed(i) * tam_binX)
      end if
    end do
    a1 = dsqrt(norma/var)   !hacer varianza y norma 1
    a2 = 1/(norma*a1)
    mixed_parameter(:) = mixed_parameter(:)*a1
    tam_binX = tam_binX*a1
    histogram_reweighted_mixed(:) = histogram_reweighted_mixed(:)*a2
    
   !! test propiedades    Poner en true para activar el test
    if ( .true. ) then 
      norma=0
      media_test=0  !calcula media_test para no sobreescribir media. media sera utilizado por call binder
      var=0
      do i=1,n_bin_mixed
        if (  histogram_reweighted_mixed(i) > 0 ) then
          norma = norma + histogram_reweighted_mixed(i)*tam_binX
          media_test = media_test + ( mixed_parameter(i) * histogram_reweighted_mixed(i) * tam_binX)
          var = var + ( mixed_parameter(i) * mixed_parameter(i) * histogram_reweighted_mixed(i) * tam_binX)
        end if
      end do
      var = var - media_test*media_test

      if (dabs(media_test) > 1E-10) then
         print*, "Error: La media no es cero. <X> =",media_test
         call exit(-1)
      end if
      if (dabs(norma - 1) > 1E-10) then
         print*, "Error: La norma no es uno. N =",norma
         call exit(-1)
      end if
      if (dabs(var- 1) > 1E-10) then
         print*, "Error: La varianza no es uno. V =",var
         call exit(-1)
      end if

    end if

  end subroutine
  
  subroutine Calcula_Histograma_Rotado_normX(alpha)   !resultados: histogram_reweighted_mixed()
 			         	  !            mixed_parameter()
						          ! DEBE ESTAR DEFINIDO histrogram_reweighted_grid(i,w)
    real*8,intent(in) :: alpha
    real*8 :: media_test, media_rx
    real*8,dimension(:),allocatable :: rx_tam_bin
    real*8:: b, alh, bt, pi
    
    Xmax = -1E12
!    Ymax = -1E12
    Xmin = 1E12
!    Ymin = 1E12
    
    do i=1,n_bin2D  !buscar maximo y minimo
       do w=1,n_bin2D
         if (histogram_reweighted_grid(i,w) > 0) then
            E = Emin + (i-0.5)*tam_binE2D
            Dens = Densmin + (w-0.5)*tam_binDens2D
            call Rotate(E,Dens,X,Y,alpha)
!            call Mix(E,Dens,X,alpha)

            if (Xmax < X) then 
                Xmax = X
            end if
            if (Xmin > X) then 
                Xmin = X
            end if
 !           if (Ymax < Y) then
!		Ymax = Y
!	    end if
!	    if (Ymin > Y) then
!	        Ymin = Y
!	    end if
         end if
       end do
    end do
    
    rangeX = Xmax - Xmin
    tam_binX = rangeX/n_bin_mixed
   ! rangeY = Ymax - Ymin
   ! tam_binY = rangeY/n_bin_mixed
   
    ! Bins vacios a cada lado. No se puede cambiar n_bin_mixed (alojamiento de vectores)
    rangeX = rangeX + 4*tam_binX
    Xmin = Xmin - 2*tam_binX
    Xmax = Xmax + 2*tam_binX
    tam_binX = rangeX/n_bin_mixed
    
    histogram_reweighted_mixed(:) = 0  !valores que no se tendran en cuenta
  !  histogram_reweighted_rotated(:,:) = 0

    do i=1,n_bin2D
      do w=1,n_bin2D
      !  if (memo_zero(i,w) == 0) then
        if (histogram_reweighted_grid(i,w) > 0) then
          E = Emin + (i-0.5)*tam_binE2D
          Dens = Densmin + (w-0.5)*tam_binDens2D
          call Rotate(E,Dens,X,Y,alpha)
 !         call Mix(E,Dens,X,alpha)
          indiceX = int((X-Xmin)/tam_binX)+1
!          indiceY = int((Y-Ymin)/tam_binY)+1

!Histograma en la direccion X (parametro de orden)
          if(indiceX >= 1 .and. indiceX <= n_bin_mixed) then
                 histogram_reweighted_mixed(indiceX) = histogram_reweighted_mixed(indiceX) + histogram_reweighted_grid(i,w)
          else
            if ( dabs(X - Xmin) < 1E-10 ) then
                 histogram_reweighted_mixed(1) = histogram_reweighted_mixed(1) + histogram_reweighted_grid(i,w)
            else
              print*,"Error: Puntos fuera del histograma", T1, P1, S, indiceX, n_bin_mixed
              print*,"(E, Emin, Dens, Densmin, X, Xmin):",E, Emin, Dens, Densmin,X,Xmin
              call exit(1)
            end if
          end if
!Histograma rotado en 2D
     !     if(indiceX >=1 .and. indiceX <= n_bin_mixed .and. indiceY >= 1 .and. indiceY <= n_bin_mixed) then
      !      	histogram_reweighted_rotated(indiceX,indiceY) = &
       !               histogram_reweighted_rotated(indiceX,indiceY) + histogram_reweighted_grid(i,w)
        !  else
         !     print*,"Error: Puntos fuera del histograma", T1, P1, S, indiceX, n_bin_mixed
          !    call exit(1)
       !   end if
        end if
      end do
    end do

!Normalizo el histograma en 2D
 !   norma = 0
 !   do i=1,n_bin_mixed
  !     do w=1,n_bin_mixed
   !      norma = norma + histogram_reweighted_rotated(i,w)
    !   end do
  !  end do
   ! histogram_reweighted_rotated(:,:) = histogram_reweighted_rotated(:,:)/(norma*tam_binX*tam_binY)

! Histograma sobre el parameto de orden
    do i=1,n_bin_mixed
       mixed_parameter(i) = Xmin+(i-0.5)*tam_binX 
    end do
    
!    delta_x = (mixed_parameter(1)+mixed_parameter(n_bin_mixed))/2.0   
!    mixed_parameter(:) = mixed_parameter(:) - delta_x     !traslado la distribución para centrarla en 0

    norma = 0   
    do i=1,n_bin_mixed   !hacer norma 1
      norma = norma + histogram_reweighted_mixed(i)
    end do
    histogram_reweighted_mixed(:) = histogram_reweighted_mixed(:)/(norma*tam_binX)

    media = 0
    do i=1,n_bin_mixed   !hacer media 0
       media = media + (mixed_parameter(i)*histogram_reweighted_mixed(i) * tam_binX)
    end do
    mixed_parameter(:) = mixed_parameter(:) - media   !traslado la distribución para que la media sea 0

    norma = 0
    var = 0
    do i=1,n_bin_mixed
      if (  histogram_reweighted_mixed(i) > 0 ) then
        norma = norma + histogram_reweighted_mixed(i)*tam_binX
        var = var + ( mixed_parameter(i) * mixed_parameter(i) * histogram_reweighted_mixed(i) * tam_binX)
      end if
    end do
    a1 = dsqrt(norma/var)   !hacer varianza y norma 1
    a2 = 1/(norma*a1)
    mixed_parameter(:) = mixed_parameter(:)*a1
    tam_binX = tam_binX*a1
    histogram_reweighted_mixed(:) = histogram_reweighted_mixed(:)*a2
    
   !! test propiedades    Poner en true para activar el test
    if ( .true. ) then 
      norma=0
      media_test=0  !calcula media_test para no sobreescribir media. media sera utilizado por call binder
      var=0
      do i=1,n_bin_mixed
        if (  histogram_reweighted_mixed(i) > 0 ) then
          norma = norma + histogram_reweighted_mixed(i)*tam_binX
          media_test = media_test + ( mixed_parameter(i) * histogram_reweighted_mixed(i) * tam_binX)
          var = var + ( mixed_parameter(i) * mixed_parameter(i) * histogram_reweighted_mixed(i) * tam_binX)
        end if
      end do
      var = var - media_test*media_test

      if (dabs(media_test) > 1E-10) then
         print*, "Error: La media no es cero. <X> =",media_test
         call exit(-1)
      end if
      if (dabs(norma - 1) > 1E-10) then
         print*, "Error: La norma no es uno. N =",norma
         call exit(-1)
      end if
      if (dabs(var- 1) > 1E-10) then
         print*, "Error: La varianza no es uno. V =",var
         call exit(-1)
      end if

    end if
    
    !	Renornalizo x para ajustar las colas...
    !  ver cuaderno 10/11/22
    
    b=1.53     ! L=16
    !b=1.61     ! L=12     
    alh=1.0    !alpha
    bt=1.0     !beta
    pi=3.14159
    mixed_parameter(:) = alh*b*b*sin(pi*mixed_parameter(:)/b)/(pi*pi) + bt*mixed_parameter(:)
    
    !nueva norma y varianza 1
    
    !ahora el tamaño del bin cambia a lo largo de x
    
    allocate(rx_tam_bin(n_bin_mixed))
    rx_tam_bin(1) = mixed_parameter(2)-mixed_parameter(1)
    do i=2,n_bin_mixed-1   
      rx_tam_bin(i) = (mixed_parameter(i+1)-mixed_parameter(i-1))/2.0
    end do
    rx_tam_bin(n_bin_mixed) = mixed_parameter(n_bin_mixed)-mixed_parameter(n_bin_mixed-1)
    
    norma = 0   
    do i=1,n_bin_mixed   !hacer norma 1
      norma = norma + histogram_reweighted_mixed(i)*rx_tam_bin(i)
    end do
    
    histogram_reweighted_mixed(:) = histogram_reweighted_mixed(:)/norma

    media_rx = 0
    do i=1,n_bin_mixed   !hacer media 0
       media_rx = media_rx + mixed_parameter(i)*histogram_reweighted_mixed(i)*rx_tam_bin(i)
    end do
    
    mixed_parameter(:) = mixed_parameter(:) - media_rx   !traslado la distribución para que la media sea 0

    norma = 0
    var = 0
    do i=1,n_bin_mixed
      if (  histogram_reweighted_mixed(i) > 0 ) then
        norma = norma + histogram_reweighted_mixed(i)*rx_tam_bin(i)
        var = var + ( mixed_parameter(i) * mixed_parameter(i) * histogram_reweighted_mixed(i) * rx_tam_bin(i))
      end if
    end do
    a1 = dsqrt(norma/var)   !hacer varianza y norma 1
    a2 = 1/(norma*a1)
    mixed_parameter(:) = mixed_parameter(:)*a1
    rx_tam_bin(:) = rx_tam_bin(:)*a1
    histogram_reweighted_mixed(:) = histogram_reweighted_mixed(:)*a2
    
    !! test propiedades    Poner en true para activar el test
    if ( .true. ) then 
      norma=0
      media_test=0  !calcula media_test para no sobreescribir media. media sera utilizado por call binder
      var=0
      do i=1,n_bin_mixed
        if (  histogram_reweighted_mixed(i) > 0 ) then
          norma = norma + histogram_reweighted_mixed(i)*rx_tam_bin(i)
          media_test = media_test + ( mixed_parameter(i) * histogram_reweighted_mixed(i) * rx_tam_bin(i))
          var = var + ( mixed_parameter(i) * mixed_parameter(i) * histogram_reweighted_mixed(i) * rx_tam_bin(i))
        end if
      end do
      var = var - media_test*media_test

      if (dabs(media_test) > 1E-10) then
         print*, "Error (RX): La media no es cero. <X> =",media_test
         call exit(-1)
      end if
      if (dabs(norma - 1) > 1E-10) then
         print*, "Error (RX): La norma no es uno. N =",norma
         call exit(-1)
      end if
      if (dabs(var- 1) > 1E-10) then
         print*, "Error (RX): La varianza no es uno. V =",var
         call exit(-1)
      end if

    end if

  end subroutine


  subroutine Err_Calcula_Histograma_Rotado(alpha,error_hist,error_X)    !calcula err_hist_reweighted() , err_hist_original()
 				     !        mixed_parameter_hist_reweighted(),  mixed_parameter_original()
 				     !        devuelve error_hist y error_X (ejes Y, X respectivamente)
    real*8,intent(in) :: alpha
    real*8, dimension(:),allocatable,intent(out) :: error_hist, error_X
    real*8, dimension(:), allocatable :: err_hist_original, err_hist_reweighted, hist_input
    real*8, dimension(:), allocatable :: mixed_parameter_original, mixed_parameter_hist_reweighted, mixed_parameter_input
    real*8 :: tam_binX_original, tam_binX_hist_reweighted, dist, min_dist
    integer :: min_idx
    
    allocate(err_hist_original(n_bin),err_hist_reweighted(n_bin_mixed))
    allocate(mixed_parameter_original(n_bin),mixed_parameter_hist_reweighted(n_bin_mixed))
    allocate(error_hist(n_bin_mixed), error_X(n_bin_mixed))
    
    allocate(mixed_parameter_input(n_bin_mixed),hist_input(n_bin_mixed))

!!!!!!!!!!!!!!!!!!!!
!!! HISTOGRAMA REPESADO

    Xmax = -1E12
    Xmin = 1E12

    do i=1,n_bin2D  !buscar maximo y minimo
       do w=1,n_bin2D
        if (histogram_reweighted_grid(i,w) > 0) then
          E = Emin + (i-0.5)*tam_binE2D
          Dens = Densmin + (w-0.5)*tam_binDens2D
          call Rotate(E,Dens,X,Y,alpha)
 !         call Mix(E,Dens,X,alpha)
          if (Xmax < X) then 
              Xmax = X
          end if
          if (Xmin > X) then 
              Xmin = X
          end if
        end if
       end do
    end do

    rangeX = Xmax - Xmin
    tam_binX = rangeX/n_bin_mixed
    ! Bins vacios a cada lado. No se puede cambiar n_bin_mixed (alojamiento de vectores)
    rangeX = rangeX + 4*tam_binX
    Xmin = Xmin - 2*tam_binX
    Xmax = Xmax + 2*tam_binX
    tam_binX = rangeX/n_bin_mixed

    err_hist_reweighted(:) = 0  !valores que no se tendran en cuenta

    do i=1,n_bin2D
      E = Emin + (i-0.5)*tam_binE2D
      do w=1,n_bin2D
        Dens = Densmin + (w-0.5)*tam_binDens2D
        call Rotate(E,Dens,X,Y,alpha)
 !       call Mix(E,Dens,X,alpha)
        indiceX = int((X-Xmin)/tam_binX)+1
 
        if (memo_zero(i,w) == 0) then
          if(indiceX >= 1 .and. indiceX <= n_bin_mixed) then
              err_hist_reweighted(indiceX) = err_hist_reweighted(indiceX) + histogram_reweighted_grid(i,w)
          end if
        end if
        
      end do
    end do

    do i=1,n_bin_mixed
       mixed_parameter_hist_reweighted(i) = Xmin+(i-0.5)*tam_binX 
    end do
    
! histograma repesado: media 0, norma y varianza 1

    norma = 0
    do i=1,n_bin_mixed
       norma = norma + err_hist_reweighted(i) 
    end do
    err_hist_reweighted(:) = err_hist_reweighted(:)/(norma * tam_binX)

    media = 0
    do i=1,n_bin_mixed
       media = media + (mixed_parameter_hist_reweighted(i)*err_hist_reweighted(i) * tam_binX)
    end do
    mixed_parameter_hist_reweighted(:) = mixed_parameter_hist_reweighted(:) - media   
    
    norma = 0
    var = 0
    do i=1,n_bin_mixed
      if (  err_hist_reweighted(i) > 0 ) then
        norma = norma + err_hist_reweighted(i)*tam_binX
        var = var + ( mixed_parameter_hist_reweighted(i) * mixed_parameter_hist_reweighted(i) * err_hist_reweighted(i) * tam_binX)
      end if
    end do
    a1 = dsqrt(norma/var)   !hacer varianza y norma 1
    a2 = 1/(norma*a1)
    mixed_parameter_hist_reweighted(:) = mixed_parameter_hist_reweighted(:)*a1
    tam_binX_hist_reweighted = tam_binX*a1
    err_hist_reweighted(:) = err_hist_reweighted(:)*a2
    
!    print*,"HISTOGRAMA REPESADO"
    
 !   do i=1,n_bin_mixed
 !      print*, mixed_parameter_hist_reweighted(i), err_hist_reweighted(i)
 !   end do
   
 !   print*, " "
    
    !! test propiedades    Poner en true para activar el test
    if ( .true. ) then 
      norma=0
      media=0
      var=0
      do i=1,n_bin_mixed
        if (  err_hist_reweighted(i) > 0 ) then
          norma = norma + err_hist_reweighted(i)*tam_binX_hist_reweighted
          media = media + ( mixed_parameter_hist_reweighted(i) * err_hist_reweighted(i) * tam_binX_hist_reweighted)
          var = var + ( mixed_parameter_hist_reweighted(i) * mixed_parameter_hist_reweighted(i) * err_hist_reweighted(i))
        end if
      end do
      var = var * tam_binX_hist_reweighted
      var = var - media*media

      if (dabs(media) > 1E-10) then
         print*, "Error: (A) La media no es cero. <X> =",media
         call exit(-1)
      end if
      if (dabs(norma - 1) > 1E-10) then
         print*, "Error: (A) La norma no es uno. N =",norma
         call exit(-1)
      end if
      if (dabs(var- 1) > 1E-10) then
         print*, "Error: (A) La varianza no es uno. V =",var
         call exit(-1)
      end if

    end if

!!!!!!!!!!!!!!!!!!!!    
! HISTOGRAMA INPUT     (integracion de hist(ene,rho))
 if ( .false. ) then
    Xmax = -1E12
    Xmin = 1E12

    do i=1,n_bin2D  !buscar maximo y minimo
       do w=1,n_bin2D
        if (datos_originales(err_k,i,w) > 0) then
          E = Emin + (i-0.5)*tam_binE2D
          Dens = Densmin + (w-0.5)*tam_binDens2D
          call Rotate(E,Dens,X,Y,alpha)
 !         call Mix(E,Dens,X,alpha)
          if (Xmax < X) then 
              Xmax = X
          end if
          if (Xmin > X) then 
              Xmin = X
          end if
        end if
       end do
    end do

    rangeX = Xmax - Xmin
    tam_binX = rangeX/n_bin_mixed
    
    ! Bins vacios a cada lado. No se puede cambiar n_bin_mixed (alojamiento de vectores)
    rangeX = rangeX + 4*tam_binX
    Xmin = Xmin - 2*tam_binX
    Xmax = Xmax + 2*tam_binX
    tam_binX = rangeX/n_bin_mixed

    hist_input(:) = 0  !valores que no se tendran en cuenta

    do i=1,n_bin2D
      E = Emin + (i-0.5)*tam_binE2D
      do w=1,n_bin2D
        Dens = Densmin + (w-0.5)*tam_binDens2D
        call Rotate(E,Dens,X,Y,alpha)
   !     call Mix(E,Dens,X,alpha)
        indiceX = int((X-Xmin)/tam_binX)+1
 
        if ( indiceX >= 1 .and. indiceX <= n_bin_mixed ) then
          hist_input(indiceX) = hist_input(indiceX) + datos_originales(err_k,i,w)
        end if
        
      end do
    end do

    do i=1,n_bin_mixed
       mixed_parameter_input(i) = Xmin+(i-0.5)*tam_binX 
    end do

! histograma input: media 0, norma y varianza 1

    norma = 0
    do i=1,n_bin_mixed
       norma = norma + hist_input(i) 
    end do
    hist_input(:) = hist_input(:)/(norma * tam_binX)

    media = 0
    do i=1,n_bin_mixed
       media = media + (mixed_parameter_input(i)*hist_input(i) * tam_binX)
    end do
    mixed_parameter_input(:) = mixed_parameter_input(:) - media   
    
    norma = 0
    var = 0
    do i=1,n_bin_mixed
      if (  hist_input(i) > 0 ) then
        norma = norma + hist_input(i)*tam_binX
        var = var + ( mixed_parameter_input(i) * mixed_parameter_input(i) * hist_input(i) * tam_binX)
      end if
    end do
    a1 = dsqrt(norma/var)   !hacer varianza y norma 1
    a2 = 1/(norma*a1)
    mixed_parameter_input(:) = mixed_parameter_input(:)*a1
    tam_binX_hist_reweighted = tam_binX*a1
    hist_input(:) = hist_input(:)*a2
    
 !   print*,"HISTOGRAMA INPUT"
    
 !   do i=1,n_bin_mixed
 !      print*, mixed_parameter_input(i), hist_input(i)
 !   end do
    
 !   print*, " "
    
    !! test propiedades    Poner en true para activar el test
    if ( .true. ) then 
      norma=0
      media=0
      var=0
      do i=1,n_bin_mixed
        if (  hist_input(i) > 0 ) then
          norma = norma + hist_input(i)*tam_binX_hist_reweighted
          media = media + ( mixed_parameter_input(i) * hist_input(i) * tam_binX_hist_reweighted)
          var = var + ( mixed_parameter_input(i) * mixed_parameter_input(i) * hist_input(i))
        end if
      end do
      var = var * tam_binX_hist_reweighted
      var = var - media*media

      if (dabs(media) > 1E-10) then
         print*, "Error: (C) La media no es cero. <X> =",media
         call exit(-1)
      end if
      if (dabs(norma - 1) > 1E-10) then
         print*, "Error: (C) La norma no es uno. N =",norma
         call exit(-1)
      end if
      if (dabs(var- 1) > 1E-10) then
         print*, "Error: (C) La varianza no es uno. V =",var
         call exit(-1)
      end if

    end if
  end if
!!!!!!!!!!!!!!!!!!!!
! HISTOGRAMA ORIGINAL  (calculo directo desde input file, no integrando hist(ene,rho))

    Xmax = -1E12
    Xmin = 1E12
    
    if ( L<10 ) then
       if ( Phat(err_k) >= 0 ) then
          write(filename,'("data_L",I1,"_T",F9.7,"_P",F6.4)') L,That(err_k),Phat(err_k)
       else
          write(filename,'("data_L",I1,"_T",F9.7,"_P",F6.3)') L,That(err_k),Phat(err_k)
       end if
    else
       if ( Phat(err_k) >= 0 ) then 
          write(filename,'("data_L",I2,"_T",F9.7,"_P",F6.4)') L,That(err_k),Phat(err_k)
       else
          write(filename,'("data_L",I2,"_T",F9.7,"_P",F6.3)') L,That(err_k),Phat(err_k)
       end if
    end if
    open(unit=10,file=filename,status="old")
    do j=1,medidas(err_k)
       read(10,*) E, Dens
       Dens = rhos*Dens+rhoi  ! Dens in kg/m^3
       E = Hs*E+Hi*That(i)- 1.8*P(i)/Dens ! Energy in kJ/mol
       call Rotate(E,Dens,X,Y,alpha)
  !     call Mix(E,Dens,X,alpha)
       
       if ( X > Xmax ) then
          Xmax = X
       end if
       
       if ( X < Xmin ) then
          Xmin = X
       end if
    end do
    
    close(10)
    
    rangeX = Xmax - Xmin
    tam_binX_original = rangeX/(n_bin)

    ! Bins vacios a cada lado. No se puede cambiar n_bin (alojamiento de vectores)
    rangeX = rangeX + 4*tam_binX_original
    Xmin = Xmin - 2*tam_binX_original
    Xmax = Xmax + 2*tam_binX_original
    tam_binX_original = rangeX/n_bin
    
    err_hist_original(:) = 0  !valores que no se tendran en cuenta
    
    open(unit=10,file=filename,status="old")
    norma = medidas(err_k)
    do j=1,medidas(err_k)
    
       read(10,*) E, Dens
       Dens = rhos*Dens+rhoi  ! Dens in kg/m^3
       E = Hs*E+Hi*That(i)- 1.8*P(i)/Dens ! Energy in kJ/mol
       
       call Rotate(E,Dens,X,Y,alpha)
   !   call Mix(E,Dens,X,alpha)
       
       indiceX = int((X-Xmin)/tam_binX_original)+1
       if(indiceX >= 1 .and. indiceX <= n_bin) then
            err_hist_original(indiceX) = err_hist_original(indiceX) + 1
       else
            print*, "Warning (B): data out of histogram"
            norma = norma - 1
       end if
       
    end do
    
    close(10)

    do i=1,n_bin
       mixed_parameter_original(i) = Xmin+(i-0.5)*tam_binX_original 
    end do

! histograma original: norma y varianza 1

    err_hist_original(:) = err_hist_original(:)/(norma * tam_binX_original)
    
!    do i=1,n_bin
!       print*,mixed_parameter_original(i),  err_hist_original(i)
!    end do
    
    media = 0
    do i=1,n_bin
       media = media + (mixed_parameter_original(i)*err_hist_original(i) * tam_binX_original)
    end do
    mixed_parameter_original(:) = mixed_parameter_original(:) - media   !traslado la distribución para que la media sea 0

    norma = 0
    var = 0
    do i=1,n_bin
      if (  err_hist_original(i) > 0 ) then
        norma = norma + err_hist_original(i)*tam_binX_original
        var = var + ( mixed_parameter_original(i) * mixed_parameter_original(i) * err_hist_original(i) * tam_binX_original)
      end if
    end do
    a1 = dsqrt(norma/var)   !hacer varianza y norma 1
    a2 = 1/(norma*a1)
    mixed_parameter_original(:) = mixed_parameter_original(:)*a1
    tam_binX_original = tam_binX_original*a1
    err_hist_original(:) = err_hist_original(:)*a2
    
  !  print*,"HISTOGRAMA ORIGINAL"
    
  !  do i=1,n_bin
  !     print*, mixed_parameter_original(i), err_hist_original(i)
  !  end do
    
  !  print*, " "
    
    !! test propiedades    Poner en true para activar el test
    if ( .true. ) then 
      norma=0
      media=0
      var=0
      do i=1,n_bin
        if (  err_hist_original(i) > 0 ) then
          norma = norma + err_hist_original(i)*tam_binX_original
          media = media + ( mixed_parameter_original(i) * err_hist_original(i) * tam_binX_original)
          var = var + ( mixed_parameter_original(i) * mixed_parameter_original(i) * err_hist_original(i) )
        end if
      end do
      var = var * tam_binX_original
      var = var - media*media

      if (dabs(media) > 1E-10) then
         print*, "Error: (B) La media no es cero. <X> =",media
         call exit(-1)
      end if
      if (dabs(norma - 1) > 1E-10) then
         print*, "Error: (B) La norma no es uno. N =",norma
         call exit(-1)
      end if
      if (dabs(var- 1) > 1E-10) then
         print*, "Error: (B) La varianza no es uno. V =",var
         call exit(-1)
      end if

    end if
    
!!!! Cálculo de errores

     error_hist(:)=0
     error_X(:)=0
     do i=1,n_bin_mixed
     
       min_dist=1E12
       do j=1,n_bin
         dist = dabs(mixed_parameter_original(j)-mixed_parameter_hist_reweighted(i))
         if ( dist < min_dist ) then
            min_dist = dist
            min_idx = j
         end if
       end do
     
       error_hist(i) = dabs(err_hist_original(min_idx)-err_hist_reweighted(i))
       error_X(i) = min_dist

     end do

  end subroutine


  subroutine ZeroTest ( zero )
     real*8, intent(out) :: zero
     real*8 :: m, hist, hist_min
     integer :: i,index_hist_min
     real*8 :: cero
     
     cero = 0.0 !(local) minimum of Ising 3D
     hist_min = 1E12
     
     do i=1,n_bin_mixed
        m = mixed_parameter(i)
        hist = histogram_reweighted_mixed(i)
        
        if ( m > -1 .and. m < 1) then !not count tails of distribution
           if (hist < hist_min) then
               hist_min = hist
               index_hist_min = i
           end if
        end if

     end do
     
     if ( index_hist_min >= 1 .and. index_hist_min <= n_bin_mixed) then
        zero = (histogram_reweighted_mixed(index_hist_min)-Ising3D(cero))*(histogram_reweighted_mixed(index_hist_min)-Ising3D(cero))
        zero = zero*(6.0/0.8)*(6.0/0.8)   
        !multiplico por 6.0/0.8 para compensar el efecto visual de la figura.
        !eje X: -3 , 3   eje Y: 0, 0.8
        zero = zero + mixed_parameter(index_hist_min)*mixed_parameter(index_hist_min)
     else
        zero = 1E3
     end if
     
  end subroutine

  subroutine HeightTest ( height )
     real*8, intent(out) :: height
     real*8 :: m, hist, hist_max_1, hist_max_2
     integer :: i,index_hist_max_1, index_hist_max_2
     real*8 :: hist1, hist2, m1, m2, x1, x2
     
     hist_max_1 = -10
     hist_max_2 = -10

     index_hist_max_1 = -10
     index_hist_max_2 = -10

     do i=1,n_bin_mixed
        m = mixed_parameter(i)
        hist = histogram_reweighted_mixed(i)
        
        if (m < 0) then
           if (hist > hist_max_1) then
               hist_max_1 = hist
               index_hist_max_1 = i
           end if
        else
           if (hist > hist_max_2) then
               hist_max_2 = hist
               index_hist_max_2 = i
           end if
        end if

     end do

     x1=-1.14  !possitions of the maxima of Ising 3D
     x2=1.14
     if ( (index_hist_max_1 >=1 .and. index_hist_max_1 <= n_bin_mixed)               &
          .and. (index_hist_max_2 >=1 .and. index_hist_max_2 <= n_bin_mixed) ) then
          
          hist1=histogram_reweighted_mixed(index_hist_max_1)
          m1=mixed_parameter(index_hist_max_1)
          hist2=histogram_reweighted_mixed(index_hist_max_2)
          m2=mixed_parameter(index_hist_max_2)

                       !multiplico por 6.0/0.8 para compensar el efecto visual de la figura.
                       !eje X: -3 , 3   eje Y: 0, 0.8
          height = (hist1-Ising3D(x1))*(hist1-Ising3D(x1))*(6.0/0.8)*(6.0/0.8)   
          height = height + (hist2-Ising3D(x2))*(hist2-Ising3D(x2))*(6.0/0.8)*(6.0/0.8)
          height = height + (m1-x1)*(m1-x1)
          height = height + (m2-x2)*(m2-x2)
     else
          height = 1E3
     end if
  end subroutine

  subroutine ConfrontIsing ( kl_2D, kl_3D, liu_2D, liu_3D)
     real*8, intent(out) :: kl_2D, kl_3D, liu_2D, liu_3D
     integer :: i, n_liu
     real*8 :: m, hist
     real*8,dimension(n_bin_mixed) :: Ising2D_KL

  !   call Calculate_Ising_KL(Ising2D_KL)
    
     kl_2D = 0
     kl_3D = 0
     liu_2D = 0
     liu_3D = 0
     n_liu = 0
     
     do i=1,n_bin_mixed
         m = mixed_parameter(i)
    !     if ( dabs(m) <= 2.5 ) then
           hist = histogram_reweighted_mixed(i)
           !print*, Ising2D(m), hist, log(Ising2D(m)/hist)
           if (.not. hist == 0) then
             if ( Ising2D(m) > 0 ) then   !avoid log(0)
               kl_2D = kl_2D + (dlog(Ising2D(m))-dlog(hist))*Ising2D(m)
             end if
             if ( Ising3D(m) > 0 ) then
               kl_3D = kl_3D + (dlog(Ising3D(m))-dlog(hist))*Ising3D(m)
             end if
               liu_2D = liu_2D + sqrt(hist)*dabs(hist-Ising2D(m)) 
               liu_3D = liu_3D + sqrt(hist)*dabs(hist-Ising3D(m))
               n_liu = n_liu + 1
           end if
    !     end if
     end do

!     kl_2D = 0
!     do i=1,n_bin_mixed
!         hist = histogram_reweighted_mixed(i)
!         if (.not. hist == 0) then
!           if ( Ising2D_KL(i) > 0 ) then   !avoid log(0) and division by 0
!             kl_2D = kl_2D + dlog(Ising2D_KL(i)/hist)*Ising2D_KL(i)
!           end if
!         end if
!     end do

 !    kl_2D = dabs(kl_2D)
 !    kl_3D = dabs(kl_3D)

 !    if( kl_2D < 0) then
 !	if (dabs(kl_2D) + kl_2D > 1E-04) then
 !          print*, "ALERTA!!!"
 !       end if
 !    end if

     liu_2D = liu_2D/( (peak2D - Ising2D( dfloat(0) ) )*n_liu ) !n_bin_mixed ) 
     liu_3D = liu_3D/( (peak3D - Ising3D( dfloat(0) ) )*n_liu ) !n_bin_mixed ) 

  end subroutine

  subroutine Calculate_Ising_KL(Ising2D_KL)
     integer :: i
     real*8 :: a1,a2, norma,media,var, m
     real*8,dimension(n_bin_mixed),intent(out) :: Ising2D_KL

     do i=1,n_bin_mixed
       m = mixed_parameter(i)
       Ising2D_KL(i) = Ising2D(m)
 !      Ising3D_KL(i) = Ising3D(m)
     end do

     norma = 0
     media = 0
     var = 0
     do i=1,n_bin_mixed
       norma = norma + Ising2D_KL(i)*tam_binX
       media = media + mixed_parameter(i)*Ising2D_KL(i)*tam_binX  !0 si la distrib es simetrica
       var = var + mixed_parameter(i)*mixed_parameter(i)*Ising2D_KL(i)*tam_binX
    end do
    var = var - media*media
    a1 = sqrt(norma/var)
    a2 = 1/(norma*a1)

    mixed_parameter(:) = mixed_parameter(:)*a1
    histogram_reweighted_mixed(:) = histogram_reweighted_mixed(:)*a2

  end subroutine

  subroutine Transalate_Histogram (delta)
      real*8 :: delta

      mixed_parameter(:) = mixed_parameter(:) + delta
  end subroutine

  function Ising2D(z)
      real*8 :: a,b,c,d
      real*8,intent(in) :: z
      real*8 :: Ising2D

      a=0.058137
      b =-2.68949
      c=-0.11235
      d=2.884162

      Ising2D=dexp(-(a*(z**16) + b*z*z +c*dabs(z) +d))
      return

! Liu Panagiotopoulos JCP 2010

  end function
	
  function Ising3D(z)
	real*8:: a,c,M0,A0
        real*8,intent(in) :: z
        real*8 :: Ising3D   
   
	a=0.196
 	c=0.807
 	M0=1.141059404
	A0=0.436591191

 ! i valori di a e c sono riportati nell'articolo di Tsypin, PRE 62, 2000. Non sono asintotici però dai grafici sembra 
 ! che abbiano raggiunto un plateau. Per quanto rigurada M0 dipende dalla size, e scala linearmente con L^(-0.8)
 ! secondo la legge lineare 2.19847*x+0.093522 che ho trovato io fittando i dati dell'articolo 
 ! Per rendere la distribuzione con norma e varianza unitaria per ogni valore di M0 avrò un fattore A0 davanti
 ! all'esponenziale. Con la scelta dei parametri fatta soddisfo la condizione di norma e varianza unitari
	
	Ising3D=A0*dexp(-(z*z/(M0*M0)-1)*(z*z/(M0*M0)-1)*(a*z*z/(M0*M0)+c))
	
	return
  end function 

  subroutine Print_Histogram (T1, P1, alpha1,B1)  !T in K and P in bar
    real*8, intent(in) :: T1,P1,alpha1,B1
    real*8 :: s

    call Calcula_Histograma(T1,P1,log_P,log_Z)

    norma = 0
    do i=1,n_bin2D
       do w=1,n_bin2D
          if (memo_zero(i,w) == 0) then
             hist = exp(log_P(i,w)-log_Z)
             histogram_reweighted_grid(i,w) = hist
             norma = norma + hist
          end if
       end do
    end do
    norma = norma*tam_binE2D*tam_binDens2D
 
    do i=1,n_bin2D
      do w=1,n_bin2D
        if (memo_zero(i,w) == 0) then
          histogram_reweighted_grid(i,w) = histogram_reweighted_grid(i,w)/norma
        end if
      end do
    end do

    print*, "alpha =", alpha1
  !  call Calcula_Histograma_Rotado(alpha1)
    call Calcula_Histograma_Rotado_normX(alpha) 
    
 !   s = tan(-alpha1)
    s=alpha1

    if ( L<10 ) then

      if (P1 < -1000.0) then
        write(filename,'("Histogram_Reweighted_OP_L",I1,"_T",F9.5,"_P",F7.1,"_A",F8.6,"_B",F5.3)') L,T1,P1,s,B1
      else if (P1 < 0 .or. P1 > 1000.0) then
        write(filename,'("Histogram_Reweighted_OP_L",I1,"_T",F9.5,"_P",F6.1,"_A",F8.6,"_B",F5.3)') L,T1,P1,s,B1
      else  
        write(filename,'("Histogram_Reweighted_OP_L",I1,"_T",F9.5,"_P",F5.1,"_A",F8.6"_B",F5.3)') L,T1,P1,s,B1
      end if
 
    else
 
      if (P1 < -1000.0) then
        write(filename,'("Histogram_Reweighted_OP_L",I2,"_T",F9.5,"_P",F7.1,"_A",F8.6,"_B",F5.3)') L,T1,P1,s,B1
      else if (P1 < 0 .or. P1 > 1000.0) then
        write(filename,'("Histogram_Reweighted_OP_L",I2,"_T",F9.5,"_P",F6.1,"_A",F8.6,"_B",F5.3)') L,T1,P1,s,B1
      else  
        write(filename,'("Histogram_Reweighted_OP_L",I2,"_T",F9.5,"_P",F5.1,"_A",F8.6"_B",F5.3)') L,T1,P1,s,B1
      end if
 
    end if
    
    open(unit=14,file=filename)
    write(14,*) "#Data simulations (Temperature internal units): ", That(:)
    do i=1,n_bin_mixed
  !    if (.not. histogram_reweighted_mixed(i) == 0) then
        write(14,*) mixed_parameter(i)*B1, histogram_reweighted_mixed(i), error_X(i), error_hist(i)
  !    end if
    end do
    close(14)

    print*, "Mixed Param (min) = ", mixed_parameter(1), " ; Mixed_param (max) = ", mixed_parameter(n_bin_mixed)
    print*, "Normalization factor=",a1

  end subroutine

  subroutine Print_Histogram_2D_Rotado (T1, P1, alpha1)
    real*8, intent(in) :: T1,P1,alpha1

    call Calcula_Histograma(T1,P1,log_P,log_Z)

    norma = 0
    do i=1,n_bin
       do w=1,n_bin
          if (memo_zero(i,w) == 0) then
             hist = exp(log_P(i,w)-log_Z)
             histogram_reweighted_grid(i,w) = hist
             norma = norma + hist
          end if
       end do
    end do
    norma = norma*tam_binE*tam_binDens
 
    do i=1,n_bin
      do w=1,n_bin
        if (memo_zero(i,w) == 0) then
          histogram_reweighted_grid(i,w) = histogram_reweighted_grid(i,w)/norma
        end if
      end do
    end do

    s = alpha1

    if (L < 10 ) then
      if (s < 0) then
        write(filename,'("Histogram_Reweighted_Rotado_L",I1,"_T",F9.5,"_P",F5.1,"_S",F9.6)') L,T1,P1,s
      else if (s < 10) then
        write(filename,'("Histogram_Reweighted_Rotado_L",I1,"_T",F9.5,"_P",F5.1,"_S",F8.6)') L,T1,P1,s
      else 
        write(filename,'("Histogram_Reweighted_Rotado_L",I1,"_T",F9.5,"_P",F5.1,"_S",F9.6)') L,T1,P1,s
      end if
    else
      if (s < 0) then
        write(filename,'("Histogram_Reweighted_Rotado_L",I2,"_T",F9.5,"_P",F5.1,"_S",F9.6)') L,T1,P1,s
      else if (s < 10) then
        write(filename,'("Histogram_Reweighted_Rotado_L",I2,"_T",F9.5,"_P",F5.1,"_S",F8.6)') L,T1,P1,s
      else 
        write(filename,'("Histogram_Reweighted_Rotado_L",I2,"_T",F9.5,"_P",F5.1,"_S",F9.6)') L,T1,P1,s
      end if    
    end if
    
    open(unit=14,file=filename)
    write(14,*) "#Data simulations (Temperature): ", T(:)

    do i=1,n_bin
      E = Emin + (i-0.5)*tam_binE
      do w=1,n_bin
        Dens = Densmin + (w-0.5)*tam_binDens
        if (memo_zero(i,w) == 0) then
          call Rotate(E,Dens,X,Y,alpha)
          write(14,*) X, Y, histogram_reweighted_grid(i,w)
        else
    !      write(14,*) E, Dens, 0
        end if
      end do
   !   write(14,*) " "
    end do
    close(14)

  end subroutine


  subroutine Print_Histogram_2D (T1, P1)
    real*8, intent(in) :: T1,P1
    real*8 :: min_free

    call Calcula_Histograma(T1,P1,log_P,log_Z)

    norma = 0
    do i=1,n_bin2D
       do w=1,n_bin2D
          if (memo_zero(i,w) == 0) then
             hist = exp(log_P(i,w)-log_Z)
             histogram_reweighted_grid(i,w) = hist
             norma = norma + hist
          end if
       end do
    end do
    norma = norma*tam_binE2D*tam_binDens2D
 
    do i=1,n_bin2D
      do w=1,n_bin2D
        if (memo_zero(i,w) == 0) then
          histogram_reweighted_grid(i,w) = histogram_reweighted_grid(i,w)/norma
        end if
      end do
    end do

 !   if (L < 10 ) then
 !      write(filename,'("Histogram_Reweighted_ene_rho_L",I1,"_T",F9.5,"_P",F5.1)') L,T1,P1
 !   else
 !      write(filename,'("Histogram_Reweighted_ene_rho_L",I2,"_T",F9.5,"_P",F5.1)') L,T1,P1    
 !   end if
    
 !   open(unit=14,file=filename)
 !   write(14,*) "#Data simulations (Temperature): ", T(:)

 !   do i=1,n_bin2D
 !     E = Emin + (i-0.5)*tam_binE2D
 !     do w=1,n_bin2D
 !       Dens = Densmin + (w-0.5)*tam_binDens2D
 !       if (memo_zero(i,w) == 0) then
 !         write(14,*) E, Dens, histogram_reweighted_grid(i,w)
 !       else
 !   !      write(14,*) E, Dens, 0
 !       end if
 !     end do
 !  !   write(14,*) " "
 !   end do
 !   close(14)

    if ( L < 10 ) then
      if ( P1 < 0 .or. P1 > 1000.0 ) then
        write(filename,'("Histogram_Reweighted_Free_Energy_L",I1,"_T",F9.5,"_P",F6.1)') L,T1,P1
      else
        write(filename,'("Histogram_Reweighted_Free_Energy_L",I1,"_T",F9.5,"_P",F5.1)') L,T1,P1
      end if
    else
      if ( P1 < 0 .or. P1 > 1000.0 ) then
        write(filename,'("Histogram_Reweighted_Free_Energy_L",I2,"_T",F9.5,"_P",F6.1)') L,T1,P1
      else
        write(filename,'("Histogram_Reweighted_Free_Energy_L",I2,"_T",F9.5,"_P",F5.1)') L,T1,P1
      end if
    end if
    open(unit=14,file=filename)
    write(14,*) "#Data simulations (Temperature): ", T(:)

    min_free = 1E12
    do i=1,n_bin2D
      E = Emin + (i-0.5)*tam_binE2D
      do w=1,n_bin2D
        Dens = Densmin + (w-0.5)*tam_binDens2D
        if (memo_zero(i,w) == 0) then
          if ( min_free > -log(histogram_reweighted_grid(i,w)) ) then
             min_free = -log(histogram_reweighted_grid(i,w));
          end if
        end if
      end do
    end do
    
    do i=1,n_bin2D
      E = Emin + (i-0.5)*tam_binE2D
      do w=1,n_bin2D
        Dens = Densmin + (w-0.5)*tam_binDens2D
        if (memo_zero(i,w) == 0) then
          write(14,*) E, Dens, -log(histogram_reweighted_grid(i,w))-min_free
     !   else
     !     write(14,*) E, Dens, 0
        end if
      end do
    !  write(14,*) " "
    end do
    close(14)

!    print*, "Emin =", Emin , " ; Emax =", Emax
!    print*, "Densmin =", Densmin, " ; Densmax =", Densmax

  end subroutine


  subroutine Check(array)
     real*8,dimension(n_bin_mixed),intent(in) :: array
     real*8 :: norma,media,var

     norma = 0
     media = 0
     var = 0
     do i=1,n_bin_mixed
       norma = norma + array(i)*tam_binX
       media = media + mixed_parameter(i)*array(i)*tam_binX  !0 si la distrib es simetrica
       var = var + mixed_parameter(i)*mixed_parameter(i)*array(i)*tam_binX
    end do
    var = var - media*media

    print*, "Media=",media, "Varianza=",var, "Norma=", norma

  end subroutine


  subroutine Calcula_His_tograma_Rotado_backup(alpha)    !resultados: histogram_reweighted_mixed()
 						 !            mixed_parameter()
    real*8,intent(in) :: alpha

    do i=1,n_bin
       do w=1,n_bin
          E = Emin + (i-0.5)*tam_binE
          Dens = Densmin + (w-0.5)*tam_binDens
          call Rotate(E,Dens,X,Y,alpha)
  !        call Mix(E,Dens,X,alpha)
          if (i==1 .and. w==1) then
            Xmax = X
            Xmin = X
          else
            if (Xmax < X) then 
                Xmax = X
            end if
            if (Xmin > X) then 
                Xmin = X
            end if
          end if
       end do
    end do

    rangeX = Xmax - Xmin
    tam_binX = rangeX/n_bin_mixed

    histogram_reweighted_mixed(:) = 0
    do i=1,n_bin
      do w=1,n_bin
        if (memo_zero(i,w) == 0) then
          E = Emin + (i-0.5)*tam_binE
          Dens = Densmin + (w-0.5)*tam_binDens
          call Rotate(E,Dens,X,Y,alpha)
  !        call Mix(E,Dens,X,alpha)
          indiceX = int((X-Xmin)/tam_binX)+1
          if(indiceX >= 1 .and. indiceX <= n_bin_mixed) then
              histogram_reweighted_mixed(indiceX) = histogram_reweighted_mixed(indiceX) + histogram_reweighted_grid(i,w)
          end if
        end if
      end do
    end do

    do i=1,n_bin_mixed
       mixed_parameter(i) = Xmin+(i-0.5)*tam_binX 
    end do

    delta_x = (mixed_parameter(1)+mixed_parameter(n_bin_mixed))/2.0
    
    mixed_parameter(:) = mixed_parameter(:) - delta_x     !traslado la distribución para centrarla en 0

    norma = 0
    media = 0
    var = 0
    do i=1,n_bin_mixed
      norma = norma + histogram_reweighted_mixed(i)*tam_binX
      media = media + mixed_parameter(i)*histogram_reweighted_mixed(i)*tam_binX  !0 si la distrib es simetrica
      var = var + mixed_parameter(i)*mixed_parameter(i)*histogram_reweighted_mixed(i)*tam_binX
    end do
    var = var - media*media
    a1 = dsqrt(norma/var)
    a2 = 1/(norma*a1)
    mixed_parameter(:) = mixed_parameter(:)*a1
    histogram_reweighted_mixed(:) = histogram_reweighted_mixed(:)*a2

    call Check(histogram_reweighted_mixed)

   !  histogram_reweighted_mixed(:) = histogram_reweighted_mixed(:)/norma
   !  mixed_parameter(:) = mixed_parameter(:)/sqrt(var)

  end subroutine
  
  subroutine binder(res)
     integer :: i
     real*8 :: m2,m4,u2,u4,u
     real*8,intent(out) :: res
     real*8 :: norma
     real*8, dimension(n_bin_mixed) :: order_parameter    
     
  !   print*, media, a1
     order_parameter(:) = (mixed_parameter(:)/a1)+media ! recalculate M=rho+s*e 
     
     norma=0
     do i=1,n_bin_mixed
        norma = norma + histogram_reweighted_mixed(i)
     end do
     
     u=0
     u2=0
     u4=0
     do i=1,n_bin_mixed
       m2 = order_parameter(i)*order_parameter(i)
       m4 = m2*m2
       u = u+order_parameter(i)*histogram_reweighted_mixed(i)/norma
       u2 = u2+m2*histogram_reweighted_mixed(i)/norma
       u4 = u4+m4*histogram_reweighted_mixed(i)/norma
    
     end do
 !  print*, "M ->", u2, u4, norma, u
     res = 1.0 - u4/(3.0*u2*u2)
  end subroutine
  
  subroutine binder_energy_density(res_ene,res_rho)
     integer :: i
     real*8 :: m2,m4,u2,u4
     real*8,intent(out) :: res_ene, res_rho
     real*8,dimension(n_bin2D) :: hist_ene, hist_rho
     real*8 :: ene, rho, norma
     
     hist_ene(:) = 0
     hist_rho(:) = 0
     norma = 0
     
     do i=1,n_bin2D
       do w=1,n_bin2D
         if (memo_zero(i,w) == 0) then
           hist_ene(i) = hist_ene(i) + histogram_reweighted_grid(i,w)
           hist_rho(w) = hist_rho(w) + histogram_reweighted_grid(i,w)
           norma = norma + histogram_reweighted_grid(i,w)
         end if
       end do
     end do
     
     hist_ene(:) = hist_ene(:) / norma
     hist_rho(:) = hist_rho(:) / norma
     
     ! Binder of energy
     u2=0
     u4=0
     do i=1,n_bin2D
       ene = Emin + (i-0.5)*tam_binE2D
       m2 = ene*ene
       m4 = m2*m2
       u2 = u2+m2*hist_ene(i)
       u4 = u4+m4*hist_ene(i)
      
     end do

     res_ene = 1.0 - u4/(3.0*u2*u2)
     
     ! Binder of density
     u2=0
     u4=0
     do i=1,n_bin2D
       rho = Densmin + (i-0.5)*tam_binDens2D
       m2 = rho*rho
       m4 = m2*m2
       u2 = u2+m2*hist_rho(i)
       u4 = u4+m4*hist_rho(i)
 
     end do

     res_rho = 1.0 - u4/(3.0*u2*u2)
     
  end subroutine
  
end program
