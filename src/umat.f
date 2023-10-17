!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2022 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!
      implicit none
!
!     Note that the following fields are not supported so far: 
!       sse,spd,scd,rpl,ddsddt,drplde,drpldt,predef, dpred,pnewdt,
!       celent,layer,kspt
!
!     Input stresses are Cauchy stresses (Eulerian)
!     Input strains are Logarithmic (ie True/Hencky) strains (Eulerian)
!
      character*80 cmname
!
      integer ndi,nshr,ntens,nstatv,nprops,noel,npt,layer,kspt,
     &  kstep,kinc
!
      real*8 stress(ntens),statev(nstatv),
     &  ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     &  stran(ntens),dstran(ntens),time(2),celent,
     &  props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     &  sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,predef,dpred,
     &  pnewdt
!
C     Material model for a 3D linear, viscoelastic material at small 
C     strains. The standard linear solid model is implemented, which is
C     equivalent to the Generalized Maxwell model with the summation 
C     only going to one. The derivations for the stress update algorithm
C     are based on the work from Kaliske and Rothert, "Formulation and
C     implementation of three-dimensional viscoelasticity at small and
C     finite strains" (1997). 
C
C                         Rheological Analogue Schematic
C          
C              |                 Linear elastic branch                |
C              |                   /\    /\    /\                     |
C              |------------------/  \  /  \  /  \  /-----------------|
C              |                      \/    \/    \/                  |
C              |                                                      |
C  	   <------ |                                                      | ------> 
C	           |      Linear elastic spring   Linear viscocity damper | 
C              |           /\    /\                 |-----            |
C              |----------/  \  /  \  /-------------|  |--------------|
C              |              \/    \/              |-----            |
C              |                                                      |  
!
C=============== Start New Variable Definitions Here ===================
      double precision, parameter :: SQRT2=1.414213562373095048801688724
!
      integer :: m, n
!      
      double precision :: CtermL, dose_prop, gamma, lameLambda,   
     &  lameShear, log_flag, poiss, prev_log_true, target_log_time, 
     &  tau, temp_prop, term_exp, tratio, trationeg, young
!
      double precision :: deps(6), dsign0(6), hn0(6), hnp1(6), 
     &  htermL(6), htermR(6), sign0(6), sign0p1(6), sigOld(6) 
!
      double precision :: Calgtan(6,6), Cel(6,6), CelInv(6,6),
     &   CmatVoi(6,6)
!
C================== Initialize Material Properties =====================  
      temp_prop = PROPS(1) ! [K] Temperature
      dose_prop = PROPS(2) ! [Gy] Dosage
      log_flag = PROPS(3)
      target_log_time = PROPS(4)
!
c     Calculate gamma [unitless] from surface fit
      call gamma_func_poly71(temp_prop, dose_prop, gamma)
!
c     Calculate tau [s], the relaxation time, from surface fit
      call tau_func_poly51(temp_prop, dose_prop, tau)
!
c     Calculate Young's modulus [Pa] from linear fit
      call young_func_poly1(dose_prop, young)
!
c     Poisson ratio [Unitless] from experiments (constant)
      poiss = 0.228250439894909
!      
C     Lame's constants; Lame's 1st parameter (lambda), and Lame's 2nd 
C     parameter also known as shear modulus (mu)
      lameShear = young/(2.D0*(1.D0+poiss))
      lameLambda = young*poiss/((1.D0+poiss)*(1.D0-2.D0*poiss))

C     Put stress (at time = n) into local variables.
C     Using SQRT(2.0) due to preference for Mandel notation.
      sigOld(1) = stress(1) ! [11]
      sigOld(2) = stress(2) ! [22]
      sigOld(3) = stress(3) ! [33]
      sigOld(4) = stress(4)*SQRT2 ! [12]
      sigOld(5) = stress(5)*SQRT2 ! [13]
      sigOld(6) = stress(6)*SQRT2 ! [23]
!
C     Put strain increment into local variables; SQRT(2.0) --> Mandel
C     UMAT shear strains are stored as engineering strains. Dividing by 
C     2 to convert them back into tensor components.
      deps(1) = dstran(1) ! [11]
      deps(2) = dstran(2) ! [22]
      deps(3) = dstran(3) ! [33]
      deps(4) = (dstran(4)/2.D0)*SQRT2 ! [12]
      deps(5) = (dstran(5)/2.D0)*SQRT2 ! [13]
      deps(6) = (dstran(6)/2.D0)*SQRT2 ! [23]
!
C     Grab the state variable, which is the elastic stress tensor of
C     the outermost spring element at time = n. 
      sign0(:) = statev(1:6)
!
C     Grab the state variable, which is the internal stress tensor of
C     the maxwell spring-dashpot element at time = n. 
      hn0(:) = statev(7:12)
!
C     Build the elastic stiffness matrix (in Mandel notation)
      Cel = 0.D0 ! Elastic stiffness matrix
      CelInv = 0.D0 ! Inverse in case it is needed later
      call buildMandelCMat(lameShear, lameLambda, Cel, CelInv)
!
C========================== Perform Updates ============================
C     Calculate increment in stress (dsign0), and update the elastic 
C     component of stress     
      call matdotvec(Cel, deps, 6, dsign0)
!
      sign0p1(1:6) = sign0(1:6) + dsign0(1:6) 
!
C     Update the internal stress tensor
      tratio = dtime/tau
      trationeg = -tratio
      term_exp = DEXP(trationeg)
!
      htermL(1:6) = term_exp*hn0(1:6)
      htermR(1:6) = (gamma*(1.D0 - term_exp)/tratio)*dsign0(1:6)
!   
      hnp1(1:6) = htermL(1:6) + htermR(1:6) ! Vector addition
!
C=============== Return Stresses and Internal Variables ================
C     Return stresses and take out of Mandel notation
      stress(1) = sign0p1(1) + hnp1(1)
      stress(2) = sign0p1(2) + hnp1(2)
      stress(3) = sign0p1(3) + hnp1(3)
      stress(4) = (sign0p1(4) + hnp1(4))/SQRT2
      stress(5) = (sign0p1(5) + hnp1(5))/SQRT2
      stress(6) = (sign0p1(6) + hnp1(6))/SQRT2
!
C     Leave internal stresses in Mandel notation for convenience
      statev(1:6) = sign0p1(1:6)
      statev(7:12) = hnp1(1:6)
!
C================== Return Consistent Tangent Modulus ==================
C     CalculiX expects the consistent tangent modulus in Voigt notation
      call mandelMat2VoigtMat(Cel, CmatVoi)
!
      CtermL = 1.D0 + gamma*(1.D0 - term_exp)/tratio
      Calgtan(1:6,1:6) = CtermL*CmatVoi(1:6,1:6)
!
      ddsdde(1:6,1:6) = Calgtan(1:6,1:6)
!
!
C======================== Log UMAT Information =========================
C     Log UMAT information, but only do so once.
      prev_log_true = statev(16)  
!$omp critical
      if (log_flag .gt. 0.0) then
        if (time(2) .ge. target_log_time) then
          if (time(2) .lt. (target_log_time + 1.1*dtime)) then
            if (prev_log_true .lt. 0.5) then

              open(unit=777, file='umat_log.txt', status='REPLACE', 
     &          action='WRITE', access='APPEND')

              write(777, *) 'Time (s),Temperature (K),Dosage (Gy),',
     &          'Youngs Modulus (Pa),Gamma,Relaxation Time (s),',
     &          'Poisson Ratio'

              write(777, 888) time(2), ',', temp_prop, ',', dose_prop, 
     &          ',', young, ',', gamma, ',', tau, ',', poiss

              prev_log_true = 1.0
              statev(16) = prev_log_true
              
            end if

            close(unit=777)

          end if
        end if
      end if
!
  888 format (F16.7, A, F16.7, A, F16.7, A, F16.7, A, F16.7, A,
     &  F16.7, A, F16.7)  
!
!$omp end critical
      end subroutine umat
!
!
      subroutine gamma_func_poly71(temp, dose, gamma)
      implicit none
      double precision, intent(in) :: temp, dose
      double precision, intent(inout) :: gamma
c     Calculates the material parameter, gamma, for the implemented
c     SLS viscoelastic model. Inputs are temperature (temp) and 
c     radiation dosage (dose). Temperatures should be in Kelvin and
c     range between 253 K and 333 K. Dosage should be in Grays and
c     range between 0 Gy and 16200 Gy. 
      double precision :: b10, b11, b12, b13, b14, b15,
     & a70, a61, a60, a51, a50, a41, a40, a31, a30, a21, a20, a11, a10,
     & a01, a00,
     & f8, f7, f6, f5, f4, f3, f2, f1, f0
!
      b10 = 287.6095684094453
      b11 = 1.0
      b12 = 1.0
      b13 = -32.21249762356403
      b14 = 1.0
      b15 = 21.78455023878911
      a70 = -1.594650344099281e-13
      a61 = 3.907025667117231e-11
      a60 = 2.570537624024945e-08
      a51 = -6.972097287568258e-08
      a50 = -3.411580223656593e-05
      a41 = 5.170611417078181e-05
      a40 = 0.016417639504081677
      a31 = -0.020396444964129942
      a30 = -2.8636085457493277
      a21 = 4.513287864875477
      a20 = -190.87085143382808
      a11 = -531.1417972716597
      a10 = 122295.43001152425
      a01 = 25970.90147915391
      a00 = -10978805.046977852
!
      f8 = b10*tanh(b11*temp + b12) + b13*tanh(b14*dose + b15)
      f7 = a70*(temp**7) + a61*(temp**6)*dose
      f6 = a60*(temp**6) + a51*(temp**5)*dose
      f5 = a50*(temp**5) + a41*(temp**4)*dose
      f4 = a40*(temp**4) + a31*(temp**3)*dose
      f3 = a30*(temp**3) + a21*(temp**2)*dose
      f2 = a20*(temp**2) + a11*temp*dose
      f1 = a10*temp + a01*dose
      f0 = a00
!
      gamma = f8 + f7 + f6 + f5 + f4 + f3 + f2 + f1 + f0 ! [Unitless]
      end subroutine gamma_func_poly71
!      
!
      subroutine tau_func_poly51(temp, dose, tau)
      implicit none
      double precision, intent(in) :: temp, dose
      double precision, intent(inout) :: tau
c     Calculates the material parameter, tau, for the implemented
c     SLS viscoelastic model. Tau is the relaxation time in seconds.
c     Inputs are temperature (temp) and radiation dosage (dose). 
c     Temperatures should be in Kelvin and range between 253 K and 
c     333 K. Dosage should be in Grays and range between 0 Gy and 
c     16200 Gy. 
      double precision :: a50, a41, a40, a31, a30, a21, a20, a11, a10,
     & a01, a00,
     & f5, f4, f3, f2, f1, f0
!
      a50 = -4.227010936820378e-09
      a41 = 7.79458980940624e-12
      a40 = 6.581726975789612e-06
      a31 = -9.56299640007316e-09
      a30 = -0.00409652524634327
      a21 = 4.402198819577896e-06
      a20 = 1.2743264615947765
      a11 = -0.0009016025756702695
      a10 = -198.1665980531809
      a01 = 0.06935013859524257
      a00 = 12326.614896755611
!
      f5 = a50*(temp**5) + a41*(temp**4)*dose
      f4 = a40*(temp**4) + a31*(temp**3)*dose
      f3 = a30*(temp**3) + a21*(temp**2)*dose
      f2 = a20*(temp**2) + a11*temp*dose
      f1 = a10*temp + a01*dose
      f0 = a00
!
      tau = f5 + f4 + f3 + f2 + f1 + f0 ! [s]
      end subroutine tau_func_poly51
!
!
      subroutine young_func_poly1(dose, young)
      implicit none
      double precision, intent(in) :: dose
      double precision, intent(inout) :: young
c     Comments
      double precision :: a1, a0
!
      a1 = -7.3111661435e-01
      a0 = 2.4251564601e+04
!
      young = a1*dose + a0 ! [Pa]
      end subroutine young_func_poly1
!
!
      subroutine buildMandelCMat(lameShear, lameLambda, cOut, cInvOut)
      implicit none
      double precision, intent(in) :: lameShear, lameLambda
      double precision, intent(inout) :: cOut(6,6), cInvOut(6,6)
c     Builds the elastic stiffness matrix into, cOut, as well as its 
c     inverse, cInvOut. NOTE, the matrix assumes TENSOR strains to 
c     stresses! A factor of 2 has been added since engineering 
c     shear strains are not used. In other words, Mandel notation is 
c     adopted.
      integer :: m, n, np
      double precision :: young, poiss
      np = 3
      young = (lameShear*(3.D0*lameLambda+2.D0*lameShear))/
     &        (lameLambda+lameShear)
      poiss = lameLambda/(2.D0*(lameLambda+lameShear))

      cOut = 0.D0
      do m=1,np
        do n=1,np
          if (m .ne. n) cOut(m,n) = lameLambda
        end do
        cOut(m,m) = 2.D0*lameShear + lameLambda
        cOut(m+3,m+3) = 2.D0*lameShear !NOTE: There's a 2 here
      end do

      cInvOut = 0.D0
      do m=1,np
        do n=1,np
          if (m .ne. n) cInvOut(m,n) = -poiss/young
        end do
        cInvOut(m,m) = 1.D0/young
        cInvOut(m+3,m+3) = (1.D0/young)*(1.D0 + poiss) !NOTE: Lack of 2 
      end do
      end subroutine buildMandelCMat
!
!
      subroutine buildVoigtCMat(lameShear, lameLambda, cOut, cInvOut)
      implicit none
      double precision, intent(in) :: lameShear, lameLambda
      double precision, intent(inout) :: cOut(6,6), cInvOut(6,6)
c     Builds the elastic stiffness matrix into, cOut, as well as its 
c     inverse, cInvOut. NOTE, the matrix assumes ENGINEERING strains to 
c     stresses! A factor of 2 has been added since engineering 
c     shear strains are not used. In other words, Voigt notation is 
c     adopted.
      integer :: m, n, np
      double precision :: young, poiss
      np = 3
      young = (lameShear*(3.D0*lameLambda+2.D0*lameShear))/
     &        (lameLambda+lameShear)
      poiss = lameLambda/(2.D0*(lameLambda+lameShear))

      cOut = 0.D0
      do m=1,np
        do n=1,np
          if (m .ne. n) cOut(m,n) = lameLambda
        end do
        cOut(m,m) = 2.D0*lameShear + lameLambda
        cOut(m+3,m+3) = lameShear !NOTE: Lack of 2 here
      end do

      cInvOut = 0.D0
      do m=1,np
        do n=1,np
          if (m .ne. n) cInvOut(m,n) = -poiss/young
        end do
        cInvOut(m,m) = 1.D0/young
        cInvOut(m+3,m+3) = (2.D0/young)*(1.D0 + poiss) !NOTE: The 2
      end do
      end subroutine buildVoigtCMat
!
!
      subroutine matdotvec(Amat, bvec, np, cvec)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: Amat(np,np), bvec(np)
      double precision, intent(inout) :: cvec(np)
c     Calculates A [dot] b = c, or index notation, c_i = A_ij * b_j
c      integer :: i, j, m
c
c      m = np
      cvec = 0.
c     do j=1,m
c       do i=1,m
c         cvec(i) = cvec(i) + Amat(i,j)*bvec(j)
c       end do
c     end do
c
      cvec = MATMUL(Amat, bvec)
      end subroutine matdotvec
!
!
      subroutine mandelMat2VoigtMat(matIn, matOut)
      implicit none
c     Takes a matrix in Mandel notation and converts it to Voigt (for 
c     what is usually done with the kinematic terms).    
      double precision, parameter :: SQRT2=1.414213562373095048801688724
      double precision, intent(in) :: matIn(6,6) 
      double precision, intent(inout) :: matOut(6,6)
      integer :: m, n
      double precision :: matTemp(6,6)
c      
      matOut = 0.D0
      matTemp = 0.D0
      do n=1,3
        do m=1,3
          matTemp(m,n) = matIn(m,n)
        end do
      end do

      do n=4,6 
        do m=1,3
          matTemp(m,n) = matIn(m,n)/SQRT2
        end do
      end do

      do n=1,3 
        do m=4,6
          matTemp(m,n) = matIn(m,n)/SQRT2
        end do
      end do

      do n=4,6
        do m=4,6
          matTemp(m,n) = matIn(m,n)/2.D0
        end do
      end do

      matOut(:,:) = matTemp(1:6,1:6)
      end subroutine mandelMat2VoigtMat
!
!
      subroutine umat_BACKUP(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!
      implicit none
!
!     Note that the following fields are not supported so far: 
!       sse,spd,scd,rpl,ddsddt,drplde,drpldt,predef, dpred,pnewdt,
!       celent,layer,kspt
!
!     Input stresses are Cauchy stresses (Eulerian)
!     Input strains are Logarithmic (ie True/Hencky) strains (Eulerian)
!
      character*80 cmname
!
      integer ndi,nshr,ntens,nstatv,nprops,noel,npt,layer,kspt,
     &  kstep,kinc
!
      real*8 stress(ntens),statev(nstatv),
     &  ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     &  stran(ntens),dstran(ntens),time(2),celent,
     &  props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     &  sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,predef,dpred,
     &  pnewdt
!
C     Material model for a 3D linear, viscoelastic material at small 
C     strains. The standard linear solid model is implemented, which is
C     equivalent to the Generalized Maxwell model with the summation 
C     only going to one. The derivations for the stress update algorithm
C     are based on the work from Kaliske and Rothert, "Formulation and
C     implementation of three-dimensional viscoelasticity at small and
C     finite strains" (1997).   
!
C=============== Start New Variable Definitions Here ===================
      double precision, parameter :: SQRT2=1.414213562373095048801688724
!
      integer :: m, n
!      
      double precision :: CtermL, gamma, lameLambda, lameShear, poiss,  
     &  tau, term_exp, tratio, trationeg, young
!
      double precision :: deps(6), dsign0(6), hn0(6), hnp1(6), 
     &  htermL(6), htermR(6), sign0(6), sign0p1(6), sigOld(6) 
!
      double precision :: Calgtan(6,6), Cel(6,6), CelInv(6,6),
     &   CmatVoi(6,6)
!
C================== Initialize Material Properties =====================  
C     Material properties for the elastic contribution
      young = PROPS(1) ! [Pa] Young's modulus 
      poiss = PROPS(2) ! Poisson's ratio
!
C     Material properties for the internal stress variables h_j
      gamma = PROPS(3) ! Normalized elastic contribution
      tau = PROPS(4)   ! [s] Relaxation time of the Maxwell element     
!      
C     Lame's constants; Lame's 1st parameter (lambda), and Lame's 2nd 
C     parameter also known as shear modulus (mu)
      lameShear = young/(2.D0*(1.D0+poiss))
      lameLambda = young*poiss/((1.D0+poiss)*(1.D0-2.D0*poiss))

C     Put stress (at time = n) into local variables.
C     Using SQRT(2.0) due to preference for Mandel notation.
      sigOld(1) = stress(1) ! [11]
      sigOld(2) = stress(2) ! [22]
      sigOld(3) = stress(3) ! [33]
      sigOld(4) = stress(4)*SQRT2 ! [12]
      sigOld(5) = stress(5)*SQRT2 ! [13]
      sigOld(6) = stress(6)*SQRT2 ! [23]
!
C     Put strain increment into local variables; SQRT(2.0) --> Mandel
C     UMAT shear strains are stored as engineering strains. Dividing by 
C     2 to convert them back into tensor components.
      deps(1) = dstran(1) ! [11]
      deps(2) = dstran(2) ! [22]
      deps(3) = dstran(3) ! [33]
      deps(4) = (dstran(4)/2.D0)*SQRT2 ! [12]
      deps(5) = (dstran(5)/2.D0)*SQRT2 ! [13]
      deps(6) = (dstran(6)/2.D0)*SQRT2 ! [23]
!
C     Grab the state variable, which is the elastic stress tensor of
C     the outermost spring element at time = n. 
      sign0(:) = statev(1:6)
!
C     Grab the state variable, which is the internal stress tensor of
C     the maxwell spring-dashpot element at time = n. 
      hn0(:) = statev(7:12)
!
C     Build the elastic stiffness matrix (in Mandel notation)
      Cel = 0.D0 ! Elastic stiffness matrix
      CelInv = 0.D0 ! Inverse in case it is needed later
      call buildMandelCMat(lameShear, lameLambda, Cel, CelInv)
!
C========================== Perform Updates ============================
C     Calculate increment in stress (dsign0), and update the elastic 
C     component of stress     
      call matdotvec(Cel, deps, 6, dsign0)
!
      sign0p1(1:6) = sign0(1:6) + dsign0(1:6) 
!
C     Update the internal stress tensor
      tratio = dtime/tau
      trationeg = -tratio
      term_exp = DEXP(trationeg)
!
      htermL(1:6) = term_exp*hn0(1:6)
      htermR(1:6) = (gamma*(1.D0 - term_exp)/tratio)*dsign0(1:6)
!   
      hnp1(1:6) = htermL(1:6) + htermR(1:6) ! Vector addition
!
C=============== Return Stresses and Internal Variables ================
C     Return stresses and take out of Mandel notation
      stress(1) = sign0p1(1) + hnp1(1)
      stress(2) = sign0p1(2) + hnp1(2)
      stress(3) = sign0p1(3) + hnp1(3)
      stress(4) = (sign0p1(4) + hnp1(4))/SQRT2
      stress(5) = (sign0p1(5) + hnp1(5))/SQRT2
      stress(6) = (sign0p1(6) + hnp1(6))/SQRT2
!
C     Leave internal stresses in Mandel notation for convenience
      statev(1:6) = sign0p1(1:6)
      statev(7:12) = hnp1(1:6)
!
C================== Return Consistent Tangent Modulus ==================
C     CalculiX expects the consistent tangent modulus in Voigt notation
      call mandelMat2VoigtMat(Cel, CmatVoi)
!
      CtermL = 1.D0 + gamma*(1.D0 - term_exp)/tratio
      Calgtan(1:6,1:6) = CtermL*CmatVoi(1:6,1:6)
!
      ddsdde(1:6,1:6) = Calgtan(1:6,1:6)
!
      end subroutine umat_BACKUP
