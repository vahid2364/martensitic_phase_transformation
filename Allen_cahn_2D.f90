!! filename = CH.f90
!! Vahid Attari
!! Created: 30 Feb. 2016
!! Modified: ....
!! Arroyave Research Group, Department of Materials Science & Engineering, Texas A&M University
!!
!! Acknowledgements:  Based on Cahn-Hilliard 1965 paper
!!
!! Purpose:
!!   - Phase Field Modeling with dynamic coupling to thermodyanmic and kinetic databases
!!     to self consistantly model the Spinodal Composition Phenomenon
!!
!! General Algorithm function:
!!
!!   1. Retrieve parameter data from file "parameters.dat"
!!   2. Assess thermodynamics of the associated system
!!   3. Reads initial phase distribution from "phase.dat" file
!!   4. Calculate Phase Evolution with time integration
!!      -  Nucleate Phases
!!      -  Resolve boundary conditions
!!         -- Periodic boundaries in all directions
!!      -  Solve differential equations via 9-stencil finite difference
!!      -  Update phase information and concentration data
!!
!! Compilation instructions: >> make
!!    - Manual: >>  ifort -mkl $MKLROOT/include/mkl_dfti.f90 -o a.out 'filename'.f90
!!
!! Execution: >> ./a.out
!!
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!====================================================================================

!   This code simulates the early stage of Spinodal Decomposition...
!   Solving 2D Allen-Cahn Eq using pseudo-spectral with Implicit/Explicit
!   u_t= epsilon(u_{xx}+u_{yy}) + u - u^3
!   where u-u^3 is treated explicitly and epsilon(u_{xx} + u_{yy}) is treated implicitly
!   BC = Periodic
!   IC=v=sin(2*pi*x)+0.001*cos(16*pi*x);

!   Allen- Cahn solver
!   with periodic boundary conditions.

module mod_geometry

    CHARACTER (*), parameter :: mic_type = 'random' !2_grain !4_grain !Read_input_file !inclusion ! sphere ! random

    integer, parameter:: kind=8

    integer, parameter::nex=256, ney=256, nez=nex

    integer, parameter::IG=nex
    integer, parameter::JG=ney
    integer, parameter::ZG=nex

    REAL(8) , PARAMETER :: Lx = 50.0D-9;    !........................................Length of x domain
    REAL(8) , PARAMETER :: Ly = 50.0D-9;    !........................................Length of y domain
    REAL(8) , PARAMETER :: Lz = 50.0D-9;    !........................................Length of z domain (if used)

    real(8) , parameter :: h  = Lx/(nex-1);

    real(8) , parameter :: dx = h
    real(8) , parameter :: dy = h

    INTEGER, PARAMETER :: case_num  = 1;        ! 1=circle , 2=random , 3=sphere
    INTEGER, PARAMETER :: dimen     = 2;        !........................................cartesian dimensions used (i.e., 1-D, 2-D, 3-D)
    !INTEGER, PARAMETER :: centers   = 5;       !........................................

    REAL(8), PARAMETER :: max_radius = 10;      !........................................Maximum Nuclei Radius
    REAL(8), PARAMETER :: grad_factor= 0;       !........................................
    INTEGER, PARAMETER :: sep_dis    = 5;       !........................................

    real, dimension(nex,ney) :: lap_phi,f
    real :: phitot
    real :: Betta_C,Betta_Max,Landa_C,Landa_Max,R_Max,betta
    integer, save :: l1, m1

end module mod_geometry

!==============================================================

MODULE mtls_property

    USE mod_geometry

    real(8) :: X11 = 0.5D0

    real(8), parameter :: T = 720.0D0 + 273.0D0
    real(8), parameter :: R = 8.314D0
    real(8), parameter :: VM = 1.0 !1.205883D-5            				  ! m^3/mol

    ! Kinetic
    real(8), parameter :: L = 2596.5/VM             				  ! m^2/N.s
    real(8), parameter :: beta = 4.50D-10/VM;        				  ! (N)
    real(8), parameter :: delta_f = 10.0D9/VM;         				  !J/m^3

    ! Landau parameters
    real(8) :: A = 0.14D0
    real(8) :: B = 12.42
    real(8) :: C = 12.28


    
    REAL(8), PARAMETER :: PI  = 4.0D0*ATAN(1.0D0) ! PI VALUE 3.141592......

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE mtls_property

!==============================================================

MODULE elastic_property

    USE mod_geometry, only : kind,nex,ney
    USE mtls_property, only : T

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL(8), PARAMETER :: theta = 0.0;        !...............................................rotation
    REAL(8), PARAMETER :: A_b = 1.0;          !......................................Constant that determines the energy barrier between the m and p phases

    !!
    !! Crystalographic pts...
    !!

!    REAL(4), DIMENSION(3,3) :: eps_T1  = eps_T*[ 0.40, 0.0, 0.0, 0.0, 0.40, 0.0, 0.0, 0.0, 0.0 ]    				  ! AlN wrt to TiN   
!    REAL(4), DIMENSION(3,3) :: eps_T2  = eps_T*[ 0.80, 0.0, 0.0, 0.0, 0.80, 0.0, 0.0, 0.0, 0.0 ]    				  ! ZrN wrt to TiN   
!    REAL(4), DIMENSION(3,3) :: eps_T3  = eps_T*[ 0.00, 0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0, 0.0 ]   					  ! TiN wrt to TiN   

    !!
    !! Physical pts...
    !!
    
    REAL(8) , PARAMETER :: V_m  = 11.85e-6   ! kg/m^3 ! AlN (volume)
    REAL(8) , PARAMETER :: V_p1 = 12.57e-6   ! kg/m^3 ! TiN (volume) 
    REAL(8) , PARAMETER :: V_p2 = 14.40e-6   ! kg/m^3 ! ZrN (volume)

    REAL(8) , PARAMETER :: Nv_m  = 1.0 ! kg/m^3 ! TiN (volume)
    REAL(8) , PARAMETER :: Nv_p1 = 1.0 ! kg/m^3 ! AlN (volume)
    REAL(8) , PARAMETER :: Nv_p2 = 1.0 ! kg/m^3 ! ZrN (volume)
    
    !!
    !! Elastic Modulus pts...
    !!

    !! M phase TiN (Fm3m cubic)
    REAL(4) , PARAMETER ::    C44_p1   = 279.0D9;               ! GPa
    REAL(4) , PARAMETER ::    C12_p1   = 250.0D9;               ! GPa
    REAL(4) , PARAMETER ::    C11_p1   = 384.0D9;               ! GPa
    REAL(4) , PARAMETER ::    C22_p1   = 384.0D9;               ! GPa
    REAL(4) , PARAMETER ::    G_p1     = C44_p1             												! GPa

    !! p phase AlN (F43m cubic)
    REAL(4) , PARAMETER ::    C44_m   = C44_p1
    REAL(4) , PARAMETER ::    C12_m   = C12_p1        	    !.............Circularly averaged Voigt constant C11
    REAL(4) , PARAMETER ::    C11_m   = C11_p1               !.............Circularly averaged Voigt constant C12
    REAL(4) , PARAMETER ::    C22_m   = C22_p1               !.............Circularly averaged Voigt constant C12
    REAL(4) , PARAMETER ::    G_m     = C44_m            												! GPa

    !! p phase ZrN ( cubic)
    REAL(4) , PARAMETER ::    C44_p2   = 118.0D9        
    REAL(4) , PARAMETER ::    C12_p2   = 107.0D9       	            !.............Circularly averaged Voigt constant C11
    REAL(4) , PARAMETER ::    C11_p2   = 524.0D9                      !.............Circularly averaged Voigt constant C12
    REAL(4) , PARAMETER ::    C22_p2   = 524.0D9                      !.............Circularly averaged Voigt constant C12

    REAL(4) , PARAMETER ::    C44      = 141.0D9        !..................................
    REAL(4) , PARAMETER ::    C11      = 488.0D9        !..........................Circularly averaged Voigt constant C11
    REAL(4) , PARAMETER ::    C22      = 488.0D9        !..........................Circularly averaged Voigt constant C11
    REAL(4) , PARAMETER ::    C12      = 137.0D9        !..........................Circularly averaged Voigt constant C12

    REAL(kind) , PARAMETER :: A_z   = 2.0*C44  /(C11  +C12  );               !...............................................Average Zener anisotropy parameter
    REAL(kind) , PARAMETER :: A_zm  = 2.0*C44_m/(C11_m+C12_m);               !...............................................Average Zener anisotropy parameter
    REAL(kind) , PARAMETER :: A_zp1 = 2.0*C44_p1/(C11_p1+C12_p1);            !...............................................Average Zener anisotropy parameter
    REAL(kind) , PARAMETER :: A_zp2 = 2.0*C44_p2/(C11_p2+C12_p2);            !...............................................Average Zener anisotropy parameter
    REAL(kind) , PARAMETER :: Ap_m = C44_m /(C11_m+2*C12_m)
    REAL(kind) , PARAMETER :: Ap_p1= C44_p1/(C11_p1+2*C12_p1)
    REAL(kind) , PARAMETER :: Ap_p2= C44_p2/(C11_p2+2*C12_p2)
	REAL(kind) , PARAMETER :: B_m = (C11_m+C12_m);!..................................
	REAL(kind) , PARAMETER :: B_p1= (C11_p1+C12_p1);!..................................
	REAL(kind) , PARAMETER :: B_p2= (C11_p1+C12_p1);!..................................
	REAL(kind) , PARAMETER :: B_v = (C11+C12);		!..................................
	REAL(kind) , PARAMETER :: delta1 = C44_p1/C44_m;!..................................
	REAL(kind) , PARAMETER :: delta2 = C44_p2/C44_m;!..................................


    !! Voigt Notation ...
    INTEGER,DIMENSION(3,3)               :: voigt = (/ 1,6,5,6,2,4,5,4,3 /)

    !! ARRAY PARAMETERS ...
    REAL(kind), DIMENSION(1:nex,1:ney)   :: Beta_C1_prime,Beta_C2_prime!,!Beta_C1,Beta_C2,
    REAL(kind), DIMENSION(1:nex,1:ney,1:6,1:6) :: C_c
    REAL(kind), DIMENSION(1:nex,1:ney)   :: Alpha_C1,Alpha_C2,Alpha_C1_prime,Alpha_C2_prime,Alpha_Beta_C1,Alpha_Beta_C2

END MODULE elastic_property


MODULE solver_opts

    USE mod_geometry !, only : dx
    USE mtls_property!, only : Mobility

    !REAL(8), PARAMETER :: dt= 6.195D-3*(dx)**2.0/(Mobility)   !1.0D+55*((dx)**2.0/2.0*Mobility)
    REAL(8) :: dt = 2.0D-12;

    INTEGER, PARAMETER :: nth_order  = 20;                !....................................The order of approximation of the displacement field
    INTEGER, PARAMETER :: order      = nth_order + 1;     !....................................The order of approximation of the displacement field adjusted for python
    INTEGER, PARAMETER :: els_update = 25;                !.................................... 
    REAL(8), PARAMETER :: elas_err   = 1.0D-7;            !....................................The order of approximation of the displacement field
    
    REAL(8) :: TIME

END MODULE solver_opts


MODULE wrt_opts

    CHARACTER ch*9
    CHARACTER(LEN=100) :: VAL
    CHARACTER(len=255) :: cwd,fileplace

    INTEGER :: wrt_cycle = 100     !100000
    INTEGER :: NNN3      = 1500000

    integer            :: IPC
    integer            :: itimes
    CHARACTER(len=9)   :: num

CONTAINS

    SUBROUTINE mkdir

        !!***** MAKE DIRs *****
        CALL getcwd(cwd)
        WRITE(fileplace,*) ADJUSTL(TRIM(cwd))

        call system('rm -r results')
        call system('rm -r images')
        call system('rm -r microstructure')
        call system('rm -r elasticity')
        call system('rm -r *.ppm')
    
        call system('mkdir images')
        call system('mkdir results')
        call system('mkdir microstructure')
        call system('mkdir elasticity')
        call system('mkdir elasticity/strain')
        call system('mkdir elasticity/stress')
        call system('mkdir elasticity/displacement')

    END SUBROUTINE mkdir

    SUBROUTINE wrt_disk

        !!***** MAKE DIRs *****
        CALL getcwd(cwd)
        WRITE(fileplace,*) ADJUSTL(TRIM(cwd))

    END SUBROUTINE wrt_disk

END MODULE wrt_opts

Module pixelfrt

	USE mod_geometry, only : nex, ney, nez

       integer ihpixf, jvpixf
       parameter(ihpixf = nex, jvpixf = ney) ! pixel and data size
       
	Contains
	
! --------------------------------------------
! Sample FORTRAN 77 program to make PPM / BMP
!
!  1998-Apr-03 (long long ago)
!      :
!  2005-Jul-22 (added P3)
!      :
!  2006 Dec-22 (added BMP)
!      :
!  2016 Jan-15 (last change)
!
!  Image array rgb(3,*,*) is filled in subroutine mkbitmap()
!  and
!  saved in PPM or BMP format in subroutine pixout().
!
!                                   K. Hayashi
! --------------------------------------------
!

    Subroutine pixelout(phi,nframe)
       implicit none
       real(8), dimension(1:nex,1:ney,1:2)    :: PHI    
       character*1 rgb(3,ihpixf,jvpixf) ! RGB image array (integer)
       real*8      dat(3,ihpixf,jvpixf) ! three, 2D data array (float or double)
       integer nframe, nf2
		
		nf2     = nframe
		dat(1,:,:)     = phi(:,:,1) !reshape( phi(:,:,1), (/ , 3 /) )
		dat(2,:,:)     = phi(:,:,1) 
		dat(3,:,:)     = phi(:,:,1) 

       call mkbitmap(dat,rgb)
       call pixout(rgb,nf2)		
		
       !do nframe = 1, 50
       !  nf2 = nframe
       !  call mkdata(dat,nf2)
       !  call mkbitmap(dat,rgb)
       !  call pixout(rgb,nf2)
       !enddo
       
       !do nframe = 1, 50
       !  nf2 = nframe
       !  call mkbitmap(rgb,nf2)
       !  call pixout(rgb,nf2)
       !enddo

       
    end Subroutine pixelout


! --------------------------------------------
!
! Notes
! o With a parameter ipixout set at 1, 2 or others,
!   this subroutine will generate PPM-P6(binary), PPM-P3(text),
!   or BMP(24bit depth without color table).
!
! o Some parts follow DEC-FORTRAN convention that had been defacto-standard long ago.
!   Some compilers today may not accept if "ipixout" is set other than 2.
!
! o g77 (ver. 3.3.3) works for all three choices.
! o Intel compiler (ver. 9 or later) works for all three choices.
!
! --------------------------------------------
!
       subroutine pixout(rgb,nframe)
       implicit none
! interface arg.
       integer ihpixf, jvpixf
       parameter(ihpixf = nex, jvpixf = ney) ! pixel size, eacg must be multiple of 4, if BMP is chosen as output format.
       character*1 rgb(3,ihpixf,jvpixf)      ! RGB data array
       integer nframe
! local
       character*19 fnameout
       integer i, j, k
       integer itmp, icnt
       character*14 frmtstr
       character*54 headmsw
       character*4  byt4
       character*2  byt2
! choices
       integer ipixout
       parameter(ipixout = 1) ! 1 / 2 / other= PPM6, PPM3, BMP(24bit)

       if (ipixout .EQ. 1) then

! PPM P6

         write(fnameout,'(''images/smpl'',i3.3,''.ppm'')') nframe ! name of PPM file
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now writing PPM (P6) file : ', fnameout
! header
         write(2,'(''P6'', 2(1x,i4),'' 255 '',$)') ihpixf, jvpixf        ! some compiler may not accept this line.
     
! image data
         itmp = ihpixf * jvpixf * 3
         write(frmtstr,'(''('',i8.8,''A,$)'')') itmp     ! make output "format"
         write(2,fmt=frmtstr) (((rgb(k,i,j),k=1,3),i=1,ihpixf),j=jvpixf,1,-1)                             ! some compiler may not accept this line.      ! here, j (vertical address) runs from top to bottom.
         close(2)

       else if (ipixout .EQ. 2) then

! PPM P3 ! rather "safer" choice for many Fortran compiler(s).

         write(fnameout,'(''smpl'',i3.3,''.ppm'')') nframe ! name of PPM file
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now writing PPM (P3) file : ', fnameout
! header
         write(2,'(A)') 'P3'
         write(2,'(2(1x,i4),'' 255 '')')  ihpixf, jvpixf
         icnt = 0
! image data
         do j = jvpixf, 1, -1                              ! here, j (vertical address) runs from top to bottom.
         do i = 1, ihpixf, 1
         do k = 1, 3
           itmp = ichar(rgb(k,i,j))
           icnt = icnt + 4
           if (icnt .LT. 60) then
             write(2,fmt='(1x,i3,$)') itmp                 ! "$" is not standard.
           else
             write(2,fmt='(1x,i3)') itmp
             icnt = 0
           endif
         enddo
         enddo
         enddo
         write(2,'(A)') ' '
         close(2)

       else

! BMP (24bit depth)... this part works only when width is multiple of 4.

         itmp = mod(ihpixf, 4)
         if (itmp .NE. 0) then
           write(*,*) 'width must be multiple of 4'
           stop
         endif

         write(fnameout,'(''images/smpl'',i3.3,''.bmp'')') nframe ! name of BMP file
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now writing BMP(24bit) file : ', fnameout
! header 1 (file header ; 1--14 byte)
         headmsw( 1: 2) = 'BM'             ! declaring this is BMP file
         itmp = 54 + ihpixf * jvpixf * 3 ! total file size = header + data
         call num2bit4(itmp,byt4)
         headmsw( 3: 6) = byt4(1:4)
         itmp = 0                        ! may be 0
         call num2bit2(itmp,byt2)
         headmsw( 7: 8) = byt2(1:2)
         itmp = 0                        ! may be 0
         call num2bit2(itmp,byt2)
         headmsw( 9:10) = byt2(1:2)
         itmp = 54                       ! must be 54 : total length of header
         call num2bit4(itmp,byt4)
         headmsw(11:14) = byt4(1:4)
! header 2 (bit-map header ; 13--54 byte)
         itmp = 40                       ! must be 40 : length of bit-map header
         call num2bit4(itmp,byt4)
         headmsw(15:18) = byt4(1:4)
         itmp = ihpixf                   ! width
         call num2bit4(itmp,byt4)
         headmsw(19:22) = byt4(1:4)
         itmp = jvpixf                   ! height
         call num2bit4(itmp,byt4)
         headmsw(23:26) = byt4(1:4)
         itmp = 1                        ! must be 1
         call num2bit2(itmp,byt2)
         headmsw(27:28) = byt2(1:2)
         itmp = 24                       ! must be 24 : color depth in bit.
         call num2bit2(itmp,byt2)
         headmsw(29:30) = byt2(1:2)
         itmp = 0                        ! may be 0 : compression method index
         call num2bit4(itmp,byt4)
         headmsw(31:34) = byt4(1:4)
         itmp = 0                        ! may be 0 : file size if compressed
         call num2bit4(itmp,byt4)
         headmsw(35:38) = byt4(1:4)
         itmp = 0                        ! arbit. : pixel per meter, horizontal
         call num2bit4(itmp,byt4)
         headmsw(39:42) = byt4(1:4)
         itmp = 0                        ! arbit. : pixel per meter, vertical
         call num2bit4(itmp,byt4)
         headmsw(43:46) = byt4(1:4)
         itmp = 0                        ! may be 0 here : num. of color used
         call num2bit4(itmp,byt4)
         headmsw(47:50) = byt4(1:4)
         itmp = 0                        ! may be 0 here : num. of important color
         call num2bit4(itmp,byt4)
         headmsw(51:54) = byt4(1:4)

! writing header part
         write(2,'(a54,$)') headmsw(1:54)
! image data
         itmp = ihpixf * jvpixf * 3
         write(frmtstr,'(''('',i8.8,''A,$)'')') itmp
         write(2,fmt=frmtstr)    (((rgb(k,i,j),k=3,1,-1),i=1,ihpixf),j=1,jvpixf) ! writing in BGR order, not RGB.
         close(2)

       endif

       return
       end subroutine pixout

! --------------------------------------
! convert integer values to 4 8-bit characters
! --------------------------------------

       subroutine num2bit4(inum,byt4)
       implicit none
       integer inum
       character*4 byt4
       integer itmp1, itmp2
       itmp1 = inum
       itmp2 = itmp1 / 256**3
       byt4(4:4) = char(itmp2)
       itmp1 =-itmp2 * 256**3 +itmp1
       itmp2 = itmp1 / 256**2
       byt4(3:3) = char(itmp2)
       itmp1 =-itmp2 * 256**2 +itmp1
       itmp2 = itmp1 / 256
       byt4(2:2) = char(itmp2)
       itmp1 =-itmp2 * 256    +itmp1
       byt4(1:1) = char(itmp1)
       return
       end subroutine num2bit4

! --------------------------------------
! convert integer values to 2 8-bit characters
! --------------------------------------

       subroutine num2bit2(inum,byt2)
       implicit none
       integer inum
       character*2 byt2
       integer itmp1, itmp2
       itmp1 = inum
       itmp2 = itmp1 / 256
       byt2(2:2) = char(itmp2)
       itmp1 =-itmp2 * 256 + itmp1
       byt2(1:1) = char(itmp1)
       return
       end subroutine num2bit2


! --------------------------------------------
!   making dummy data
! --------------------------------------------

       subroutine mkdata(dat,nframe)
       implicit none
       integer nframe
       integer ihpixf, jvpixf
       parameter(ihpixf = 128, jvpixf = 128) ! pixel and data size
       real*8      dat(3,ihpixf,jvpixf)
! choices for generating three two-dim. float data
       integer ichoice
       parameter(ichoice = 1)
! local
       integer i, j, itmp
! some parameters for making dummy data
       real*8  ofst
       parameter(ofst = 0.7D+00)
       real*8   aa, bb, cc, rr, xx, yy, tt
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
       real*8  v1, v2, v3

       do 100 j = 1, jvpixf
       do 100 i = 1, ihpixf
         if (ichoice .EQ. 0) then
           itmp = i*3*nframe + j*2
           itmp = mod(itmp,256)    ! assuming 8-bit color depth, rangeing from 0 to 255
           v1 = dfloat(itmp) / 256.0D+00
           itmp = i*1*nframe + j*3
           itmp = mod(itmp,256)
           v2 = dfloat(itmp) / 256.0D+00
           itmp = i*5*nframe + j*7
           itmp = mod(itmp,256)
           v3 = dfloat(itmp) / 256.0D+00
         else
! v1-ball
           tt = dfloat(nframe) / 25.0D+00 !                  time/period
           xx = dfloat(i) / dfloat(ihpixf) - 0.33D+00  !     center x
           yy = dfloat(j) / dfloat(jvpixf) - 0.25D+00  !     center y
           rr = dsqrt(xx**2 + yy**2 + 1.0D-30)
           aa = rr / 0.25D+00 !                              half-width
           bb =(tt - rr) * pi
           cc = dexp(-aa**2) * (dcos(bb))**2
           if (cc .LT. ofst) then
             v1 = cc / ofst
             v2 = 0.0D+00
             v3 = 0.0D+00
           else
             v1 = 1.0D+00
             v2 =(cc - ofst) / (1.0D+00 - ofst)
             v3 =(cc - ofst) / (1.0D+00 - ofst)
           endif
! green-ball
           tt = dfloat(nframe) / 50.0D+00
           xx = dfloat(i) / dfloat(ihpixf) - 0.40D+00
           yy = dfloat(j) / dfloat(jvpixf) - 0.65D+00
           rr = dsqrt(xx**2 + yy**2 + 1.0D-30)
           aa = rr / 0.40D+00
           bb =(tt - rr) * pi
           cc = dexp(-aa**2) * (dcos(bb))**2
           if (cc .LT. ofst) then ! here and hereafter, additive rgb color is simply added.
             v2 = v2 + cc / ofst
           else
             v1 = v1 +(cc - ofst) / (1.0D+00 - ofst)
             v2 = v2 + 1.0D+00
             v3 = v3 +(cc - ofst) / (1.0D+00 - ofst)
           endif
! blue-ball
           tt = dfloat(nframe) / 12.5D+00
           xx = dfloat(i) / dfloat(ihpixf) - 0.75D+00
           yy = dfloat(j) / dfloat(jvpixf) - 0.70D+00
           rr = dsqrt(xx**2 + yy**2 + 1.0D-30)
           aa = rr / 0.30D+00
           bb =(tt - rr) * pi
           cc = dexp(-aa**2) * (dcos(bb))**2
           if (cc .LT. ofst) then
             v3 = v3 + cc / ofst
           else
             v1 = v1 +(cc - ofst) / (1.0D+00 - ofst)
             v2 = v2 +(cc - ofst) / (1.0D+00 - ofst)
             v3 = v3 + 1.0D+00
           endif
! yellow-ball
           tt = dfloat(nframe) / 16.66666666666D+00
           xx = dfloat(i) / dfloat(ihpixf) - 0.75D+00
           yy = dfloat(j) / dfloat(jvpixf) - 0.30D+00
           rr = dsqrt(xx**2 + yy**2 + 1.0D-30)
           aa = rr / 0.25D+00
           bb =(tt - rr) * pi
           cc = dexp(-aa**2) * (dcos(bb))**2
           if (cc .LT. ofst) then
             v1 = v1 + cc / ofst
             v2 = v2 + cc / ofst
           else
             v1 = v1 + 1.0D+00
             v2 = v2 + 1.0D+00
             v3 = v3 +(cc - ofst) / (1.0D+00 - ofst)
           endif
         endif
         dat(1,i,j) = v1
         dat(2,i,j) = v2
         dat(3,i,j) = v3
 100   continue

       return
       end subroutine mkdata

! --------------------------------------------
!   fill rgb data array with something
! --------------------------------------------
       subroutine mkbitmap(dat,rgb)
       implicit none
       integer ihpixf, jvpixf
       parameter(ihpixf = nex, jvpixf = ney) ! pixel and data size
       character*1 rgb(3,ihpixf,jvpixf) !      RGB pixel data array
       real*8      dat(3,ihpixf,jvpixf)
! local
       integer i, j
       real*8  red, gre, blu
       integer ired, igre, iblu

       do 100 j = 1, jvpixf
       do 100 i = 1, ihpixf
         red = dat(1,i,j)
         gre = dat(2,i,j)
         blu = dat(3,i,j)
         ired = int(red * 255.0D+00)
         igre = int(gre * 255.0D+00)
         iblu = int(blu * 255.0D+00)
         if (ired .GT. 255) ired = 255
         if (igre .GT. 255) igre = 255
         if (iblu .GT. 255) iblu = 255
         if (ired .LT.   0) ired =   0
         if (igre .LT.   0) igre =   0
         if (iblu .LT.   0) iblu =   0
         rgb(1,i,j) = char(ired)
         rgb(2,i,j) = char(igre)
         rgb(3,i,j) = char(iblu)
100    continue

!!
! Make a white dot that will appear at a point,
!             1/3 horizontal size (width) from the leftmost
!         and 1/4 vertical (height) from the bottom.
!!
!       rgb(1,ihpixf/3,jvpixf/4) = char(255)
!       rgb(2,ihpixf/3,jvpixf/4) = char(255)
!       rgb(3,ihpixf/3,jvpixf/4) = char(255)
!!
       return
       end subroutine mkbitmap


! --------------------------------------------
! end of this file, thank you.
! --------------------------------------------	
	
END Module pixelfrt

!==============================================================

MODULE math_opts

    implicit none

CONTAINS

    subroutine meshgrid_3D(xgv, ygv, zgv, X, Y, Z)

        use mod_geometry, only: nex,ney
        implicit none
        real(8),intent(in)   :: xgv(:), ygv(:), zgv(:)
        real(8),intent(out)  :: X(:,:,:), Y(:,:,:), Z(:,:,:)
        integer           :: sX, sY, sZ, i

        sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)

        do i=1,sZ
            X(:,:,i) = spread( xgv, 1, sY )
            Y(:,:,i) = spread( ygv, 2, sX )
        enddo ! i
        do i=1,sX
            Z(i,:,:) = spread( zgv, 1, sY)
        enddo ! i
    end subroutine meshgrid_3D

       !!!
       !!!
       !!!
       !!!

    subroutine meshgrid_2D(xgv, ygv, arrayx, arrayy)

        use mod_geometry, only: nex,ney
        IMPLICIT NONE
        real(8),intent(in)   :: xgv(:), ygv(:)
        real(8),intent(out)  :: arrayx(:,:), arrayy(:,:)!, Z(:,:,:)
        integer              :: sX, sY, sZ, i

        sX = size(xgv) ; sY = size(ygv) ;        !sZ = size(zgv)

        do i=1,sX
            arrayx = spread( xgv, 2, sY );
            arrayy = spread( ygv, 1, sX );
        enddo ! i

    !        do i=1,sX(i)
    !            Z(i,:,:) = spread( zgv, 1, sY)
    !        enddo ! i

    end subroutine meshgrid_2D

       !!!
       !!!
       !!!
       !!!

    REAL(8) FUNCTION FindDet(matrix, n)
        IMPLICIT NONE
        REAL(8), DIMENSION(n,n) :: matrix
        INTEGER, INTENT(IN) :: n
        REAL(8) :: m, temp
        INTEGER :: i, j, k, l
        LOGICAL :: DetExists = .TRUE.
        l = 1
        !Convert to upper triangular form
        DO k = 1, n-1
            IF (matrix(k,k) == 0) THEN
                DetExists = .FALSE.
                DO i = k+1, n
                    IF (matrix(i,k) /= 0) THEN
                        DO j = 1, n
                            temp = matrix(i,j)
                            matrix(i,j)= matrix(k,j)
                            matrix(k,j) = temp
                        END DO
                        DetExists = .TRUE.
                        l=-l
                        EXIT
                    ENDIF
                END DO
                IF (DetExists .EQV. .FALSE.) THEN
                    FindDet = 0
                    return
                END IF
            ENDIF
            DO j = k+1, n
                m = matrix(j,k)/matrix(k,k)
                DO i = k+1, n
                    matrix(j,i) = matrix(j,i) - m*matrix(k,i)
                END DO
            END DO
        END DO

        !Calculate determinant by finding product of diagonal elements
        FindDet = l
        DO i = 1, n
            FindDet = FindDet * matrix(i,i)
        END DO

    END FUNCTION FindDet

        !!!
        !!!
        !!!
        !!!

    FUNCTION cross(a1, b1)

        IMPLICIT NONE
        REAL(8), DIMENSION(3) :: cross
        REAL(8), DIMENSION(3), INTENT(IN):: a1, b1

        cross(1) = a1(2) * b1(3) - a1(3) * b1(2)
        cross(2) = a1(3) * b1(1) - a1(1) * b1(3)
        cross(3) = a1(1) * b1(2) - a1(2) * b1(1)

    END FUNCTION cross

        !!!
        !!!
        !!!
        !!!

    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    !FUNCTION inverse(A,n) result(Ainv)
    function inverse(A) result(Ainv)

        INCLUDE 'mkl_rci.fi'
        !INCLUDE 'lapack.f90'

        real(8), dimension(:,:), intent(in) :: A
        real(8), dimension(size(A,1),size(A,2)) :: Ainv

        real(8), dimension(size(A,1))  :: work   ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n,info

        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv(:,:) = A(:,:)
        n = size(A,1)

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)

        if (info /= 0) then
            stop 'Matrix inversion failed!'
        end if

    end function inverse

            !!!
            !!!
            !!!
            !!!

    REAL(8) function mean(array_in,xLength,yLength)

        integer :: xLength,yLength
        REAL(8), DIMENSION(xLength,yLength) :: array_in

        mean = SUM(SUM(array_in,DIM=1),DIM=1)/(size(array_in,1)*size(array_in,2))

    end function



    !!!
    !!!
    !!!
    !!!

    !! ! NOT IN PLACE:

    subroutine twoD_fft( in_out_2D , FFT_type ) ! NOT IN PLACE

        ! Fortran example.
        ! 2D complex to complex, and real to conjugate-even
        USE mod_geometry
        Use MKL_DFTI

        CHARACTER(LEN=8) :: FFT_type

        Complex(8) , DIMENSION(ig,jg)   :: in_out_2D
        Complex(8) , DIMENSION(ig,jg)   :: in_out_temp_2D
        Complex(8) , DIMENSION(1:ig*jg) :: in_out_temp_1D
        Equivalence (in_out_temp_2D, in_out_temp_1D)

        type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
        Integer :: Status, L(2)
        Real(8) :: fscale = 1.0D0 !1.0D0/real(ig)
        Real(8) :: bscale = 1.0D0/real(ig*jg)

        !!
        integer nth, len(2)
        ! 4 OMP threads, each does 2D FFT 50x100 points
        parameter (nth = 4, len = (/ig, jg/))
        complex(8)  :: x(len(2)*len(1), nth)
        type(dfti_descriptor), pointer :: myFFT
        integer th, myStatus

        in_out_temp_2D = in_out_2D

        !*****
        L(1) = ig
        L(2) = jg

        SELECT CASE (FFT_type)

            CASE ("Forward ")

                !! FORWARD C - C (OPENMP)
                !! FORWARD C - C
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, L)
                if (status .ne. 0) then
                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                        print *, 'Error: ', DftiErrorMessage(status)
                    endif
                endif
                !Status = DftiSetValue(My_Desc1_Handle,DFTI_PLACEMENT, DFTI_NOT_INPLACE )
                Status = DftiCommitDescriptor( My_Desc1_Handle)
                if (status .ne. 0) then
                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                        print *, 'Error: ', DftiErrorMessage(status)
                    endif
                endif
                Status = DftiComputeForward( My_Desc1_Handle, in_out_temp_1D)
                if (status .ne. 0) then
                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                        print *, 'Error: ', DftiErrorMessage(status)
                    endif
                endif
                Status = DftiFreeDescriptor( My_Desc1_Handle)
                if (status .ne. 0) then
                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                        print *, 'Error: ', DftiErrorMessage(status)
                    endif
                endif

                in_out_2D = in_out_temp_2D*fscale

            CASE ("Backward")

                !! BACKWARD C - C
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                Status = DftiCreateDescriptor( My_Desc2_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, L)
                if (status .ne. 0) then
                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                        print *, 'Error: ', DftiErrorMessage(status)
                    endif
                endif
                        !Status = DftiSetValue(My_Desc2_Handle,DFTI_PLACEMENT, DFTI_NOT_INPLACE )
                Status = DftiCommitDescriptor( My_Desc2_Handle)
                Status = DftiComputeBackward( My_Desc2_Handle, in_out_temp_1D )
                if (status .ne. 0) then
                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                        print *, 'Error: ', DftiErrorMessage(status)
                    endif
                endif
                Status = DftiFreeDescriptor( My_Desc2_Handle)
                if (status .ne. 0) then
                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                        print *, 'Error: ', DftiErrorMessage(status)
                    endif
                endif

                in_out_2D = in_out_temp_2D*bscale

            CASE DEFAULT

                WRITE(*,*) ' WRONG FFT TYPE. Please check the FFT type ...'
                STOP

        END SELECT

    end subroutine twoD_fft

    !! ! NOT IN PLACE:

!    subroutine threeD_fft( in_out_3D , FFT_type ) ! NOT IN PLACE
!
!        ! Fortran example.
!        ! 2D complex to complex, and real to conjugate-even
!        USE mod_geometry
!        Use MKL_DFTI
!
!        CHARACTER(LEN=8) :: FFT_type
!
!        Complex(8) , DIMENSION(ig,jg,zg)   :: in_out_3D
!        Complex(8) , DIMENSION(ig,jg,zg)   :: in_out_temp_3D
!        Complex(8) , DIMENSION(1:ig*jg*zg) :: in_out_temp_1D
!        Equivalence (in_out_temp_3D, in_out_temp_1D)
!
!        type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
!        Integer :: Status, L(3)
!        Real(8) :: fscale = 1.0D0 !1.0D0/real(ig)
!        Real(8) :: bscale = 1.0D0/real(ig*jg*zg)
!
!        !! Code begins here .......
!
!        in_out_temp_3D = in_out_3D
!
!        !*****
!        L(1) = ig
!        L(2) = jg
!        L(3) = zg
!
!        SELECT CASE (FFT_type)
!
!            CASE ("Forward ")
!
!                !! FORWARD C - C (OPENMP)
!                !! FORWARD C - C
!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, L)
!                if (status .ne. 0) then
!                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
!                        print *, 'Error: ', DftiErrorMessage(status)
!                    endif
!                endif
!                !Status = DftiSetValue(My_Desc1_Handle,DFTI_PLACEMENT, DFTI_NOT_INPLACE )
!                Status = DftiCommitDescriptor( My_Desc1_Handle)
!                if (status .ne. 0) then
!                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
!                        print *, 'Error: ', DftiErrorMessage(status)
!                    endif
!                endif
!                Status = DftiComputeForward( My_Desc1_Handle, in_out_temp_1D)
!                if (status .ne. 0) then
!                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
!                        print *, 'Error: ', DftiErrorMessage(status)
!                    endif
!                endif
!                Status = DftiFreeDescriptor( My_Desc1_Handle)
!                if (status .ne. 0) then
!                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
!                        print *, 'Error: ', DftiErrorMessage(status)
!                    endif
!                endif
!
!                in_out_3D = in_out_temp_3D*fscale
!
!            CASE ("Backward")
!
!                !! BACKWARD C - C
!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                Status = DftiCreateDescriptor( My_Desc2_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, L)
!                if (status .ne. 0) then
!                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
!                        print *, 'Error: ', DftiErrorMessage(status)
!                    endif
!                endif
!                        !Status = DftiSetValue(My_Desc2_Handle,DFTI_PLACEMENT, DFTI_NOT_INPLACE )
!                Status = DftiCommitDescriptor( My_Desc2_Handle)
!                Status = DftiComputeBackward( My_Desc2_Handle, in_out_temp_1D )
!                if (status .ne. 0) then
!                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
!                        print *, 'Error: ', DftiErrorMessage(status)
!                    endif
!                endif
!                Status = DftiFreeDescriptor( My_Desc2_Handle)
!                if (status .ne. 0) then
!                    if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
!                        print *, 'Error: ', DftiErrorMessage(status)
!                    endif
!                endif
!
!                in_out_3D = in_out_temp_3D*bscale
!
!            CASE DEFAULT
!
!                WRITE(*,*) ' WRONG FFT TYPE. Please check the FFT type ...'
!                STOP
!
!        END SELECT
!
!    end subroutine threeD_fft


END MODULE math_opts



!! ---------------------------------
MODULE elastic_potential

    USE mod_geometry
    USE elastic_property
    USE solver_opts
    USE wrt_opts
    USE math_opts

    IMPLICIT NONE

    REAL(kind), DIMENSION(1:6,1:6)             :: C_EFFV   ! RVE

        !! Eigenstrain tensors ...
        REAL(kind), DIMENSION(1:9,1:2)             :: epslon_00
        REAL(kind), DIMENSION(1:3,1:3)             :: eps_00_1, eps_00_2

CONTAINS

    subroutine elastic_pot(C1,C2,mu_el,sigma_el, epslon_el)

        REAL(kind), INTENT(in) , DIMENSION(1:nex,1:ney)      :: C1,C2
        REAL(kind), INTENT(out), DIMENSION(1:nex,1:ney,1:2)  :: mu_el

        !! ****************************** PHASE FIELD CODE ************************************

        REAL(kind), DIMENSION(1:nex,1:ney)                :: f_xyz
        !EQUIVALENCE(phi11,f_xyz)
        !!
        INTEGER :: x,y,z,m,n,o
        INTEGER :: i,j

        !! Reciprocal space arrays
        REAL(kind), DIMENSION(1:3)                 :: g
        REAL(kind), ALLOCATABLE, DIMENSION(:,:,:)  :: K

        !! Elastic modulus tensors ...
        REAL(kind), DIMENSION(1:12)                :: C_vals,C_p1_vals,C_p2_vals,C_m_vals
        REAL(kind), DIMENSION(1:6,1:6)             :: Cv,C_mv,C_pv1,C_pv2,delta_Cv1,delta_Cv2
        REAL(kind), DIMENSION(1:6,1:6)             :: Cpv      ! rotation...
        REAL(kind), DIMENSION(1:6,1:6)             :: Sv       ! stiffness

        REAL(kind), ALLOCATABLE, DIMENSION(:,:,:)  :: sigma_T                   ! Eigenstrain tensor
        REAL(kind), ALLOCATABLE, DIMENSION(:,:,:)  :: epslon_T                  ! Eigenstrain tensor
        REAL(kind), ALLOCATABLE, DIMENSION(:,:)    :: lattice_constant

        REAL(kind), ALLOCATABLE, DIMENSION(:,:,:,:):: G_il

        !! dis_field ...
        REAL(kind), DIMENSION(1:nex,1:ney,1:dimen) :: disp_field,disp_field_temp
        REAL(kind), DIMENSION(1:nex,1:ney,1:9)     :: Jac_disp_field

        REAL(kind), DIMENSION(1:6)                 :: E_ij                      ! Homogeneous strain tensor; mean strain tensor of the simulation cell
        REAL(kind), dimension(1:6)                 :: sigma_a

        !! Strain values
        REAL(kind), ALLOCATABLE, DIMENSION(:,:,:)  :: epslon_st, sigma_0_st, sigma_st, epslon_0, sigma_0

        REAL(kind), dimension(1:6)                 :: sigma_avg,sigma_0_avg

        REAL(kind), DIMENSION(1:nex,1:ney,1:6)     :: epslon_el, sigma_el
        !
        REAL(kind), dimension(1:order)             :: err_x,err_y
        REAL(kind), dimension(1:nex)               :: x1,y1
        
        REAL(kind), DIMENSION(1:nex,1:ney)   	   :: Beta_C1,Beta_C2,Beta_C3
		
		!! Green lagrange tensors
		REAL(kind), DIMENSION(1:nex,1:ney,1:6)     :: e
		REAL(kind), DIMENSION(1:nex,1:ney,1:9) 	   :: Jac_e_green_field
		REAL(kind), DIMENSION(1:nex,1:ney)  	   :: U_tot
		REAL(kind), DIMENSION(1:nex,1:ney)  	   :: e_2D_aspect_ratio
		REAL(kind), DIMENSION(1:nex,1:ney)  	   :: e_3D_aspect_ratio
		
        !! Temp. arrays ...
        REAL(kind) :: phi1MAX,phi1max1,phi1tot1,delta_C
        integer :: II

        !! Code begins here ...
        
        ALLOCATE(K(1:nex,1:ney,1:3));
        ALLOCATE(sigma_T(1:nex,1:ney,1:6));
        ALLOCATE(epslon_T(1:nex,1:ney,1:6));
        ALLOCATE(lattice_constant(1:nex,1:ney));
        ALLOCATE(G_il(1:nex,1:ney,1:dimen,1:dimen)) 
		ALLOCATE(epslon_st(1:nex,1:ney,1:6),sigma_0_st(1:nex,1:ney,1:6),sigma_st(1:nex,1:ney,1:6),epslon_0(1:nex,1:ney,1:6),sigma_0(1:nex,1:ney,1:6))

        !**************************************************************************
        !**************************************************************************
        !********************  Domain config. *************************************
        !**************************************************************************
        !**************************************************************************

        !f_xyz(1:nex,1:ney) = C1(1:nex,1:ney)

        !**************************************************************************
        !**************************************************************************
        !********************  Reciprocal space  **********************************
        !**************************************************************************
        !**************************************************************************

        call step0_reciprocal_lat(g ,K);


        !    !**************************************************************************
        !    !**************************************************************************
        !    !*********************          START             *************************
        !    !**************************************************************************
        !    !**************************************************************************

        !**************************************************************************
        !******  Interpolation functions for elastic modulus/eigenstrain  *********
        !**************************************************************************

        Beta_C1        = (C1**3)*(10.0D0 - 15.0*C1 + 6.0D0*(C1**2));
        Alpha_C1       = (C1**3)*(10.0D0 - 15.0*C1 + 6.0D0*(C1**2)) - 1.0D0/3.0D0;
        Beta_C1_prime  = (3.0D0*(C1)**2)*(10.0D0 - 15.0D0*C1 + 6.0D0*((C1)**2))+((C1)**3)*(12.0D0*C1-15.0D0);        !....Interpolation function for eigenstrain, a scalar function
        Alpha_C1_prime = (3.0D0*(C1)**2)*(10.0D0 - 15.0D0*C1 + 6.0D0*((C1)**2))+((C1)**3)*(12.0D0*C1-15.0D0);        !...Interpolation function for elastic modulus, a scalar function
        Alpha_Beta_C1  = Alpha_C1*Beta_C1;

        Beta_C2        = (C2**3)*(10.0D0 - 15.0*C2 + 6.0D0*(C2**2));
        Alpha_C2       = (C2**3)*(10.0D0 - 15.0*C2 + 6.0D0*(C2**2)) - 1.0D0/3.0D0;
        Beta_C2_prime  = (3.0D0*(C2)**2)*(10.0D0 - 15.0D0*C2 + 6.0D0*((C2)**2))+((C2)**3)*(12.0D0*C2-15.0D0);        !....Interpolation function for eigenstrain, a scalar function
        Alpha_C2_prime = (3.0D0*(C2)**2)*(10.0D0 - 15.0D0*C2 + 6.0D0*((C2)**2))+((C2)**3)*(12.0D0*C2-15.0D0);        !...Interpolation function for elastic modulus, a scalar function
        Alpha_Beta_C2  = Alpha_C2*Beta_C2;

        !Beta_C3        = ((1.0-C1-C2)**3)*(10.0D0 - 15.0*(1.0-C1-C2) + 6.0D0*((1.0-C1-C2)**2));

        !**************************************************************************
        !**************************************************************************
        !********************  Eigenstrain Tensor  ********************************
        !**************************************************************************
        !**************************************************************************
        !!

        epslon_00(:,1) = [ +0.1D0, -0.1D0 , 0.0D0, 0.0D0 ,0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 ]
        epslon_00(:,2) = [ -0.1D0, +0.1D0 , 0.0D0, 0.0D0 ,0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 ]

        call eg_st_tensor(eps_00_1  ,epslon_00(:,1) );          !.......
        call eg_st_tensor(eps_00_2  ,epslon_00(:,2) );          !.......
        
!        write(*,*) eps_00_1
!        write(*,*) eps_00_2
!        pause


        !**************************************************************************
        !**************************************************************************
        !********************  Elastic Modulus Tensor  ****************************
        !**************************************************************************
        !**************************************************************************
        !!

        C_vals    = (/ C11,C22,C11,C12,C12,C12,C12,C12,C12,2*C44,2*C44,2*C44 /);
        C_p1_vals = (/ C11_p1,C22_p1,C11_p1,C12_p1,C12_p1,C12_p1,C12_p1,C12_p1,C12_p1,C44_p1,C44_p1,C44_p1 /); ! AlN
        C_p2_vals = (/ C11_p2,C22_p2,C11_p2,C12_p2,C12_p2,C12_p2,C12_p2,C12_p2,C12_p2,C44_p2,C44_p2,C44_p2 /); ! ZrN
        C_m_vals  = (/ C11_m,C22_m,C11_m,C12_m,C12_m,C12_m,C12_m,C12_m,C12_m,C44_m,C44_m,C44_m /);

        call elmod_tensor(Cv  ,C_vals);            !.........................Elastic modulus tensor; composition (and hence, position) dependent
        call elmod_tensor(C_mv ,C_m_vals);         !.........................Elastic modulus tensor of the m phase
        call elmod_tensor(C_pv1,C_p1_vals);        !.........................Elastic modulus tensor of the p phase
        call elmod_tensor(C_pv2,C_p2_vals);        !.........................Elastic modulus tensor of the p phase

        C_effv   = (C_pv1 + C_mv)/2.0D0;   				!.................................

!        delta_C = maxval(C1) - X11;

!        IF(abs(delta_C) < 0.15) THEN

!            IF(itimes==0) THEN
!               OPEN(unit=50,FILE=TRIM(fileplace)//'/elasticity/0_Delta_C.dat')
!            ENDIF
!			C_effv    = Cv;										    !.................................
!			delta_Cv1 = 0.0!C_pv1 - C_mv                 !.................................
!			delta_Cv2 = 0.0!C_pv2 - C_mv                    !................................
!            !call rotation(Cpv, C_pv);
!            delta_C  = MAXVAL(C1) - X11
!            epslon_T = Beta_C1*0.001
!            write(50,*) itimes,delta_C,maxval(epslon_T)
  
!        ELSE
!            C_effv   = (C_pv1 + C_pv2 + C_mv)/3.0D0;   				!.................................
!			delta_Cv1= C_pv1 - C_mv  !C_mv- C_pv1;                  !.................................
!			delta_Cv2= C_pv2 - C_mv; !C_mv- C_pv2                   !.................................
!			!call rotation(Cpv, C_pv);
!			call step2_calc_epslon_T(epslon_T,lattice_constant,Beta_C1,Beta_C2,Beta_C3)
!        ENDIF
        
        !**************************************************************************
        !**************************************************************************
        !***************************  Rotation  ***********************************
        !**************************************************************************
        !**************************************************************************

        !call rotation(Cpv, C_pv1);

        C_c = 0.0D0;
        DO i = 1,dimen
            DO j = i,dimen
                DO m = 1,dimen
                    DO n = m,dimen
                        C_c(:,:,voigt(i,j),voigt(m,n)) = C_effv(voigt(i,j),voigt(m,n)) + Alpha_C1*delta_Cv1(voigt(i,j),voigt(m,n)) + Alpha_C2*delta_Cv2(voigt(i,j),voigt(m,n));
                    enddo
                enddo
            enddo
        enddo


        !A_lit = 1;
        !f_0 = ((A_lit**2)/4)*(f_xyz**2)*(1 - f_xyz)**2;
        !h_f = ((A_lit**2)/2)*(f_xyz)*(1 - f_xyz)**2 - ((A_lit**2)/2)*(f_xyz**2)*(1 - f_xyz);
        !F_s = fft2(f_xyz);!...................................................................................FFT of f(x,y) with Fourier nodes
        !F_s_2prime = (((1j).*(g(1).*K(:,:,1) + g(2).*K(:,:,2))).^2).*F_s;!....................................Fourier space 2nd derivative of f(x,y)
        !Laplacian_c = real(ifft2(F_s_2prime));!...............................................................Laplacian of the phase field f(x,y)
        !mu_ch = h_f -2.*Kappa.*Laplacian_c;!..................................................................Chemical potential contribution of the phase field

        !=========================================================================================
        !=========================================================================================
        !=========================================================================================
        !=========================Mechanical Equilibrium Algorithm================================
        !=========================================================================================
        !=========================================================================================
        !=========================================================================================

        !WRITE(*,*) ' *** STARTING MECHANICAL EQ. calculations ...'

        !=========================================================================================
        !=================================Acoustic tensor=========================================
        !=========================================================================================
        call step1_calc_G_il(G_il, K , g);
        !WRITE(*,*)
        !WRITE(*,*) ' *** Acoustic tensor (G_il) (inverse of C_ijkl*K_j*K_i) is calculated ...'
        !WRITE(*,*)
        
        !=========================================================================================
        !=================================Node Analysis===========================================
        !=========================================================================================
        call step2_calc_epslon_T(epslon_T,Beta_C1,Beta_C2,Beta_C3)
        !WRITE(*,*)
        !WRITE(*,*) ' *** Eigenstrain tensor is calculated ...'
        !WRITE(*,*)

        !=========================================================================================
        !=================================Eigenstrain tensor======================================
        !=========================================================================================
        call step3_calc_sigma_T(sigma_T,epslon_T,C1)
        !WRITE(*,*)
        !WRITE(*,*) ' *** Eigenstress tensor is calculated ...'
        !WRITE(*,*)
        
        !=========================================================================================
        !=================================Acoustic tensor=========================================
        !=========================================================================================
        call step4_calc_S(Sv, C_effv, delta_Cv1, delta_Cv2);
        !WRITE(*,*)
        !WRITE(*,*) ' *** Sv tensor (Sv) (Compliance tensor) is calculated ...'
        !WRITE(*,*)

        !**************************************************************************
        !**************************************************************************
        !*********************          OUTPUT            *************************
        !**************************************************************************
        !**************************************************************************

        !**************************************************************************
        !**************************************************************************
        !*********************          OUTPUT            *************************
        !**************************************************************************
        !**************************************************************************

        IF(itimes == 0) THEN

            !***** OPEN DISK *****
            call wrt_disk
            !
            OPEN(unit= 0,FILE=TRIM(fileplace)//'/elasticity/0_KxKy.dat')
            OPEN(unit=16,FILE=TRIM(fileplace)//'/elasticity/9_error.dat')

            WRITE(0,770) IG, JG
            DO i=1,IG; DO j=1,JG;
                WRITE(0 ,'(I3,2X,I3,3X,9F17.7)') i,j, K(i,j,1), K(i,j,2), K(i,j,3)
            ENDDO; ENDDO

        ENDIF

        IF(itimes == 0) THEN

            OPEN(unit=1,FILE=TRIM(fileplace)//'/elasticity/strain/1_epslon_T.dat')
551         FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"f_xyz1"',2X,'"f_xyz2"',2X,'"lattice_constant"',2X,'"epslon_T_1"',2X,'"epslon_T_2"',2X,'"epslon_T_6"')
            WRITE(1,551)
            WRITE(1,770) IG, JG
            DO i=1,IG; DO j=1,JG;
                WRITE(1 ,'(I3,2X,I3,2X,9E17.7)') i,j, C1(i,j), C2(i,j), lattice_constant(i,j), epslon_T(i,j,1), epslon_T(i,j,2), epslon_T(i,j,6)
            ENDDO; ENDDO

        ENDIF

        !=========================================================================================
        !================================= Solver        =========================================
        !=========================================================================================

        disp_field     = 0.0D0
        Jac_disp_field = 0.0D0
        disp_field_temp= 0.0D0
        E_ij           = 0.0D0
        sigma_a        = 0.0D0

        DO o = 1,order

            err_x(o) = 0.0D0
            err_y(o) = 0.0D0

            sigma_el   = 0.0D0
            sigma_el   = 0.0D0

            sigma_0_st = 0.0D0
            sigma_st   = 0.0D0
            epslon_0   = 0.0D0
            epslon_st  = 0.0D0

            call step5_calc_disp_field(disp_field, Jac_disp_field, C1, C2, sigma_T, epslon_T, G_il, Sv, Beta_C1, Beta_C2, delta_Cv1, delta_Cv2, E_ij, K, g);
            call step6_calc_Jac       (Jac_disp_field, disp_field, K, g);
            call step7_calc_stars     (epslon_st, sigma_st, sigma_0_st, epslon_0, Jac_disp_field, Beta_C1, Beta_C2, C_effv, delta_Cv1, delta_Cv2, epslon_T);
            call step8_calc_E_ij      (E_ij, sigma_0_st, sigma_st, sigma_a, Sv);
            call step9_calc_el        (epslon_el, sigma_el, E_ij, epslon_st, epslon_0, C_effv, delta_Cv1, delta_Cv2);
            !call step10_calc_mu_el    (mu_el_1,mu_el_2, C1, C2, epslon_0, epslon_T, epslon_st, E_ij, epslon_el, sigma_el, delta_Cv1, delta_Cv2)

            mu_el(:,:,1) = (sigma_el(:,:,1)*eps_00_1(1,1)+sigma_el(:,:,2)*eps_00_1(2,2) )/VM;
            mu_el(:,:,2) = (sigma_el(:,:,1)*eps_00_2(1,1)+sigma_el(:,:,2)*eps_00_2(2,2) )/VM;

            DO i=1,IG;DO j=1,JG
                err_x(o) = err_x(o) + ((disp_field(i,j,1)-disp_field_temp(i,j,1))*(disp_field(i,j,1)-disp_field_temp(i,j,1)))**0.5  ;
                err_y(o) = err_y(o) + ((disp_field(i,j,2)-disp_field_temp(i,j,2))*(disp_field(i,j,2)-disp_field_temp(i,j,2)))**0.5  ;
            ENDDO;ENDDO
            disp_field_temp = disp_field; ! disp_field_temp: n -- disp_field: n+1

!            WRITE(*,*) 'order of approx.:',o
!            WRITE(*,*) 'E_ij:'
!            WRITE(*,'(6D15.6)') E_ij
!            WRITE(*,*) 'err is :',err_x(o),err_y(o)
!            WRITE(*,*) 'max S11:',maxval(sigma_el(:,:,1)),'min S11:',minval(sigma_el(:,:,1))
!            WRITE(*,*) 'max S22:',maxval(sigma_el(:,:,2)),'min S22:',minval(sigma_el(:,:,2))

!            WRITE(13 ,'(I9,2X,7E18.7)') o, err_x(o), err_y(o)

            if ( err_y(o)  < elas_err) then
                !write(*,*) 'error tol. is satisfied...'
                exit
            endif

        enddo
        
        e(:,:,1) = (epslon_el(:,:,1) + epslon_el(:,:,2) + epslon_el(:,:,3))/sqrt(3.0);
        e(:,:,2) = (epslon_el(:,:,1) - epslon_el(:,:,2))/sqrt(2.0);
        e(:,:,3) = (epslon_el(:,:,1) + epslon_el(:,:,2) - 2.0D0*epslon_el(:,:,3))/sqrt(6.0);
        e(:,:,4) =  sqrt(2.0)*epslon_el(:,:,4) 
        e(:,:,5) =  sqrt(2.0)*epslon_el(:,:,5) 
        e(:,:,6) =  sqrt(2.0)*epslon_el(:,:,6)         

        e_2D_aspect_ratio = e(:,:,6)/e(:,:,2)
        e_3D_aspect_ratio = e(:,:,3)/e(:,:,2)
           
		call step12_calc_Jac(U_tot, Jac_e_green_field, e, K, g);

        !**************************************************************************
        !**************************************************************************
        !*********************     WRITE OUTPUT FILES      ************************
        !**************************************************************************
        !**************************************************************************

        IF (itimes == 0) then

            OPEN(unit=100,FILE=TRIM(fileplace)//'/elasticity/Elastic_modulus.dat')

            write(100,*) 'G_m  =',G_m,'G_p1  =',C44_p1,'G_p2  =',C44_p2
            write(100,*) 'delta1=',delta1
            write(100,*) 'delta2=',delta2
            WRITE(100,*)
            WRITE(100,*) 'A_z =',A_z
            WRITE(100,*) 'Az_p1 =',A_zp1
            WRITE(100,*) 'Az_p2 =',A_zp2
            WRITE(100,*) 'Az_m =',A_zm
            WRITE(100,*) 'Ap_p1=',Ap_p1
            WRITE(100,*) 'Ap_p2=',Ap_p2
            WRITE(100,*) 'Ap_m =',Ap_m
            WRITE(100,*) 'beta1=',C44_p1/C44_m
            WRITE(100,*) 'beta2=',C44_p2/C44_m
            WRITE(100,*) 'B_m  =',B_m
            WRITE(100,*) 'B_p1 =',B_p1
            WRITE(100,*) 'B_p1 =',B_p2
            WRITE(100,*) 'B_v  =',B_v
            WRITE(100,*)
            WRITE(100,*) '******* Elastic tensors *************************'
            write(100,*) 'C_matrix='
            write(100,'(2X,6E14.3)') C_mv
            WRITE(100,*)
            write(100,*) 'C_precipitate1='
            write(100,'(2X,6E14.3)') C_pv1
            WRITE(100,*)
            write(100,*) 'C_precipitate2='
            write(100,'(2X,6E14.3)') C_pv2
            WRITE(100,*)
            write(100,*) 'Cv='
            write(100,'(2X,6E14.3)') Cv
            WRITE(100,*)
            WRITE(100,*) '******* Effective E tensors *************************'
            write(100,*) 'C_Effective='
            write(100,'(2X,6E14.3)') C_effv
            WRITE(100,*)
            WRITE(100,*)
            WRITE(100,*) '******* delta_Cv tensor *****************************'
            write(100,*) 'delta_Cv1='
            write(100,'(2X,6E14.3)') delta_Cv1
            WRITE(100,*)
            write(100,*) 'delta_Cv2='
            write(100,'(2X,6E14.3)') delta_Cv2
            WRITE(100,*)
            WRITE(100,*) '******* Cv tensor in the center of domian *****************************'
            write(100,*) 'Cv='
            write(100,'(2X,6E14.3)') Alpha_C1(nex/2,ney/2)*delta_Cv1 + Alpha_C2(nex/2,ney/2)*delta_Cv2 + C_mv
            WRITE(100,*)
            WRITE(100,*)
            WRITE(100,*) '******* Cv tensor in the corner of domian *****************************'
            write(100,*) 'Cv='
            write(100,'(2X,6E14.3)') Alpha_C1(nex,ney)*delta_Cv1 + Alpha_C2(nex,ney)*delta_Cv2 + C_mv
            WRITE(100,*)
            CLOSE(UNIT=100)
        ENDIF

        !**************************************************************************
        !**************************************************************************
        !*********************     WRITE OUTPUT FILES      ************************
        !**************************************************************************
        !**************************************************************************
        !**************************************************************************
        !**************************************************************************
        !*********************     WRITE OUTPUT FILES      ************************
        !**************************************************************************
        !**************************************************************************

!        x1 = [ ( real(i), i=(-int(nex/2)),-1 ) , 0.0 ,  ( real(i), i=1,int(nex/2)) ]
!        y1 = [ ( real(i), i=(-int(nex/2)),-1 ) , 0.0 ,  ( real(i), i=1,int(nex/2)) ]

        IF (MOD(itimes,wrt_cycle) == 0) THEN

            call wrt_disk
            !! *********** WRITTING OUTPUT ***************
            WRITE(num,'(I9.9)') itimes
            OPEN(unit=2,FILE=TRIM(fileplace)//'/elasticity/stress/1_sigma_T_'//TRIM(num)//'.dat')
            OPEN(unit=3,FILE=TRIM(fileplace)//'/elasticity/0_interpolation_'//TRIM(num)//'.dat')
            OPEN(unit=4,FILE=TRIM(fileplace)//'/elasticity/displacement/3_disp_field_'//TRIM(num)//'.dat')
            OPEN(unit=9,FILE=TRIM(fileplace)//'/elasticity/strain/6_epslon_0_'//TRIM(num)//'.dat')
            OPEN(unit=10,FILE=TRIM(fileplace)//'/elasticity/strain/8_epslon_el_'//TRIM(num)//'.dat')
            OPEN(unit=11,FILE=TRIM(fileplace)//'/elasticity/stress/8_sigma_el_'//TRIM(num)//'.dat')
            !OPEN(unit=12,FILE=TRIM(fileplace)//'/elasticity/7_epslon_tot_'//TRIM(num)//'.dat')
            OPEN(unit=13,FILE=TRIM(fileplace)//'/elasticity/strain/9_e_green_'//TRIM(num)//'.dat')
            OPEN(unit=14,FILE=TRIM(fileplace)//'/elasticity/strain/9_jac_e_green_'//TRIM(num)//'.dat')
            OPEN(unit=15,FILE=TRIM(fileplace)//'/elasticity/strain/9_U_tot_'//TRIM(num)//'.dat')
            OPEN(unit=16,FILE=TRIM(fileplace)//'/elasticity/9_error.dat')

            !OPEN(unit=2,FILE=TRIM(fileplace)//'/elasticity/1_sigma_T.dat')
            !OPEN(unit=4,FILE=TRIM(fileplace)//'/elasticity/3_disp_field.dat')
            !OPEN(unit=5,FILE=TRIM(fileplace)//'/elasticity/displacement/4_Jac_disp_field_'//TRIM(num)//'.dat')
            !OPEN(unit=6,FILE=TRIM(fileplace)//'/elasticity/5_epslon_st.dat')
            !OPEN(unit=7,FILE=TRIM(fileplace)//'/elasticity/5_sigma_st.dat')
            !OPEN(unit=8,FILE=TRIM(fileplace)//'/elasticity/6_sigma_0_st.dat')
            !OPEN(unit=9,FILE=TRIM(fileplace)//'/elasticity/6_epslon_0.dat')
            !OPEN(unit=12,FILE=TRIM(fileplace)//'/elasticity/7_epslon_tot.dat')
            !OPEN(unit=11,FILE=TRIM(fileplace)//'/elasticity/8_sigma_el.dat')

            !OPEN(unit=12,FILE=TRIM(fileplace)//'/elasticity/f_el_alpha.dat')

553         FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"C1"',2X,'"C2"',2X,'"Beta_C1"',2X,'"Beta_C2"',2X,'"Beta_C3"',2X,'"Alpha_C1"',2X,'"Alpha_C2"',2X, & 
									'"Alpha_C1*delta_Cv1"',2X,'"Alpha_C2*delta_Cv2"',2X,'"C_c"',2X)
            WRITE(3,553)             
554         FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"C1"',2X,'"C2"',2X,'"e11"',2X,'"e22"',2X,'"e12"',2X)
            WRITE(10,554)
555         FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"C1"',2X,'"C2"',2X,'"s11"',2X,'"s22"',2X,'"s12"',2X)
            WRITE(11,555)

            WRITE(2,770) IG, JG
            WRITE(3,770) IG, JG
            WRITE(4,770) IG, JG
            !WRITE(5,770) IG, JG
            !WRITE(6,770) IG, JG
            !WRITE(7,770) IG, JG
            !WRITE(8,770) IG, JG
            WRITE(9,770) IG, JG
            WRITE(10,770) IG, JG
            WRITE(11,770) IG, JG
            !WRITE(12,770) IG, JG
            !WRITE(13,770) IG, JG
            WRITE(14,770) IG, JG
            WRITE(15,770) IG, JG

            DO i=1,IG; DO j=1,JG;

                WRITE(2 ,'(I3,2X,I3,2X,10E17.7)') i,j, C1(i,j), C2(i,j), sigma_T(i,j,1),sigma_T(i,j,2),sigma_T(i,j,6)
                WRITE(3 ,'(I3,2X,I3,2X,10E17.7)') i,j, C1(i,j), C2(i,j), Beta_C1(i,j),Beta_C2(i,j),Beta_C3(i,j), Alpha_C1(i,j), Alpha_C2(i,j), & 
					C_effv(1,1)+Alpha_C1(i,j)*delta_Cv1(1,1), C_effv(1,1)+Alpha_C2(i,j)*delta_Cv2(1,1), C_c(i,j,1,1)
                
                !WRITE(3 ,'(I3,2X,I3,2X,9F12.5)') i,j, phi1MAX, f_xyz(i,j),G_il(i,j,1,1),G_il(i,j,1,2),G_il(i,j,2,1),G_il(i,j,2,2)
                WRITE(4 ,'(I3,2X,I3,2X,9E17.7)') i,j, C1(i,j),C2(i,j),disp_field(i,j,1),disp_field(i,j,2),disp_field(i,j,3)
                !WRITE(5 ,'(I3,2X,I3,2X,9E17.7)') i,j, C1(i,j),Jac_disp_field(i,j,1),Jac_disp_field(i,j,2),Jac_disp_field(i,j,4),Jac_disp_field(i,j,5),Jac_disp_field(i,j,6)
                !WRITE(6 ,'(I3,2X,I3,2X,9E17.7)') i,j, f_xyz(i,j),epslon_st(i,j,1), epslon_st(i,j,2), epslon_st(i,j,3), epslon_st(i,j,4), epslon_st(i,j,5), epslon_st(i,j,6)
                !WRITE(7 ,'(I3,2X,I3,2X,9E17.7)') i,j, f_xyz(i,j),sigma_st(i,j,1), sigma_st(i,j,2), sigma_st(i,j,3), sigma_st(i,j,4), sigma_st(i,j,5), sigma_st(i,j,6)
                !WRITE(8 ,'(I3,2X,I3,2X,9E17.7)') i,j, f_xyz(i,j),sigma_0_st(i,j,1), sigma_0_st(i,j,2), sigma_0_st(i,j,3), sigma_0_st(i,j,4), sigma_0_st(i,j,5), sigma_0_st(i,j,6)
                WRITE(9 ,'(I4,2X,I4,2X,9E17.7)') i,j, C1(i,j),C2(i,j),epslon_0(i,j,1) , epslon_0(i,j,2), epslon_0(i,j,6)
                WRITE(10,'(I4,2X,I4,2X,9E17.7)') i,j, C1(i,j),C2(i,j),epslon_el(i,j,1), epslon_el(i,j,2), epslon_el(i,j,6)
                WRITE(11,'(I4,2X,I4,2X,9E17.7)') i,j, C1(i,j),C2(i,j),sigma_el(i,j,1) ,sigma_el(i,j,2),sigma_el(i,j,6)
                !WRITE(12,'(I3,2X,I3,2X,9E17.7)') i,j, f_xyz(i,j),epslon_tot(i,j,1),epslon_tot(i,j,2),epslon_tot(i,j,3),epslon_tot(i,j,4),epslon_tot(i,j,5),epslon_tot(i,j,6)
                WRITE(13,'(I4,2X,I4,2X,9E17.7)') i,j, C1(i,j),C2(i,j),e(i,j,1) ,e(i,j,2),e(i,j,3),e(i,j,4),e(i,j,5),e(i,j,6)
                WRITE(14,'(I4,2X,I4,2X,9E17.7)') i,j, C1(i,j),C2(i,j),Jac_e_green_field(i,j,1) ,Jac_e_green_field(i,j,2),Jac_e_green_field(i,j,3),Jac_e_green_field(i,j,4),Jac_e_green_field(i,j,5),Jac_e_green_field(i,j,6), &
																      Jac_e_green_field(i,j,7) ,Jac_e_green_field(i,j,8),Jac_e_green_field(i,j,9)
				WRITE(15,'(I4,2X,I4,2X,9E17.7)') i,j, C1(i,j),C2(i,j),U_tot(i,j),e_2D_aspect_ratio(i,j),e_3D_aspect_ratio(i,j)

               !WRITE(12,'(I3,2X,I3,2X,9E17.7)') i,j, phi1MAX, f_xyz(i,j),f_el_alpha(i,j),mu_el_alpha(i,j)

            ENDDO; ENDDO;

            WRITE(16 ,'(I9,2X,7E12.5)') itimes, err_x(order), err_y(order)

        ENDIF


770     FORMAT('ZONE',2X,'I=',I5,2X,'J=',I5,2X)
771     FORMAT(I3,2X,I3,2X,9E15.4)

        DEALLOCATE(K,sigma_T,epslon_T,lattice_constant,G_il);
        DEALLOCATE(epslon_st, sigma_0_st, sigma_st, epslon_0, sigma_0);

        !WRITE(*,*) 'ELASTICITY CALCULATIONS FINISHED! ...'
        !STOP

        RETURN

    END subroutine elastic_pot

    !!!
    !!!
    !!!
    !!!

    SUBROUTINE step0_reciprocal_lat(g, K1)

        IMPLICIT NONE

        !output
        REAL(kind), INTENT(out), DIMENSION(1:nex,1:ney,1:3) :: K1
        REAL(kind), INTENT(out), DIMENSION(1:3)             :: g

        !Temp.
        INTEGER :: i
        REAL(kind), DIMENSION(1:3)                            :: g1,g2,g3
        REAL(kind), DIMENSION(1:nex) :: k_1Dx
        REAL(kind), DIMENSION(1:ney) :: k_1Dy
        REAL(kind), DIMENSION(1:nez) :: k_1Dz
        REAL(kind), DIMENSION(1:3)   :: x_hat,y_hat,z_hat
        REAL(kind), DIMENSION(1:nex,1:ney) :: X,Y
        REAL(kind), DIMENSION(1:nex,1:ney,1:nez) :: XX,YY,ZZ

        k_1Dx(1:nex) = 0.0D0
        k_1Dy(1:ney) = 0.0D0
        k_1Dz(1:nez) = 0.0D0

        K1(1:nex,1:ney,1:3)    = 0.0D0

        !! initing the k space

        if (dimen == 1) THEN

        elseif (dimen == 2) THEN

!            k_1Dx(1:(nex))     = [ ( real(i), i=(-int(nex/2)),-1 ) , 0.0 ,  ( real(i), i=1,int(nex/2)) ]
!            k_1Dy(1:(ney))     = [ ( real(i), i=(-int(ney/2)),-1 ) , 0.0 ,  ( real(i), i=1,int(ney/2)) ]
!            k_1Dz = k_1Dx;

            k_1Dx(1:(nex))     = [ ( real(i), i=0,(nex/2-1)) , 0.0 , ( real(i), i=(-nex/2+1),-1 ) ]
            k_1Dy(1:(nex))     = [ ( real(i), i=0,(nex/2-1)) , 0.0 , ( real(i), i=(-nex/2+1),-1 ) ]
            k_1Dz = k_1Dx;

            x_hat = (/1.0, 0.0, 0.0/);
            y_hat = (/0.0, 1.0, 0.0/);
            z_hat = (/0.0, 0.0, 1.0/);

            g1 = 2.0*pi*cross((Ly*y_hat),(Lz*z_hat))/(DOT_PRODUCT(Lx*x_hat,cross((Ly*y_hat),(Lz*z_hat))));
            g2 = 2.0*pi*cross((Lz*z_hat),(Lx*x_hat))/(DOT_PRODUCT(Ly*y_hat,cross((Lz*z_hat),(Lx*x_hat))));
            g3 = 2.0*pi*cross((Lx*x_hat),(Ly*y_hat))/(DOT_PRODUCT(Lz*z_hat,cross((Lx*x_hat),(Ly*y_hat))));
            g = g1 + g2 + g3;

            call meshgrid_2D(k_1Dx, k_1Dy, X, Y)

            K1(1:nex,1:ney,1) = X(1:nex,1:ney);
            K1(1:nex,1:ney,2) = Y(1:nex,1:ney);

        elseif (dimen == 3) THEN

        endif

        RETURN

    END SUBROUTINE step0_reciprocal_lat

        !!!
        !!!
        !!!
        !!!

    !    subroutine meshgrid(xgv, ygv, zgv, X, Y, Z)
    !        implicit none
    !        real(kind),intent(in)   :: xgv(:), ygv(:), zgv(:)
    !        real(kind),intent(out)  :: X(:,:,:), Y(:,:,:), Z(:,:,:)
    !        integer           :: sX, sY, sZ, i
    !
    !        sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)
    !
    !        do i=1,sZ
    !            X(:,:,i) = spread( xgv, 1, sY )
    !            Y(:,:,i) = spread( ygv, 2, sX )
    !        enddo ! i
    !        do i=1,sX
    !            Z(i,:,:) = spread( zgv, 1, sY)
    !        enddo ! i
    !    end subroutine

    !!!
    !!!
    !!!
    !!!

    SUBROUTINE eg_st_tensor(eps_0_ten,eps_0_vec)

        IMPLICIT NONE

        INTEGER i,j,k,l

        REAL(kind), INTENT(IN), DIMENSION(1:9)     :: eps_0_vec
        REAL(kind), INTENT(OUT),DIMENSION(1:3,1:3) :: eps_0_ten

        !**************************************************************************
        !**************************************************************************
        !********************  Elastic Modulus Tensor  ****************************
        !**************************************************************************
        !**************************************************************************

        eps_0_ten(1:3,1:3) = 0.0D0;

        eps_0_ten(1,1) = eps_0_vec(1);
        eps_0_ten(2,2) = eps_0_vec(2);
        eps_0_ten(3,3) = eps_0_vec(3);
        eps_0_ten(1,2) = eps_0_vec(4);
        eps_0_ten(1,3) = eps_0_vec(5);
        eps_0_ten(2,1) = eps_0_vec(6);
        eps_0_ten(2,3) = eps_0_vec(7);
        eps_0_ten(3,1) = eps_0_vec(8);
        eps_0_ten(3,2) = eps_0_vec(9);

        RETURN

    END SUBROUTINE eg_st_tensor
    
    !!!
    !!!
    !!!
    !!!

    SUBROUTINE elmod_tensor(Cv,c)

        IMPLICIT NONE

        INTEGER i,j,k,l

        REAL(kind), INTENT(IN), DIMENSION(1:12)    :: c
        REAL(kind), INTENT(OUT),DIMENSION(1:6,1:6) :: Cv

        !**************************************************************************
        !**************************************************************************
        !********************  Elastic Modulus Tensor  ****************************
        !**************************************************************************
        !**************************************************************************

        Cv(1:6,1:6) = 0.0D0;

        Cv(1,1) = c(1);
        Cv(2,2) = c(2);
        Cv(3,3) = c(3);
        Cv(1,2) = c(4);
        Cv(1,3) = c(5);
        Cv(2,1) = c(6);
        Cv(2,3) = c(7);
        Cv(3,1) = c(kind);
        Cv(3,2) = c(9);
        Cv(4,4) = c(10);
        Cv(5,5) = c(11);
        Cv(6,6) = c(12);

        RETURN

    END SUBROUTINE elmod_tensor
    !
    !    !!!
    !    !!!
    !    !!!
    !    !!!
    !
    !    SUBROUTINE calc_lattice_constant(lattice_constant,a,phi11,C1)
    !
    !        IMPLICIT NONE
    !
    !        !! output
    !        REAL(kind), INTENT(OUT),DIMENSION(1:NPP,1:nex,1:ney)     :: a
    !        REAL(kind), INTENT(OUT),DIMENSION(1:nex,1:ney)           :: lattice_constant
    !        REAL(kind), DIMENSION(1:NPP,1:nex,1:ney)                 :: Delta_a
    !        REAL(kind), DIMENSION(1:NPP,1:nex,1:ney)                 :: r
    !
    !        !! input
    !        REAL(kind), INTENT(IN) ,DIMENSION(1:NPP,-2:IG+2,-2:JG+2) :: phi11
    !        REAL(kind), INTENT(IN) ,DIMENSION(-2:IG+2,-2:JG+2)       :: C1
    !
    !        !! temp. arrays
    !        INTEGER :: i,j,NN,N,KK
    !        REAL(kind), DIMENSION(1:nex,1:ney)        :: h
    !        REAL(kind) :: Delta_a_temp
    !        REAL(kind) :: SUMphi1I(-2:IG+2,-2:JG+2)
    !        REAL(kind) :: SUMphi1K(-2:IG+2,-2:JG+2)
    !        !! Code begins here ....................
    !
    !        !**************************************************************************
    !        !**************************************************************************
    !        !*************************** lattice_constant *****************************
    !        !*************************** **** radius **** *****************************
    !        !**************************************************************************
    !        !**************************************************************************
    !
    !        DO N=NPPK,NPP-3
    !            SUMphi1I(-2:IG+2,-2:JG+2) = SUMphi1I(-2:IG+2,-2:JG+2) + phi11(N,-2:IG+2,-2:JG+2)
    !        ENDDO
    !        DO N=2,NPPK-1
    !            SUMphi1K(-2:IG+2,-2:JG+2) = SUMphi1K(-2:IG+2,-2:JG+2) + phi11(N,-2:IG+2,-2:JG+2)
    !        ENDDO
    !
    !        a(1:NPP,1:nex,1:ney)         = 0.0D0
    !        r(1:NPP,1:nex,1:ney)         = 0.0D0
    !        lattice_constant(1:nex,1:ney)= 0.0D0
    !
    !        r(1,:,:)            = phi11(1         ,1:nex,1:ney)*r_cu
    !        r(2:NPPK-1,:,:)     = phi11(2:NPPK-1  ,1:nex,1:ney)*r_cu3
    !        r(NPPK:NPP-1,:,:)   = phi11(NPPK:NPP-1,1:nex,1:ney)*r_cu6
    !        r(NPP,:,:)          = phi11(NPP       ,1:nex,1:ney)*r_sn
    !
    !        a(1,:,:)            = phi11(1         ,1:nex,1:ney)*a_cu
    !        a(2:NPPK-1,:,:)     = phi11(2:NPPK-1  ,1:nex,1:ney)*a_cu3
    !        a(NPPK:NPP-1,:,:)   = phi11(NPPK:NPP-1,1:nex,1:ney)*a_cu6
    !        a(NPP,:,:)          = phi11(NPP       ,1:nex,1:ney)*a_sn
    !
    !        Do N = 1,NPP
    !            Do i = 1,nex
    !                Do j = 1,ney
    !
    !                    h(I,J) = phi11(N,I,J)**3  * ( 6*phi11(N,I,J)**2 - 15*phi11(N,I,J) + 10 )
    !                    lattice_constant(i,j) = lattice_constant(i,j) + r(N,i,j)*phi11(N,i,j)
    !
    !                    !lattice_constant(i,j) = lattice_constant(i,j) + phi11(1,i,j)*r_cu + SUMphi1I(i,j)*r_cu6 + SUMphi1K(i,j)*r_cu3 + phi11(NPP,i,j)*r_sn !
    !
    !                ENDDO
    !            ENDDO
    !        ENDDO
    !
    !        RETURN
    !
    !    END SUBROUTINE calc_lattice_constant
    !
    !    !!!
    !    !!!
    !    !!!
    !    !!!
    !
    SUBROUTINE rotation(Cvp, Cv)

        IMPLICIT NONE

        ! output
        real(kind), INTENT(OUT), DIMENSION(1:6,1:6) :: Cvp

        ! input
        real(kind), INTENT(IN) , DIMENSION(1:6,1:6) :: Cv

        ! temp.
        integer :: i,j,k,l,m,n,p,q
        real(kind), DIMENSION(1:3,1:3)              :: R
        real(8), DIMENSION(1:3,1:3)                 :: R_temp_8
        real(8) 								    :: theta_r
        
        real(kind), DIMENSION(1:3,1:3,1:3,1:3)      :: C,Cp,C_p


        theta_r = DBLE(pi)*(DBLE(theta)/180);
	    R_temp_8 = reshape( (/ dcos(theta_r),dcos(pi/2.0-theta_r),0.0D0,dcos(pi/2.0+theta_r),dcos(theta_r),0.0D0,0.0D0,0.0D0,1.0D0 /), (/3,3/) );		
		IF(kind==4) THEN
			R = SNGL(R_temp_8) 
		ELSEIF(kind==4) THEN
			R = R_temp_8;
		ENDIF
         
        !   RESHAPE( x, (/2, 2/))
        Do i = 1,3
            Do j = 1,3
                Do k = 1,3
                    Do l = 1,3
                        C(i,j,k,l) = Cv(voigt(i,j),voigt(k,l));
                    enddo
                enddo
            enddo
        enddo

        C_p(1:3,1:3,1:3,1:3) = 0.0;

        Do i = 1,3
            Do j = 1,3
                Do k = 1,3
                    Do l = 1,3
                        Do m = 1,3
                            Do n = 1,3
                                Do p = 1,3
                                    Do q = 1,3
                                        C_p(i,j,k,l) = C_p(i,j,k,l) + R(i,m)*R(j,n)*R(k,p)*R(l,q)*C(m,n,p,q);
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        Do i = 1,3
            Do j = i,3
                Do k = 1,3
                    Do l = k,3
                        Cvp(voigt(i,j),voigt(k,l)) = C_p(i,j,k,l);
                    enddo
                enddo
            enddo
        enddo

        RETURN

    END SUBROUTINE rotation

        !!!
        !!!
        !!!
        !!!
    SUBROUTINE step1_calc_G_il(G_il, K1, g1)

        IMPLICIT NONE

        INTEGER :: i,j,k,l,m,n,o

        ! Output
        REAL(kind), INTENT(OUT), DIMENSION(1:nex,1:ney,1:dimen,1:dimen) :: G_il

        ! Input
        REAL(kind), INTENT(IN) , DIMENSION(1:3)           :: g1
        REAL(kind), INTENT(IN) , DIMENSION(1:nex,1:ney,3) :: K1

        ! Temp. array...
        REAL(kind), ALLOCATABLE, DIMENSION(:,:,:,:) :: G_temp
        REAL(kind), DIMENSION(1:dimen,1:dimen) :: G_t_inv
        REAL(kind), DIMENSION(1:dimen,1:dimen) :: G_tempm

        !! Code begins here ...

        ALLOCATE(G_temp(1:nex,1:ney,1:dimen,1:dimen))

        G_temp(1:nex,1:ney,1:dimen,1:dimen) = 0.0D0  !................Acoustic tensor prior to taking inverse
        G_il  (1:nex,1:ney,1:dimen,1:dimen) = 0.0D0  !................Acoustic tensor prior to taking inverse

        !G_temp = zeros(nex,ney,6,6);
        do i = 1,dimen
            do k = 1,dimen
                do j = 1,dimen
                    do l = 1,dimen
                        G_temp(:,:,i,l) = G_temp(:,:,i,l) + C_effv(voigt(i,j),voigt(k,l))*K1(:,:,j)*g1(j) * K1(:,:,k)*g1(k);
                    enddo
                enddo
            enddo
        enddo

        DO m = 1,ney
            DO n = 1,nex
                !Inverse Acoustic tensor

                G_tempm(1:dimen,1:dimen) = 0.0D0

                DO i = 1,dimen
                    DO j = 1,dimen
                        G_tempm(i,j) = G_temp(m,n,i,j);                    !2-D Acoustic Tensor
                    enddo
                enddo
                !Mark for comment
                Do i = 1,dimen
                    Do j = 1,dimen
                        if ((i == j) .and. (G_tempm(i,j) == 0.0D0)) then
                            G_tempm(i,j) = 1.0D0;
                        endif
                    enddo
                enddo

                !!!!!!!!!!!!!!!!!!!!!!!!!
                G_t_inv = inverse(G_tempm);
                !!!!!!!!!!!!!!!!!!!!!!!!!

                do i = 1,dimen
                    do j = 1,dimen
                        G_il(m,n,i,j) = G_t_inv(i,j);                    ! inverse of C_ijkl*K_j*K_i
                    enddo
                enddo
            enddo
        enddo

        DEALLOCATE(G_temp)

        RETURN

    END SUBROUTINE step1_calc_G_il
    !
    !        !!!
    !        !!!
    !        !!!
    !        !!!
    ! 
    SUBROUTINE step2_calc_epslon_T(epslon_T,Beta_C1,Beta_C2,Beta_C3)

        !! output
        REAL(kind), INTENT(out), DIMENSION(1:nex,1:ney,1:6) :: epslon_T !..... Strength of the eigenstrain; related to epslon_0 through ij epslon_0 = epslon_T (Kronecker delta)_ij,

        !! input
        REAL(kind), INTENT(in),  DIMENSION(1:nex,1:ney) :: Beta_C1,Beta_C2,Beta_C3
        
        !C1,C2,

        !! Temp. arrays
        INTEGER :: i,j
        real(kind) :: r_p,R_m
                


        epslon_T = 0.0D0

        ! epslon_T is the Eigenstrain tensor
        DO i = 1,dimen
            DO j = i,dimen
                  epslon_T(:,:,voigt(i,j)) = epslon_T(:,:,voigt(i,j))  + ( eps_00_1( i,j )*Beta_C1 + eps_00_2( i,j )*Beta_C2 );
            enddo
        enddo        
        
        
        


!        SELECT CASE (mic_type)

!            CASE ("read")

!                epslon_T = eps_T1*Beta_C1 + eps_T2*Beta_C2 + eps_T3*Beta_C3; ! eps_T1:AlN eps_T2:ZrN eps_T3:TiN 

!            CASE ("inclusion")

!                !epslon_T = eps_T1*C1 + eps_T2*C2 + eps_T3*(1.0 - C1 - C2);

!            CASE ("sphere")

!                !epslon_T = eps_T1*C1 + eps_T2*C2 + eps_T3*(1.0 - C1 - C2);

!            CASE ("random")

!                epslon_T = eps_T1*Beta_C1 + eps_T2*Beta_C2 + eps_T3*Beta_C3; ! eps_T1:AlN eps_T2:ZrN eps_T3:TiN 

!            CASE DEFAULT

!                WRITE(*,*) ' STH. went wrong... Please check the error code ...'
!                STOP

!        END SELECT

        RETURN

    END SUBROUTINE step2_calc_epslon_T

        !!!
        !!!
        !!!
        !!!

    SUBROUTINE step3_calc_sigma_T(sigma_T,epslon_T,f_xyz)

        IMPLICIT NONE
        INTEGER :: i,j,m,n

        !! input
        REAL(kind), INTENT(IN), DIMENSION(1:nex,1:ney)    :: f_xyz
        REAL(kind), INTENT(IN), DIMENSION(1:nex,1:ney,1:6)  :: epslon_T

        !! output
        REAL(kind), INTENT(out), DIMENSION(1:nex,1:ney,1:6) :: sigma_T                       ! Eigenstrain tensor

        !! temp.
        REAL(kind), DIMENSION(1:dimen,1:dimen) :: KRONIJ

        !! Code begins here ...

        sigma_T(1:nex,1:ney,1:6) = 0.0D0    ! Eigenstrain tensor

        ! Kronecker delta function
        DO i=1,dimen; DO j=1,dimen
            IF(i==j) THEN
                KRONIJ(i,j) = 1.0D0
            ELSE
                KRONIJ(i,j) = 0.0D0
            ENDIF
        ENDDO; ENDDO

        ! sigma_T is the Eigenstrain tensor
        DO i = 1,dimen
            DO j = i,dimen
                DO m = 1,dimen
                    DO n = m,dimen
                        !KRONIJ(m,n) = INT(( FLOAT((m+n)-ABS(m-n)) )/( FLOAT((m+n)+ABS(m-n))))
                        sigma_T(:,:,voigt(i,j)) = sigma_T(:,:,voigt(i,j)) + C_effv(voigt(i,j),voigt(m,n))*KRONIJ(m,n) * (epslon_T(:,:,voigt(m,n)) );
                    enddo
                enddo
            enddo
        enddo

        RETURN

    END SUBROUTINE step3_calc_sigma_T
    !
    !    !!!
    !    !!!
    !    !!!
    !    !!!
    !
    subroutine step4_calc_S(Sv,C_effv, delta_Cv1,delta_Cv2)

        IMPLICIT NONE
        INTEGER :: i,j
        REAL(kind), DIMENSION(1:6,1:6)            :: S_temp,S_temp1,Sv,C_effv,delta_Cv1,delta_Cv2,dummy
        REAL(kind), ALLOCATABLE, DIMENSION(:,:)   :: S_temp2

        ALLOCATE(S_temp2(1:nex,1:ney))

        S_temp (1:6,1:6) = 0.0D0;
        S_temp1(1:6,1:6) = 0.0D0;

        DO i = 1,6
            DO j = 1,6

                S_temp2(:,:) = C_effv(i,j) + Alpha_C1*delta_Cv1(i,j) + Alpha_C2*delta_Cv2(i,j);
                S_temp(i,j)  = SUM(SUM(S_temp2,DIM=1),DIM=1)/(size(S_temp2,1)*size(S_temp2,2))

            enddo
        enddo

        S_temp1 = S_temp

        if (FindDet(S_temp1,6) == 0.0D0) then
            do i = 1,6
                if (S_temp(i,i) == 0.0D0) then
                    S_temp = 1.0D0;
                endif
            enddo
        endif

        Sv = inverse(S_temp);

        DEALLOCATE(S_temp2)

        RETURN

    end subroutine step4_calc_S

    !!!
    !!!
    !!!
    !!!

    SUBROUTINE step5_calc_disp_field(disp_field, Jac_disp_field, f_xyz1, f_xyz2, sigma_T, epslon_T, G_il, Sv, Beta_C1, Beta_C2, delta_Cv1, delta_Cv2, E_ij, K, g )

        IMPLICIT NONE

        !! output
        REAL(kind), INTENT(inout) , DIMENSION(1:nex,1:ney,1:dimen) :: disp_field

        !! input
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,1:9)       :: Jac_disp_field
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney)           :: f_xyz1, f_xyz2
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,6)         :: sigma_T        ! Eigenstrain tensor
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney)           :: epslon_T       ! Eigenstrain tensor
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,1:dimen,1:dimen) :: G_il
        REAL(kind), INTENT(in)  , DIMENSION(1:6,1:6)                     :: Sv
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney)           :: Beta_C1,Beta_C2
        REAL(kind), INTENT(in)  , DIMENSION(1:6,1:6)               :: delta_Cv1,delta_Cv2
        REAL(kind), INTENT(in)  , DIMENSION(1:3)                   :: g
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,1:3)       :: K
        REAL(kind), INTENT(in)  , DIMENSION(1:6)                   :: E_ij;          ! Homogeneous strain tensor; mean strain tensor of the simulation cell

        !! temp. arrays
        INTEGER ::i,j,l,m,n,e,f
        Complex(kind), ALLOCATABLE,DIMENSION(:,:,:)         :: u_l,F_disp_field
        complex(kind), ALLOCATABLE,DIMENSION(:,:)           :: F_Sigma_Beta , F_space
        complex(kind), ALLOCATABLE,DIMENSION(:)             :: F_Sigma_Beta_vector1
        complex(kind), ALLOCATABLE,DIMENSION(:)             :: F_Sigma_Beta_vector2
        Complex(kind), ALLOCATABLE,DIMENSION(:,:)           :: F_disp_field_tmp1
        complex(kind), ALLOCATABLE,DIMENSION(:)             :: F_disp_field_tmp_vector1

        !! KRONIJ
        REAL(kind), DIMENSION(1:dimen,1:dimen)              :: KRONIJ

        !! Code begins here ...
        !! Code begins here ...
        !! Code begins here ...

        ALLOCATE(u_l(1:nex,1:ney,1:dimen))      ; ALLOCATE(F_disp_field(1:nex,1:ney,1:dimen));
        ALLOCATE(F_Sigma_Beta(1:nex,1:ney))     ; ALLOCATE(F_space(1:nex,1:ney));
        ALLOCATE(F_Sigma_Beta_vector1(1:nex*ney)) ;
        ALLOCATE(F_Sigma_Beta_vector2(1:nex*ney));
        ALLOCATE(F_disp_field_tmp1(1:nex,1:ney));
        ALLOCATE(F_disp_field_tmp_vector1(1:nex*ney));

        u_l         (1:nex,1:ney,1:dimen) = 0.0D0    !.........................Nodal Displacement Field
        F_disp_field(1:nex,1:ney,1:dimen) = 0.0D0    !.........................Fourier Space Displacement Field

        ! Kronecker delta function
        DO i=1,dimen; DO j=1,dimen
            IF(i==j) THEN
                KRONIJ(i,j) = 1.0D0
            ELSE
                KRONIJ(i,j) = 0.0D0
            ENDIF
        ENDDO; ENDDO

        If (order == 1) then
            !Displacement field
            DO l = 1,dimen
                DO i = 1,dimen
                    DO j = 1,dimen

                        F_Sigma_Beta = sigma_T(:,:,voigt(i,j)); !Beta_C1*sigma_T(:,:,voigt(i,j));                        ! correct

                        !! CALL FFT2F FORWARD fft
                        call twoD_fft( F_Sigma_Beta, 'Forward ' )

                        u_l(:,:,l) = u_l(:,:,l) + DCMPLX(0.0D0,-1.0D0)*G_il(:,:,i,l)*K(:,:,j)*g(j)*F_Sigma_Beta ;

                    enddo
                enddo
                F_disp_field(:,:,l) = u_l(:,:,l);            ! should be complex arrays
            enddo

        else             !Higher Order Approximation of the Displacement Field

            u_l(1:ney,1:nex,1:dimen) = 0.0D0

            DO l = 1,dimen
                DO i = 1,dimen
                    DO j = 1,dimen

                        F_Sigma_Beta = sigma_T(:,:,voigt(i,j));  !Beta_C1*sigma_T(:,:,voigt(i,j));

                        call twoD_fft( F_Sigma_Beta, 'Forward ' ) ! NOT IN PLACE

                        u_l(:,:,l) = u_l(:,:,l) + DCMPLX(0,-1)*G_il(:,:,i,l)*K(:,:,j)*g(j)*F_Sigma_Beta    !(reshape(F_Sigma_Beta_vector2, (/ nex,ney /) ));

                        DO m = 1,dimen
                            DO n = 1,dimen

                                F_space = (- delta_Cv1(voigt(i,j),voigt(m,n))*E_ij(voigt(m,n))*Alpha_C1 - delta_Cv2(voigt(i,j),voigt(m,n))*E_ij(voigt(m,n))*Alpha_C2          &
                                + delta_Cv1(voigt(i,j),voigt(m,n))*KRONIJ(m,n)*Alpha_C1*epslon_T*Beta_C1 + delta_Cv2(voigt(i,j),voigt(m,n))*KRONIJ(m,n)*Alpha_C2*epslon_T*Beta_C2    &
                                - delta_Cv1(voigt(i,j),voigt(m,n))*Alpha_C1*(Jac_disp_field(:,:,((m-1)*3+n))) - delta_Cv2(voigt(i,j),voigt(m,n))*Alpha_C2*(Jac_disp_field(:,:,((m-1)*3+n))) );
                                
                                !! FFT2 of F_space
                                !F_space_vector2 = FFT2F(F_space_vector1,nex,ney)
                                call twoD_fft( F_space , 'Forward ' ) ! NOT IN PLACE

                                u_l(:,:,l) = u_l(:,:,l) + DCMPLX(0,-1)*G_il(:,:,i,l)*K(:,:,j)*g(j)*F_space  ! reshape(F_space_vector2, (/ nex,ney /) );

                            enddo
                        enddo
                    enddo
                enddo
                F_disp_field(:,:,l) = u_l(:,:,l);            ! complex array
            enddo

        endif


        DO i = 1,dimen

            F_disp_field_tmp1 = F_disp_field(:,:,i)

            !! Inverse FFT2 of the F_disp_field(:,:,i)
            !F_disp_field_tmp_vector2 = FFT2B(F_disp_field_tmp_vector1,nex,ney)
            call twoD_fft( F_disp_field_tmp1 , 'Backward' ) ! NOT IN PLACE

            disp_field(:,:,i) = real(F_disp_field_tmp1)          !.....Cartesian space displacement field

        enddo

        DEALLOCATE(u_l,F_disp_field,F_Sigma_Beta , F_space,F_Sigma_Beta_vector1,F_Sigma_Beta_vector2,F_disp_field_tmp1,F_disp_field_tmp_vector1)

        RETURN

    END SUBROUTINE step5_calc_disp_field

    !    !!!
    !    !!!
    !    !!!
    !    !!!
    !
    SUBROUTINE step6_calc_Jac(Jac_disp_field, disp_field, K, g);

        IMPLICIT NONE

        !! output
        REAL(kind) , INTENT(out), DIMENSION(1:nex,1:ney,1:9)     :: Jac_disp_field

        !! input
        REAL(kind), INTENT(in)  , DIMENSION(1:ney,1:ney,1:dimen) :: disp_field
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,3)       :: K
        REAL(kind), INTENT(in)  , DIMENSION(1:3)                 :: g

        !! temp. arrays
        COMPLEX(kind) , ALLOCATABLE, DIMENSION(:,:)                    :: DF2D_1,DF2D_2,DF2D_3
        COMPLEX(kind) , ALLOCATABLE, DIMENSION(:,:)                    :: U2D_1x,U2D_2x,U2D_3x,U2D_1y,U2D_2y,U2D_3y,U2D_1z,U2D_2z,U2D_3z

        !! Code begings here ...

        ALLOCATE(DF2D_1(1:nex,1:ney))
        ALLOCATE(DF2D_2(1:nex,1:ney))
        ALLOCATE(DF2D_3(1:nex,1:ney))

        ALLOCATE(U2D_1x(1:nex,1:ney))
        ALLOCATE(U2D_2x(1:nex,1:ney))
        ALLOCATE(U2D_3x(1:nex,1:ney))
        ALLOCATE(U2D_1y(1:nex,1:ney))
        ALLOCATE(U2D_2y(1:nex,1:ney))
        ALLOCATE(U2D_3y(1:nex,1:ney))
        ALLOCATE(U2D_1z(1:nex,1:ney))
        ALLOCATE(U2D_2z(1:nex,1:ney))
        ALLOCATE(U2D_3z(1:nex,1:ney))

        DF2D_1 = disp_field(:,:,1)
        DF2D_2 = disp_field(:,:,2)
        DF2D_3 = disp_field(:,:,3)

        !! CALL FFT2F FORWARD fft calculator of F_Sigma_Beta in 4 steps
        call twoD_fft( DF2D_1, 'Forward ' ) ! NOT IN PLACE
        call twoD_fft( DF2D_2, 'Forward ' ) ! NOT IN PLACE
        call twoD_fft( DF2D_3, 'Forward ' ) ! NOT IN PLACE

        U2D_1x = DCMPLX(0.0D0,1.0D0)*(g(1)*K(:,:,1)) * DF2D_1 !u1 in the x
        U2D_2x = DCMPLX(0.0D0,1.0D0)*(g(1)*K(:,:,1)) * DF2D_2 !u2 in the x
        U2D_3x = DCMPLX(0.0D0,1.0D0)*(g(1)*K(:,:,1)) * DF2D_3 !u3 in the x

        U2D_1y = DCMPLX(0.0D0,1.0D0)*(g(2)*K(:,:,2)) * DF2D_1 !u1 in the y
        U2D_2y = DCMPLX(0.0D0,1.0D0)*(g(2)*K(:,:,2)) * DF2D_2 !u2 in the y
        U2D_3y = DCMPLX(0.0D0,1.0D0)*(g(2)*K(:,:,2)) * DF2D_3 !u3 in the y

        U2D_1z = DCMPLX(0.0D0,1.0D0)*(g(3)*K(:,:,3)) * DF2D_1 !u1 in the z
        U2D_2z = DCMPLX(0.0D0,1.0D0)*(g(3)*K(:,:,3)) * DF2D_2 !u2 in the z
        U2D_3z = DCMPLX(0.0D0,1.0D0)*(g(3)*K(:,:,3)) * DF2D_3 !u2 in the z

        !!!!!!!
        !! CALL FFT2B BACKWARD fft calculator in 4 steps
        call twoD_fft( U2D_1x, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( U2D_2x, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( U2D_3x, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( U2D_1y, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( U2D_2y, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( U2D_3y, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( U2D_1z, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( U2D_2z, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( U2D_3z, 'Backward' ) ! NOT IN PLACE

        Jac_disp_field(1:nex,1:ney,1:9) = 0.0D0

        Jac_disp_field(:,:,1) = real( U2D_1x );        !u1 in the x
        Jac_disp_field(:,:,2) = real( U2D_1y );        !u1 in the y
        Jac_disp_field(:,:,3) = real( U2D_1z );        !u1 in the z

        Jac_disp_field(:,:,4) = real( U2D_2x );        !u2 in the x
        Jac_disp_field(:,:,5) = real( U2D_2y );        !u2 in the y
        Jac_disp_field(:,:,6) = real( U2D_2z );        !u2 in the y

        Jac_disp_field(:,:,7) = real( U2D_3x );        !u3 in the x
        Jac_disp_field(:,:,8) = real( U2D_3y );        !u3 in the y
        Jac_disp_field(:,:,9) = real( U2D_3z );        !u3 in the y

        DEALLOCATE(DF2D_1,DF2D_2,DF2D_3)
        DEALLOCATE(U2D_1x,U2D_2x,U2D_3x,U2D_1y,U2D_2y,U2D_3y,U2D_1z,U2D_2z,U2D_3z)


        !Jac_disp_field(:,:,1) = real(ifft2(((1j).*(g(1).*K(:,:,1))).*fft2(disp_field(:,:,1)))); !u1 in the x
        !Jac_disp_field(:,:,4) = real(ifft2(((1j).*(g(1).*K(:,:,1))).*fft2(disp_field(:,:,2)))); !u2 in the x
        !Jac_disp_field(:,:,2) = real(ifft2(((1j).*(g(2).*K(:,:,2))).*fft2(disp_field(:,:,1)))); !u1 in the y
        !Jac_disp_field(:,:,5) = real(ifft2(((1j).*(g(2).*K(:,:,2))).*fft2(disp_field(:,:,2)))); !u2 in the y

        RETURN

    END SUBROUTINE step6_calc_Jac

    !!!
    !!!
    !!!
    !!!

    SUBROUTINE step7_calc_stars(epslon_st, sigma_st, sigma_0_st, epslon_0, Jac_disp_field, Beta_C1, Beta_C2, C_effv, delta_Cv1, delta_Cv2, epslon_T)

        IMPLICIT NONE

        !! output
        REAL(kind), INTENT(out) , dimension(1:nex,1:ney,1:6) :: epslon_st          ! ......... Period strain
        REAL(kind), INTENT(out) , dimension(1:nex,1:ney,1:6) :: sigma_0_st         ! ......... Period stress
        REAL(kind), INTENT(out) , dimension(1:nex,1:ney,1:6) :: sigma_st           ! ......... Period strain
        REAL(kind), INTENT(out) , dimension(1:nex,1:ney,1:6) :: epslon_0           ! ......... Eigenstrain

        !! input
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,9)   :: Jac_disp_field
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney)     :: Beta_C1,Beta_C2
        REAL(kind), INTENT(in)  , DIMENSION(1:6,1:6)         :: C_effv !
        REAL(kind), INTENT(in)  , DIMENSION(1:6,1:6)         :: delta_Cv1,delta_Cv2
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,6)     :: epslon_T

        !! temp. arrays
        integer :: i,j,k,l
        REAL(kind), ALLOCATABLE, dimension(:,:,:) ::  sigma_el,sigma_0
        REAL(kind), DIMENSION(1:dimen,1:dimen)  :: KRONIJ
        REAL(kind), DIMENSION(1:6)              :: eps_f

        !! Code begins here ....................

        ALLOCATE(sigma_el(1:nex,1:ney,1:6))
        ALLOCATE(sigma_0 (1:nex,1:ney,1:6))

        ! Kronecker delta function
        DO i=1,dimen; DO j=1,dimen
            IF(i==j) THEN
                KRONIJ(i,j) = 1.0D0
            ELSE
                KRONIJ(i,j) = 0.0D0
            ENDIF
        ENDDO; ENDDO

        epslon_st(1:nex,1:ney,1:6) = 0.0D0;        !............. (compatible total strain)
        Do i = 1,dimen
            Do j = i,dimen
                epslon_st(:,:,voigt(i,j)) = 0.5D0*Jac_disp_field(:,:,((i-1)*3+j))+ 0.5D0*Jac_disp_field(:,:,((j-1)*3+i));
            enddo
        enddo

        !eps_f = [1; 1; 1; 2; 2; 2];
        eps_f = (/ 1,1,1,2,2,2 /);
        epslon_0   = 0.0D0;        !............. eigenstrain (position dependent dilatational eigenstrain field)
        sigma_el   = 0.0D0;        !.............
        sigma_0    = 0.0D0;        !.............
        sigma_st   = 0.0D0;        !.............
        sigma_0_st = 0.0D0;        !.............

        !! eigenstrain
        Do i = 1,dimen
            Do j = i,dimen
                epslon_0(:,:,voigt(i,j)) = (epslon_T(:,:,voigt(i,j)))*KRONIJ(i,j);
            enddo
        enddo

        !!
        Do i = 1,dimen
            Do j = i,dimen
                Do k = 1,dimen
                    Do l = k,dimen
                        sigma_0_st(:,:,voigt(i,j)) = sigma_0_st(:,:,voigt(i,j)) + &
                        (C_effv(voigt(i,j),voigt(k,l)) + Alpha_C1*delta_Cv1(voigt(i,j),voigt(k,l)) + Alpha_C2*delta_Cv2(voigt(i,j),voigt(k,l)))*epslon_0(:,:,voigt(k,l))*eps_f(voigt(k,l));

                        sigma_st(:,:,voigt(i,j))   = sigma_st(:,:,voigt(i,j))   + & 
                        (C_effv(voigt(i,j),voigt(k,l)) + Alpha_C1*delta_Cv1(voigt(i,j),voigt(k,l)) + Alpha_C2*delta_Cv2(voigt(i,j),voigt(k,l)) )*epslon_st(:,:,voigt(k,l))*eps_f(voigt(k,l));
                    enddo
                enddo
            enddo
        enddo

        DEALLOCATE(sigma_el, sigma_0)

    END SUBROUTINE step7_calc_stars

    !!!
    !!!
    !!!
    !!!

    SUBROUTINE step8_calc_E_ij(E_ij, sigma_0_st, sigma_st, sigma_a, Sv)

        !! output
        REAL(kind), INTENT(out) , DIMENSION(1:6)           :: E_ij          !........................ Homogeneous strain tensor; mean strain tensor of the simulation cell

        !! input
        REAL(kind), INTENT(in)  , dimension(1:nex,1:ney,1:6) :: sigma_0_st, sigma_st
        REAL(kind), INTENT(in)  , dimension(1:6)     :: sigma_a
        REAL(kind), INTENT(in)  , DIMENSION(1:6,1:6) :: Sv

        !! temp. arrays
        integer :: i,j,k,l
        REAL(kind), dimension(1:6) :: sigma_avg,sigma_0_avg

        sigma_avg(1:6)   = 0.0D0;
        sigma_0_avg(1:6) = 0.0D0;

        Do i = 1,dimen
            Do j = i,dimen

                sigma_0_avg(voigt(i,j)) = mean(sigma_0_st(:,:,voigt(i,j)),nex,ney);
                sigma_avg(voigt(i,j))   = mean(sigma_st  (:,:,voigt(i,j)),nex,ney);

            enddo
        enddo

        E_ij(1:6) = 0.0D0;

        Do i = 1,dimen
            Do j = i,dimen
                Do k = 1,dimen
                    Do l = k,dimen
                        E_ij(voigt(i,j)) = E_ij(voigt(i,j)) + Sv(voigt(i,j),voigt(k,l))*(sigma_a(voigt(k,l)) + sigma_0_avg(voigt(k,l)) - sigma_avg(voigt(k,l)));
                    enddo
                enddo
            enddo
        enddo

    END SUBROUTINE step8_calc_E_ij

    !!!
    !!!
    !!!
    !!!

    SUBROUTINE step9_calc_el(epslon_el, sigma_el, E_ij, epslon_st, epslon_0, C_effv, delta_Cv1,delta_Cv2)

        IMPLICIT NONE

        !! output
        REAL(kind), INTENT(out) , DIMENSION(1:nex,1:ney,1:6) :: epslon_el, sigma_el

        !! input
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,1:6) :: epslon_st             !........... Period strain
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,1:6) :: epslon_0              !........... Eigenstrain
        REAL(kind), INTENT(in)  , DIMENSION(1:6)             :: E_ij                  !........... Homogeneous strain tensor; mean strain tensor of the simulation cell
        REAL(kind), INTENT(in)  , DIMENSION(1:6,1:6)         :: C_effv
        REAL(kind), INTENT(in)  , DIMENSION(1:6,1:6)         :: delta_Cv1,delta_Cv2

        !! temp. arrays
        integer :: i,j,k,l
        REAL(kind), DIMENSION(1:nex,1:ney,1:6) :: epslon
        REAL(kind), DIMENSION(1:6)             :: eps_f

        !! Code begins here ...

        sigma_el(1:nex,1:ney,6) = 0.0D0;        !................Elastic stress tensor
        epslon  (1:nex,1:ney,6) = 0.0D0;        !................Total strain

        Do i = 1,dimen
            Do j = i,dimen
                epslon(:,:,voigt(i,j)) = E_ij(voigt(i,j)) + epslon_st(:,:,voigt(i,j));
            enddo
        enddo


        epslon_el = epslon - epslon_0;

        !eps_f = [1; 1; 1; 2; 2; 2];
        eps_f = (/ 1,1,1,2,2,2 /);


        Do i = 1,dimen
            Do j = i,dimen
                Do k = 1,dimen
                    Do l = k,dimen
                        sigma_el(:,:,voigt(i,j)) = sigma_el(:,:,voigt(i,j)) + & 
                        ( C_effv(voigt(i,j),voigt(k,l)) + Alpha_C1*delta_Cv1(voigt(i,j),voigt(k,l)) + Alpha_C2*delta_Cv2(voigt(i,j),voigt(k,l)) )*epslon_el(:,:,voigt(k,l))*eps_f(voigt(k,l));
                    enddo
                enddo
            enddo
        enddo

        RETURN

    END SUBROUTINE step9_calc_el

    !    !!!
    !    !!!
    !    !!!
    !    !!!

    SUBROUTINE step10_calc_mu_el(mu_el_1,mu_el_2,f_xyz1,f_xyz2,epslon_0,epslon_T,epslon_st,E_ij,epslon_el,sigma_el,delta_Cv1,delta_Cv2)

        !! output
        REAL(kind), INTENT(OUT), DIMENSION(1:nex,1:ney)   :: mu_el_1,mu_el_2

        !! input
        REAL(kind), INTENT(IN) , DIMENSION(1:nex,1:ney)   :: f_xyz1,f_xyz2
        REAL(kind), INTENT(IN) , DIMENSION(1:nex,1:ney,6) :: epslon_0                        ! Eigenstrain tensor
        REAL(kind), INTENT(IN) , DIMENSION(1:nex,1:ney,6) :: epslon_st
        REAL(kind), INTENT(IN) , DIMENSION(1:nex,1:ney)   :: epslon_T
        REAL(kind), INTENT(IN) , DIMENSION(1:nex,1:ney,6) :: epslon_el
        REAL(kind), INTENT(IN) , DIMENSION(1:nex,1:ney,6) :: sigma_el
        REAL(kind), INTENT(in) , DIMENSION(1:6,1:6)       :: delta_Cv1,delta_Cv2
        REAL(kind), INTENT(IN) , DIMENSION(1:6)           :: E_ij;                           ! mean strain tensor of the simulation cell

        !! temp. arrays
        integer :: i,j,k,l,p
        REAL(kind), ALLOCATABLE, DIMENSION(:,:,:)   	  :: epsilon_star_alpha,epsilon_alpha
        REAL(kind), DIMENSION(1:dimen,1:dimen) 		  	  :: KRONIJ
        REAL(kind), DIMENSION(1:6)             		 	  :: eps_f

        !! Code begins here ....................
        
        ALLOCATE(epsilon_star_alpha(1:nex,1:ney,6))
        ALLOCATE(epsilon_alpha(1:nex,1:ney,6))

        ! Kronecker delta function
        DO i=1,dimen; DO j=1,dimen
            IF(i==j) THEN
                KRONIJ(i,j) = 1.0D0
            ELSE
                KRONIJ(i,j) = 0.0D0
            ENDIF
        ENDDO; ENDDO

        epsilon_star_alpha = epslon_0
        epsilon_alpha      = epslon_el
        !f_el_alpha         = 0.0D0
        mu_el_1	           = 0.0D0
        mu_el_2	           = 0.0D0

        eps_f = (/ 1,1,1,2,2,2 /);

        Do i = 1,dimen
            Do j = i,dimen
                Do k = 1,dimen
                    Do l = k,dimen

                        !f_el_alpha(:,:,1) = f_el_alpha(:,:,1) + f_xyz(:,:)*( epsilon_alpha(:,:,voigt(i,j)) - epsilon_star_alpha(:,:,voigt(i,j) ) ) * C_effv(voigt(i,j),voigt(k,l))   &
                        !*( epsilon_alpha(:,:,voigt(k,l)) - epsilon_star_alpha(:,:,voigt(k,l) ) )

                        !f_el_alpha(:,:)  = f_el_alpha(:,:) + 0.5D0*(epslon_el(:,:,voigt(i,j))* sigma_el(:,:,voigt(i,j)) )/((Nv_p*(f_xyz1) + Nv_m*(1-f_xyz1))*(dx*dy))

!                        mu_el_1(:,:) = mu_el_1(:,:) + ( (0.5*Alpha_C1_prime*delta_Cv1(voigt(i,j),voigt(k,l))*(E_ij(voigt(i,j)) &
!                        + epslon_st(:,:,voigt(i,j)) - epslon_0(:,:,voigt(i,j)))*(E_ij(voigt(k,l))   &
!                        + epslon_st(:,:,voigt(k,l)) - epslon_0(:,:,voigt(k,l)))*eps_f(voigt(i,j))*eps_f(voigt(k,l))) &
!                        - Beta_C1_prime*(epslon_T*KRONIJ(i,j)*C_c(:,:,voigt(i,j),voigt(k,l))*(E_ij(voigt(k,l)) &
!                        + epslon_st(:,:,voigt(k,l)) - epslon_0(:,:,voigt(k,l)))*eps_f(voigt(k,l))) & 
!                        + (0.5*Alpha_C2_prime*delta_Cv1(voigt(i,j),voigt(k,l))*(E_ij(voigt(i,j)) &
!                        + epslon_st(:,:,voigt(i,j)) - epslon_0(:,:,voigt(i,j)))*(E_ij(voigt(k,l))   &
!                        + epslon_st(:,:,voigt(k,l)) - epslon_0(:,:,voigt(k,l)))*eps_f(voigt(i,j))*eps_f(voigt(k,l))) &
!                        - Beta_C2_prime*(epslon_T*KRONIJ(i,j)*C_c(:,:,voigt(i,j),voigt(k,l))*(E_ij(voigt(k,l)) &
!                        + epslon_st(:,:,voigt(k,l)) - epslon_0(:,:,voigt(k,l)))*eps_f(voigt(k,l))) ) &
!                        /((Nv_p*(f_xyz1) + Nv_m*(1-f_xyz1))*(dx**dimen));

                        mu_el_1(:,:) = mu_el_1(:,:) + ( (0.5*Alpha_C1_prime*delta_Cv1(voigt(i,j),voigt(k,l))*(E_ij(voigt(i,j)) &
                        + epslon_st(:,:,voigt(i,j)) - epslon_0(:,:,voigt(i,j)))*(E_ij(voigt(k,l))   &
                        + epslon_st(:,:,voigt(k,l)) - epslon_0(:,:,voigt(k,l)))*eps_f(voigt(i,j))*eps_f(voigt(k,l))) &
                        - Beta_C1_prime*(epslon_T*KRONIJ(i,j)*C_c(:,:,voigt(i,j),voigt(k,l))*(E_ij(voigt(k,l)) &
                        + epslon_st(:,:,voigt(k,l)) - epslon_0(:,:,voigt(k,l)))*eps_f(voigt(k,l))) ) & 
                        *(V_p1*(f_xyz1) + V_p2*(f_xyz2) + V_m*(1.0-f_xyz1-f_xyz2))!/(dx**dimen);
                        
                        
                        !/((Nv_p1*(f_xyz1) + Nv_m*(1-f_xyz1))*(dx**dimen));
                        
                        !*(V_p1*(f_xyz) + V_m*(1.0-f_xyz))!/(dx**2);

                        !/((Nv_p1*(f_xyz) + Nv_m*(1-f_xyz))*(dx**dimen));

                    Enddo
                Enddo
            Enddo
        Enddo
        
        Do i = 1,dimen
            Do j = i,dimen
                Do k = 1,dimen
                    Do l = k,dimen

                        !f_el_alpha(:,:,1) = f_el_alpha(:,:,1) + f_xyz(:,:)*( epsilon_alpha(:,:,voigt(i,j)) - epsilon_star_alpha(:,:,voigt(i,j) ) ) * C_effv(voigt(i,j),voigt(k,l))   &
                        !*( epsilon_alpha(:,:,voigt(k,l)) - epsilon_star_alpha(:,:,voigt(k,l) ) )

                        !f_el_alpha(:,:)  = f_el_alpha(:,:) + 0.5D0*(epslon_el(:,:,voigt(i,j))* sigma_el(:,:,voigt(i,j)) )/((Nv_p*(f_xyz1) + Nv_m*(1-f_xyz1))*(dx*dy))

                        mu_el_2(:,:) = mu_el_2(:,:) + ( (0.5*Alpha_C2_prime*delta_Cv2(voigt(i,j),voigt(k,l))*(E_ij(voigt(i,j)) &
                        + epslon_st(:,:,voigt(i,j)) - epslon_0(:,:,voigt(i,j)))*(E_ij(voigt(k,l))   &
                        + epslon_st(:,:,voigt(k,l)) - epslon_0(:,:,voigt(k,l)))*eps_f(voigt(i,j))*eps_f(voigt(k,l))) &
                        - Beta_C2_prime*( epslon_T*KRONIJ(i,j)*C_c(:,:,voigt(i,j),voigt(k,l))*(E_ij(voigt(k,l)) &
                        + epslon_st(:,:,voigt(k,l)) - epslon_0(:,:,voigt(k,l)))*eps_f(voigt(k,l)) ) ) &
                        *(V_p1*(f_xyz1) + V_p2*(f_xyz2) + V_m*(1.0-f_xyz1-f_xyz2))!/(dx**dimen);

                        !*(V_p1*(f_xyz) + V_m*(1.0-f_xyz))!/(dx**2);

                        !/((Nv_p1*(f_xyz) + Nv_m*(1-f_xyz))*(dx**dimen));

                    Enddo
                Enddo
            Enddo
        Enddo

        DEALLOCATE(epsilon_star_alpha,epsilon_alpha)

    END SUBROUTINE step10_calc_mu_el

    !!!
    !!!
    !!!
    !!!

    SUBROUTINE step12_calc_Jac(U_tot, Jac_e_green_field, e, K, g);

        IMPLICIT NONE

        !! output
        REAL(kind) , INTENT(out), DIMENSION(1:nex,1:ney,1:9)     :: Jac_e_green_field
        REAL(kind) , INTENT(out), DIMENSION(1:nex,1:ney)	     :: U_tot

        !! input
        REAL(kind), INTENT(in)  , DIMENSION(1:ney,1:ney,1:dimen) :: e
        REAL(kind), INTENT(in)  , DIMENSION(1:nex,1:ney,3)       :: K
        REAL(kind), INTENT(in)  , DIMENSION(1:3)                 :: g

        !! temp. arrays
        COMPLEX(kind) , ALLOCATABLE, DIMENSION(:,:)              :: e_1,e_2,e_3,e_4,e_5,e_6
        COMPLEX(kind) , ALLOCATABLE, DIMENSION(:,:)              :: e_11,e_12,e_13,e_21,e_22,e_23,e_31,e_32,e_33
        REAL(kind), DIMENSION(1:nex,1:ney)   				     :: PHI,PSI
		REAL(kind)												 :: A,B,C,gg,h

        !! Code begings here ...

        ALLOCATE(e_1(1:nex,1:ney))
        ALLOCATE(e_2(1:nex,1:ney))
        ALLOCATE(e_3(1:nex,1:ney))
        ALLOCATE(e_4(1:nex,1:ney))
        ALLOCATE(e_5(1:nex,1:ney))
        ALLOCATE(e_6(1:nex,1:ney))

        ALLOCATE(e_11(1:nex,1:ney))
        ALLOCATE(e_12(1:nex,1:ney))
        ALLOCATE(e_13(1:nex,1:ney))
        ALLOCATE(e_21(1:nex,1:ney))
        ALLOCATE(e_22(1:nex,1:ney))
        ALLOCATE(e_23(1:nex,1:ney))
        ALLOCATE(e_31(1:nex,1:ney))
        ALLOCATE(e_32(1:nex,1:ney))
        ALLOCATE(e_33(1:nex,1:ney))

        e_1 = e(:,:,1)
        e_2 = e(:,:,2)
        e_3 = e(:,:,3)

        !! CALL FFT2F FORWARD fft calculator of F_Sigma_Beta in 4 steps
        call twoD_fft( e_1, 'Forward ' ) ! NOT IN PLACE
        call twoD_fft( e_2, 'Forward ' ) ! NOT IN PLACE
        call twoD_fft( e_3, 'Forward ' ) ! NOT IN PLACE

        e_11 = DCMPLX(0.0D0,1.0D0)*(g(1)*K(:,:,1)) * e_1 !u1 in the x
        e_12 = DCMPLX(0.0D0,1.0D0)*(g(2)*K(:,:,2)) * e_1 !u2 in the x
        e_13 = DCMPLX(0.0D0,1.0D0)*(g(3)*K(:,:,3)) * e_1 !u3 in the x

        e_21 = DCMPLX(0.0D0,1.0D0)*(g(1)*K(:,:,1)) * e_2 !u1 in the y
        e_22 = DCMPLX(0.0D0,1.0D0)*(g(2)*K(:,:,2)) * e_2 !u2 in the y
        e_23 = DCMPLX(0.0D0,1.0D0)*(g(3)*K(:,:,3)) * e_2 !u3 in the y

        e_31 = DCMPLX(0.0D0,1.0D0)*(g(1)*K(:,:,1)) * e_3 !u1 in the z
        e_32 = DCMPLX(0.0D0,1.0D0)*(g(2)*K(:,:,2)) * e_3 !u2 in the z
        e_33 = DCMPLX(0.0D0,1.0D0)*(g(3)*K(:,:,3)) * e_3 !u2 in the z

        !!!!!!!
        !! CALL FFT2B BACKWARD fft calculator in 4 steps
        call twoD_fft( e_11, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( e_12, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( e_13, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( e_21, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( e_22, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( e_23, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( e_31, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( e_32, 'Backward' ) ! NOT IN PLACE
        call twoD_fft( e_33, 'Backward' ) ! NOT IN PLACE

        Jac_e_green_field(1:nex,1:ney,1:9) = 0.0D0

        Jac_e_green_field(:,:,1) = real( e_11 );        !u1 in the x
        Jac_e_green_field(:,:,2) = real( e_12 );        !u1 in the y
        Jac_e_green_field(:,:,3) = real( e_13 );        !u1 in the z

        Jac_e_green_field(:,:,4) = real( e_21 );        !u2 in the x
        Jac_e_green_field(:,:,5) = real( e_22 );        !u2 in the y
        Jac_e_green_field(:,:,6) = real( e_23 );        !u2 in the y

        Jac_e_green_field(:,:,7) = real( e_31 );        !u3 in the x
        Jac_e_green_field(:,:,8) = real( e_32 );        !u3 in the y
        Jac_e_green_field(:,:,9) = real( e_33 );        !u3 in the y



        A = 1.0D0;
        B = 1.0D0;
        C = 1.0D0;
        gg= 1.0D0;
        h = 1.0D0;
        
        PHI = A*(e(:,:,2)**2 + e(:,:,3)**2) + B*e(:,:,3)*( e(:,:,3)**2 - 3.0*e(:,:,2)**2) + C*( e(:,:,2)**2 + e(:,:,3)**2 )**2
		PSI = gg*( e_21**2 + e_22**2 + (e_31**2 + e_32**2)/3.0 + 2.0*(e_21*e_31-e_22*e_32)/sqrt(3.0) ) + h*( e_31**2 + e_32**2 - sqrt(3.0)*(e_21*e_31 - e_22*e_32) );
		
		U_tot = PHI + PSI; 

        DEALLOCATE(e_1,e_2,e_3)
        DEALLOCATE(e_11,e_12,e_13,e_21,e_22,e_23,e_31,e_32,e_33)



        !Jac_disp_field(:,:,1) = real(ifft2(((1j).*(g(1).*K(:,:,1))).*fft2(disp_field(:,:,1)))); !u1 in the x
        !Jac_disp_field(:,:,4) = real(ifft2(((1j).*(g(1).*K(:,:,1))).*fft2(disp_field(:,:,2)))); !u2 in the x
        !Jac_disp_field(:,:,2) = real(ifft2(((1j).*(g(2).*K(:,:,2))).*fft2(disp_field(:,:,1)))); !u1 in the y
        !Jac_disp_field(:,:,5) = real(ifft2(((1j).*(g(2).*K(:,:,2))).*fft2(disp_field(:,:,2)))); !u2 in the y

        RETURN

    END SUBROUTINE step12_calc_Jac

    !!!
    !!!
    !!!
    !!!

END MODULE elastic_potential

!! ---------------------------------









!! ---------------------------------

MODULE chemical_potential

    use mod_geometry
    use mtls_property
    use solver_opts
    use wrt_opts
    use math_opts

    implicit none

CONTAINS

    subroutine chemical_pot(phi, chem_pot, int_pot, kx, ky, k2x, k2y)

        ! output
        real(8), intent(out), dimension(nex,ney,1:2) :: chem_pot,int_pot

        ! input
        real(8), intent(in), dimension(nex,ney,1:2)  :: phi
        REAL(8), intent(in), DIMENSION(1:nex,1:ney)  :: kx,ky,k2x,k2y


        complex(8), dimension(1:nex,1:ney)       :: phi_tmp
        real(8), dimension(1:nex,1:ney)          :: lap_phi

        chem_pot = 0.0D0
        int_pot  = 0.0D0

        ! Chemical potential
        chem_pot(:,:,1) = - delta_f*( A*phi(:,:,1) - B*phi(:,:,1)**2 + C*phi(:,:,1)*(phi(:,:,1)**2 + phi(:,:,2)**2) )
        chem_pot(:,:,2) = - delta_f*( A*phi(:,:,2) - B*phi(:,:,2)**2 + C*phi(:,:,2)*(phi(:,:,1)**2 + phi(:,:,2)**2) )

        !chem_pot = chem_pot

!       ! Interfacial potential 1
        phi_tmp = phi(:,:,1)
        call twoD_fft(phi_tmp,'Forward ')
        !phi_tmp = ((DCMPLX(0.0,1.0)*(g(1)*K2x(:,:) + g(2)*K2y(:,:)))**2)*phi_tmp;        !.....................Fourier space 2nd derivative of f(x,y)
        phi_tmp = ( k2x(:,:) + k2y(:,:) )*phi_tmp;        !.....................Fourier space 2nd derivative of f(x,y)
        call twoD_fft(phi_tmp,'Backward')
        lap_phi = real(phi_tmp)
        int_pot(:,:,1)  = beta*lap_phi(:,:)

!        ! Interfacial potential 2
        phi_tmp = phi(:,:,2)
        call twoD_fft(phi_tmp,'Forward ')
        !phi_tmp = ((DCMPLX(0.0,1.0)*(g(1)*K2x(:,:) + g(2)*K2y(:,:)))**2)*phi_tmp;        !.....................Fourier space 2nd derivative of f(x,y)
        phi_tmp = ( k2x(:,:) + k2y(:,:) )*phi_tmp;        !.....................Fourier space 2nd derivative of f(x,y)
        call twoD_fft(phi_tmp,'Backward')
        lap_phi = real(phi_tmp)
        int_pot(:,:,2)  = beta*lap_phi(:,:)
        






       ! Interfacial potential 1
!        phi_tmp = phi(:,:,1)
!        ! step 1
!        call twoD_fft(phi_tmp,'Forward ')
!        phi_tmp = ( kx(:,:) )*phi_tmp;        !.....................Fourier space 2nd derivative of f(x,y)
!        call twoD_fft(phi_tmp,'Backward')
!        ! step 2
!        call twoD_fft(phi_tmp,'Forward ')
!        phi_tmp = ( ky(:,:) )*phi_tmp;        !.....................Fourier space 2nd derivative of f(x,y)
!        call twoD_fft(phi_tmp,'Backward')
!        lap_phi = real(phi_tmp)
!        int_pot(:,:,1)  = beta*lap_phi(:,:)
!
!
       ! Interfacial potential 2
!        phi_tmp = phi(:,:,2)
!        ! step 1
!        call twoD_fft(phi_tmp,'Forward ')
!        phi_tmp = ( kx(:,:) )*phi_tmp;        !.....................Fourier space 2nd derivative of f(x,y)
!        call twoD_fft(phi_tmp,'Backward')
!        ! step 2
!        call twoD_fft(phi_tmp,'Forward ')
!        phi_tmp = ( ky(:,:) )*phi_tmp;        !.....................Fourier space 2nd derivative of f(x,y)
!        call twoD_fft(phi_tmp,'Backward')
!        lap_phi = real(phi_tmp)
!        int_pot(:,:,2)  = beta*lap_phi(:,:)



        RETURN

    end subroutine chemical_pot

     !=======================================================================================================================


END MODULE chemical_potential

!! ---------------------------------

program main

    USE omp_lib
    use mod_geometry
    use mtls_property
    use solver_opts
    use wrt_opts
    use chemical_potential
    use elastic_potential
    use pixelfrt

    implicit none
    integer::i,j

    real(8), dimension(1:nex,1:ney,1:2)    :: PHI
    real(8), dimension(1:nex,1:ney,1:2)    :: chem_pot
    real(8), dimension(1:nex,1:ney,1:2)    :: elas_pot
    real(8), dimension(1:nex,1:ney,1:2)    :: int_pot
    

    REAL(8), DIMENSION(1:nex,1:ney,1:6):: epslon_el, sigma_el

    !! Reciprocal space arrays
    REAL(8), DIMENSION(1:3)            :: g
    REAL(8), DIMENSION(1:nex,1:ney)    :: Kxx,Kyy,k2x,k2y
    REAL(8), DIMENSION(1:nex)          :: X,Y
    real(8) :: time_begin,time_end
    real(8) :: area_si,Vf_si,err,Vf_new,dt_temp,Vf_old_sn,Vf_old_si
    INTEGER :: Reason

    !! code begins here ..............
    CALL CPU_TIME(time_begin)

!# of parameters    List of parameters for Mg2Sn-Mg2Si system   symbol in matlab    Lower bound Upper bound
!1  Initial Composition N/A 0.3 0.5
!2  Molar volume (phase 1) (Vm) Vm_cryst_Mg2Sn  see doc see doc
!3  Molar volume (phase 2) (Vm) Vm_cryst_Mg2Si  see mat or doc file see mat or doc file
!4  Lattice mismatch    N/A -0.02   0.02
!5,6,7  Elastic constants 3 for phase 1 (C11,C12 and C44)   C11_mg2si_calc, ... see mat or doc file see mat or doc file
!8,9,10 Elastic constants 3 for phase 2 (C11,C12 and C44)   C11_mg2sn_calc, ... see mat or doc file see mat or doc file
!11 Interface mobility
!   M   see mat or doc file see mat or doc file
!12 Gradient energy coefficient
!   kappa   see mat or doc file see mat or doc file
!13-14-15-16-17-18
!   thermodynamic parameters 6 total    N/A see mat or doc file see mat or doc file
    
    call mkdir
    call wrt_disk
    call system('rm -r *.mod')
    OPEN(unit=53,FILE=TRIM(fileplace)//'/results/K.dat')
    OPEN(unit=54,FILE=TRIM(fileplace)//'/results/parameter.dat')
    OPEN(unit=55,FILE=TRIM(fileplace)//'/results/time_map.dat')
    OPEN(unit=1002,FILE=TRIM(fileplace)//'/results/screen.dat')

    !! Code begins here ... 
    

!    !***** INITIALIZATION ***********************************************
!    WRITE(*,*) '******* INPUT FILE READ RESULTS ******************'
!    open(unit=1000, FILE=TRIM(fileplace)//'/input/var.txt',ACTION="READ")
!    READ(1000,'(18E14.6)',IOSTAT=Reason) X11,V_m,V_p,eps_T,C11_m,C12_m,C44_m,C11_p,C12_p,C44_p,mobility,kappa,L0_Si_Sn,L1_Si_Sn,L2_Si_Sn,L0_Si_Sn_liq,L1_Si_Sn_liq,L2_Si_Sn_liq
!    IF(Reason.gt.0) write(*,*) 'var.txt : sth wrong'
!    IF(Reason.lt.0) write(*,*) 'var.txt : end of file'
!    IF(Reason.eq.0) write(*,*) 'var.txt : SUCCESSFUL READING...'
!    CLOSE(unit=1000, STATUS='KEEP');
!    IF(Reason.ne.0) stop
    !WRITE(*,'(18E14.6)') X11,V_m,V_p,eps_T,C11_m,C12_m,C44_m,C11_p,C12_p,C44_p,mobility,kappa,L0_Si_Sn_liq,L1_Si_Sn_liq,L2_Si_Sn_liq,L0_Si_Sn,L1_Si_Sn,L2_Si_Sn
        
    !! **************************
    !! **************************   
    !! **************************
    
    IPC    = 0
    itimes = 0
    TIME   = 0.0D0
    chem_pot = 0.0D0
    elas_pot = 0.0D0
    !Vf_sn    = 0.0D0
    Vf_si    = 0.0D0
    Vf_old_si= 0.0D0

    !! ************* geometry
    !! ************* geometry
    !! ************* geometry
    
    call init_microstructure(phi)
    CALL pixelout(phi,itimes)
    
    call recip_lat(g , Kxx, Kyy , k2x, k2y);

    !! **************************
    !! **************************
    !! **************************

    !! **************************
    !! **************************
    !! **************************

        !! *********** WRITTING OUTPUT ***************
        WRITE(53,*) 'ZONE ','I=',nex,'J= ',ney
        DO i=1,nex
            DO j=1,ney
                    WRITE(53,771) i,j,k2x(i,j),k2y(i,j)
            END DO
        END DO
        CLOSE(53)
771     Format(I3,1X,i3,2X,5(1X,E16.7))
        !! *********************************************

    write(55,*) 'itimes     IPC      TIME     MAXVAL(phi)    MINVAL(phi)'
    write(55,'(2I6,2X,5E14.5)') itimes,IPC,TIME,MAXVAL(phi),MINVAL(phi)

    !! ************* CALCULATING
    do while (IPC < 200 )

        call chemical_pot(phi, chem_pot, int_pot, Kxx, Kyy ,k2x, k2y)
        call elastic_pot (phi(:,:,1), phi(:,:,2), elas_pot, sigma_el, epslon_el)
        call evolution   (phi, int_pot, chem_pot, elas_pot, sigma_el, k2x,k2y,g)

        IF ( MOD(itimes,wrt_cycle).EQ.0 ) THEN

            write(55,'(2I6,2X,5E14.5)') itimes,IPC,TIME,MAXVAL(phi),MINVAL(phi),Vf_si
            CALL printdata(phi,chem_pot,elas_pot,int_pot,sigma_el, epslon_el)

            !!! ****** !!!

            IPC  = IPC + 1 ;
                      
        ENDIF

        itimes = itimes + 1;
        TIME   = TIME + dt;

    end do

    IPC  = IPC  + 1 ;
    write(55,'(2I6,2X,5E14.5)') itimes,IPC,TIME,MAXVAL(phi),MINVAL(phi),1.0-Vf_si,Vf_si
    CALL printdata(phi,chem_pot,elas_pot,int_pot,sigma_el, epslon_el)

    CALL CPU_TIME(time_end)

    print *, 'SIMULATION FINISHED... GO HOME!!!'
    WRITE (*,*) 'Time of operation was ', time_end - time_begin, ' seconds'
    stop

contains

    !=================================================================

    subroutine init_microstructure(phi)

        use mod_geometry
        use IFPORT

        implicit none
        integer::i,j,K,L,Reason
        real::ranum

        REAL(8), INTENT(OUT), DIMENSION(1:nex,1:ney,1:2) :: phi
        REAL(8), DIMENSION(1:nex,1:ney) :: num,noise
        INTEGER :: type,centers

        phitot = 0.0D0

        !***** INITIALIZATION ***********************************************
        !***** INITIALIZATION ***********************************************
        !***** INITIALIZATION ***********************************************

        !***** Random domain ***********************************************
        call seed(500) ! initialize
        DO i=1,nex,1
            DO j=1,ney,1

                num(I,J) = random(0)

                num(I,J) = (2.0D0*num(I,J) - 1.0D0)         ! Generate random number between -1 and +1
                NOISE(I,J) = (0.0050)*num(I,J)
                phi(i,j,1) = X11 + NOISE(I,J)

            ENDDO
        ENDDO
        call seed(200) ! initialize
        DO i=1,nex,1
            DO j=1,ney,1

                num(I,J) = random(0)
 
                num(I,J) = (2.0D0*num(I,J) - 1.0D0)         ! Generate random number between -1 and +1
                NOISE(I,J) = (0.0050)*num(I,J)
                phi(i,j,2) = X11 + NOISE(I,J)
            ENDDO
        ENDDO
        
!        !***** circle domain ******************************************
!        type    = 1;
!        centers = 1;
!        call ICs(phi(:,:,1), type, centers);
!        phi(:,:,2) = 1.0D0 - phi(:,:,1)

!        !***** sinusoidal domain ***************************************
!        type    = 3;
!        centers = 0;
!        call ICs(phi(:,:,1), type, centers);
!        !phi(:,:,2) = 1.0D0 - phi(:,:,1)
       
        
        
        !***** Write domain ********************************************       
        OPEN(unit=1,FILE=TRIM(fileplace)//'/microstructure/initial_microstructure.dat')
        WRITE(1,*) 'ZONE ','I=',nex,'J= ',ney
        DO i=1,nex
            DO j=1,ney
                    WRITE(1,'(2X,I3,2X,I3,2X,2E12.5)') i,j,phi(i,j,1),phi(i,j,2)
            END DO
        END DO
        CLOSE(1)

!        phitot = sum(sum(sum(phi,DIM=1),Dim=1),DIM=1)/size(phi,1)

153     format((F5.2))

        phitot = sum(sum(sum(phi,DIM=1),DIM=1),DIM=1)/size(phi,1)**2

    end subroutine init_microstructure

    !=======================================================================================================================

    SUBROUTINE ICs(f_xyz, type, centers)

        USE IFPORT
        IMPLICIT NONE
        INTEGER :: i,j,c
        REAL(8) :: cx,cy,cz,TX,TY,BX,BY,rand_number,radius
        INTEGER :: type,centers

        REAL(8), DIMENSION(1:nex,1:ney)     :: f_xyz
        REAL(8), DIMENSION(1:nex,1:ney)     :: func_xyz
        REAL(8), DIMENSION(1:nex,1:ney)     :: XX,YY
        REAL(8), DIMENSION(1:nex)           :: x1
        REAL(8), DIMENSION(1:ney)           :: y1
        REAL(8), DIMENSION(1:nex)           :: z1

!        x1 = [ ( real(i), i=(-int(nex/2)),-1 ) , 0.0 ,  ( real(i), i=1,int(nex/2)) ] !...........................define the domain discretization
!        y1 = [ ( real(i), i=(-int(nex/2)),-1 ) , 0.0 ,  ( real(i), i=1,int(nex/2)) ] !...........................define the domain discretization
!        z1 = [ ( real(i), i=(-int(nex/2)),-1 ) , 0.0 ,  ( real(i), i=1,int(nex/2)) ] !...........................define the domain discretization

!            x1(1:(nex))     = [ ( real(i), i=0,(nex/2-1)) , 0.0 , ( real(i), i=(-nex/2+1),-1 ) ]
!            y1(1:(nex))     = [ ( real(i), i=0,(nex/2-1)) , 0.0 , ( real(i), i=(-nex/2+1),-1 ) ]
!            z1(1:(nex))     = [ ( real(i), i=0,(nex/2-1)) , 0.0 , ( real(i), i=(-nex/2+1),-1 ) ]

        x1 = (/((i + 0.5D0),i=-nex/2,+nex/2-1)/)  !...........................define the domain discretization
        y1 = (/((i + 0.5D0),i=-ney/2,+ney/2-1)/)  !...........................define the domain discretization
        z1 = (/((i + 0.5D0),i=-nez/2,+nez/2-1)/)  !...........................define the domain discretization

        if (type == 1) then !cylinder

            call meshgrid_2D(x1, y1, XX, YY)     !......................................X, Y, Z field values
            f_xyz(1:nex,1:ney) = 0.0D0;            !......................................1-D array of f(x,y)

            If (centers == 1) then ! 1 circle

                WRITE(*,*) 'Initiating cylinder precipiate in the center...'

                DO i=1,nex
                    DO j=1,ney
                        if (DSQRT(XX(i,j)**2 + YY(i,j)**2) < max_radius) then
                            f_xyz(i,j) = 1.0;
                        elseif (DSQRT(XX(i,j)**2 + YY(i,j)**2) < (max_radius + grad_factor*dx) ) then
                            f_xyz(i,j) = (grad_factor*dx  - (DSQRT(XX(i,j)**2 + YY(i,j)**2) - real(max_radius)))/(grad_factor*dx);
                        else
                            f_xyz(i,j) = 0.0D0;
                        endif
                    enddo
                enddo

            elseif (centers == 2) then;                ! 2 circle

                DO j = 1,ney
                    Do i = 1,nex
                        if ( (sqrt((XX(i,j)+max_radius+sep_dis)**2 + YY(i,j)**2) < max_radius) .or. (sqrt((XX(i,j)-max_radius-sep_dis)**2 + YY(i,j)**2) < max_radius)) then
                            f_xyz(i,j) = 1.0;
                        elseif (sqrt((XX(i,j)+max_radius+sep_dis)**2 + YY(i,j)**2) < max_radius + grad_factor*dx) then
                            f_xyz(i,j) = (grad_factor*dx  - (sqrt((XX(i,j)+max_radius+sep_dis)**2 + YY(i,j)**2) - max_radius))/(grad_factor*dx);
                        elseif (sqrt((XX(i,j)-max_radius-sep_dis)**2 + YY(i,j)**2) < max_radius + grad_factor*dx) then
                            f_xyz(i,j) = (grad_factor*dx  - (sqrt((XX(i,j)-max_radius-sep_dis)**2 + YY(i,j)**2) - max_radius))/(grad_factor*dx);
                        else
                            f_xyz(i,j) = 0.0;
                        endif
                    enddo
                enddo

                phitot = sum(sum(f_xyz,DIM=1),DIM=1)/size(f_xyz,1)**2


            else

                !call srand(500);
                CALL RANDOM_SEED()

                Do c = 1,centers
                    radius = rand(0)*max_radius;
                    cx = 1.9*(rand(0)-0.5)*(Lx/2. - radius*(2./3.) - 5.*grad_factor*dx);
                    cy = 1.9*(rand(0)-0.5)*(Ly/2. - radius*(2./3.) - 5.*grad_factor*dy);
                    bx = ceiling(nex/2.) + ceiling(real(cx/dx)) - ceiling(real(radius/dx)) - grad_factor;
                    tx = ceiling(nex/2.) + ceiling(real(cx/dx)) + ceiling(real(radius/dx)) + grad_factor;
                    by = ceiling(ney/2.) + ceiling(real(cy/dy)) - ceiling(real(radius/dy)) - grad_factor;
                    ty = ceiling(ney/2.) + ceiling(real(cy/dy)) + ceiling(real(radius/dy)) + grad_factor;

                    !write(*,*) radius
                    !write(*,*) bx,tx,by,ty
                    !pause

                    Do i = by,ty
                        Do j = bx,tx
                            if (sqrt((XX(i,j)-cx)**2 + (YY(i,j)-cy)**2) < radius) then
                                func_xyz(i,j) = 1.0;
                            elseif (sqrt((XX(i,j)-cx)**2 + (YY(i,j)-cy)**2) < radius + grad_factor*dx) then
                                if ( func_xyz(i,j) > (grad_factor*dx  - (sqrt((XX(i,j)-cx)**2 + (YY(i,j)-cy)**2) - radius))/(grad_factor*dx) ) then
                                    func_xyz(i,j) = func_xyz(i,j);
                                else
                                    func_xyz(i,j) = (grad_factor*dx  - (sqrt((XX(i,j)-cx)**2 + (YY(i,j)-cy)**2) - radius))/(grad_factor*dx);
                                endif
                            else
                               !func_xyz(i,j) = 0.1;
                            endif

                            f_xyz(i,j) = func_xyz(i,j)

                        enddo
                    enddo

                enddo

            endif


        elseif (type == 2) then !sphere

            call meshgrid_2D(x1, y1, XX, YY)               !......................................X, Y, Z field values
            func_xyz(1:nex,1:ney) = 0.0D0;            !......................................1-D array of f(x,y,z)

            Do i = 1,nex
                Do j = 1,ney

                    if (dsqrt(XX(i,j)**2 + YY(i,j)**2) < max_radius) then
                        func_xyz(i,j) = 1.0D0;
                    elseif ( (dsqrt(XX(i,j)**2 + YY(i,j)**2)) < max_radius + grad_factor*dx) then
                        func_xyz(i,j) = (grad_factor*dx  - (sqrt(XX(i,j)**2 + YY(i,j)**2 ) - max_radius))/(grad_factor*dx);
                    else
                        func_xyz(i,j) = 0.0D0;
                    endif

                    f_xyz(i,j) = func_xyz(i,j)

                enddo
            enddo

        elseif (type == 3) then !sphere

	        x1 = (/((i ),i=1,nex)/)  !...........................define the domain discretization
	        y1 = (/((i ),i=1,ney)/)  !...........................define the domain discretization
	        z1 = (/((i ),i=1,nez)/)  !...........................define the domain discretization
	
	            call meshgrid_2D(x1, y1, XX, YY)     !......................................X, Y, Z field values
	            f_xyz(1:nex,1:ney) = 0.0D0;          !......................................1-D array of f(x,y)
	
	            Do i = 1,nex
	                Do j = 1,ney
	
						f_xyz(i,j)=DSIN(2.0*pi*XX(i,j)/nex)+0.001*DCOS(16.0*pi*XX(i,j)/nex);
	
	                enddo
	            enddo
	            
        endif

        RETURN

    END SUBROUTINE ICs

    !=================================================================

    subroutine evolution(phi,int_pot,chem_pot,elas_pot,sigma_el,k2x,k2y,g)
        use mod_geometry
        use mtls_property
        use wrt_opts
        use solver_opts, only : dt
        use math_opts

        ! output
        real(8),  dimension(1:nex,1:ney,1:2) :: phi

        ! input
        real(8), dimension(1:nex,1:ney,1:2)   :: chem_pot
        real(8), dimension(1:nex,1:ney,1:2)   :: int_pot
        real(8), dimension(1:nex,1:ney,1:2)   :: elas_pot
        real(8), dimension(1:nex,1:ney)       :: k2x,k2y
        real(8), dimension(1:3)               :: g

        ! temp.
        complex(8), dimension(1:nex,1:ney,1:2) :: phi_temp
        complex(8), dimension(1:nex,1:ney,1:2) :: dphi_temp

        complex(8), dimension(1:nex,1:ney,1:2) :: h_f
        complex(8), dimension(1:nex,1:ney,1:2) :: elas_pot_temp

        real(8), dimension(1:nex,1:ney) :: lap_chem_pot
        real(8), dimension(1:nex,1:ney) :: lap_elas_pot

        REAL(8), DIMENSION(1:nex,1:ney,1:6):: epslon_el, sigma_el


        !! code begins here .............

        phi_temp(:,:,1) 	 = phi(:,:,1)
        phi_temp(:,:,2)      = phi(:,:,2)
        h_f(:,:,1)           = chem_pot(:,:,1)
        h_f(:,:,2)           = chem_pot(:,:,2)
        elas_pot_temp(:,:,1) = elas_pot(:,:,1)
        elas_pot_temp(:,:,2) = elas_pot(:,:,2)

        call twoD_fft(phi_temp(:,:,1),'Forward ')
        call twoD_fft(phi_temp(:,:,2),'Forward ')
        call twoD_fft(h_f(:,:,1),'Forward ')
        call twoD_fft(h_f(:,:,2),'Forward ')
        call twoD_fft(elas_pot_temp(:,:,1),'Forward ')
        call twoD_fft(elas_pot_temp(:,:,2),'Forward ')

        !!!!!!!!!!!

        phitot  = 0.0D0

        !! Explicit- Implicit method	
		!! Implicit/Explicit timestepping
!		phi_temp(:,:,1)=(phi_temp(:,:,1)*(1.0+1.0/dt) - h_f(:,:,1))/ (-(k2x+k2y)*beta+1.0/dt); 
		phi_temp(:,:,1)=( phi_temp(:,:,1)/dt +  h_f(:,:,1) + elas_pot_temp(:,:,1) ) / (-(k2x+k2y)*beta+1.0/dt) !             *(1.0+1.0/dt) - h_f(:,:,1))/ (-(k2x+k2y)*beta+1.0/dt); 
		phi_temp(:,:,2)=( phi_temp(:,:,2)/dt +  h_f(:,:,2) + elas_pot_temp(:,:,2) ) / (-(k2x+k2y)*beta+1.0/dt) !             *(1.0+1.0/dt) - h_f(:,:,1))/ (-(k2x+k2y)*beta+1.0/dt); 

		!! converts to real space in x-direction
        call twoD_fft(phi_temp(:,:,1),'Backward')
        call twoD_fft(phi_temp(:,:,2),'Backward')
        phi(:,:,1) = real( phi_temp(:,:,1) )
        phi(:,:,2) = real( phi_temp(:,:,2) )
        		
        phitot = sum(sum(phi(:,:,1),DIM=1),DIM=1)/(size(phi(:,:,1),1)**2)

        if(isnan(phitot)) then
            write(*,*) 'NaN found...'
            write(*,*) 'itimes: ',itimes
            stop! 'NaN found'
        endif


    end subroutine evolution
    !==================================================================

     subroutine printdata(phi,chem_pot,elas_pot,int_pot, sigma_el, epslon_el)
        use mod_geometry
        use mtls_property
        use wrt_opts

        IMPLICIT NONE

        ! input
        real(8), intent(in), dimension(1:nex,1:ney,1:2) :: phi
        ! input
        real(8), dimension(1:nex,1:ney,1:2) :: chem_pot
        real(8), dimension(1:nex,1:ney,1:2) :: elas_pot
        real(8), dimension(1:nex,1:ney,1:2) :: int_pot
        REAL(8), DIMENSION(1:nex,1:ney,1:6):: epslon_el, sigma_el

        INTEGER :: i,j,k

!        IF(itimes == 0) then
!        
!        
!
!            WRITE(*,*)
!            WRITE(*,'(A20,1X,5F6.2)') 'Alloy composition =',X11,1.0-X11
!            WRITE(*,*) 'Grid spacing                     =',dx
!            WRITE(*,*) '# of grids                       =',nex
!            WRITE(*,*) 'K-Space Total Length             =',g(1)*nex
!            WRITE(*,*) 'R-Space Total Length             =',Lx
!            WRITE(*,*) 'a_m:',r_m,'a_p:',r_p
!            WRITE(*,*) 'eps_T                            =',eps_T
!            WRITE(*,*)
!
!            OPEN(unit=1002,FILE=TRIM(fileplace)//'/results/screen_output.dat')
!
!            WRITE(*,*)
!            WRITE(*,*) ' ****** Composition Parameters ************ '
!            WRITE(*,*) '1: ', X11
!
!            WRITE(*,*) ' ****** Physical Parameters *************** '
!            WRITE(*,*) '2: ', V_m
!            WRITE(*,*) '3: ', V_p
!
!            WRITE(*,*) ' ****** Elastic Parameters **************** '
!            WRITE(*,*) '4: ', eps_T
!            WRITE(*,*) ' 5-6-7: Mg2Si  ---- 8-9-10:Mg2Sn '
!            WRITE(*,*) '5: ', C11_m
!            WRITE(*,*) '6: ', C12_m
!            WRITE(*,*) '7: ', C44_m
!            WRITE(*,*) '8: ', C11_p
!            WRITE(*,*) '9: ', C12_p
!            WRITE(*,*) '10:', C44_p
!
!            WRITE(*,*) ' ****** Kinetic Parameters **************** '
!            WRITE(*,*) '11:', mobility
!            WRITE(*,*) '12:', kappa
!
!            WRITE(*,*) ' ****** Thermodynamic Parameters ********** '
!            WRITE(*,*) '13:', L0_Si_Sn_liq
!            WRITE(*,*) '14:', L1_Si_Sn_liq
!            WRITE(*,*) '15:', L2_Si_Sn_liq
!            WRITE(*,*) '16:', L0_Si_Sn
!            WRITE(*,*) '17:', L1_Si_Sn
!            WRITE(*,*) '18:', L2_Si_Sn
!            WRITE(*,*)
!
!        ENDIF

        !! *********** WRITTING OUTPUT ***************
        WRITE(num,'(I9.9)') itimes
        OPEN(unit=52,FILE=TRIM(fileplace)//'/microstructure/phi'//TRIM(num)//'.dat')
155     FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"phi1"',2X,'"phi2"',2X,'"chem_pot1"',2X,'"chem_pot2"',2X,'"int_pot1"',2X,'"int_pot2"',2X,'"elas_pot1"',2X,'"elas_pot2"')
        WRITE(52,155)
        WRITE(52,770) IG, JG
        DO i=1,nex
            DO j=1,ney
                    WRITE(52,771) i,j,phi(i,j,1),phi(i,j,2),chem_pot(i,j,1),chem_pot(i,j,2),int_pot(i,j,1),int_pot(i,j,2),elas_pot(i,j,1),elas_pot(i,j,2)
            END DO
        END DO

        CLOSE(52)
770     FORMAT('ZONE',2X,'I=',I5,2X,'J=',I5,2X)       
771     Format(I3,1X,i3,2X,10(1X,E16.7))
        !! *********************************************
		
        !! *********** WRITTING IMAGE FILE *************
		CALL pixelout(phi,ipc)
        !! *********** WRITTING IMAGE FILE *************

        !! ************* DISPLAY OUTPUT ****************
        WRITE(*,*) '******Calculating...******'
        WRITE(*,'(1X,A4,I9,1X,A4,E12.5,1X,A4,E12.5,1X,A8,I9)') 'IPC=',IPC,'Dx=',dx,'DT=',dt,'itimes=',itimes
        WRITE(*,*) 'TIME=',TIME
        WRITE(*,*) 'PHITOT=',phitot
        WRITE(*,*) 'PHIMAX=',MAXVAL(PHI),'PHIMIN=',MINVAL(PHI)
        WRITE(*,*) '******Material Pro...******'
        !WRITE(*,*) 'Nv_m =',Nv_m ,'Nv_p=',Nv_p
        WRITE(*,*) 'beta=',beta,'L   =',L
        WRITE(*,*) '******Chemical Pro...******'
        WRITE(*,*) 'MAX_chem1     =',MAXVAL(chem_pot(:,:,1)) ,'MIN_chem1     =',MINVAL(chem_pot(:,:,1))
        WRITE(*,*) '******Interface Pro...******'
        WRITE(*,*) 'MAX_int1      =',MAXVAL(int_pot(:,:,1))  ,'MIN_int1      =',MINVAL(int_pot(:,:,1))
        WRITE(*,*) '******Elastic Pro...******'
        WRITE(*,*) 'MAX_elas1     =',MAXVAL(elas_pot(:,:,1)) ,'MIN_elas1     =',MINVAL(elas_pot(:,:,1))
        WRITE(*,*) '******Elastic stress and strain...******'
        WRITE(*,*) 'MAX_eps_el_11=',MAXVAL(epslon_el(:,:,1)),'MIN_eps_el_11=',MINVAL(epslon_el(:,:,1))
        WRITE(*,*) 'MAX_sig_el_11=',MAXVAL(sigma_el (:,:,1)),'MIN_sig_el_11=',MINVAL(sigma_el(:,:,1))
        WRITE(*,*) 'MAX_eps_el_22=',MAXVAL(epslon_el(:,:,2)),'MIN_eps_el_22=',MINVAL(epslon_el(:,:,2))
        WRITE(*,*) 'MAX_sig_el_22=',MAXVAL(sigma_el (:,:,2)),'MIN_sig_el_22=',MINVAL(sigma_el(:,:,2))
        WRITE(*,*) 'dx   =',dx
        WRITE(*,*)
        !! *********************************************

        !! ************* DISPLAY OUTPUT ****************
        WRITE(1002,*) '******Calculating...******'
        WRITE(1002,*) 'IPC =',IPC,'DT=',dt,'itimes=',itimes
        WRITE(1002,*) 'TIME=',TIME
        WRITE(1002,*) 'PHITOT=',phitot
        WRITE(1002,*) 'PHIMAX=',MAXVAL(PHI),'PHIMIN=',MINVAL(PHI)
        WRITE(1002,*) '******Material Pro...******'
        !WRITE(1002,*) 'Nv_m =',Nv_m ,'Nv_p=',Nv_p
        WRITE(1002,*) 'beta=',beta,'L   =',L
!        WRITE(1002,*) 'MAX_chem   =',MAXVAL(chem_pot) ,'MIN_chem     =',MINVAL(chem_pot)
!        WRITE(1002,*) 'MAX_chem   =',MAXVAL(chem_pot) ,'MIN_chem     =',MINVAL(chem_pot)
!        WRITE(1002,*) 'MAX_elas   =',MAXVAL(elas_pot) ,'MIN_elas     =',MINVAL(elas_pot)
!        WRITE(1002,*) 'MAX_elas   =',MAXVAL(elas_pot) ,'MIN_elas     =',MINVAL(elas_pot)
        WRITE(1002,*) '******Elastic stress and strain...******'
        WRITE(1002,*) 'MAX_eps_el_11=',MAXVAL(epslon_el(:,:,1)),'MIN_eps_el_11=',MINVAL(epslon_el(:,:,1))
        WRITE(1002,*) 'MAX_sig_el_11=',MAXVAL(sigma_el (:,:,1)),'MIN_sig_el_11=',MINVAL(sigma_el(:,:,1))
        WRITE(1002,*) 'dx   =',dx
        WRITE(1002,*)

        RETURN

    end subroutine printdata

    !! ------

    SUBROUTINE recip_lat(g, Kx, Ky, k2x, k2y)

        use math_opts
        IMPLICIT NONE

        !output
        REAL(8), INTENT(out), DIMENSION(1:nex,1:ney) :: Kx,Ky,k2x,k2y
        REAL(8), INTENT(out), DIMENSION(1:3)         :: g

        !Temp.
        INTEGER :: i
        REAL(8), DIMENSION(1:3)                            :: g1,g2,g3
        REAL(8), DIMENSION(1:nex) :: k_1Dx
        REAL(8), DIMENSION(1:ney) :: k_1Dy
        REAL(8), DIMENSION(1:nex) :: k_1Dz
        REAL(8), DIMENSION(1:3)   :: x_hat,y_hat,z_hat
        REAL(8), DIMENSION(1:nex,1:ney) :: X,Y
        REAL(8), DIMENSION(1:nex,1:ney,1:nez) :: XX,YY,ZZ

        k_1Dx(1:nex) = 0.0D0
        k_1Dy(1:ney) = 0.0D0
        k_1Dz(1:nez) = 0.0D0

        if (dimen == 1) THEN

        !! initing the k space
        k_1Dx(1:(nex/2))     = (/ ( real(i), i=0,(nex/2-1)) /)
        k_1Dx((nex/2+1):nex) = (/ ( real(i), i=(-nex/2),-1 ) /)

            !x_hat(1) = 1.0D0;
            !g(1)     = (2.*pi*x_hat(1))/Lx;
            !K1(1:nex,1,1) = k_1Dx;

        elseif (dimen == 2) THEN

            !! initing the k space
            k_1Dx(1:(nex))     = [ ( real(i), i=0,(nex/2-1)) , 0.0 , ( real(i), i=(-nex/2)+1,-1 ) ]
            k_1Dy(1:(ney))     = [ ( real(i), i=0,(ney/2-1)) , 0.0 , ( real(i), i=(-ney/2)+1,-1 ) ]
            k_1Dz(:)           = k_1Dx(:);

            x_hat = (/1.0, 0.0, 0.0/);
            y_hat = (/0.0, 1.0, 0.0/);
            z_hat = (/0.0, 0.0, 1.0/);

            g1 = 2.0*pi*cross((Ly*y_hat),(Lz*z_hat))/(DOT_PRODUCT(Lx*x_hat,cross((Ly*y_hat),(Lz*z_hat))));
            g2 = 2.0*pi*cross((Lz*z_hat),(Lx*x_hat))/(DOT_PRODUCT(Ly*y_hat,cross((Lz*z_hat),(Lx*x_hat))));
            g3 = 2.0*pi*cross((Lx*x_hat),(Ly*y_hat))/(DOT_PRODUCT(Lz*z_hat,cross((Lx*x_hat),(Ly*y_hat))));
            g = g1 + g2 + g3;

            call meshgrid_2D(k_1Dx, k_1Dy, X, Y)

            Kx(1:nex,1:ney) = g(1)*X(1:nex,1:ney);
            Ky(1:nex,1:ney) = g(2)*Y(1:nex,1:ney);

            k2x(1:nex,1:ney) = (DCMPLX(0.0,1.0)*Kx(1:nex,1:ney))**2
            k2y(1:nex,1:ney) = (DCMPLX(0.0,1.0)*Ky(1:nex,1:ney))**2

            !k2x(1:nex,1:ney) = (2.0*pi*DCMPLX(0.0,1.0)*Kx(1:nex,1:ney))**2
            !k2y(1:nex,1:ney) = (2.0*pi*DCMPLX(0.0,1.0)*Ky(1:nex,1:ney))**2

        elseif (dimen == 3) THEN


        endif

        RETURN

    END SUBROUTINE recip_lat

end program main

