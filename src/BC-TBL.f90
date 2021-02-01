module tbl

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  real(mytype),save, allocatable, dimension(:) :: ttd1a,ttd1b,ttd1c
  real(mytype),save, allocatable, dimension(:) :: ttd2a,ttd2b,ttd2c
  real(mytype),save, allocatable, dimension(:) :: ttd3a,ttd3b,ttd3c


  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_tbl, boundary_conditions_tbl, postprocess_tbl, tbl_flrt, momentum_forcing_tbl


contains

  subroutine init_tbl (ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,fh,ierror,ii,is,it,code
    integer (kind=MPI_OFFSET_KIND) :: disp

    integer, dimension (:), allocatable :: seed

    allocate(ttd1a(ysize(2)),ttd1b(ysize(2)),ttd1c(ysize(2)))
    allocate(ttd2a(ysize(2)),ttd2b(ysize(2)),ttd2c(ysize(2)))
    allocate(ttd3a(ysize(2)),ttd3b(ysize(2)),ttd3c(ysize(2)))

    if (iscalar==1) then

       phi1(:,:,:,:) = 0.25 !change as much as you want
          if ((nclyS1.eq.2).and.(xstart(2).eq.1)) then
             !! Generate a hot patch on bottom boundary
             phi1(:,1,:,:) = one
          endif
          if ((nclySn.eq.2).and.(xend(2).eq.ny)) THEN
             phi1(:,xsize(2),:,:) = 0.25
          endif

    endif
    ux1=zero;uy1=zero;uz1=zero

    !a blasius profile is created in ecoule and then duplicated for the all domain
    if(nclx1==2)call blasius !spatial framework
    if(nclx1==0)call perturb_init_tbl(ux1,uy1,uz1) !temporal framework 

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_tbl
  !********************************************************************
  subroutine boundary_conditions_tbl (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype) :: x, y, z, um
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx

    integer :: i, j, k, is

    !INFLOW with an update of bxx1, byy1 and bzz1 at the inlet

    if (nclx1 == 2) THEN

    call blasius()
    !INLET FOR SCALAR, TO BE CONSISTENT WITH INITIAL CONDITION
    if (iscalar.eq.1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(1,:,:,:)=0.25
             if ((xstart(2).eq.1)) then
                phi(:,1,:,:) = one
             endif
             if ((xend(2).eq.ny)) THEN
                phi(:,xsize(2),:,:) = 0.25
             endif
          enddo
       enddo
    endif

    !OUTFLOW based on a 1D convection equation

    udx=one/dx
    udy=one/dy
    udz=one/dz
    uddx=half/dx
    uddy=half/dy
    uddz=half/dz

    do k=1,xsize(3)
       do j=1,xsize(2)

          cx=ux(nx,j,k)*gdt(itr)*udx

          if (cx.LT.0.0) cx=0.0
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
          if (iscalar.eq.1) phi(nx,:,:,:) =  phi(nx,:,:,:) - cx*(phi(nx,:,:,:)-phi(nx-1,:,:,:))
          enddo
    enddo

   endif !if (nclx1 == 2)


    !! Bottom Boundary
    if (ncly1 == 2) THEN
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byx1(i, k) = zero
          byy1(i, k) = zero
          byz1(i, k) = zero
        enddo
      enddo
    endif
    !! Top Boundary
    if (nclyn == 2) then
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byxn(i, k) = ux(i, xsize(2) - 1, k)
             byyn(i, k) = uy(i, xsize(2) - 1, k)
             byzn(i, k) = uz(i, xsize(2) - 1, k)
          enddo
       enddo
    !update of the flow rate (what is coming in the domain is getting out)
    call tbl_flrt(ux,uy,uz)
    endif

    !SCALAR ! 
    if (itimescheme.ne.7) then
     if (iscalar.ne.0) then
          if ((nclyS1.eq.2).and.(xstart(2).eq.1)) then
             !! Generate a hot patch on bottom boundary
             phi(1,1,:,:) = one
          endif
          if ((nclySn.eq.2).and.(xend(2).eq.ny)) THEN
             phi(1,xsize(2),:,:) = phi(1,xsize(2)-1,:,:)
          endif
     endif
    endif

    return
  end subroutine boundary_conditions_tbl

  !********************************************************************

  !********************************************************************
!
subroutine tbl_flrt (ux1,uy1,uz1)
  !
  !********************************************************************

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE MPI

  implicit none
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2

  integer :: j,i,k,code
  real(mytype) :: can,ut1,ut2,ut3,ut4,utt1,utt2,utt3,utt4,udif

  ux1(1,:,:)=bxx1(:,:)
  ux1(nx,:,:)=bxxn(:,:)

  call transpose_x_to_y(ux1,ux2)
  call transpose_x_to_y(uy1,uy2)
  ! Flow rate at the inlet
  ut1=zero;utt1=zero
  if (ystart(1)==1) then !! CPUs at the inlet
    do k=1,ysize(3)
      do j=1,ysize(2)-1
        ut1=ut1+(yp(j+1)-yp(j))*(ux2(1,j+1,k)-half*(ux2(1,j+1,k)-ux2(1,j,k)))
      enddo
    enddo
    ! ut1=ut1/real(ysize(3),mytype)
  endif
  call MPI_ALLREDUCE(ut1,utt1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  utt1=utt1/real(nz,mytype) !! Volume flow rate per unit spanwise dist
  ! Flow rate at the outlet
  ut2=zero;utt2=zero
  if (yend(1)==nx) then !! CPUs at the outlet
    do k=1,ysize(3)
      do j=1,ysize(2)-1
        ut2=ut2+(yp(j+1)-yp(j))*(ux2(ysize(1),j+1,k)-half*(ux2(ysize(1),j+1,k)-ux2(ysize(1),j,k)))
      enddo
    enddo
    ! ut2=ut2/real(ysize(3),mytype)
  endif

  call MPI_ALLREDUCE(ut2,utt2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  utt2=utt2/real(nz,mytype) !! Volume flow rate per unit spanwise dist

  ! Flow rate at the top and bottom
  ut3=zero
  ut4=zero
  do k=1,ysize(3)
    do i=1,ysize(1)
      ut3=ut3+uy2(i,1,k)
      ut4=ut4+uy2(i,ny,k)
    enddo
  enddo
  call MPI_ALLREDUCE(ut3,utt3,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE(ut4,utt4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  utt3=utt3/(real(nx*nz,mytype))*xlx  !!! Volume flow rate per unit spanwise dist
  utt4=utt4/(real(nx*nz,mytype))*xlx  !!! Volume flow rate per unit spanwise dist

  !! velocity correction
  udif=(utt1-utt2+utt3-utt4)/yly
  if (nrank==0 .and. mod(itime,1)==0) then
    write(*,"(' Mass balance: L-BC, R-BC,',2f12.6)") utt1,utt2
    write(*,"(' Mass balance: B-BC, T-BC, Crr-Vel',3f11.5)") utt3,utt4,udif
  endif
  ! do k=1,xsize(3)
  !   do j=1,xsize(2)
  !     ux1(nx,i,k)=ux1(nx,i,k)+udif
  !   enddo
  ! enddo
  do k=1,xsize(3)
    do j=1,xsize(2)
      bxxn(j,k)=bxxn(j,k)+udif
    enddo
  enddo


  return
end subroutine tbl_flrt

!********************************************************************
  subroutine blasius()


    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype) :: eta_bl, f_bl, g_bl, x_bl,h_bl
    real(mytype) :: delta_int, delta_eta, eps_eta


    real(mytype) :: x, y, z
    integer :: i, j, k, is

    do k=1,xsize(3)
       do j=1,xsize(2)
          if (istret.eq.0) y=(j+xstart(2)-1-1)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)

          eta_bl=y*4.91/9.0

          !OLD POLYNOMIAL FITTING

          delta_eta=0.0
          eps_eta=0.0
          delta_int=0.2

          if (eta_bl .ge. (7.5/9.0)) then
             delta_eta=eta_bl-7.5/9.0
             eta_bl=7.5/9.0
             eps_eta=0.00015
          end if

          f_bl=1678.64209592595*eta_bl**14-11089.6925017429*eta_bl**13 &
               +31996.4350140670*eta_bl**12-52671.5249779799*eta_bl**11 &
               +54176.1691167667*eta_bl**10-35842.8204706097*eta_bl**9  &
               +15201.3088871240*eta_bl**8 -4080.17137935648*eta_bl**7  &
               +702.129634528103*eta_bl**6 -56.2063925805318*eta_bl**5  &
               -17.0181128273914*eta_bl**4 +0.819582894357566*eta_bl**3  &
               -0.0601348202321954*eta_bl**2 +2.98973991270405*eta_bl**1

          f_bl=f_bl+(1-exp(-delta_eta/delta_int))*eps_eta


          if (eta_bl .ge. (7.15/9.0)) then
             delta_int=0.8
             delta_eta=eta_bl-7.15/9.0
             eta_bl=7.15/9.0
             eps_eta=0.0005
          end if

          g_bl=4924.05284779754*eta_bl**14-34686.2970972733*eta_bl**13 &
               +108130.253843618*eta_bl**12-195823.099139525*eta_bl**11 &
               +227305.908339065*eta_bl**10-176106.001047617*eta_bl**9  &
               +92234.5885895112*eta_bl**8 -32700.3687158807*eta_bl**7  &
               +7923.51008739107*eta_bl**6 -1331.09245288739*eta_bl**5  &
               +130.109496961069*eta_bl**4 -7.64507811014497*eta_bl**3  &
               +6.94303207046209*eta_bl**2 -0.00209716712558639*eta_bl**1 ! &

          g_bl=g_bl+(1-exp(-delta_eta/delta_int))*eps_eta


          x_bl=1.0/(4.91**2*xnu)

          bxx1(j,k)=f_bl/1.0002014996204402/1.0000000359138641 !To assure 1.0 in infinity
          bxy1(j,k)=g_bl*sqrt(xnu/x_bl)/1.000546554
          bxz1(j,k)=0.0

       enddo
    enddo

    !STORE VALUE F_BL_INF G_BL_INF (ONLY ONE MORE TIME)------------------

    y=yly
    eta_bl=y*4.91/9.0  !The 9 is due to interpolation

    delta_eta=0.0
    eps_eta=0.0
    delta_int=0.2

    if (eta_bl .ge. (7.5/9.0)) then
       delta_eta=eta_bl-7.5/9.0
       eta_bl=7.5/9.0
       eps_eta=0.00015
    end if

    f_bl_inf=1678.64209592595*eta_bl**14-11089.6925017429*eta_bl**13 &
         +31996.4350140670*eta_bl**12-52671.5249779799*eta_bl**11 &
         +54176.1691167667*eta_bl**10-35842.8204706097*eta_bl**9  &
         +15201.3088871240*eta_bl**8 -4080.17137935648*eta_bl**7  &
         +702.129634528103*eta_bl**6 -56.2063925805318*eta_bl**5  &
         -17.0181128273914*eta_bl**4 +0.819582894357566*eta_bl**3  &
         -0.0601348202321954*eta_bl**2 +2.98973991270405*eta_bl**1


    f_bl_inf=f_bl_inf+(1-exp(-delta_eta/delta_int))*eps_eta
    f_bl_inf=f_bl_inf/1.0002014996204402/1.0000000359138641 !To assure 1.0 in infinity

    if (eta_bl .ge. (7.15/9.0)) then
       delta_int=0.8
       delta_eta=eta_bl-7.15/9.0
       eta_bl=7.15/9.0
       eps_eta=0.0005
    end if

    g_bl_inf=4924.05284779754*eta_bl**14-34686.2970972733*eta_bl**13 &
         +108130.253843618*eta_bl**12-195823.099139525*eta_bl**11 &
         +227305.908339065*eta_bl**10-176106.001047617*eta_bl**9  &
         +92234.5885895112*eta_bl**8 -32700.3687158807*eta_bl**7  &
         +7923.51008739107*eta_bl**6 -1331.09245288739*eta_bl**5  &
         +130.109496961069*eta_bl**4 -7.64507811014497*eta_bl**3  &
         +6.94303207046209*eta_bl**2 -0.00209716712558639*eta_bl**1


    g_bl_inf=g_bl_inf+(1-exp(-delta_eta/delta_int))*eps_eta
    g_bl_inf=g_bl_inf/1.000546554

    return
  end subroutine blasius

  !############################################################################
  subroutine postprocess_tbl(ux1,uy1,uz1,ep1) !By Felipe Schuch

    USE MPI
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    character(len=30) :: filename

       !! Write vorticity as an example of post processing
    !x-derivatives
!    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
!    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    !y-derivatives
!    call transpose_x_to_y(ux1,td2)
!    call transpose_x_to_y(uy1,te2)
!    call transpose_x_to_y(uz1,tf2)
!    call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!    call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
!    call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    !!z-derivatives
!    call transpose_y_to_z(td2,td3)
!    call transpose_y_to_z(te2,te3)
!    call transpose_y_to_z(tf2,tf3)
!    call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
!    call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
!    call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    !!all back to x-pencils
!    call transpose_z_to_y(ta3,td2)
!    call transpose_z_to_y(tb3,te2)
!    call transpose_z_to_y(tc3,tf2)
!    call transpose_y_to_x(td2,tg1)
!    call transpose_y_to_x(te2,th1)
!    call transpose_y_to_x(tf2,ti1)
!    call transpose_y_to_x(ta2,td1)
!    call transpose_y_to_x(tb2,te1)
!    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

!    di1(:,:,:)=sqrt((tf1(:,:,:)-th1(:,:,:))**2+(tg1(:,:,:)-tc1(:,:,:))**2+&
!         (tb1(:,:,:)-td1(:,:,:))**2)
!    if (iibm==2) then
!       di1(:,:,:) = (one - ep1(:,:,:)) * di1(:,:,:)
!    endif
!    uvisu=0.
!    call fine_to_coarseV(1,di1,uvisu)
!994 format('vort',I3.3)
!    write(filename, 994) itime/ioutput
!    call decomp_2d_write_one(1,uvisu,filename,2)

    return
  end subroutine postprocess_tbl

  subroutine perturb_init_tbl(ta1,tb1,tc1)
   USE decomp_2d
   USE decomp_2d_io
   USE variables
   USE param
   USE MPI
   implicit none
   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1
   real(mytype) :: y,r,r3,x,z,h,aform
   integer :: k,j,i,ijk,fh,ierror,ii,is,code
   integer (kind=MPI_OFFSET_KIND) :: disp
   real(mytype), dimension(ysize(2)) :: um
   integer, dimension (:), allocatable :: seed

  !call system_clock(count=code)
   code=0 !fixing seed!
   call random_seed(size = ii)
   call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))
   call random_number(ta1)
   call random_number(tb1)
   call random_number(tc1)
   
   aform = 0.23369497_mytype
   do k=1,xsize(3)
      z=real((k+xstart(3)-1-1),mytype)*dz
      do j=1,xsize(2)
      if (istret.eq.0) y=(j+xstart(2)-1-1)*dy
      if (istret.ne.0) y=yp(j+xstart(2)-1)
      
         um = (one-erf(aform*y))
         do i=1,xsize(1)
      
         x=real(i-1,mytype)*dx
         ta1(i,j,k)=init_noise*um(j)*(two*ta1(i,j,k)-one)+erf(y*aform)+&
         4*init_noise*(y*exp(-y)/0.3678784468818499_mytype)*&
         ( cos(z*pi/five)*cos(x*pi/five) + &
            cos((x+((one+sqrt(five))*half))*pi/five)*cos((z+((one+sqrt(five))*half))*pi/five) )

         tb1(i,j,k)=init_noise*um(j)*(two*tb1(i,j,k)-one)
         tc1(i,j,k)=init_noise*um(j)*(two*tc1(i,j,k)-one)
         enddo
      enddo
   enddo

   return
  end subroutine perturb_init_tbl
    !############################################################################
  subroutine momentum_forcing_tbl(dux1,duy1,duz1,ux1,uy1,uz1,phi1)
    !
    !*******************************************************************************
  
      USE param
      USE variables
      USE decomp_2d
  
      implicit none
  
      real(mytype),dimension(xsize(1),xsize(2),xsize(3), ntime) :: dux1, duy1, duz1
      real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
      real(mytype),dimension(xsize(1),xsize(2),xsize(3), numscalar) :: phi1
  
      integer :: j
     
      ! !! BL Forcing (Pressure gradient or geostrophic wind)
      ! if (iPressureGradient.eq.1) then
      !    dux1(:,:,:,1)=dux1(:,:,:,1)+ustar**2./dBL
      ! else if (iCoriolis.eq.1.and.iPressureGradient.eq.0) then
      !    dux1(:,:,:,1)=dux1(:,:,:,1)+CoriolisFreq*(-UG(3))
      !    duz1(:,:,:,1)=duz1(:,:,:,1)-CoriolisFreq*(-UG(1))
      ! endif
  
      ! ! Coriolis terms
      ! if (iCoriolis.eq.1) then
      !    dux1(:,:,:,1)=dux1(:,:,:,1)+CoriolisFreq*uz1(:,:,:)
      !    duz1(:,:,:,1)=duz1(:,:,:,1)-CoriolisFreq*ux1(:,:,:)
      ! endif
  
      ! ! Damping zone
      ! if (idamping.eq.1) then
      !    call damping_zone(dux1,duy1,duz1,ux1,uy1,uz1)
      ! endif
  
      ! !! Buoyancy terms
      ! if (iscalar.eq.1.and.ibuoyancy.eq.1) then
      !    duy1(:,:,:,1)=duy1(:,:,:,1)+gravv*phi1(:,:,:,1)/Tref
      ! endif
  
      return
    end subroutine momentum_forcing_tbl
  subroutine comp_thetad(ux1,uy1,thetad)
  USE MPI
  USE decomp_2d
  USE decomp_2d_io
  ! USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
  ! USE var, only : uvisu
  ! USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  ! USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
  implicit none

  integer  :: i,j,k,code,istheta
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2p,dudy2,di2
  real(mytype),dimension(ysize(2)) :: ux2m,ux2m_old,uxuy2pm
  real(mytype),dimension(ysize(2)) :: ydudy2m,dudy2m,du2dy22m,duxuy2pm
  real(mytype) :: resi,resi1,cf_ref,delta
  real(mytype) :: ttheta,thetad,theta1,theta2
  real(mytype) :: theta3,thetad1,thetad2,thetad3
  real(mytype) :: y,tau_wall,tau_wall1,ufric,temp,utau
  real(8)      :: tstart
  tstart=MPI_WTIME()
  istheta = 10
 
    if(itime.eq.ifirst)then
     if(nrank.eq.0)print *,'first dy=',yp(2)
     thetad=zero
    endif
 
    call transpose_x_to_y(ux1,ux2)
    call horizontal_avrge(ux2,ux2m)
    ttheta = sum((ux2m*(one-ux2m))*ypw)
    delta = sum((one-ux2m)*ypw)
    resi = ttheta - one
    ux2m_old = ux2m
 
    if(mod(itime,istheta).eq.0)then
    call extract_fluctuat(ux2,ux2m,ux2p) ! ux'
    call transpose_x_to_y(uy1,uy2)
    call horizontal_avrge(ux2p*uy2,uxuy2pm) ! ux'uy' in profile
    call dery  (duxuy2pm,uxuy2pm ,di2,sy,ffy ,fsy ,fwy ,ppy  ,1,ysize(2),1,0) !0*1=0 
    call dery  (dudy2m  ,ux2m    ,di2,sy,ffyp,fsyp,fwyp,ppy  ,1,ysize(2),1,1)
    call dery  (ydudy2m ,ux2m*yp ,di2,sy,ffyp,fsyp,fwyp,ppy  ,1,ysize(2),1,1)
    !ydudy2m = ydudy2m - ux2m !ydudy = (duydy+u) - u = ydudy
    call deryy (du2dy22m,ux2m    ,di2,sy,sfyp,ssyp,swyp      ,1,ysize(2),1,1)
    if (istret.ne.0) then!correct the second derivative
     do j = 1,ysize(2)
      du2dy22m(j) = du2dy22m(j)*pp2y(j)-pp4y(j)*dudy2m(j)
     enddo
    endif

    !DEBUGGGGGGGGGGGGGGGGGGGGG
    !call outpdt(ux2m    ,'ux2m')
    !call outpdt(dudy2m  ,'dudy2m')
    !call outpdt(du2dy22m,'du2dy22m')
    !call outpdt(duxuy2pm,'duxuy2pm')
    !ux2m(:) = ux2m_old(:)+dt*(thetad1*yp(:)*dudy2m(:)+xnu*du2dy22m(:)-duxuy2pm(:))
    !theta1 = sum((ux2m*(one-ux2m))*ypw)
    !call outpdt(ux2m    ,'ux2m_theta1')
    !ux2m(:) = ux2m_old(:)+dt*(thetad2*yp(:)*dudy2m(:)+xnu*du2dy22m(:)-duxuy2pm(:))
    !theta2 = sum((ux2m*(one-ux2m))*ypw)
    !call outpdt(ux2m    ,'ux2m_theta2')
 
    if(itime.gt.ifirst+3)then !skipping Euler and AB2
     !if(mod(itime,3)==0)then
      if(nrank.eq.0)then
 
      !target 1.4e-3
      i=0;
 
      !if(thetad.eq.zero)then
       thetad1=1.0e-6
       thetad2=1.0e-1
      !else
      ! thetad1=thetad -istheta*dt
      ! if(thetad1.lt.zero)thetad1=zero
      ! thetad2=thetad +istheta*dt
      !endif     
 
      do while(abs(resi).gt.1.e-11)
      i=i+1!;theta1=zero;theta2=zero;theta3=zero
 
      if(i.gt.50)then
       print *,'exiting after 50',theta3,thetad3
       exit
      endif
 
 !phi1(ijk,1,1,is)=phi1(ijk,1,1,is)+adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1,is)+cdt(itr)*phiss1(ijk,1,1,is)
 !phiss1(ijk,1,1,is)=phis1(ijk,1,1,is)
 !phis1(ijk,1,1,is)=ta1(ijk,1,1)
 
    ttd1a(:) = thetad1*ydudy2m(:)+xnu*du2dy22m(:)-duxuy2pm(:)
    ux2m(:) = ux2m_old(:)+adt(1)*ttd1a(:)+bdt(1)*ttd1b(:)+cdt(1)*ttd1c(:)
    theta1 = sum((ux2m*(one-ux2m))*ypw)
 
    ttd2a(:) = thetad2*ydudy2m(:)+xnu*du2dy22m(:)-duxuy2pm(:)
    ux2m(:) = ux2m_old(:)+adt(1)*ttd2a(:)+bdt(1)*ttd2b(:)+cdt(1)*ttd2c(:)
    theta2 = sum((ux2m*(one-ux2m))*ypw)
 
    thetad3=half*(thetad1+thetad2)
    ttd3a(:) = thetad3*ydudy2m(:)+xnu*du2dy22m(:)-duxuy2pm(:)
    ux2m(:) = ux2m_old(:)+adt(1)*ttd3a(:)+bdt(1)*ttd3b(:)+cdt(1)*ttd3c(:)
    theta3 = sum((ux2m*(one-ux2m))*ypw)
 
    if( (theta1-one)*(theta3-one) .lt. zero)then
     thetad2 = thetad3
    else
     thetad1 = thetad3
    endif
 
    resi=theta1-one
    end do
 
    if(abs(theta1-one).lt.1.e-11)then
     print *,'Converged after ',i
     thetad=thetad3
    endif
 
     print *,'Time computing theta dot (s)', real(MPI_WTIME()-tstart,4)   ! USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean

     endif !nrank.eq.0
     call MPI_BCAST(thetad,1,real_type,0,MPI_COMM_WORLD,code)
     endif !Euler + AB2 skip
     endif !mod
 
   ttd1c(:)=ttd1b(:);ttd1b(:)=ttd1a(:)
   ttd2c(:)=ttd2b(:);ttd2b(:)=ttd2a(:)
   ttd3c(:)=ttd3b(:);ttd3b(:)=ttd3a(:)
 
   call dery(dudy2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   !call dery(dwdy2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   temp = sum(dudy2(:,1,:))!**2+dwdy2(:,1,:)**2)
   call MPI_ALLREDUCE(temp,tau_wall1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
   tau_wall=tau_wall1/real(nx*nz,mytype)
   ufric=sqrt(tau_wall*xnu)
 
   if (nrank .eq. 0) then
     if(mod(itime,itest)==0)then
      print *,'dx+=',real(dx   *ufric*re,4),'<8'
      print *,'dy+=',real(yp(2)*ufric*re,4),'<1'
      print *,'dz+=',real(dz   *ufric*re,4),'<5'
     endif
     
     print *,'theta,thetad=',ttheta,thetad
     print *,'cf   ,utau  =',2*ufric**2,ufric
 
     FS = 2 ;write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
     FS = FS*14+1  !Line width
     open(67,file='s_utau',status='unknown',form='formatted',access='direct',recl=FS)
     write(67,fileformat,rec=itime+1) t,ufric,NL ;close(67)
 
     open(67,file='s_thetad',status='unknown',form='formatted',access='direct',recl=FS)
     write(67,fileformat,rec=itime+1) t,thetad,NL ;close(67)
 
     open(67,file='s_delta',status='unknown',form='formatted',access='direct',recl=FS)
     write(67,fileformat,rec=itime+1) t,delta,NL ;close(67)
 
     if(itime.eq.ifirst) then
      write(fileformat, '( "(",I4,"(E16.8),A)" )' ) ny
      open(67,file='./out/00_yp',status='unknown',form='formatted',access='direct',recl=(ny*16+1))
      write(67,fileformat,rec=1) yp(:) ;close(67)
     endif
 
   endif
   return
   end
   !############################################################################
  !############################################################################
   subroutine outp(field,filename) !outpost to disk buffered quantities
      USE param
      USE variables
      implicit none
      real(mytype),intent(IN),dimension(ysize(2),iprocessing) :: field
      character(len=*), intent(IN) :: filename
      integer :: i
      character(len=1),parameter :: NL=char(10)
      character(len=300) :: fileformat
      write(fileformat, '( "(",I4,"(E16.8),A)" )' ) ny
      open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=(ny*16+1))
      do i=1, iprocessing/itest
      !if(nrank.eq.0)print *, 'itime/itest=',itime/itest-iprocessing/itest+i
        write(67,fileformat,rec=(itime-initstat)/itest -iprocessing/itest + i) field(:,i),NL
      enddo
      close(67)
    end subroutine outp
    !############################################################################
    subroutine outpdt(field,filename)
      USE param
      USE variables
      implicit none
      real(mytype),intent(IN),dimension(ysize(2)) :: field
      character(len=*), intent(IN) :: filename
      integer :: i
      character(len=1),parameter :: NL=char(10)
      character(len=300) :: fileformat
      write(fileformat, '( "(",I4,"(E16.8),A)" )' ) ny
      open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=(ny*16+1))
      write(67,fileformat,rec=itime) field,NL
      close(67)
    end subroutine outpdt
    !############################################################################
    subroutine outpd(field,filename) !outpost to disk just computed
      USE param
      USE variables
      implicit none
      real(mytype),intent(IN),dimension(ysize(2)) :: field
      character(len=*), intent(IN) :: filename
      integer :: i
      character(len=1),parameter :: NL=char(10)
      character(len=300) :: fileformat
      write(fileformat, '( "(",I4,"(E16.8),A)" )' ) ny
      open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=(ny*16+1))
      write(67,fileformat,rec=(itime-initstat)/iprocessing) field,NL
      close(67)
    end subroutine outpd
    !############################################################################
    subroutine horizontal_avrge(field2,profile)
      USE param
      USE variables
      USE MPI
      implicit none
      real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: field2
      real(mytype),dimension(ysize(2)) :: sxz, sxz1
      real(mytype),intent(out),dimension(ysize(2)) :: profile
      integer :: j,code
      sxz=zero
      do j=1,ysize(2)
        sxz(j) = sum(field2(:,j,:))
      enddo
      call MPI_ALLREDUCE(sxz,sxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      profile = sxz1/real(nx*nz,mytype)
    end subroutine horizontal_avrge
    !############################################################################
    subroutine extract_fluctuat(field2,profile,field2p)
      USE param
      USE variables
      USE MPI
      implicit none
      real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: field2
      real(mytype),intent(in),dimension(ysize(2)) :: profile
      real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: field2p
      integer :: j
      do j=1,ysize(2)
        field2p(:,j,:) = field2(:,j,:) - profile(j)
      enddo
    end subroutine extract_fluctuat
    !############################################################################
end module tbl
