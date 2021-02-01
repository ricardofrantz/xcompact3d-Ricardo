module flow_type
  use decomp_2d, only : mytype
  integer :: scalar_reset
end module flow_type

subroutine ft_parameter(arg)
  USE param
  USE variables
  USE flow_type
  USE complex_geometry
  USE decomp_2d, only : nrank

  implicit none
  logical,intent(in) :: arg
  integer :: is
  character :: a

  nphi=1; ro=zero

  open(10,file='BC-ttbl.prm',status='old',form='formatted')
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D computational parameters
  read (10,*) a !
  read (10,*) p_row
  read (10,*) p_col
  read (10,*) a
  read (10,*) nx
  read (10,*) ny
  read (10,*) nz
  if (arg) then
    close(10)
    return
  endif
  read (10,*) xlx
  read (10,*) yly
  read (10,*) zlz
  read (10,*) re
  read (10,*) iscalar
  read (10,*) ri_b(1)
  read (10,*) noise
  read (10,*) dt
  read (10,*) a
  read (10,*) iturb
  read (10,*) ilit
  read (10,*) ifirst
  read (10,*) ilast
  read (10,*) istatsig
  read (10,*) isave
  read (10,*) imodulo
  read (10,*) iprocessing
  read (10,*) itest
  read (10,*) iskips
  read (10,*) ifprobes
  read (10,*) istheta
  read (10,*) beta
  read (10,*) jLES
  read (10,*) fpi2

  nscheme=2 !# Temporal scheme (2:AB3)
  istret=3 !y-mesh refinement (3:bottom)
  cont_phi=0
  nclx1=0;nclxn=0;nclxS1=0;nclxSn=0
  nclz1=0;nclzn=0;nclzS1=0;nclzSn=0

  ncly1=2 !BC in y=0  (2: Dirichlet)
  nclyn=2 !BC in y=Ly (2:Dirichlet)

  if(ri_b(1).eq.0.0)then
   nclyS1=1 !BC in y=0  (2: Dirichlet)
   nclySn=1 !BC in y=Ly (2: Dirichlet)
  else
   nclyS1=2 !BC in y=0  (2: Dirichlet)
   nclySn=2 !BC in y=Ly (2: Dirichlet)
  endif
  
  return
end subroutine ft_parameter
!********************************************************************
subroutine init (ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)
  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI
  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1,phis1,phiss1

  real(mytype) :: y,r,r3,x,z,h,aform!ct
  !real(mytype) :: cx0,cy0,cz0,hg,lg
  !real(mytype) :: a_a,a_0,a_1,b_0,b_1,b_2,b_3,c_c,c_0,c_1,c_2,c_3,d_0,d_1,d_2,d_3

  integer :: k,j,i,ijk,fh,ierror,ii,is,code
  integer (kind=MPI_OFFSET_KIND) :: disp
  real(mytype), dimension(ysize(2)) :: um
  integer, dimension (:), allocatable :: seed
  !real(mytype), save, allocatable, dimension(:) :: init_u, init_v, init_w
  !allocate(init_u(ysize(2)),init_v(ysize(2)),init_w(ysize(2)))

  phi1=zero
  phis1=phi1
  phiss1=phis1
  ux1=zero;uy1=zero;uz1=zero

  !call system_clock(count=code)
  code=0 !fixing seed!
  call random_seed(size = ii)
  call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))
  call random_number(ux1)
  call random_number(uy1)
  call random_number(uz1)

  aform = 0.23369497_mytype
       do k=1,xsize(3)
          z=real((k+xstart(3)-1-1),mytype)*dz
          do j=1,xsize(2)
           if (istret.eq.0) y=real((j+xstart(2)-1-1),mytype)*dy
           if (istret.ne.0) y=yp(j+xstart(2)-1)
            !um=exp(-zptwo*y*y)
             um = (one-erf(aform*y))
             do i=1,xsize(1)
             x=real(i-1,mytype)*dx
             ux1(i,j,k)=noise*um(j)*(two*ux1(i,j,k)-one)+erf(y*aform)+&
             4*noise*(y*exp(-y)/0.3678784468818499_mytype)*&
             ( cos(z*pi/five)*cos(x*pi/five) + &
               cos((x+((one+sqrt(five))*half))*pi/five)*cos((z+((one+sqrt(five))*half))*pi/five) )
             uy1(i,j,k)=noise*um(j)*(two*uy1(i,j,k)-one)
             uz1(i,j,k)=noise*um(j)*(two*uz1(i,j,k)-one)
          enddo
        enddo
      enddo

      !call transpose_x_to_y(ux1,ux2)
      !call transpose_x_to_y(uy1,uy2)
      !call transpose_x_to_y(uz1,uz2)

      !open(10,file='ux1_init.dat',status='unknown',form='formatted')
      !open(11,file='uy1_init.dat',status='unknown',form='formatted')
      !open(12,file='uz1_init.dat',status='unknown',form='formatted')

    !write(fileformat, '( "(",I4,"(E16.8),A)" )' ) ny
    !open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=(ny*16+1))
    !do i=1, iprocessing
    !  write(67,fileformat,rec=(itime-istatsig)-iprocessing+i) field(:,i),NL
    !enddo

       !do i=1,ysize(2)
        !read(10,*) init_u(i)
        !if(nrank.eq.0)write(*,*) 'i init_u =',i,init_u(i)
        !read(11,*) init_v(i)
        !if(nrank.eq.0)write(*,*) 'i init_v =',i,init_v(i)
        !read(12,*) init_w(i)
        !if(nrank.eq.0)write(*,*) 'i init_w =',i,init_w(i)
       !enddo

      !do i=1,ysize(2)
        !ux2(:,i,:) = ux2(:,i,:) + init_u(i)
        !uy2(:,i,:) = uy2(:,i,:) + init_v(i)
        !uz2(:,i,:) = uz2(:,i,:) + init_w(i)
      !enddo

      !call transpose_y_to_x(ux2,ux1)
      !call transpose_y_to_x(uy2,uy1)
      !call transpose_y_to_x(uz2,uz1)

      !call decomp_2d_read_one(1,ux1,'ux1_init.dat')
      !call decomp_2d_read_one(1,uy1,'uy1_init.dat')
      !call decomp_2d_read_one(1,uz1,'uz1_init.dat')
  !  endif


  if(iturb.eq.1)then
    call decomp_2d_read_one(1,ux1,'ux1.ic')
    call decomp_2d_read_one(1,uy1,'uy1.ic')
    call decomp_2d_read_one(1,uz1,'uz1.ic')
  endif

  !INIT FOR G AND U=MEAN FLOW + NOISE
  gx1=ux1;gy1=uy1;gz1=uz1
  hx1=gx1;hy1=gy1;hz1=gz1
  return
end subroutine init
!********************************************************************
subroutine boundary_conditions (ux1,uy1,uz1,phi1,ep1)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  USE flow_type
  USE MPI
  implicit none

  integer  :: i,j,k
  real(mytype),intent(inout),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1

  !new diverge free BC to allow vertical outflow and minimize pressure bump

  !dyp2 = (yp(ny)-yp(ny-1)/(yp(ny-1)-yp(ny-2))

  !byxn(i,k)= (-d2udx2(:,ny-1,:) -d2udz2(:,ny-1,:))*dyp2 +ux(:,ny-1,:) -ux(:,ny-2,:)
  !byyn(i,k)= (-d2vdx2(:,ny-1,:) -d2vdz2(:,ny-1,:))*dyp2 +uy(:,ny-1,:) -uy(:,ny-2,:)
  !byzn(i,k)= (-d2wdx2(:,ny-1,:) -d2wdz2(:,ny-1,:))*dyp2 +uz(:,ny-1,:) -uz(:,ny-2,:)

  !du/dy = -lambda*u !lambda <<1 Newman; lambda >> 1 Dirichlet; lambda = 1/Ly -> mixed
  !byxn(i,k)=ux(:,ny-1,:) - dy*lambda*ux(:,ny-1,:) = (1-lambda*dy)*ux(:,ny-1,:)

   if (xend(2)==ny) then
    do k=1,xsize(3)
     do i=1,xsize(1) 
      byxn(i,k)=one
      byyn(i,k)=zero
      byzn(i,k)=zero
     enddo
    enddo
   endif
  return

end subroutine boundary_conditions
#ifdef POST
!********************************************************************
module post_processing

  USE decomp_2d
  USE variables
  USE param
  USE flow_type
  implicit none
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character
  integer, save :: nprobes
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes
  real(mytype),save,allocatable,dimension(:) :: usum, vsum, wsum, uusum, uvsum, uwsum, vvsum, vwsum, wwsum
  real(mytype),save,allocatable,dimension(:,:) :: phisxz, phipphipsxz, upphipsxz, vpphipsxz, wpphipsxz, dphidysxz
  real(mytype),save,allocatable,dimension(:,:) :: enssxz, psxz, usxz, vsxz, wsxz, upupsxz, enspenspsxz
  real(mytype),save,allocatable,dimension(:,:) :: vpvpsxz, wpwpsxz, upvpsxz, vpwpsxz, upwpsxz
  real(mytype),save,allocatable,dimension(:,:) :: tkesxz, Isxz ,epssxz
  real(mytype),save,allocatable, dimension(:,:,:) :: vol1

  real(mytype),save, allocatable, dimension(:) :: ttd1a,ttd1b,ttd1c
  real(mytype),save, allocatable, dimension(:) :: ttd2a,ttd2b,ttd2c
  real(mytype),save, allocatable, dimension(:) :: ttd3a,ttd3b,ttd3c
contains

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: x, xprobes, yprobes, zprobes, dxdydz, dxdz
    integer :: i,j,k,code
    character :: a

    allocate(ttd1a(ysize(2)),ttd1b(ysize(2)),ttd1c(ysize(2)))
    allocate(ttd2a(ysize(2)),ttd2b(ysize(2)),ttd2c(ysize(2)))
    allocate(ttd3a(ysize(2)),ttd3b(ysize(2)),ttd3c(ysize(2)))

    allocate(usum(ysize(2)),vsum(ysize(2)),wsum(ysize(2)))
    allocate(uusum(ysize(2)),uvsum(ysize(2)),uwsum(ysize(2)))
    allocate(vvsum(ysize(2)),vwsum(ysize(2)),wwsum(ysize(2)))

    call alloc_x(vol1, opt_global=.true.)

    allocate(enssxz(ysize(2),      iprocessing/iskips))
    allocate(psxz(ysize(2),        iprocessing/iskips))
    allocate(usxz(ysize(2),        iprocessing/iskips))
    allocate(vsxz(ysize(2),        iprocessing/iskips))
    allocate(wsxz(ysize(2),        iprocessing/iskips))
    allocate(upupsxz(ysize(2),     iprocessing/iskips))
    allocate(enspenspsxz(ysize(2), iprocessing/iskips))
    allocate(vpvpsxz(ysize(2),     iprocessing/iskips))
    allocate(wpwpsxz(ysize(2),     iprocessing/iskips))
    allocate(upvpsxz(ysize(2),     iprocessing/iskips))
    allocate(vpwpsxz(ysize(2),     iprocessing/iskips))
    allocate(upwpsxz(ysize(2),     iprocessing/iskips))
    allocate(tkesxz(ysize(2),      iprocessing/iskips))
    allocate(Isxz(ysize(2),        iprocessing/iskips))
    allocate(epssxz(ysize(2),      iprocessing/iskips))
    allocate(phisxz(ysize(2),      iprocessing/iskips))
    allocate(phipphipsxz(ysize(2), iprocessing/iskips))
    allocate(upphipsxz(ysize(2),   iprocessing/iskips))
    allocate(vpphipsxz(ysize(2),   iprocessing/iskips))
    allocate(wpphipsxz(ysize(2),   iprocessing/iskips))
    allocate(dphidysxz(ysize(2),   iprocessing/iskips))

    !X PENCILS !Utilizar para integral volumétrica dentro do domínio físico (método de Simpson)
    vol1 = zero
    dxdydz=dx*dy*dz
    do k=xstart(3),xend(3)
      do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
          vol1(i,j,k)=dxdydz
          if (i .eq. 1 .or. i .eq. nx) vol1(i,j,k) = vol1(i,j,k) * five/twelve
          if (j .eq. 1 .or. j .eq. ny) vol1(i,j,k) = vol1(i,j,k) * five/twelve
          if (k .eq. 1 .or. k .eq. nz) vol1(i,j,k) = vol1(i,j,k) * five/twelve
          if (i .eq. 2 .or. i .eq. nx-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
          if (j .eq. 2 .or. j .eq. ny-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
          if (k .eq. 2 .or. k .eq. nz-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
        end do
      end do
    end do

    if(ifprobes.gt.0)then

      nprobes = int( (ny - 3)/ ifprobes )
      allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
      rankprobes(:)=0

      j=1
      do i=2, ny-2, ifprobes
        nyprobes(j) = i
        if (nrank .eq. 0) write(*,*) 'Probe number ',j,'ny =',nyprobes(j)
        j=j+1
      enddo
      nxprobes(:) = nx / 2
      nzprobes(:) = nz / 2

      do i=1, nprobes
        if       (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
          if    (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
            if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
              rankprobes(i)=1
            endif
          endif
        endif
      enddo

    else
      if (nrank .eq. 0) write(*,*) 'Probes OFF '
    endif

  end subroutine init_post
  !############################################################################
  subroutine comp_thetad(ux1,uy1)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  USE flow_type
  USE MPI
  implicit none

  integer  :: i,j,k,code
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2p,dudy2,di2
  real(mytype),dimension(ysize(2)) :: ux2m,ux2m_old,uxuy2pm
  real(mytype),dimension(ysize(2)) :: ydudy2m,dudy2m,du2dy22m,duxuy2pm
  real(mytype) :: resi,resi1,cf_ref,delta
  real(mytype) :: y,tau_wall,tau_wall1,ufric,temp,utau
  real(8)      :: tstart
  tstart=MPI_WTIME()

   if(itime.eq.ifirst)then
    if(nrank.eq.0)print *,'first dy=',yp(2)
    thetad=zero
   endif

   call transpose_x_to_y(ux1,ux2)
   call horizontal_avrge(ux2,ux2m)
   theta = sum((ux2m*(one-ux2m))*ypw)
   delta = sum((one-ux2m)*ypw)
   resi = theta - one
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

    print *,'Time computing theta dot (s)', real(MPI_WTIME()-tstart,4)
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
    
    print *,'theta,thetad=',theta,thetad
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
  subroutine postprocessing(ux1,uy1,uz1,phi1,ep1,pre1,diss1)

    USE param
    USE variables
    USE decomp_2d
    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, pre1, ep1
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3)) :: diss1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2 , uy2,  uz2, phi2, pre2

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1p, uy1p, uz1p, phi1p, dphi1pdx
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2p, uy2p, uz2p, phi2p, pre2p, dphi2pdx, dphi2pdy, dphi2pdz, dpre2pdy, ens2p
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) ::                   phi3p,        dphi3pdz

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: duy2pdy
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,temp1,tke3d1,ens1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,di2,temp2,tke3d2,ens2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,temp3

    real(mytype),dimension(ysize(2)) :: sxz, sxz1
    real(mytype),dimension(ysize(2)) :: dup2dxsxz,dvp2dxsxz,dwp2dxsxz
    real(mytype),dimension(ysize(2)) :: dup2dysxz,dvp2dysxz,dwp2dysxz
    real(mytype),dimension(ysize(2)) :: dup2dzsxz,dvp2dzsxz,dwp2dzsxz
    real(mytype),dimension(ysize(2)) :: vdpre2pdysxz,pduy2pdysxz
    real(mytype),dimension(ysize(2)) :: pre_sxz, dvpppsxz, vpppsxz, vpdppsxz, ppppsxz !pressure transport
    real(mytype),dimension(ysize(2)) :: kurt_x,kurt_y,kurt_z,kurt_phi,skew_x,skew_y,skew_z,skew_phi
    real(mytype),dimension(ysize(2)) :: dupupvpsxz, dvpwpwpsxz, dvpvpvpsxz !turbulent transport
    real(mytype),dimension(ysize(2)) :: dusxz,dwsxz !shear production
    real(mytype),dimension(ysize(2)) :: d2upupsxz, d2vpvpsxz, d2wpwpsxz, dtkesxzdy, d2tkesxzdy !visocus diffusion
    real(mytype),dimension(ysize(2)) :: dvpphipphipsxz !buoyancy turbulent transport
    real(mytype),dimension(ysize(2)) :: d2phipphipsxz, dphipphipsxz !buoyancy diffusion
    real(mytype),dimension(ysize(2)) :: dphisxz !buoyancy production
    real(mytype),dimension(ysize(2)) :: b_eps !buoyancy dissipation
    real(mytype),dimension(ysize(2)) :: dusxzdy,dwsxzdy!dudy and dwdy

    real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu
    real(mytype) :: mp(nphi)
    real(8) :: tstart
    integer :: i,j,k,is,code,b
    character(len=1),parameter :: NL=char(10)
    character(len=60) :: filename
    character(len=300) :: fileformat

    call comp_thetad(ux1,uy1)

    if (itime.gt.istatsig.AND.mod(itime,iskips).eq.0)then
    tstart=MPI_WTIME()

    b = mod(itime,iprocessing)/iskips
    if (mod(itime,iprocessing).eq.0) b = iprocessing/iskips
    !if(nrank.eq.0)print *,'b=',b

    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    !y-derivatives
    call transpose_x_to_y(ux1,td2)
    call transpose_x_to_y(uy1,te2)
    call transpose_x_to_y(uz1,tf2)

    call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    !!z-derivatives
    call transpose_y_to_z(td2,td3)
    call transpose_y_to_z(te2,te3)
    call transpose_y_to_z(tf2,tf3)
    call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)

    ens1=(tf1-th1)**2+(tg1-tc1)**2+(tb1-td1)**2
    call transpose_x_to_y(ens1,ens2)
    call horizontal_avrge(ens2,enssxz(:,b)) ! enstrophy
    call extract_fluctuat(ens2,enssxz(:,b),ens2p) ! enstrophy'
    call horizontal_avrge(ens2p**2,enspenspsxz(:,b)) ! enstrophy' enstrophy'

    call transpose_x_to_y(ux1,ux2)
    call horizontal_avrge(ux2,usxz(:,b)) ! ux
    call extract_fluctuat(ux2,usxz(:,b),ux2p) ! ux'
    call transpose_y_to_x(ux2p,ux1p)
    
    call transpose_x_to_y(uy1,uy2)
    call horizontal_avrge(uy2,vsxz(:,b)) ! uy
    call extract_fluctuat(uy2,vsxz(:,b),uy2p) ! uy'
    !uy2p = uy2 !damian case
    call transpose_y_to_x(uy2p,uy1p)

    call transpose_x_to_y(uz1,uz2)
    call horizontal_avrge(uz2,wsxz(:,b)) ! uz
    call extract_fluctuat(uz2,wsxz(:,b),uz2p) ! uz'
    call transpose_y_to_x(uz2p,uz1p)

    call horizontal_avrge(ux2p**2,upupsxz(:,b)) ! ux' ux'
    call horizontal_avrge(uy2p**2,vpvpsxz(:,b)) ! uy' uy'
    call horizontal_avrge(uz2p**2,wpwpsxz(:,b)) ! uz' uz'

    call horizontal_avrge(ux2p*uy2p,upvpsxz(:,b)) ! ux' uy'
    call horizontal_avrge(uy2p*uz2p,vpwpsxz(:,b)) ! uy' uz'
    call horizontal_avrge(ux2p*uz2p,upwpsxz(:,b)) ! ux' uz'

    tke3d2 = half*(ux2p**2 + uy2p**2 + uz2p**2)
    call horizontal_avrge(tke3d2,tkesxz(:,b)) ! tke
    call transpose_y_to_x(tke3d2,tke3d1)
    call tke(tke3d1*vol1)
    !if (mod(itime,imodulo).eq.0) then
    ! write(filename,"('./data/tke',I4.4)") itime/imodulo
    ! call decomp_2d_write_one(1,tke3d1uvisu,filename,2)
    !endif

    !call horizontal_avrge(sqrt((two*tke3d2)/(ux2**2+uy2**2+uz2**2)),tkesxz(:,b)) ! turbulent intensity !! sqrt(w*tke)/|U|

    if(iscalar==1)then
     call transpose_x_to_y(phi1(:,:,:,1),phi2)
     call horizontal_avrge(phi2,phisxz(:,b)) ! phi
     call extract_fluctuat(phi2,phisxz(:,b),phi2p) ! phi'
     !phi2p = phi2 !damian case
     call transpose_y_to_x(phi2p,phi1p)
     call horizontal_avrge(phi2p**2,phipphipsxz(:,b)) ! phi' phi'
     call horizontal_avrge(ux2p*phi2p,upphipsxz(:,b)) ! ux' phi'
     call horizontal_avrge(uy2p*phi2p,vpphipsxz(:,b)) ! uy' phi'
     call horizontal_avrge(uz2p*phi2p,wpphipsxz(:,b)) ! uz' phi'
     call dery (dphidysxz(:,b),phisxz(:,b),di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1)!1
    endif

    call derx (ta1,ux1p,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tb1,uy1p,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (tc1,uz1p,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    !y-derivatives
    call dery (ta2,ux2p,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (tb2,uy2p,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (tc2,uz2p,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    !!z-derivatives
    call transpose_y_to_z(ux2p,td3)
    call transpose_y_to_z(uy2p,te3)
    call transpose_y_to_z(uz2p,tf3)
    call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

    !!!!!!!!!!!!!!!!!!!!!! diss sij sij - viscous dissipation rate

    call dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1) !sij sij
    call transpose_x_to_y(diss1,temp2)
    call horizontal_avrge(temp2,epssxz(:,b)) ! diss
    if (mod(itime,imodulo).eq.0) then
     write(filename,"('./data/dissp',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,diss1,filename,2)
    endif

    if (nrank.eq.0) print *,'Time computing statistics (s)', real(MPI_WTIME()-tstart,4)

    if (mod(itime,iprocessing).eq.0 .AND. nrank.eq.0 .AND. itime.ge.iprocessing) then
      tstart=MPI_WTIME()
      call outp(enssxz,     './out/02_ens')
      call outp(usxz,       './out/03_ux')
      call outp(vsxz,       './out/04_uy')
      call outp(wsxz,       './out/05_uz')
      call outp(enspenspsxz,'./out/06_ensp2')
      call outp(upupsxz,    './out/07_uxp2')
      call outp(vpvpsxz,    './out/08_uyp2')
      call outp(wpwpsxz,    './out/09_uzp2')
      call outp(upvpsxz,    './out/10_uxp_uyp')
      call outp(vpwpsxz,    './out/11_uyp_uzp')
      call outp(upwpsxz,    './out/12_uxp_uzp')
      call outp(tkesxz,     './out/13_tke')
      !call outp(Isxz,       './out/14_tint')
      call outp(-epssxz,     './out/15_diss')
 
     if(iscalar==1)then

       call outp(phisxz,     './out/16_phi')
       call outp(phipphipsxz,'./out/17_phip2')
       call outp(upphipsxz,  './out/18_phip_uxp')
       call outp(vpphipsxz,  './out/19_phip_uyp')
       call outp(wpphipsxz,  './out/20_phip_uzp')
       call outp(dphidysxz,  './out/21_brunt_vaisala')

      endif

      print *,'Time writing statistics (s)', real(MPI_WTIME()-tstart,4)

    endif

    if (mod(itime,iprocessing).eq.0) then

      tstart=MPI_WTIME()

      call transpose_x_to_y(ta1,temp2) !ta1=du/dx
      call horizontal_avrge(temp2**2,dup2dxsxz)

      call transpose_x_to_y(tb1,temp2) !tb1=dv/dx
      call horizontal_avrge(temp2**2,dvp2dxsxz)

      call transpose_x_to_y(tc1,temp2) !tc1=dw/dx
      call horizontal_avrge(temp2**2,dwp2dxsxz)

      call transpose_x_to_y(td1,temp2) !td1=du/dy
      call horizontal_avrge(temp2**2,dup2dysxz)

      call transpose_x_to_y(te1,temp2) !te1=dv/dy
      call horizontal_avrge(temp2**2,dvp2dysxz)

      call transpose_x_to_y(tf1,temp2) !tf1=dw/dy
      call horizontal_avrge(temp2**2,dwp2dysxz)

      call transpose_x_to_y(tg1,temp2) !tg1=du/dz
      call horizontal_avrge(temp2**2,dup2dzsxz)

      call transpose_x_to_y(th1,temp2) !th1=dv/dz
      call horizontal_avrge(temp2**2,dvp2dzsxz)

      call transpose_x_to_y(ti1,temp2) !ti1=dw/dz
      call horizontal_avrge(temp2**2,dwp2dzsxz)

    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

      call transpose_x_to_y(pre1,pre2) ! p
      call horizontal_avrge(pre2,pre_sxz)
      call extract_fluctuat(pre2,pre_sxz,pre2p) ! p'
      call horizontal_avrge(pre2p**2,ppppsxz) ! p' p'

      !pressure transport rate !d(w' p')/dy
      call horizontal_avrge(uy2p*pre2p,vpppsxz) ! p' v'
      call dery (dvpppsxz,vpppsxz,di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1)

      !turbulent transport
      call horizontal_avrge(ux2p*ux2p*uy2p,sxz1) ! u' u' v'
      call dery (dupupvpsxz,sxz1,di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1) !1*1*0=1
      call horizontal_avrge(uy2p*uz2p*uz2p,sxz1) ! v' w' w'
      call dery (dvpwpwpsxz,sxz1,di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1) !1*1*0=1
      call horizontal_avrge(uy2p*uy2p*uy2p,sxz1) ! v' v' v'
      call dery (dvpvpvpsxz,sxz1,di2,sy,ffy,fsy,fwy,ppy,1,ysize(2),1,0) !0*0*0=0

      !shear production rate
      call dery (dusxz,usxz(:,b),di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1) !1
      call dery (dwsxz,wsxz(:,b),di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1) !1

      !viscous diffusion rate !ux/y -> 1, uy/y -> 0, uz/y -> 1
      d2tkesxzdy=zero
      dtkesxzdy=zero !auxiliary field - ignore the name
      call dery (dtkesxzdy,tkesxz(:,b),di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1) !0*0=1
      call deryy (d2tkesxzdy,tkesxz(:,b),di2,sy,sfyp,ssyp,swyp,1,ysize(2),1,1) !0*0=1
      if (istret.ne.0) then
        do j = 1,ysize(2)
          d2tkesxzdy(j) = d2tkesxzdy(j)*pp2y(j)-dtkesxzdy(j)*pp4y(j)
        enddo
      endif

      if(iscalar==1)then

       !buoyancy production
       call dery (dphisxz,phisxz(:,b),di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1) !1

       !buoyancy turbulent transport
       call horizontal_avrge(uy2p*phi2p*phi2p,sxz1) ! vi phi' phi'
       call dery (dvpphipphipsxz,sxz1,di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1) !1*1*0=1

       !buoyancy diffusion rate !ux/y -> 1, uy/y -> 0, uz/y -> 1
       d2phipphipsxz=zero
       call deryy (d2phipphipsxz,phipphipsxz(:,b),di2,sy,sfyp,ssyp,swyp,1,ysize(2),1,1) !0*0=1
       if (istret.ne.0) then
         dphipphipsxz=zero
         call dery (dphipphipsxz,phipphipsxz(:,b),di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1) !0*0=1
         do j = 1,ysize(2)
           d2phipphipsxz(j) = d2phipphipsxz(j)*pp2y(j)-dphipphipsxz(j)*pp4y(j)
         enddo
 
      endif

      !buoyancy dissipation rate
      call derx (dphi1pdx,phi1p,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
      call transpose_x_to_y(dphi1pdx,dphi2pdx)
      call dery (dphi2pdy,phi2p,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
      call transpose_y_to_z(phi2p,phi3p)
      call derz (dphi3pdz,phi3p,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
      call transpose_z_to_y(dphi3pdz,dphi2pdz)
      call horizontal_avrge(dphi2pdx**2+dphi2pdy**2+dphi2pdz**2,b_eps)

      endif

      !dudy and dwdy
      call dery (dusxzdy,usxz(:,b),di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1)
      call dery (dwsxzdy,wsxz(:,b),di2,sy,ffyp,fsyp,fwyp,ppy,1,ysize(2),1,1)

      !higher order moments
      call horizontal_avrge(ux2p**three,skew_x) 
      call horizontal_avrge(ux2p**four,kurt_x)       
      call horizontal_avrge(uy2p**three,skew_y)
      call horizontal_avrge(uy2p**four,kurt_y)
      call horizontal_avrge(uz2p**three,skew_z)
      call horizontal_avrge(uz2p**four,kurt_z)

      if(iscalar==1)then
       call horizontal_avrge(phi2p**three,skew_phi)
       call horizontal_avrge(phi2p**four,kurt_phi)
      endif
 
      if (mod(itime,iprocessing).eq.0 .AND. nrank.eq.0 .AND. itime.ge.iprocessing) then

        call outpd(skew_x/(upupsxz(:,b)**1.5),'./out/34_skew_x')!mean(y^3)/mean(y^2)^1.5
        call outpd(skew_y/(vpvpsxz(:,b)**1.5),'./out/34_skew_y')
        call outpd(skew_z/(wpwpsxz(:,b)**1.5),'./out/34_skew_z')
        if(iscalar==1)call outpd(skew_phi/(wpwpsxz(:,b)**1.5),'./out/34_skew_phi')

        call outpd(kurt_x/(upupsxz(:,b)**2),'./out/35_kurt_x')!mean(y^4)/mean(y^2)^2
        call outpd(kurt_y/(vpvpsxz(:,b)**2),'./out/35_kurt_y')
        call outpd(kurt_z/(wpwpsxz(:,b)**2),'./out/35_kurt_z')
        if(iscalar==1)call outpd(kurt_phi/(phipphipsxz(:,b)**2),'./out/35_kurt_phi')

        call outpd(two*dup2dxsxz/dwp2dxsxz,'./out/20_K1')
        call outpd(two*dup2dxsxz/dvp2dxsxz,'./out/20_K2')
        call outpd(two*dup2dxsxz/dup2dzsxz,'./out/20_K3')
        call outpd(two*dup2dxsxz/dup2dysxz,'./out/20_K4')

        call outpd(-two*xnu*(dup2dxsxz+dup2dysxz+dup2dzsxz),'./out/31_eps_x')
        call outpd(-two*xnu*(dvp2dxsxz+dvp2dysxz+dvp2dzsxz),'./out/31_eps_y')
        call outpd(-two*xnu*(dwp2dxsxz+dwp2dysxz+dwp2dzsxz),'./out/31_eps_z')

        call outpd(dusxzdy,'./out/32_dudy')
        call outpd(dwsxzdy,'./out/32_dwdy')

        call outpd(pre_sxz ,'./out/21_pre')
        call outpd(ppppsxz ,'./out/22_pre2p')

        call outpd(-dvpppsxz,                               './out/28_velp_grad')
        call outpd(-half*(dupupvpsxz+dvpwpwpsxz+dvpvpvpsxz),'./out/28_t_diff')
        call outpd(-(upvpsxz(:,b)*dusxz+vpwpsxz(:,b)*dwsxz),'./out/28_prod')
        call outpd(xnu*d2tkesxzdy,                          './out/28_vis_diff')
        call outpd(vsxz(:,b)*dtkesxzdy,                     './out/28_meanf_transp')

        if(iscalar==1)then
         call outpd(ri_b(1)*vpphipsxz(:,b),             './out/30_b_flux')
         call outpd(-two*ri_b(1)*vpphipsxz(:,b)*dphisxz,'./out/30_b_prod')
         call outpd(-ri_b(1)*dvpphipphipsxz,            './out/30_b_ttransp')
         call outpd((ri_b(1)*xnu/nsc(1))*d2phipphipsxz, './out/30_b_diff')
         call outpd((-two*ri_b(1)*xnu/nsc(1))*b_eps ,   './out/30_b_diss')
        endif
      endif!(mod(itime,iprocessing).eq.0 .AND. nrank.eq.0 .AND. itime.ge.iprocessing)
      mp=zero
      call budget(ux1,uy1,uz1,phi1,vol1)
      call suspended(phi1,vol1,mp)
       if (nrank .eq. 0) then
         FS = 1+nphi !Number of columns
         write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
         FS = FS*14+1  !Line width
         open(67,file='./out/0_suspended',status='unknown',form='formatted',access='direct',recl=FS)
         write(67,fileformat,rec=itime/iprocessing+1) t,mp,NL
         close(67)
         print *,'Time computing budgets (s)', real(MPI_WTIME()-tstart,4)
       endif

    endif !(mod(itime,iprocessing).eq.0)
    endif!if (itime.gt.istatsig.AND.mod(itime,iskips).eq.0)
  end subroutine postprocessing
  !############################################################################
  subroutine tke(field)
    USE param
    USE variables
    USE MPI
    implicit none
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: field
    integer :: i,FS,code
    real(mytype) :: tkee,temp
    character(len=1),parameter :: NL=char(10)
    character(len=300) :: fileformat

    temp = sum(field)
    call MPI_REDUCE(temp,tkee,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    tkee=tkee/real(nx*ny*nz,mytype)

    if (nrank .eq. 0) then
      FS = 1+1 !Number of columns
      write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
      FS = FS*14+1  !Line width
      open(67,file='./out/0_tke_global',status='unknown',form='formatted',access='direct',recl=FS)
      write(67,fileformat,rec=itime/iskips+1) t,tkee,NL
      close(67)
      if(mod(itime,itest)==0)print *,'Global TKE=', real(tkee,4)
    endif

  end subroutine tke
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
    do i=1, iprocessing/iskips
    !if(nrank.eq.0)print *, 'itime/iskips=',itime/iskips-iprocessing/iskips+i
      write(67,fileformat,rec=(itime-istatsig)/iskips -iprocessing/iskips + i) field(:,i),NL
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
    write(67,fileformat,rec=(itime-istatsig)/iprocessing) field,NL
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
  subroutine write_probes(ux1,uy1,uz1,pre1,diss1,phi1) !By Felipe Schuch
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1, pre1, diss1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),nphi) :: phi1
    integer :: i
    character(len=30) :: filename
    !FS = 1+3+2+nphi !Number of columns
    FS = 1+3+nphi !Number of columns
    write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
    FS = FS*14+1  !Line width
    do i=1, nprobes
      if (rankprobes(i) .eq. 1) then
        write(filename,"('./probes/probe_ny_',I4.4)") nyprobes(i)
        open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=FS)
        write(67,fileformat,rec=itime) t,&                         !1
        ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
        uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
        uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
        !pre1(nxprobes(i),nyprobes(i),nzprobes(i)),&           !5
        !diss1(nxprobes(i),nyprobes(i),nzprobes(i)),&          !6
        phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !nphi
        NL                                                    !+1
        close(67)
      endif
    enddo
  end subroutine write_probes
  !############################################################################
  subroutine dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1)
    USE param
    USE variables
    USE decomp_2d
    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1
    real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A
    integer :: k,j,i,m,l
    !INSTANTANEOUS DISSIPATION RATE
    diss1=0._mytype
    A(:,:,:,:,:)=0._mytype
    A(1,1,:,:,:)=ta1(:,:,:)!du/dx=ta1
    A(2,1,:,:,:)=tb1(:,:,:)!dv/dx=tb1
    A(3,1,:,:,:)=tc1(:,:,:)!dw/dx=tc1
    A(1,2,:,:,:)=td1(:,:,:)!du/dy=td1
    A(2,2,:,:,:)=te1(:,:,:)!dv/dy=te1
    A(3,2,:,:,:)=tf1(:,:,:)!dw/dy=tf1
    A(1,3,:,:,:)=tg1(:,:,:)!du/dz=tg1
    A(2,3,:,:,:)=th1(:,:,:)!dv/dz=th1
    A(3,3,:,:,:)=ti1(:,:,:)!dw/dz=ti1
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          do m=1,3
            do l=1,3
              diss1(i,j,k)=diss1(i,j,k)+two*xnu*half*half*(A(l,m,i,j,k)+A(m,l,i,j,k))**2
            enddo
          enddo
        enddo
      enddo
    enddo
    return
  end subroutine dissipation
  !############################################################################
   subroutine suspended(phi1,vol1,mp1)
  
     USE decomp_2d_io
     USE MPI
  
     implicit none
  
     real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
     real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: vol1
     real(mytype),intent(out) :: mp1(1:nphi)
   
     real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1
     real(mytype) :: mp(1:nphi)
     integer :: is,code
  
     mp=zero; mp1=zero
  
     do is=1, nphi
       temp1 = phi1(:,:,:,is)*vol1(:,:,:)
       mp(is)= sum(temp1)
     end do
  
     call MPI_REDUCE(mp,mp1,nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  
     return
   end subroutine suspended
   !############################################################################
   subroutine budget(ux1,uy1,uz1,phi1,vol1)
  
     USE decomp_2d
     USE decomp_2d_io
     USE MPI
  
     implicit none
  
     real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,vol1
     real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
  
     real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: diss1,dphiy1, dphixx1, dphiyy1, dphizz1, ddphi1, di1, temp1
     real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi2, dphiy2, dphixx2, dphiyy2, dphizz2, ddphi2, di2, vol2, temp2
     real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: phi3, dphiy3, dphixx3, dphiyy3, dphizz3, ddphi3, di3, temp3
  
     real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1
     real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2
     real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3
     real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A
  
     real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu
  
     real(8) :: ek,ek1,dek,dek1,ep,ep1,dep,dep1,xvol
     integer :: ijk,i,j,k,l,m,is,code
     character(len=30) :: filename
  
     ek=zero;ek1=zero;dek=zero;dek1=zero;ep=zero;ep1=zero;dep=zero;dep1=zero;diss1=zero
  
     !x-derivatives
     call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
     call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
     call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
     !y-derivatives
     call transpose_x_to_y(ux1,td2)
     call transpose_x_to_y(uy1,te2)
     call transpose_x_to_y(uz1,tf2)
     call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
     call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
     call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
     !!z-derivatives
     call transpose_y_to_z(td2,td3)
     call transpose_y_to_z(te2,te3)
     call transpose_y_to_z(tf2,tf3)
     call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
     call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
     call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
     !!all back to x-pencils
     call transpose_z_to_y(ta3,td2)
     call transpose_z_to_y(tb3,te2)
     call transpose_z_to_y(tc3,tf2)
     call transpose_y_to_x(td2,tg1)
     call transpose_y_to_x(te2,th1)
     call transpose_y_to_x(tf2,ti1)
     call transpose_y_to_x(ta2,td1)
     call transpose_y_to_x(tb2,te1)
     call transpose_y_to_x(tc2,tf1)
  
     A(:,:,:,:,:)=zero
     A(1,1,:,:,:)=ta1(:,:,:)
     A(2,1,:,:,:)=tb1(:,:,:)
     A(3,1,:,:,:)=tc1(:,:,:)
     A(1,2,:,:,:)=td1(:,:,:)
     A(2,2,:,:,:)=te1(:,:,:)
     A(3,2,:,:,:)=tf1(:,:,:)
     A(1,3,:,:,:)=tg1(:,:,:)
     A(2,3,:,:,:)=th1(:,:,:)
     A(3,3,:,:,:)=ti1(:,:,:)
  
     do k=1,xsize(3)
       do j=1,xsize(2)
         do i=1,xsize(1)
           do m=1,3
             do l=1,3
               diss1(i,j,k)=diss1(i,j,k)+two*xnu*half*half*(A(l,m,i,j,k)+A(m,l,i,j,k))**two
             enddo
           enddo
         enddo
       enddo
     enddo
  
     do ijk=1,xsize(1)*xsize(2)*xsize(3)
       xvol=real(vol1(ijk,1,1),8)
       ek = ek + half * xvol * (ux1(ijk,1,1)**two+uy1(ijk,1,1)**two+uz1(ijk,1,1)**two)
       dek = dek + xvol * diss1(ijk,1,1)
     enddo
  
     call transpose_x_to_y(vol1,vol2)
  
     if (ivirt==2) then
       ilag=0
     endif
     do is=1, nphi
       if (ri_b(is) .eq. 0.) cycle
       call derxxS (dphixx1,phi1(:,:,:,is),di1,sx,sfxpS,ssxpS,swxpS,xsize(1),xsize(2),xsize(3),1)
  
       call transpose_x_to_y(dphixx1,dphixx2)
  
       call transpose_x_to_y(phi1(:,:,:,is),phi2)
  
       call deryS (dphiy2,phi2,di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)  
       call deryyS (dphiyy2,phi2,di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1)
       if (istret.ne.0) then
       !call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
       !call dery (te2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
        do k = 1,ysize(3)
         do j = 1,ysize(2)
          do i = 1,ysize(1)
              dphiyy2(i,j,k) = dphiyy2(i,j,k)*pp2y(j)-pp4y(j)*dphiy2(i,j,k)
              !td2(i,j,k) = td2(i,j,k)*pp2y(j)-pp4y(j)*te2(i,j,k)
          enddo
         enddo 
        enddo
       endif

       call transpose_y_to_z(phi2,phi3)
  
       call derzzS (dphizz3,phi3,di3,sz,sfzpS,sszpS,swzpS,zsize(1),zsize(2),zsize(3),1)
  
       call transpose_z_to_y(dphizz3,dphizz2)
  
       do ijk=1,ysize(1)*ysize(2)*ysize(3)
         ddphi2(ijk,1,1)=dphixx2(ijk,1,1)+dphiyy2(ijk,1,1)+dphizz2(ijk,1,1)
       enddo
  
      do k=1,ysize(3)
         do j=1,ysize(2)
           do i=1,ysize(1)
             xvol=real(vol2(i,j,k),8)
             ep=ep + xvol * (phi2(i,j,k)*(j-1)*dy)
             dep=dep - xvol * (ddphi2(i,j,k)*xnu/nsc(is))*(j-1)*dy
           enddo
         enddo
       enddo
     enddo
     if (ivirt==2) then
       ilag=1
     endif
  
     call MPI_REDUCE(ek,ek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(dek,dek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(ep,ep1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(dep,dep1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
  
     if (nrank .eq. 0) then
       open(67,file='./out/0_budget',status='unknown',form='formatted',access='direct',recl=71) !71=5*14+1
       write(67,"(5E14.6,A)",rec=itime/iprocessing+1) t,ek1,dek1,ep1,dep1,NL
       close(67)
     end if
  
   end subroutine budget

end module post_processing
#endif
