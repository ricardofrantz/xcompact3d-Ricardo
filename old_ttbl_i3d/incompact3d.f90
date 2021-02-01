PROGRAM incompact3d

  USE decomp_2d
  USE decomp_2d_poisson
  use decomp_2d_io
  USE variables
  USE param
  USE var
  USE MPI
  USE derivX
  USE derivZ
#ifdef POST
  USE post_processing
#endif
#ifdef FORCES
  USE forces
#endif
  implicit none

  integer :: code,nlock,i,j,k,ii,bcx,bcy,bcz,fh,ierror
  real(mytype) :: x,y,z,timeleft
  real(8) :: tstart,t1,trank,tranksum,ttotal,tremaining,telapsed

  TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

  call ft_parameter(.true.)

  call MPI_INIT(code)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true

  if (nrank==0) call system('mkdir data out probes')

  call parameter()

  call init_variables()


  call schemes()

#ifdef ELES
  call filter()
#endif

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

#ifdef ELES
    call init_explicit_les()
#endif

  if (ilit==0) call init(ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)

  if (ilit==1) call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
       px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,0)

  !bypasss sauvel 0120000 without scalar 
  !if (ilit==1.and.scalar_reset.eq.1) then
  ! phi1=zero
  ! phis1=phi1 
  ! phiss1=phi1
  !endif

#ifdef IBM
  if(ilag.eq.1)  call genepsi3d(ep1)
  if(ivirt.eq.1) call body(ux1,uy1,uz1,ep1,0)
#endif

#ifdef POST
  call init_post(ep1)
#endif

#ifdef FORCES
  call init_forces(ep1)
  if (ilit==1) call restart_forces(0)
#endif

  if (ilit==0) then
     itime=0
#ifdef VISU
     call VISU_INSTA (ux1,uy1,uz1,phi1,ep1,diss1,.false.)
#endif
#ifdef POST
     call postprocessing(ux1,uy1,uz1,phi1,ep1,pre1,diss1)
#endif
  endif

  call test_speed_min_max(ux1,uy1,uz1)
  if (iscalar==1) call test_scalar_min_max(phi1)

#ifdef STATS
  u1sum=zero;v1sum=zero;w1sum=zero
  u2sum=zero;v2sum=zero;w2sum=zero
  u3sum=zero;v3sum=zero;w3sum=zero
  u4sum=zero;v4sum=zero;w4sum=zero
  uvsum=zero;uwsum=zero;vwsum=zero
  upsum=zero;vpsum=zero;wpsum=zero
  ppsum=zero; psum=zero; dudx=zero
  u1sum_tik=zero;u1sum_tak=zero
  tik1=zero;tik2=zero;tak1=zero;tak2=zero
#endif

  tstart=zero;t1=zero;trank=zero;tranksum=zero;ttotal=zero
  call cpu_time(tstart)

  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)
 
  !dtini=dt
  !t=ifirst*dt ! account for restart time  

  do itime=ifirst,ilast
     t=itime*dt !original
     
     if (nrank==0)  print *,'--------------------------------------------------'
     !if(dtini.eq.0)then
     !if(itime.eq.ifirst)then
     ! if (nrank.eq.0)print *,"   Variable time step activated: computing initial dt"
     ! call compute_cfl(ux1,uy1,uz1)
     ! dt=0.3_mytype*dtmax !dt=0 purely viscous condition
     ! if(ilit.eq.1)dt=0.05_mytype*dtmax
     ! dt1=zero
     ! dt2=zero
     ! if (nrank.eq.0)write(*,"(' initial dt 0,1,2=',F9.4,F9.4,F9.4)") dt,dt1,dt2
     ! if (nrank==0)  print *,'--------------------------------------------------'
     !elseif(itime.eq.ifirst+1)then
     ! if (nrank.eq.0)print *," Variable time step activated: computing second dt"
     ! call compute_cfl(ux1,uy1,uz1)
     !  dt=0.7_mytype*dtmax !dt=0 purely viscous condition
     ! if (nrank.eq.0)write(*,"(' second dt 0,1,2=',F9.4,F9.4,F9.4)") dt,dt1,dt2
     !endif
     !else
     !if(itime.eq.ifirst)then
     !dt=dtini*half
     !elseif(itime.eq.ifirst+1)then
     !dt=dtini 
     !endif
     !endif
    
     !t=t+dt
     !dt2=dt1
     !dt1=dt
     !call settime
     !call compute_cfl(ux1,uy1,uz1)

     call cpu_time(t1)
     if (nrank==0)write(*,"(' Time step =',i7,'/',i7,', Time unit =',F9.4,', dt =',F9.4)") itime,ilast,t,dt
  
     do itr=1,iadvance_time

        call boundary_conditions(ux1,uy1,uz1,phi1,ep1)

        !X-->Y-->Z-->Y-->X
        call convdiff(ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
             ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
             ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phi1,ep1,nut1)

        if (iscalar==1) then
           call scalar(ux1,uy1,uz1,phi1,phis1,phiss1,di1,tg1,th1,ti1,td1,&
                uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,ep1,nut1)
        endif

        !X PENCILS
        call intt(ux1,uy1,uz1,gx1,gy1,gz1,hx1,hy1,hz1,ta1,tb1,tc1)
        call pre_correc(ux1,uy1,uz1,ep1)

#ifdef IBM
        if (ivirt==1) then !solid body old school
           !we are in X-pencil
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
           call body(ux1,uy1,uz1,ep1,1)
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
        endif
#endif

        !X-->Y-->Z
        call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
             td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
             nxmsize,nymsize,nzmsize,ph1,ph3,ph4,1)

        !POISSON Z-->Z
        call poisson(pp3)

        !Z-->Y-->X
        call gradp(px1,py1,pz1,di1,td2,tf2,ta2,tb2,tc2,di2,&
             ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)

        !X PENCILS
        call corgp(ux1,ux2,uy1,uz1,px1,py1,pz1)

     if (mod(itime,itest)==0) then
        !does not matter --> output=DIV U=0 (in dv3)
        call divergence(ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
             td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,dv3,&
             nxmsize,nymsize,nzmsize,ph1,ph3,ph4,2)

        call test_speed_min_max(ux1,uy1,uz1)
        if (iscalar==1) call test_scalar_min_max(phi1)
        call compute_cfl(ux1,uy1,uz1)
     endif

     enddo

#ifdef STATS
     call CONVERGENCE_STATISTIC(ux1,ep1,u1sum_tik,u1sum_tak,tsum)
     call CONVERGENCE_STATISTIC2(ux1,ep1,tik1,tik2,tak1,tak2)
#endif
#ifdef FORCES
     call force(ux1,uy1,uz1,ux03,ux13,uy03,uy13,ep1,epcv3,pp3,nzmsize,phG,ph2,ph3)
     if (mod(itime,isave).eq.0) call restart_forces(1)
#endif
#ifdef POST
     if (mod(itime,iprocessing).eq.0)call VISU_PRE(pp3,ta1,tb1,di1,ta2,tb2,di2,ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu,pre1)
     !if (itime.gt.istatsig.AND.mod(itime,iskips).eq.0) !moved inside postprocessing
     call postprocessing(ux1,uy1,uz1,phi1,ep1,pre1,diss1)
     if(ifprobes.gt.0)call write_probes(ux1,uy1,uz1,pre1,diss1,phi1)

#endif

!if(itime.gt.ifirst+2)then
! if(mod(itime,itest)==0)call compute_cfl(ux1,uy1,uz1)
! !dt=dt*1.0001
! !call setab3
!endif

    if (mod(itime,itest)==0) then

     call cpu_time(trank)
     if (nrank==0) print *,'Time per this time step (s):',real(trank-t1)
     if (nrank==0) print *,'Snapshot current/final ',int(itime/imodulo),int(ilast/imodulo)

    endif

#ifdef VISU
     if (mod(itime,imodulo).eq.0) then
        call VISU_INSTA (ux1,uy1,uz1,phi1,ep1,.false.)

        if (save_pre.eq.1.OR.save_prem.eq.1) call VISU_PRE (pp3,ta1,tb1,di1,&
             ta2,tb2,di2,ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu,pre1)
     endif
#endif

#ifdef STATS
     if (t.ge.callstat) then
        call cpu_time(t1)
        call OVERALL_STATISTIC(ux1,uy1,uz1,phi1,pre1,diss1,ta1,&
             u1sum,v1sum,w1sum,u2sum,v2sum,w2sum,&
             u3sum,v3sum,w3sum,u4sum,v4sum,w4sum,&
             uvsum,vwsum,vwsum,disssum,tsum,&
             psum,ppsum,upsum,vpsum,wpsum)

        if (save_dudx.eq.1) call EXTRA_STAT (ux1,uy1,uz1,uvisu,tsum,dudxsum,utmapsum)
        call cpu_time(trank)
        if (nrank==0) print *,'Time in statistics (s):',real(trank-t1,4)

     endif
#endif

#ifdef VISU
     if (mod(itime,isave).eq.0) call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
          px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,1)
#endif

     if (mod(itime,itest)==0) then

       call cpu_time(trank)

       telapsed = (trank-tstart)/thirtysixthousand
       tremaining  = telapsed*(ilast-itime)/(itime-ifirst)

       if (nrank==0) then
          write(*,"(' Remaining time:',I8,' h ',I2,' min')") int(tremaining), int((tremaining-int(tremaining))*sixty)
          write(*,"(' Elapsed time:  ',I8,' h ',I2,' min')") int(telapsed), int((telapsed-int(telapsed))*sixty)
       endif

     endif

  enddo

  call cpu_time(trank); ttotal=trank-tstart

  if (nrank==0) then
     print *,'==========================================================='
     print *,''
     print *,'Good job! Incompact3d finished successfully!'
     print *,''
     print *,'2DECOMP with p_row*p_col=',p_row,p_col
     print *,''
     print *,'nx*ny*nz=',nx*ny*nz
     print *,'nx,ny,nz=',nx,ny,nz
     print *,'dx,dy,dz=',dx,dy,dz
     print *,''
     print *,'Averaged time per step (s):',real(ttotal/(ilast-ifirst),4)
     print *,'Total wallclock (s):',real(ttotal,4)
     print *,'Total wallclock (m):',real(ttotal/sixty,4)
     print *,'Total wallclock (h):',real(ttotal/thirtysixthousand,4)
     print *,'Total wallclock (d):',real(ttotal*1.1574e-5,4)
     print *,''
  endif

  call decomp_2d_finalize
  CALL MPI_FINALIZE(code)

end PROGRAM incompact3d
