MODULE mc_route_module

! muskingum-cunge routing

USE nrtype
USE dataTypes,     ONLY: STRFLX          ! fluxes in each reach
USE dataTypes,     ONLY: STRSTA          ! state in each reach
USE dataTypes,     ONLY: RCHTOPO         ! Network topology
USE dataTypes,     ONLY: RCHPRP          ! Reach parameter
USE dataTypes,     ONLY: mcRCH           ! MC specific state data structure
USE public_var,    ONLY: iulog           ! i/o logical unit number
USE public_var,    ONLY: realMissing     ! missing value for real number
USE public_var,    ONLY: integerMissing  ! missing value for integer number
USE public_var,    ONLY: dt              ! simulation time step [sec]
USE public_var,    ONLY: qmodOption      ! qmod option (use 1==direct insertion)
USE public_var,    ONLY: is_flux_wm      ! logical water management components fluxes should be read
USE globalData,    ONLY: idxMC           ! routing method index for muskingum method
USE water_balance, ONLY: comp_reach_wb   ! compute water balance error
USE base_route,    ONLY: base_route_rch  ! base (abstract) reach routing method class

implicit none

private
public::mc_route_rch

type, extends(base_route_rch) :: mc_route_rch
 CONTAINS
   procedure, pass :: route => mc_rch
end type mc_route_rch

CONTAINS

 ! *********************************************************************
 ! subroutine: perform muskingum-cunge routing for one segment
 ! *********************************************************************
 SUBROUTINE mc_rch(this,           & ! mc_route_rch object to bound this procedure
                   iEns, segIndex, & ! input: index of runoff ensemble to be processed
                   ixDesire,       & ! input: reachID to be checked by on-screen pringing
                   T0,T1,          & ! input: start and end of the time step
                   NETOPO_in,      & ! input: reach topology data structure
                   RPARAM_in,      & ! input: reach parameter data structure
                   RCHSTA_out,     & ! inout: reach state data structure
                   RCHFLX_out,     & ! inout: reach flux data structure
                   ierr, message)    ! output: error control

 implicit none
 ! Argument variables
 class(mc_route_rch)                       :: this              ! mc_route_rch object to bound this procedure
 integer(i4b),  intent(in)                 :: iEns              ! runoff ensemble to be routed
 integer(i4b),  intent(in)                 :: segIndex          ! segment where routing is performed
 integer(i4b),  intent(in)                 :: ixDesire          ! index of the reach for verbose output
 real(dp),      intent(in)                 :: T0,T1             ! start and end of the time step (seconds)
 type(RCHTOPO), intent(in),    allocatable :: NETOPO_in(:)      ! River Network topology
 type(RCHPRP),  intent(inout), allocatable :: RPARAM_in(:)      ! River reach parameter
 type(STRSTA),  intent(inout)              :: RCHSTA_out(:,:)   ! reach state data
 type(STRFLX),  intent(inout)              :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b),  intent(out)                :: ierr              ! error code
 character(*),  intent(out)                :: message           ! error message
 ! Local variables
 logical(lgt)                              :: verbose           ! check details of variables
 logical(lgt)                              :: isHW              ! headwater basin?
 integer(i4b)                              :: nUps              ! number of upstream segment
 integer(i4b)                              :: iUps              ! upstream reach index
 integer(i4b)                              :: iRch_ups          ! index of upstream reach in NETOPO
 real(dp)                                  :: Qlat              ! lateral flow into channel [m3/s]
 real(dp)                                  :: Qabs              ! maximum allowable water abstraction rate [m3/s]
 real(dp)                                  :: q_upstream        ! total discharge at top of the reach [m3/s]
 real(dp)                                  :: q_upstream_mod    ! total discharge at top of the reach after water abstraction [m3/s]
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ierr=0; message='mc_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   verbose = .true.
 end if

 ! get discharge coming from upstream
 nUps = count(NETOPO_in(segIndex)%goodBas) ! reminder: goodBas is reach with >0 total contributory area
 isHW = .true.
 q_upstream = 0.0_dp
 if (nUps>0) then
   isHW = .false.
   do iUps = 1,nUps
     if (.not. NETOPO_in(segIndex)%goodBas(iUps)) cycle ! skip upstream reach which does not any flow due to zero total contributory areas
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     if (qmodOption==1 .and. RCHFLX_out(iens,iRch_ups)%Qobs/=realMissing) then
       RCHFLX_out(iens, iRch_ups)%ROUTE(idxMC)%REACH_Q = RCHFLX_out(iens,iRch_ups)%Qobs
     end if
     q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxMC)%REACH_Q
   end do
 endif

 q_upstream_mod  = q_upstream
 Qlat = RCHFLX_out(iens,segIndex)%BASIN_QR(1)
 Qabs = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! initial water abstraction (positive) or injection (negative)
 RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! initialize actual water abstraction

 ! update volume at previous time step
 RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(0) = RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1)

 ! Water management - water injection or abstraction (irrigation or industrial/domestic water usage)
 ! For water abstraction, water is extracted from the following priorities:
 ! 1. existing storage(REACH_VOL(0), 2. upstream inflow , 3 lateral flow (BASIN_QR)
 if((RCHFLX_out(iens,segIndex)%REACH_WM_FLUX /= realMissing).and.(is_flux_wm)) then
   if (Qabs > 0) then ! positive == abstraction
     if (RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1)/dt > Qabs) then ! take out all abstraction from strorage
       RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1) - Qabs*dt
     else ! if inital abstraction is greater than volume
       Qabs = Qabs - RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1)/dt ! get residual Qabs after extracting from strorage
       RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1) = 0._dp ! voluem gets 0
       if (q_upstream > Qabs) then ! then take out all residual abstraction from upstream inflow
         q_upstream_mod = q_upstream - Qabs
       else ! if residual abstraction is still greater than lateral flow
         Qabs = Qabs - q_upstream ! get residual abstraction after extracting upstream inflow and storage.
         q_upstream_mod = 0._dp ! upstream inflow gets 0 (all is gone to abstracted flow).
         if (Qlat > Qabs) then ! then take residual abstraction out from lateral flow
           Qlat = Qlat - Qabs
         else ! if residual abstraction is greater than upstream inflow
           Qabs = Qabs - Qlat ! take out residual abstraction from lateral flow
           Qlat = 0._dp ! lateral flow gets 0 (all are gone to abstraction)
           RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX - Qabs
         end if
       end if
     end if
   else ! negative == injection
     Qlat = Qlat - Qabs
   endif
 endif

 if(verbose)then
   write(iulog,'(2A)') new_line('a'), '** Check muskingum-cunge routing **'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,1X,I12,1X,G12.5)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps), &
             RCHFLX_out(iens, iRch_ups)%ROUTE(idxMC)%REACH_Q
     enddo
   end if
   write(iulog,'(A,1X,G12.5)') ' RCHFLX_out(iEns,segIndex)%BASIN_QR(1)=',RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
 endif

 ! solve muskingum-cunge alogorithm
 call muskingum_cunge(RPARAM_in(segIndex),                     & ! input: parameter at segIndex reach
                      T0,T1,                                   & ! input: start and end of the time step
                      q_upstream_mod,                          & ! input: total discharge at top of the reach being processed
                      qlat,                                    & ! input: lateral flow [m3/s]
                      isHW,                                    & ! input: is this headwater basin?
                      RCHSTA_out(iens,segIndex)%MC_ROUTE,      & ! inout:
                      RCHFLX_out(iens,segIndex),               & ! inout: updated fluxes at reach
                      verbose,                                 & ! input: reach index to be examined
                      ierr, cmessage)                            ! output: error control
 if(ierr/=0)then
   write(message, '(A,1X,I10,1X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage); return
 endif

 if(verbose)then
   write(iulog,'(A,1X,G12.5)') ' RCHFLX_out(iens,segIndex)%REACH_Q=', RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_Q
 endif

 if (RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1) < 0) then
   write(iulog,'(A,1X,G12.5,1X,A,1X,I9)') ' ---- NEGATIVE VOLUME = ', RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1), &
         'at ', NETOPO_in(segIndex)%REACHID
 end if

 call comp_reach_wb(NETOPO_in(segIndex)%REACHID, idxMC, q_upstream, RCHFLX_out(iens,segIndex), verbose, lakeFlag=.false.)

 END SUBROUTINE mc_rch

 ! *********************************************************************
 ! subroutine: solve muskingum equation
 ! *********************************************************************
 SUBROUTINE muskingum_cunge(rch_param,     & ! input: river parameter data structure
                            T0,T1,         & ! input: start and end of the time step
                            q_upstream,    & ! input: discharge from upstream
                            q_lat,         & ! input: lateral discharge into chaneel [m3/s]
                            isHW,          & ! input: is this headwater basin?
                            rstate,        & ! inout: reach state at a reach
                            rflux,         & ! inout: reach flux at a reach
                            verbose,       & ! input: reach index to be examined
                            ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Perform muskingum-cunge routing
 !
 !
 !
 !
 ! *
 ! *
 ! *
 !
 ! state array:
 ! (time:0:1, loc:0:1) 0-previous time step/inlet, 1-current time step/outlet.
 ! Q or A(1,2,3,4): 1: (t=0,x=0), 2: (t=0,x=1), 3: (t=1,x=0), 4: (t=1,x=1)

 ! -- EBK 06/26/2023 -- comment out isnan check, doesn't seem to be needed
 ! Use of shr_infnan_isnan will require changes to the standalone build, and
 ! this version is required to work on all compilers.
 !use shr_infnan_mod, only : isnan => shr_infnan_isnan
 implicit none
 ! Argument variables
 type(RCHPRP), intent(in)                 :: rch_param    ! River reach parameter
 real(dp),     intent(in)                 :: T0,T1        ! start and end of the time step (seconds)
 real(dp),     intent(in)                 :: q_upstream   ! total discharge at top of the reach being processed
 real(dp),     intent(in)                 :: q_lat        ! lateral discharge into chaneel [m3/s]
 logical(lgt), intent(in)                 :: isHW         ! is this headwater basin?
 type(mcRCH),  intent(inout)              :: rstate       ! curent reach states
 type(STRFLX), intent(inout)              :: rflux        ! current Reach fluxes
 logical(lgt), intent(in)                 :: verbose      ! reach index to be examined
 integer(i4b), intent(out)                :: ierr         ! error code
 character(*), intent(out)                :: message      ! error message
 ! Local variables
 real(dp)                                 :: alpha        ! sqrt(slope)(/mannings N* width)
 real(dp)                                 :: beta         ! constant, 5/3
 real(dp)                                 :: theta        ! dT/dX
 real(dp)                                 :: X            ! X-factor in descreterized kinematic wave
 real(dp)                                 :: dx           ! length of segment [m]
 real(dp)                                 :: Q(0:1,0:1)   ! discharge at computational molecule
 real(dp)                                 :: Qbar         ! 3-point average discharge [m3/s]
 real(dp)                                 :: Abar         ! 3-point average flow area [m2]
 real(dp)                                 :: Vbar         ! 3-point average velocity [m/s]
 real(dp)                                 :: Ybar         ! 3-point average flow depth [m]
 real(dp)                                 :: B            ! flow top width [m]
 real(dp)                                 :: ck           ! kinematic wave celerity [m/s]
 real(dp)                                 :: Cn           ! Courant number [-]
 real(dp)                                 :: dTsub        ! time inteval for sut time-step [sec]
 real(dp)                                 :: C0,C1,C2     ! muskingum parameters
 real(dp), allocatable                    :: QoutLocal(:) ! out discharge [m3/s] at sub time step
 real(dp), allocatable                    :: QinLocal(:)  ! in discharge [m3/s] at sub time step
 integer(i4b)                             :: ix           ! loop index
 integer(i4b)                             :: ntSub        ! number of sub time-step
 character(len=strLen)                    :: cmessage     ! error message from subroutine
 real(dp), parameter                      :: Y = 0.5      ! muskingum parameter Y (this is fixed)

 ierr=0; message='muskingum-cunge/'

 Q(0,0) = rstate%molecule%Q(1) ! inflow at previous time step (t-1)
 Q(0,1) = rstate%molecule%Q(2) ! outflow at previous time step (t-1)
 Q(1,1) = realMissing

 if (.not. isHW) then

   ! Get the reach parameters
   ! A = (Q/alpha)**(1/beta)
   ! Q = alpha*A**beta
   alpha = sqrt(rch_param%R_SLOPE)/(rch_param%R_MAN_N*rch_param%R_WIDTH**(2._dp/3._dp))
   beta  = 5._dp/3._dp
   dx = rch_param%RLENGTH
   theta = dt/dx     ! [s/m]

   ! compute total flow rate and flow area at upstream end at current time step
   Q(1,0) = q_upstream

   if (verbose) then
     write(iulog,'(A,1X,G12.5)') ' length [m]        =',rch_param%RLENGTH
     write(iulog,'(A,1X,G12.5)') ' slope [-]         =',rch_param%R_SLOPE
     write(iulog,'(A,1X,G12.5)') ' channel width [m] =',rch_param%R_WIDTH
     write(iulog,'(A,1X,G12.5)') ' manning coef [-]  =',rch_param%R_MAN_N
     write(iulog,'(A)')          ' Initial 3 point discharge [m3/s]: '
     write(iulog,'(3(A,1X,G12.5))') ' Qin(t-1) Q(0,0)=',Q(0,0),' Qin(t) Q(1,0)=',Q(1,0),' Qout(t-1) Q(0,1)=',Q(0,1)
   end if

   ! first, using 3-point average in computational molecule, check Cournat number is less than 1, otherwise subcycle within one time step
   Qbar = (Q(0,0)+Q(1,0)+Q(0,1))/3.0  ! average discharge [m3/s]
   Abar = (Qbar/alpha)**(1/beta)      ! average flow area [m2] (from manning equation)
   Vbar = Qbar/Abar                   ! average velocity [m/s]
   ck   = beta*Vbar                   ! kinematic wave celerity [m/s]
   Cn   = ck*theta                    ! Courant number [-]

   ! time-step adjustment so Courant number is less than 1
   ntSub = 1
   dTsub = dt
   if (Cn>1.0_dp) then
     ntSub = ceiling(dt/dx*cK)
     dTsub = dt/ntSub
   end if
   if (verbose) then
     write(iulog,'(A,1X,I3,A,1X,G12.5)') ' No. sub timestep=',nTsub,' sub time-step [sec]=',dTsub
   end if

   allocate(QoutLocal(0:ntSub), QinLocal(0:ntSub), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   QoutLocal(:) = realMissing
   QoutLocal(0) = Q(0,1)        ! outfloe last time step
   QinLocal(0)  = Q(0,0)        ! inflow at last time step
   QinLocal(1:ntSub)  = Q(1,0)  ! infllow at sub-time step in current time step

   ! solve outflow at each sub time step
   do ix = 1, nTsub

     Qbar = (QinLocal(ix)+QinLocal(ix-1)+QoutLocal(ix-1))/3.0 ! 3 point average discharge [m3/s]
     Abar = (Qbar/alpha)**(1/beta)                            ! flow area [m2] (manning equation)
     Ybar = Abar/rch_param%R_WIDTH                            ! flow depth [m] (rectangular channel)
     B = rch_param%R_WIDTH                                    ! top width at water level [m] (rectangular channel)
     ck = beta*(Qbar/Abar)                                    ! kinematic wave celerity [m/s]

     X = 0.5*(1.0 - Qbar/(B*rch_param%R_SLOPE*ck*dX))         ! X factor for descreterized kinematic wave equation
     Cn = ck*dTsub/dx                                         ! Courant number [-]

     C0 = (-X+Cn*(1-Y))/(1-X+Cn*(1-Y))
     C1 = (X+Cn*Y)/(1-X+Cn*(1-Y))
     C2 = (1-X-Cn*Y)/(1-X+Cn*(1-Y))

     QoutLocal(ix) = C0* QinLocal(ix)+ C1* QinLocal(ix-1)+ C2* QoutLocal(ix-1)
     QoutLocal(ix) = max(0.0, QoutLocal(ix))

     ! -- EBK 06/26/2023 -- comment out isnan check, doesn't seem to be needed.
     !if (isnan(QoutLocal(ix))) then
     !  ierr=10; message=trim(message)//'QoutLocal is Nan; activate vodose for this segment for diagnosis';return
     !end if

     if (verbose) then
       write(iulog,'(A,I3,1X,A,G12.5,1X,A,G12.5)') '   sub time-step= ',ix,'Courant number= ',Cn, 'Q= ',QoutLocal(ix)
     end if
   end do

   Q(1,1) = sum(QoutLocal(1:nTsub))/real(nTsub,kind=dp)

 else ! if head-water

   Q(1,0) = 0._dp
   Q(1,1) = 0._dp

   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif

 endif

 ! For very low flow condition, outflow - inflow > current storage, so limit outflow and adjust Q(1,1)
 Q(1,1) = min(rflux%ROUTE(idxMC)%REACH_VOL(1)/dt + Q(1,0)*0.999, Q(1,1))
 rflux%ROUTE(idxMC)%REACH_VOL(1) = rflux%ROUTE(idxMC)%REACH_VOL(1) + (Q(1,0)-Q(1,1))*dt

 ! add catchment flow
 rflux%ROUTE(idxMC)%REACH_Q = Q(1,1)+ q_lat

 if (verbose) then
   write(iulog,'(A,1X,G12.5)') ' Qout(t)=',Q(1,1)
 endif

 ! save inflow (index 1) and outflow (index 2) at current time step
 rstate%molecule%Q(1) = Q(1,0)
 rstate%molecule%Q(2) = Q(1,1)

 END SUBROUTINE muskingum_cunge

END MODULE mc_route_module