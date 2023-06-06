MODULE network_route

USE nrtype
! data types
USE dataTypes,   ONLY: STRFLX            ! fluxes in each reach
USE dataTypes,   ONLY: STRSTA            ! states in each reach
USE dataTypes,   ONLY: RCHTOPO           ! Network topology
USE dataTypes,   ONLY: RCHPRP            ! Reach parameter
USE dataTypes,   ONLY: subbasin_omp      ! mainstem+tributary data strucuture
USE perf_mod,    ONLY: t_startf,t_stopf  ! timing start/stop
USE model_utils, ONLY: handle_err        ! error handling

implicit none

private
public:: route_network

ABSTRACT INTERFACE
  SUBROUTINE sub_route_rch(iEns_dummy, iSeg_dummy,    & ! input: array indices
                           ixDesire_dummy,            & ! input: index of verbose reach
                           T0_dummy,T1_dummy,         & ! input: start and end of the time step
                           NETOPO_dummy,              & ! input: reach topology data structure
                           RPARAM_dummy,              & ! input: reach parameter data structure
                           RCHSTA_dummy,              & ! inout: reach state data structure
                           RCHFLX_dummy,              & ! inout: reach flux data structure
                           ierr_dummy,cmessage_dummy)   ! output: error control
    USE nrtype
    USE dataTypes, ONLY: STRFLX            ! fluxes in each reach
    USE dataTypes, ONLY: STRSTA            ! states in each reach
    USE dataTypes, ONLY: RCHTOPO           ! Network topology
    USE dataTypes, ONLY: RCHPRP            ! Reach parameter
    implicit none
    integer(i4b),       intent(in)                 :: iEns_dummy             ! ensemble member
    integer(i4b),       intent(in)                 :: iSeg_dummy             ! ensemble member
    integer(i4b),       intent(in)                 :: ixDesire_dummy         ! index of the reach for verbose output
    real(dp),           intent(in)                 :: T0_dummy,T1_dummy      ! start and end of the time step (seconds)
    type(RCHTOPO),      intent(in),    allocatable :: NETOPO_dummy(:)        ! River Network topology
    type(RCHPRP),       intent(in),    allocatable :: RPARAM_dummy(:)        ! River reach parameter
    type(STRSTA),       intent(inout)              :: RCHSTA_dummy(:,:)      ! reach state data
    type(STRFLX),       intent(inout)              :: RCHFLX_dummy(:,:)      ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
    integer(i4b),       intent(out)                :: ierr_dummy             ! error code
    character(*),       intent(out)                :: cmessage_dummy         ! error message
  END SUBROUTINE sub_route_rch
END INTERFACE

CONTAINS

 ! *********************************************************************
 ! subroutine: route kinematic waves through the river network
 ! *********************************************************************
 SUBROUTINE route_network(route_rch,            & ! input: rech routing subroutine
                          iens,                 & ! input: ensemble index
                          river_basin,          & ! input: river basin information (mainstem, tributary outlet etc.)
                          T0,T1,                & ! input: start and end of the time step
                          ixDesire,             & ! input: reachID to be checked by on-screen pringing
                          NETOPO_in,            & ! input: reach topology data structure
                          RPARAM_in,            & ! input: reach parameter data structure
                          RCHSTA_out,           & ! inout: reach state data structure
                          RCHFLX_out,           & ! inout: reach flux data structure
                          ierr,message,         & ! output: error control
                          ixSubRch)               ! optional input: subset of reach indices to be processed

   implicit none
   ! Argument variables
   procedure(sub_route_rch)                       :: route_rch
   integer(i4b),       intent(in)                 :: iEns                 ! ensemble member
   type(subbasin_omp), intent(in),    allocatable :: river_basin(:)       ! river basin information (mainstem, tributary outlet etc.)
   real(dp),           intent(in)                 :: T0,T1                ! start and end of the time step (seconds)
   integer(i4b),       intent(in)                 :: ixDesire             ! index of the reach for verbose output
   type(RCHTOPO),      intent(in),    allocatable :: NETOPO_in(:)         ! River Network topology
   type(RCHPRP),       intent(in),    allocatable :: RPARAM_in(:)         ! River reach parameter
   type(STRSTA),       intent(inout)              :: RCHSTA_out(:,:)      ! reach state data
   type(STRFLX),       intent(inout)              :: RCHFLX_out(:,:)      ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   integer(i4b),       intent(out)                :: ierr                 ! error code
   character(*),       intent(out)                :: message              ! error message
   integer(i4b),       intent(in), optional       :: ixSubRch(:)          ! subset of reach indices to be processed
   ! local variables
   character(len=strLen)                          :: cmessage             ! error message for downwind routine
   logical(lgt),                      allocatable :: doRoute(:)           ! logical to indicate which reaches are processed
   integer(i4b)                                   :: nOrder               ! number of stream order
   integer(i4b)                                   :: nTrib                ! number of tributary basins
   integer(i4b)                                   :: nSeg                 ! number of reaches in the network
   integer(i4b)                                   :: iSeg, jSeg           ! loop indices - reach
   integer(i4b)                                   :: iTrib                ! loop indices - branch
   integer(i4b)                                   :: ix                   ! loop indices stream order

   ierr=0; message='network_route/'

   nSeg = size(RCHFLX_out(iens,:))

   ! number of reach check
   if (size(NETOPO_in)/=nSeg) then
     ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
   endif

   allocate(doRoute(nSeg), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating space for [doRoute]'; return; endif

   if (present(ixSubRch))then
    doRoute(:)=.false.
    doRoute(ixSubRch) = .true. ! only subset of reaches are on
   else
    doRoute(:)=.true. ! every reach is on
   endif

   nOrder = size(river_basin)

   call t_startf('route_network') ! timing start

   ! routing through river network
   do ix = 1, nOrder
     nTrib=size(river_basin(ix)%branch)
!$OMP PARALLEL DO schedule(dynamic,1)                   & ! chunk size of 1
!$OMP          private(jSeg, iSeg)                      & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(T0,T1)                            & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(doRoute)                          & ! data array shared
!$OMP          shared(NETOPO_in)                        & ! data structure shared
!$OMP          shared(RPARAM_in)                        & ! data structure shared
!$OMP          shared(RCHSTA_out)                       & ! data structure shared
!$OMP          shared(RCHFLX_out)                       & ! data structure shared
!$OMP          shared(ix, iEns, ixDesire)               & ! indices shared
!$OMP          firstprivate(nTrib)
     do iTrib = 1,nTrib
       do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
         jSeg  = river_basin(ix)%branch(iTrib)%segIndex(iSeg)
         if (.not. doRoute(jSeg)) cycle
         call route_rch(iEns,jSeg,           & ! input: array indices
                        ixDesire,            & ! input: index of verbose reach
                        T0,T1,               & ! input: start and end of the time step
                        NETOPO_in,           & ! input: reach topology data structure
                        RPARAM_in,           & ! input: reach parameter data structure
                        RCHSTA_out,          & ! inout: reach state data structure
                        RCHFLX_out,          & ! inout: reach flux data structure
                        ierr,cmessage)         ! output: error control
         if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
       end do ! branch
     end do ! tributary
!$OMP END PARALLEL DO
   end do ! basin loop

   call t_stopf('route_network')

 END SUBROUTINE route_network

END MODULE network_route
