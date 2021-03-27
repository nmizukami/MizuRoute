MODULE write_restart_pio

! Moudle wide shared data/external routines
USE nrtype

USE var_lookup,        ONLY: ixStateDims, nStateDims
USE var_lookup,        ONLY: ixIRF, nVarsIRF
USE var_lookup,        ONLY: ixIRFbas, nVarsIRFbas
USE var_lookup,        ONLY: ixKWE, nVarsKWE
USE var_lookup,        ONLY: ixKWT, nVarsKWT
USE var_lookup,        ONLY: ixIRFbas, nVarsIRFbas

USE dataTypes,         ONLY: STRFLX            ! fluxes in each reach
USE dataTypes,         ONLY: STRSTA            ! state in each reach
USE dataTypes,         ONLY: RCHTOPO           ! Network topology
USE dataTypes,         ONLY: states
USE dataTypes,         ONLY: time

USE public_var,        ONLY: iulog             ! i/o logical unit number
USE public_var,        ONLY: integerMissing
USE public_var,        ONLY: realMissing
USE public_var,        ONLY: verySmall
USE public_var,        ONLY: rpntfil           ! ascii containing last restart file (used in coupled mode)

USE globalData,        ONLY: meta_stateDims  ! states dimension meta
USE globalData,        ONLY: meta_irf        ! IRF routing
USE globalData,        ONLY: meta_irf_bas
USE globalData,        ONLY: meta_kwt
USE globalData,        ONLY: meta_kwe

USE globalData,        ONLY: pid, nNodes
USE globalData,        ONLY: masterproc
USE globalData,        ONLY: mpicom_route
USE globalData,        ONLY: pio_netcdf_format
USE globalData,        ONLY: pio_typename
USE globalData,        ONLY: pio_numiotasks
USE globalData,        ONLY: pio_rearranger
USE globalData,        ONLY: pio_root
USE globalData,        ONLY: pio_stride

USE nr_utility_module, ONLY: arth
USE pio_utils

implicit none

type(iosystem_desc_t),save :: pioSystemState
type(file_desc_t),    save :: pioFileDescState     ! contains data identifying the file
type(io_desc_t),      save :: iodesc_state_int
type(io_desc_t),      save :: iodesc_state_double
type(io_desc_t),      save :: iodesc_wave_int
type(io_desc_t),      save :: iodesc_wave_double
type(io_desc_t),      save :: iodesc_mesh_double
type(io_desc_t),      save :: iodesc_irf_double
type(io_desc_t),      save :: iodesc_vol_double
type(io_desc_t),      save :: iodesc_irf_bas_double

integer(i4b),         save :: kTime                ! Time index in restart netCDF

integer(i4b),    parameter :: recordDim=-999       ! record dimension Indicator
integer(i4b),    parameter :: currTimeStep = 1
integer(i4b),    parameter :: nextTimeStep = 2

private

public::restart_fname
public::restart_output
public::main_restart     ! used for stand-alone

CONTAINS

 ! *********************************************************************
 ! public subroutine: restart write main routine
 ! *********************************************************************
 SUBROUTINE main_restart(ierr, message)

  implicit none
  ! output variables
  integer(i4b),   intent(out)          :: ierr             ! error code
  character(*),   intent(out)          :: message          ! error message
  ! local variables
  logical(lgt)                         :: restartAlarm     ! logical to make alarm for restart writing
  character(len=strLen)                :: cmessage         ! error message of downwind routine

  ierr=0; message='main_restart/'

  call restart_alarm(restartAlarm, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (restartAlarm) then
    call restart_output(ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

 END SUBROUTINE main_restart


 ! *********************************************************************
 ! private subroutine: restart alarming
 ! *********************************************************************
 SUBROUTINE restart_alarm(restartAlarm, ierr, message)

   USE public_var,        ONLY: calendar
   USE public_var,        ONLY: restart_write  ! restart write options
   USE public_var,        ONLY: restart_day
   USE globalData,        ONLY: restCal        ! restart Calendar time
   USE globalData,        ONLY: dropCal        ! restart drop off Calendar time
   USE globalData,        ONLY: modTime        ! previous and current model time
   ! external routine
   USE time_utils_module, ONLY: ndays_month    ! compute number of days in a month

   implicit none

   ! output
   logical(lgt),   intent(out)          :: restartAlarm     ! logical to make alarm for restart writing
   integer(i4b),   intent(out)          :: ierr             ! error code
   character(*),   intent(out)          :: message          ! error message
   ! local variables
   character(len=strLen)                :: cmessage         ! error message of downwind routine
   integer(i4b)                         :: nDays            ! number of days in a month

   ierr=0; message='restart_alarm/'

   ! adjust restart dropoff day if the dropoff day is outside number of days in particular month
   dropCal%id=restart_day
   call ndays_month(modTime(1)%iy, modTime(1)%im, calendar, nDays, ierr, cmessage)
   if (dropCal%id > nDays) then
     dropCal%id=nDays
   end if

   ! adjust dropoff day further if restart day is actually outside number of days in a particular month
   if (restCal%id > nDays) then
     dropCal%id=dropCal%id-1
   end if

   select case(trim(restart_write))
     case('Specified','specified','Last','last')
       restartAlarm = (dropCal%iy==modTime(1)%iy .and. dropCal%im==modTime(1)%im .and. dropCal%id==modTime(1)%id .and. &
                       dropCal%ih==modTime(1)%ih .and. dropCal%imin==modTime(1)%imin .and. nint(dropCal%dsec)==nint(modTime(1)%dsec))
     case('Annual','annual')
       restartAlarm = (dropCal%im==modTime(1)%im .and. dropCal%id==modTime(1)%id .and. &
                       dropCal%ih==modTime(1)%ih .and. dropCal%imin==modTime(1)%imin .and. nint(dropCal%dsec)==nint(modTime(1)%dsec))
     case('Monthly','monthly')
       restartAlarm = (dropCal%id==modTime(1)%id .and. &
                       dropCal%ih==modTime(1)%ih .and. dropCal%imin==modTime(1)%imin .and. nint(dropCal%dsec)==nint(modTime(1)%dsec))
     case('Daily','daily')
       restartAlarm = (dropCal%ih==modTime(1)%ih .and. dropCal%imin==modTime(1)%imin .and. nint(dropCal%dsec)==nint(modTime(1)%dsec))
     case('Never','never')
       restartAlarm = .false.
     case default
       ierr=20; message=trim(message)//'Current accepted <restart_write> options: L[l]ast, N[n]ever, S[s]pecified, Annual, Monthly, or Daily '; return
   end select

 END SUBROUTINE restart_alarm


 ! *********************************************************************
 ! Public subroutine: output restart netcdf
 ! *********************************************************************
 SUBROUTINE restart_output(ierr, message)

  USE public_var, ONLY: restart_dir
  USE globalData, ONLY: modTime           ! previous and current model time

  implicit none
  ! output variables
  integer(i4b),   intent(out)          :: ierr             ! error code
  character(*),   intent(out)          :: message          ! error message
  ! local variables
  character(len=strLen)                :: cmessage         ! error message of downwind routine
  character(len=strLen)                :: fnameRestart     ! name of the restart file name
  character(len=50),parameter          :: fmtYMDHMS = '(2a,I0.4,a,I0.2,a,I0.2,x,I0.2,a,I0.2,a,I0.2)'

  ierr=0; message='restart_output/'

  if (masterproc) then
    write(iulog,fmtYMDHMS) new_line('a'),'Write restart file at ', &
                           modTime(1)%iy,'-',modTime(1)%im, '-', modTime(1)%id, &
                           modTime(1)%ih,':',modTime(1)%imin,':',nint(modTime(1)%dsec)
  end if

  call restart_fname(fnameRestart, nextTimeStep, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call define_state_nc(fnameRestart, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call write_state_nc(fnameRestart, ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  open (1, file = trim(restart_dir)//trim(rpntfil), status='replace', action='write')
  write(1,*) trim(fnameRestart)
  close(1)

 END SUBROUTINE restart_output


 ! *********************************************************************
 ! public subroutine: define restart NetCDF file name
 ! *********************************************************************
 SUBROUTINE restart_fname(fnameRestart, timeStamp, ierr, message)

   USE public_var,          ONLY: restart_dir
   USE public_var,          ONLY: case_name        ! simulation name ==> output filename head
   USE public_var,          ONLY: calendar
   USE public_var,          ONLY: secprday
   USE public_var,          ONLY: dt
   USE globalData,          ONLY: modJulday        ! current model Julian day
   USE process_time_module, ONLY: conv_julian2cal  ! compute data and time from julian day

   implicit none

   ! input
   integer(i4b),   intent(in)           :: timeStamp        ! optional:
   ! output
   character(*),   intent(out)          :: fnameRestart     ! name of the restart file name
   integer(i4b),   intent(out)          :: ierr             ! error code
   character(*),   intent(out)          :: message          ! error message
   ! local variables
   character(len=strLen)                :: cmessage         ! error message of downwind routine
   real(dp)                             :: timeStampJulday  ! Julidan days corresponding to file name time stamp
   integer(i4b)                         :: sec_in_day       ! second within day
   type(time)                           :: timeStampCal     ! calendar date at next time step (for restart file name)
   character(len=50),parameter          :: fmtYMDS='(a,I0.4,a,I0.2,a,I0.2,a,I0.5,a)'

   ierr=0; message='restart_fname/'

   select case(timeStamp)
     case(currTimeStep); timeStampJulday = modJulday
     case(nextTimeStep); timeStampJulday = modJulday + dt/secprday
     case default;       ierr=20; message=trim(message)//'time stamp option in restart filename: invalid -> 1: current time Step or 2: next time step'; return
   end select

   call conv_julian2cal(timeStampJulday, calendar, timeStampCal, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   sec_in_day = timeStampCal%ih*60*60+timeStampCal%imin*60+nint(timeStampCal%dsec)

   write(fnameRestart, fmtYMDS) trim(restart_dir)//trim(case_name)//'.mizuRoute.r.', &
                                timeStampCal%iy, '-', timeStampCal%im, '-', timeStampCal%id, '-',sec_in_day,'.nc'

 END SUBROUTINE restart_fname


 ! *********************************************************************
 ! subroutine: define restart NetCDF file
 ! *********************************************************************
 SUBROUTINE define_state_nc(fname,           &  ! input: filename
                            ierr, message)      ! output: error control

 USE public_var, ONLY: time_units               ! time units (seconds, hours, or days)
 USE public_var, ONLY: calendar                 ! calendar name
 USE public_var, ONLY: routOpt
 USE public_var, ONLY: allRoutingMethods
 USE public_var, ONLY: doesBasinRoute
 USE public_var, ONLY: kinematicWave
 USE public_var, ONLY: kinematicWaveEuler
 USE public_var, ONLY: impulseResponseFunc
 USE globalData, ONLY: rch_per_proc             ! number of reaches assigned to each proc (size = num of procs+1)
 USE globalData, ONLY: nRch                     ! number of ensembles and river reaches

 implicit none
 ! input variables
 character(*),   intent(in)      :: fname            ! filename
 ! output variables
 integer(i4b),   intent(out)     :: ierr             ! error code
 character(*),   intent(out)     :: message          ! error message
 ! local variables
 integer(i4b)                    :: ix1, ix2         ! frst and last indices of global array for local array chunk
 integer(i4b)                    :: ixRch(nRch)      ! global reach index used for output
 integer(i4b)                    :: jDim             ! loop index for dimension
 integer(i4b)                    :: ixDim_common(4)  ! custom dimension ID array
 character(len=strLen)           :: cmessage         ! error message of downwind routine

 ! Initialize error control
 ierr=0; message='define_state_nc/'

 ! Initialize netCDF time index
 kTime = 0

 associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimId,     &
           dim_ens     => meta_stateDims(ixStateDims%ens)%dimId,     &
           dim_time    => meta_stateDims(ixStateDims%time)%dimId,    &
           dim_tbound  => meta_stateDims(ixStateDims%tbound)%dimId)

 ! ----------------------------------
 ! pio initialization for restart netCDF
 ! ----------------------------------
 pio_numiotasks = nNodes/pio_stride
 call pio_sys_init(pid, mpicom_route,          & ! input: MPI related parameters
                   pio_stride, pio_numiotasks, & ! input: PIO related parameters
                   pio_rearranger, pio_root,   & ! input: PIO related parameters
                   pioSystemState)               ! output: PIO system descriptors

 ! ----------------------------------
 ! Create file
 ! ----------------------------------
 call createFile(pioSystemState, trim(fname), pio_typename, pio_netcdf_format, pioFileDescState, ierr, cmessage)
 if(ierr/=0)then; message=trim(cmessage)//'cannot create state netCDF'; return; endif

 ! For common dimension/variables - seg id, time, time-bound -----------
 ixDim_common = [ixStateDims%seg, ixStateDims%ens, ixStateDims%time, ixStateDims%tbound]

 ! ----------------------------------
 ! Define dimensions
 ! ----------------------------------
 do jDim = 1,size(ixDim_common)
  associate(ixDim => ixDim_common(jDim))
  if (meta_stateDims(ixDim)%dimLength == integerMissing) then
   call set_dim_len(ixDim, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//' for '//trim(meta_stateDims(ixDim)%dimName); return; endif
  endif
  call def_dim(pioFileDescState, trim(meta_stateDims(ixDim)%dimName), meta_stateDims(ixDim)%dimLength, meta_stateDims(ixDim)%dimId)
  end associate
 end do

 ! ----------------------------------
 ! Define variable
 ! ----------------------------------
 call def_var(pioFileDescState, 'reachID', [dim_seg], ncd_int, ierr, cmessage, vdesc='reach ID',  vunit='-' )
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call def_var(pioFileDescState, 'time', [dim_time], ncd_float, ierr, cmessage, vdesc='time', vunit=trim(time_units), vcal=calendar)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call def_var(pioFileDescState, 'time_bound', [dim_tbound, dim_time], ncd_float, ierr, cmessage, vdesc='time bound at last time step', vunit='sec')
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end associate

 ! Routing specific variables --------------
 ! basin IRF
 if (doesBasinRoute == 1) then
  call define_IRFbas_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! KWT routing
 if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
  call define_KWT_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! KWE routing
 if (routOpt==kinematicWaveEuler) then
  call define_KWE_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! IRF routing
 if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
  call define_IRF_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! ----------------------------------
 ! pio initialization of decomposition
 ! ----------------------------------
 associate(nSeg     => meta_stateDims(ixStateDims%seg)%dimLength,     &
           nEns     => meta_stateDims(ixStateDims%ens)%dimLength,     &
           ntdh     => meta_stateDims(ixStateDims%tdh)%dimLength,     & ! maximum future q time steps among basins
           ntdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimLength, & ! maximum future q time steps among reaches
           nTbound  => meta_stateDims(ixStateDims%tbound)%dimLength,  & ! time bound
           nFdmesh  => meta_stateDims(ixStateDims%fdmesh)%dimLength,  & ! finite difference mesh points
           nWave    => meta_stateDims(ixStateDims%wave)%dimLength)      ! maximum waves allowed in a reach

 if (masterproc) then
   ix1 = 1
 else
   ix1 = sum(rch_per_proc(-1:pid-1))+1
 endif
 ix2 = sum(rch_per_proc(-1:pid))
 ixRch = arth(1,1,nSeg)

 ! type: float  dim: [dim_seg, dim_ens, dim_time]  -- channel runoff coming from hru
 call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                 ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nSeg,nEns],            & ! input: dimension length == global array size
                 ixRch(ix1:ix2),         & ! input:
                 iodesc_state_double)

 ! type: int  dim: [dim_seg, dim_ens, dim_time]  -- number of wave or uh future time steps
 call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                 ncd_int,                & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nSeg,nEns],            & ! input: dimension length == global array size
                 ixRch(ix1:ix2),         & ! input:
                 iodesc_state_int)

 if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
   ! type: int, dim: [dim_seg, dim_wave, dim_ens, dim_time]
   call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                   ncd_int,                & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                   [nSeg,nWave,nEns],      & ! input: dimension length == global array size
                   ixRch(ix1:ix2),         & ! input:
                   iodesc_wave_int)

   ! type: float, dim: [dim_seg, dim_wave, dim_ens, dim_time]
   call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                   ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                   [nSeg,nWave,nEns],      & ! input: dimension length == global array size
                   ixRch(ix1:ix2),         & ! input:
                   iodesc_wave_double)
 end if

 if (routOpt==kinematicWaveEuler) then
   ! type: float, dim: [dim_seg, dim_fdmesh, dim_ens, dim_time]
   call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                   ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                   [nSeg,nFdmesh,nEns],    & ! input: dimension length == global array size
                   ixRch(ix1:ix2),         & ! input:
                   iodesc_mesh_double)
 end if

 if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
   ! type: float dim: [dim_seg, dim_tdh_irf, dim_ens, dim_time]
   call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                   ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                   [nSeg,ntdh_irf,nEns],   & ! input: dimension length == global array size
                   ixRch(ix1:ix2),         & ! input:
                   iodesc_irf_double)

   ! type: float dim: [dim_seg, dim_ens, dim_time]
   call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                   ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                   [nSeg,nTbound,nEns],    & ! input: dimension length == global array size
                   ixRch(ix1:ix2),         & ! input:
                   iodesc_vol_double)

 end if
 ! type: float dim: [dim_seg, dim_tdh_irf, dim_ens, dim_time]
 call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                 ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nSeg,ntdh,nEns],       & ! input: dimension length == global array size
                 ixRch(ix1:ix2),         & ! input:
                 iodesc_irf_bas_double)

 end associate

 ! Finishing up definition -------
 call end_def(pioFileDescState, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 CONTAINS

  SUBROUTINE set_dim_len(ixDim, ierr, message1)
   ! populate state netCDF dimension size
   USE public_var, ONLY: MAXQPAR
   USE globalData, ONLY: FRAC_FUTURE     ! To get size of q future for basin IRF
   USE globalData, ONLY: nEns, nRch      ! number of ensembles and river reaches
   implicit none
   ! input
   integer(i4b), intent(in)   :: ixDim    ! ixDim
   ! output
   integer(i4b), intent(out)  :: ierr     ! error code
   character(*), intent(out)  :: message1  ! error message

   ierr=0; message1='set_dim_len/'

   select case(ixDim)
    case(ixStateDims%time);    meta_stateDims(ixStateDims%time)%dimLength    = recordDim
    case(ixStateDims%seg);     meta_stateDims(ixStateDims%seg)%dimLength     = nRch
    case(ixStateDims%ens);     meta_stateDims(ixStateDims%ens)%dimLength     = nEns
    case(ixStateDims%tbound);  meta_stateDims(ixStateDims%tbound)%dimLength  = 2
    case(ixStateDims%tdh);     meta_stateDims(ixStateDims%tdh)%dimLength     = size(FRAC_FUTURE)
    case(ixStateDims%tdh_irf); meta_stateDims(ixStateDims%tdh_irf)%dimLength = 20   !just temporarily
    case(ixStateDims%fdmesh);  meta_stateDims(ixStateDims%fdmesh)%dimLength  = 4
    case(ixStateDims%wave);    meta_stateDims(ixStateDims%wave)%dimLength    = MAXQPAR
    case default; ierr=20; message1=trim(message1)//'unable to identify dimension variable index'; return
   end select

  END SUBROUTINE set_dim_len

  SUBROUTINE define_IRFbas_state(ierr, message1)
   implicit none
   ! output
   integer(i4b), intent(out)         :: ierr          ! error code
   character(*), intent(out)         :: message1      ! error message
   ! local
   integer(i4b)                      :: iVar, ixDim   ! index loop
   integer(i4b)                      :: nDims         ! number of dimensions
   integer(i4b),allocatable          :: dim_IRFbas(:) ! dimension id array

   ierr=0; message1='define_IRFbas_state/'

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%tdh)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%tdh, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%tdh)%dimName); return; endif
   end if

   call def_dim(pioFileDescState, trim(meta_stateDims(ixStateDims%tdh)%dimName), meta_stateDims(ixStateDims%tdh)%dimLength, meta_stateDims(ixStateDims%tdh)%dimId)
   if(ierr/=0)then; ierr=20; message1=trim(message1)//'cannot define dimension'; return; endif

   do iVar=1,nVarsIRFbas

     nDims = size(meta_irf_bas(iVar)%varDim)
     if (allocated(dim_IRFbas)) then
       deallocate(dim_IRFbas)
     end if
     allocate(dim_IRFbas(nDims))
     do ixDim = 1, nDims
       dim_IRFbas(ixDim) = meta_stateDims(meta_irf_bas(iVar)%varDim(ixDim))%dimId
     end do

     call def_var(pioFileDescState, trim(meta_irf_bas(iVar)%varName), dim_IRFbas, meta_irf_bas(iVar)%varType, ierr, cmessage, vdesc=trim(meta_irf_bas(iVar)%varDesc), vunit=trim(meta_irf_bas(iVar)%varUnit))
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   end do

  END SUBROUTINE define_IRFbas_state

  SUBROUTINE define_KWE_state(ierr, message1)
   implicit none
   ! output
   integer(i4b), intent(out)         :: ierr          ! error code
   character(*), intent(out)         :: message1      ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim    ! index loop for variables
   integer(i4b)                      :: nDims         ! number of dimensions
   integer(i4b),allocatable          :: dim_KWE(:)    ! dimension Id array

   ierr=0; message1='define_KWE_state/'

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%fdmesh)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%fdmesh, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%fdmesh)%dimName); return; endif
   end if

   call def_dim(pioFileDescState, trim(meta_stateDims(ixStateDims%fdmesh)%dimName), meta_stateDims(ixStateDims%fdmesh)%dimLength, meta_stateDims(ixStateDims%fdmesh)%dimId)
   if(ierr/=0)then; ierr=20; message1=trim(message1)//'cannot define dimension'; return; endif

   do iVar=1,nVarsKWE

     nDims = size(meta_KWE(iVar)%varDim)
     if (allocated(dim_KWE)) then
       deallocate(dim_KWE)
     end if
     allocate(dim_KWE(nDims))

     do ixDim = 1, nDims
       dim_KWE(ixDim) = meta_stateDims(meta_KWE(iVar)%varDim(ixDim))%dimId
     end do

     call def_var(pioFileDescState, trim(meta_KWE(iVar)%varName), dim_KWE, meta_KWE(iVar)%varType, ierr, cmessage, vdesc=trim(meta_KWE(iVar)%varDesc), vunit=trim(meta_KWE(iVar)%varUnit))
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   end do

  END SUBROUTINE define_KWE_state

  SUBROUTINE define_KWT_state(ierr, message1)
   implicit none
   ! output
   integer(i4b), intent(out)         :: ierr        ! error code
   character(*), intent(out)         :: message1    ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim  ! index loop for variables
   integer(i4b)                      :: nDims       ! number of dimensions
   integer(i4b),allocatable          :: dim_kwt(:)  ! dimensions ID array

   ierr=0; message1='define_KWT_state/'

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%wave)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%wave, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%wave)%dimName); return; endif
   end if

   ! Define dimension needed for this routing specific state variables
   call def_dim(pioFileDescState, trim(meta_stateDims(ixStateDims%wave)%dimName), meta_stateDims(ixStateDims%wave)%dimLength, meta_stateDims(ixStateDims%wave)%dimId)
   if(ierr/=0)then; ierr=20; message1=trim(message1)//'cannot define dimension'; return; endif

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimId,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimId,     &
             dim_time    => meta_stateDims(ixStateDims%time)%dimId,    &
             dim_wave    => meta_stateDims(ixStateDims%wave)%dimId)

   call def_var(pioFileDescState, 'numWaves', [dim_seg,dim_ens,dim_time], ncd_int, ierr, cmessage, vdesc='number of waves in a reach', vunit='-')
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsKWT

     nDims = size(meta_kwt(iVar)%varDim)
     if (allocated(dim_kwt)) then
       deallocate(dim_kwt)
     endif
     allocate(dim_kwt(nDims))
     do ixDim = 1, nDims
       dim_kwt(ixDim) = meta_stateDims(meta_kwt(iVar)%varDim(ixDim))%dimId
     end do

     call def_var(pioFileDescState, trim(meta_kwt(iVar)%varName), dim_kwt, meta_kwt(iVar)%varType, ierr, cmessage, vdesc=trim(meta_kwt(iVar)%varDesc), vunit=trim(meta_kwt(iVar)%varUnit))
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   end do

   end associate

  END SUBROUTINE define_KWT_state

  SUBROUTINE define_IRF_state(ierr, message1)
   implicit none
   ! output
   integer(i4b), intent(out)         :: ierr        ! error code
   character(*), intent(out)         :: message1    ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim  ! index loop for variables
   integer(i4b)                      :: nDims       ! number of dimensions
   integer(i4b),allocatable          :: dim_irf(:)  ! dimensions combination case 4

   ierr=0; message1='define_IRF_state/'

   ! Define dimension needed for this routing specific state variables
   if (meta_stateDims(ixStateDims%tdh_irf)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%tdh_irf, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%tdh_irf)%dimName); return; endif
   endif

   call def_dim(pioFileDescState, trim(meta_stateDims(ixStateDims%tdh_irf)%dimName), meta_stateDims(ixStateDims%tdh_irf)%dimLength, meta_stateDims(ixStateDims%tdh_irf)%dimId)
   if(ierr/=0)then; ierr=20; message1=trim(message1)//'cannot define dimension'; return; endif

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimId,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimId,     &
             dim_time    => meta_stateDims(ixStateDims%time)%dimId,    &
             dim_tdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimId)

   call def_var(pioFileDescState, 'numQF', [dim_seg,dim_ens,dim_time], ncd_int, ierr, cmessage, vdesc='number of future q time steps in a reach', vunit='-')
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsIRF

     nDims = size(meta_irf(iVar)%varDim)
     if (allocated(dim_irf)) then
       deallocate(dim_irf)
     endif
     allocate(dim_irf(nDims))
     do ixDim = 1, nDims
       dim_irf(ixDim) = meta_stateDims(meta_irf(iVar)%varDim(ixDim))%dimId
     end do

     call def_var(pioFileDescState, trim(meta_irf(iVar)%varName), dim_irf, meta_irf(iVar)%varType, ierr, cmessage, vdesc=trim(meta_irf(iVar)%varDesc), vunit=trim(meta_irf(iVar)%varUnit))
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   end do

   end associate

  END SUBROUTINE define_IRF_state

 END SUBROUTINE define_state_nc


 ! *********************************************************************
 ! public subroutine: writing routing state NetCDF file
 ! *********************************************************************
 SUBROUTINE write_state_nc(fname,                &   ! Input: state netcdf name
                           ierr, message)            ! Output: error control

 USE public_var, ONLY: routOpt
 USE public_var, ONLY: doesBasinRoute
 USE public_var, ONLY: allRoutingMethods
 USE public_var, ONLY: kinematicWave
 USE public_var, ONLY: kinematicWaveEuler
 USE public_var, ONLY: impulseResponseFunc
 USE globalData, ONLY: RCHFLX_main         ! mainstem reach fluxes (ensembles, reaches)
 USE globalData, ONLY: RCHFLX_trib         ! tributary reach fluxes (ensembles, reaches)
 USE globalData, ONLY: NETOPO_main         ! mainstem reach topology
 USE globalData, ONLY: NETOPO_trib         ! tributary reach topology
 USE globalData, ONLY: RCHSTA_main         ! mainstem reach state (ensembles, reaches)
 USE globalData, ONLY: RCHSTA_trib         ! tributary reach state (ensembles, reaches)
 USE globalData, ONLY: rch_per_proc        ! number of reaches assigned to each proc (size = num of procs+1)
 USE globalData, ONLY: nRch_mainstem       ! number of mainstem reaches
 USE globalData, ONLY: reachID             ! reach ID in network
 USE globalData, ONLY: nRch                ! number of reaches in network
 USE globalData, ONLY: TSEC                ! beginning/ending of simulation time step [sec]
 USE globalData, ONLY: timeVar             ! time variable
 USE globalData, ONLY: iTime               ! time index

 implicit none

 ! input variables
 character(*), intent(in)        :: fname           ! filename
 ! output variables
 integer(i4b), intent(out)       :: ierr            ! error code
 character(*), intent(out)       :: message         ! error message
 ! local variables
 integer(i4b)                    :: iens            ! temporal
 type(states)                    :: state(0:3)      ! temporal state data structures -currently 2 river routing scheme + basin IRF routing
 type(STRFLX), allocatable       :: RCHFLX_local(:) ! reordered reach flux data structure
 type(RCHTOPO),allocatable       :: NETOPO_local(:) ! reordered topology data structure
 type(STRSTA), allocatable       :: RCHSTA_local(:) ! reordered statedata structure
 character(len=strLen)           :: cmessage        ! error message of downwind routine

 ierr=0; message='write_state_nc/'

 iens = 1
 kTime = kTime + 1

 if (masterproc) then
  associate(nRch_trib => rch_per_proc(0))
  allocate(RCHFLX_local(nRch_mainstem+nRch_trib), &
           NETOPO_local(nRch_mainstem+nRch_trib), &
           RCHSTA_local(nRch_mainstem+nRch_trib), stat=ierr)
  if (nRch_mainstem>0) then
    RCHFLX_local(1:nRch_mainstem) = RCHFLX_main(iens,1:nRch_mainstem)
    NETOPO_local(1:nRch_mainstem) = NETOPO_main(1:nRch_mainstem)
    RCHSTA_local(1:nRch_mainstem) = RCHSTA_main(iens,1:nRch_mainstem)
  end if
   if (nRch_trib>0) then
     RCHFLX_local(nRch_mainstem+1:nRch_mainstem+nRch_trib) = RCHFLX_trib(iens,:)
     NETOPO_local(nRch_mainstem+1:nRch_mainstem+nRch_trib) = NETOPO_trib(:)
     RCHSTA_local(nRch_mainstem+1:nRch_mainstem+nRch_trib) = RCHSTA_trib(iens,:)
   endif
   end associate
 else
  allocate(RCHFLX_local(rch_per_proc(pid)), &
           NETOPO_local(rch_per_proc(pid)), &
           RCHSTA_local(rch_per_proc(pid)),stat=ierr)
  RCHFLX_local = RCHFLX_trib(iens,:)
  NETOPO_local = NETOPO_trib(:)
  RCHSTA_local = RCHSTA_trib(iens,:)
 endif

 ! -- Write out to netCDF

 call openFile(pioSystemState, pioFileDescState, trim(fname),pio_typename, ncd_write, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! Miscellaneous variables - seg id, time etc
 call write_netcdf(pioFileDescState, 'reachID', reachID, [1], [nRch], ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call write_netcdf(pioFileDescState, 'time', [timeVar(iTime)], [kTime], [1], ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call write_netcdf(pioFileDescState, 'time_bound', [TSEC(0),TSEC(1)], [1,kTime], [2,1], ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 if (doesBasinRoute == 1) then
  call write_IRFbas_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
  call write_IRF_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (routOpt==kinematicWaveEuler) then
  call write_KWE_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
  call write_KWT_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 call closeFile(pioFileDescState)

 CONTAINS

  ! Basin IRF writing procedures
  SUBROUTINE write_IRFbas_state(ierr, message1)
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively

  ierr=0; message1='write_IRFbas_state/'

  associate(nSeg     => size(RCHFLX_local),                         &
            nEns     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            ntdh     => meta_stateDims(ixStateDims%tdh)%dimLength)      ! maximum future q time steps among basins

  allocate(state(0)%var(nVarsIRFbas), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRFbas

   select case(iVar)
    case(ixIRFbas%qfuture); allocate(state(0)%var(iVar)%array_3d_dp(nSeg, ntdh, nEns), stat=ierr)
    case(ixIRFbas%q);       allocate(state(0)%var(iVar)%array_2d_dp(nSeg, nEns),       stat=ierr)
    case default; ierr=20; message1=trim(message1)//'unable to identify basin routing variable index'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//'problem allocating space for basin IRF routing state '//trim(meta_irf_bas(iVar)%varName); return; endif

  end do

 ! --Convert data structures to arrays
  do iens=1,nEns
   do iSeg=1,nSeg
     do iVar=1,nVarsIRFbas
      select case(iVar)
       case(ixIRFbas%qfuture); state(0)%var(iVar)%array_3d_dp(iSeg,:,iens) = RCHFLX_local(iSeg)%QFUTURE
       case(ixIRFbas%q);       state(0)%var(iVar)%array_2d_dp(iSeg,iens)   = RCHFLX_local(iSeg)%BASIN_QR(1)
       case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF state variable index'; return
      end select
    enddo
   enddo
  enddo

  do iVar=1,nVarsIRFbas

   select case(iVar)
    case(ixIRFbas%qfuture); call write_pnetcdf_recdim(pioFileDescState, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_3d_dp, iodesc_irf_bas_double, kTime, ierr, cmessage)
    case(ixIRFbas%q);       call write_pnetcdf_recdim(pioFileDescState, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_2d_dp, iodesc_state_double, kTime, ierr, cmessage)
    case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF variable index for nc writing'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  enddo

  end associate

  END SUBROUTINE write_IRFbas_state

  ! KWE writing procedures
  SUBROUTINE write_KWE_state(ierr, message1)
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively

  ierr=0; message1='write_KWE_state/'

  associate(nSeg     => size(RCHFLX_local),                         &
            nEns     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            nMesh    => meta_stateDims(ixStateDims%fdmesh)%dimLength)     ! maximum waves allowed in a reach

  allocate(state(kinematicWaveEuler)%var(nVarsKWE), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWE
    select case(iVar)
     case(ixKWE%a, ixKWE%q)
      allocate(state(kinematicWaveEuler)%var(iVar)%array_3d_dp(nSeg, nMesh, nEns), stat=ierr)
     case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KWE routing state '//trim(meta_kwe(iVar)%varName); return; endif
  end do

  ! --Convert data structures to arrays
  do iens=1,nEns
   do iSeg=1,nSeg

    do iVar=1,nVarsKWE
     select case(iVar)
      case(ixKWE%a)
       state(kinematicWaveEuler)%var(iVar)%array_3d_dp(iSeg,1:4,iens) = RCHSTA_local(iSeg)%EKW_ROUTE%A(:)
      case(ixKWE%q)
       state(kinematicWaveEuler)%var(iVar)%array_3d_dp(iSeg,1:4,iens) = RCHSTA_local(iSeg)%EKW_ROUTE%Q(:)
      case default; ierr=20; message1=trim(message1)//'unable to identify KWE routing state variable index'; return
     end select
    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! Writing netCDF
  do iVar=1,nVarsKWE
    select case(iVar)
     case(ixKWE%a, ixKWE%q)
       call write_pnetcdf_recdim(pioFileDescState, trim(meta_kwe(iVar)%varName), state(kinematicWaveEuler)%var(iVar)%array_3d_dp, iodesc_mesh_double, kTime, ierr, cmessage)
     case default; ierr=20; message1=trim(message1)//'unable to identify KWE variable index for nc writing'; return
    end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  end do

  end associate

  END SUBROUTINE write_KWE_state


  ! KWT writing procedures
  SUBROUTINE write_KWT_state(ierr, message1)
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively
  integer(i4b), allocatable  :: RFvec(:)        ! temporal vector
  integer(i4b), allocatable  :: numWaves(:,:)   ! number of waves for each ensemble and segment

  ierr=0; message1='write_KWT_state/'

  associate(nSeg     => size(RCHFLX_local),                         &
            nEns     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            nWave    => meta_stateDims(ixStateDims%wave)%dimLength)     ! maximum waves allowed in a reach

  allocate(state(kinematicWave)%var(nVarsKWT), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  ! array to store number of wave per segment and ensemble
  allocate(numWaves(nEns,nSeg), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT

    select case(iVar)
     case(ixKWT%routed); allocate(state(kinematicWave)%var(iVar)%array_3d_int(nSeg, nWave, nEns), stat=ierr)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
      allocate(state(kinematicWave)%var(iVar)%array_3d_dp(nSeg, nWave, nEns), stat=ierr)
     case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KWT routing state '//trim(meta_kwt(iVar)%varName); return; endif
  end do

  ! --Convert data structures to arrays
  do iens=1,nEns
   do iSeg=1,nSeg

    numWaves(iens,iseg) = size(RCHSTA_local(iseg)%LKW_ROUTE%KWAVE)

    do iVar=1,nVarsKWT

     select case(iVar)
      case(ixKWT%tentry)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%TI
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%texit)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%TR
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%qwave)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%QF
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%qwave_mod)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%QM
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%routed) ! this is suppposed to be logical variable, but put it as 0 or 1 in double now
       if (allocated(RFvec)) deallocate(RFvec, stat=ierr)
       allocate(RFvec(numWaves(iens,iSeg)),stat=ierr); RFvec=0_i4b
       where (RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%RF) RFvec=1_i4b
       state(kinematicWave)%var(iVar)%array_3d_int(iSeg,1:numWaves(iens,iSeg),iens) = RFvec
       state(kinematicWave)%var(iVar)%array_3d_int(iSeg,numWaves(iens,iSeg)+1:,iens) = integerMissing
      case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
     end select
    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! Writing netCDF
  call write_pnetcdf_recdim(pioFileDescState,'numWaves', numWaves, iodesc_state_int, kTime, ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT

    select case(iVar)
     case(ixKWT%routed)
       call write_pnetcdf_recdim(pioFileDescState, trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_int, iodesc_wave_int, kTime, ierr, cmessage)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
       call write_pnetcdf_recdim(pioFileDescState, trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_dp, iodesc_wave_double, kTime, ierr, cmessage)
     case default; ierr=20; message1=trim(message1)//'unable to identify KWT variable index for nc writing'; return
    end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  end do

  end associate

  END SUBROUTINE write_KWT_state


  ! IRF writing procedures
  SUBROUTINE write_IRF_state(ierr, message1)
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively
  integer(i4b), allocatable  :: numQF(:,:)      ! number of future Q time steps for each ensemble and segment

  ierr=0; message1='write_IRF_state/'

  associate(nSeg     => size(RCHFLX_local),                         &
            nEns     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            nTbound  => meta_stateDims(ixStateDims%tbound)%dimLength, &
            ntdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimLength)    ! maximum future q time steps among reaches

  allocate(state(impulseResponseFunc)%var(nVarsIRF), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  ! array to store number of wave per segment and ensemble
  allocate(numQF(nEns,nSeg), stat=ierr)
  if(ierr/=0)then; message1=trim(message1)//'problem allocating space for numQF'; return; endif

  do iVar=1,nVarsIRF
   select case(iVar)
    case(ixIRF%qfuture); allocate(state(impulseResponseFunc)%var(iVar)%array_3d_dp(nSeg, ntdh_irf, nEns), stat=ierr)
    case(ixIRF%irfVol);  allocate(state(impulseResponseFunc)%var(iVar)%array_3d_dp(nSeg, nTbound, nEns), stat=ierr)
    case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//'problem allocating space for IRF routing state '//trim(meta_irf(iVar)%varName); return; endif
  end do

  ! --Convert data structures to arrays
  do iens=1,nEns
   do iSeg=1,nSeg

    numQF(iens,iseg) = size(NETOPO_local(iSeg)%UH)

    do iVar=1,nVarsIRF

     select case(iVar)
      case(ixIRF%qfuture)
       state(impulseResponseFunc)%var(iVar)%array_3d_dp(iSeg,1:numQF(iens,iSeg),iens) = RCHFLX_local(iSeg)%QFUTURE_IRF
       state(impulseResponseFunc)%var(iVar)%array_3d_dp(iSeg,numQF(iens,iSeg)+1:ntdh_irf,iens) = realMissing
      case(ixIRF%irfVol)
       state(impulseResponseFunc)%var(iVar)%array_3d_dp(iSeg,1:nTbound,iens) = RCHFLX_local(iSeg)%REACH_VOL(0:1)
      case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
     end select

    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! writing netcdf
  call write_pnetcdf_recdim(pioFileDescState, 'numQF', numQF, iodesc_state_int, kTime, ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRF

   select case(iVar)
    case(ixIRF%qfuture)
     call write_pnetcdf_recdim(pioFileDescState, trim(meta_irf(iVar)%varName), state(impulseResponseFunc)%var(iVar)%array_3d_dp, iodesc_irf_double, kTime, ierr, cmessage)
    case(ixIRF%irfVol)
     call write_pnetcdf_recdim(pioFileDescState, trim(meta_irf(iVar)%varName), state(impulseResponseFunc)%var(iVar)%array_3d_dp, iodesc_vol_double, kTime, ierr, cmessage)
    case default; ierr=20; message1=trim(message1)//'unable to identify IRF variable index for nc writing'; return
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end select

  end do

  end associate

  END SUBROUTINE write_IRF_state

 END SUBROUTINE write_state_nc

END MODULE write_restart_pio
