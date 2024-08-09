
module balance_mod

    implicit none

    contains

    subroutine balanceInit

        use time_evolution, only: iexit, timescale, tmax, timstep, Nstorage, &
                                  allocate_prev_variables, tmax_factor
        use grid_mod, only: mwind, rmax, rmin, setBoundaryCondition, npoib, rb
        use baseparam_mod, only: dperp
        use diag_mod, only: write_diag, write_diag_b
        use hdf5_tools, only: h5overwrite
        use h5mod, only: mode_m, mode_n
        use control_mod, only: gyro_current_study, write_gyro_current, debug_mode, &
                          ihdf5IO
        use parallelTools, only: initMPI, irank
        use wave_code_data, only: m_vals, n_vals
        use paramscan_mod, only: initialize_parameter_scan_vars, creategroupstructure
        use linear_run, only: init_background_profiles
        use plasma_parameters, only: writeInitialParameters, alloc_hold_parameters

        implicit none

        call read_config

        iexit = 0 ! 0 - don't skip, 1 - skip, 2 - stop
        mwind = 10

        write_diag = .false.
        write_diag_b = .false.

        ! if h5overwrite = true, existing data will be deleted
        ! before new one is written
        ! This is contained in hdf5_tools module
        h5overwrite = .true.
    
        if (gyro_current_study .ne. 0) then
            write_gyro_current = .true.
        else
            write_gyro_current = .false.
        end if

        call initMPI

        timescale = (rmax - rmin)**2/dperp
        tmax = timescale*tmax_factor
        timstep = tmax/Nstorage
    
        if (irank .eq. 0) then
            write(*,*) "timstep = ", timstep
        end if

        call gengrid

        call setBoundaryCondition

        CALL initialize_wave_code_interface(npoib, rb);
        CALL initialize_parameter_scan_vars

        mode_m = m_vals(1)
        mode_n = n_vals(1)
        if (debug_mode) write(*,*) 'Debug: mode_m = ', mode_m, 'mode_n = ', mode_n

        if (ihdf5IO .eq. 1) then
            CALL creategroupstructure
        end if

        call allocate_prev_variables

        call init_background_profiles

        if (irank .eq. 0) then
            CALL writeInitialParameters
            call alloc_hold_parameters
        end if

    end subroutine

end module

module baseparam_mod
    double precision, parameter :: pi=3.14159265358979d0
    double precision,parameter  :: c = 29979245800.0;
    double precision,parameter  :: e_charge=4.8032d-10
    double precision,parameter  :: e_mass=9.1094d-28
    double precision,parameter  :: p_mass=1.6726d-24
    double precision,parameter  :: ev=1.6022d-12
    double precision :: btor,rtor,dperp,Z_i,am,pertamp,omega,rsepar

    double precision :: urelax = 0.5e0 !0.5d0  !0.9d0
    double precision :: tol_max = 3.d-2 !3.d-4 !3.d-3 !3.d-2
    double precision :: factolmax = 3.d0 ! keep
    double precision :: factolred = 0.5d0 ! keep

end module baseparam_mod

module control_mod
    logical :: write_formfactors
    integer :: iwrite
    integer :: ihdf5IO ! added: Markus Markl
    logical :: paramscan ! added: Markus Markl, 03.03.2021
    logical :: timing_mode         ! added by Markus Markl 06.04.2021
    logical :: debug_mode
    logical :: suppression_mode    ! added by Markus Markl 13.04.2021
    logical :: misalign_diffusion ! trigger the calculation and addition of the diffusion due to misaligned equipotentials and flux surfaces
    logical :: diagnostics_output
    logical :: write_gyro_current
    integer :: irf
    integer :: icoll
    double precision :: relax
    integer :: zeitschritt
    integer :: readfromtimestep ! added by Markus Markl 02.06.2021. Reads the background profiles from hdf5 file in which profiles of a ql time evolution are stored.
    integer :: gyro_current_study
    character(len=1024) :: equil_path ! path to equil file containing q, psi, phi,...
    double precision :: eps
    DOUBLE PRECISION :: temperature_limit ! limits ion and electron temperatures from below, in eV
end module control_mod


module matrix_mod
    integer :: isw_rhs
    integer :: nz,nsize
    integer,          dimension(:),   allocatable :: irow,icol
    double precision, dimension(:),   allocatable :: amat,rhsvec
end module matrix_mod

module recstep_mod
    integer :: nstack
    double precision :: tol
    double precision, dimension(:),   allocatable :: tim_stack
    double precision, dimension(:),   allocatable :: timstep_arr
    double precision, dimension(:,:), allocatable :: y_stack
end module recstep_mod

module resonances_mod
    integer :: numres,iunit_res
    double precision, dimension(:), allocatable :: r_res,width_res,ampl_res
end module resonances_mod
