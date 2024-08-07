

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
