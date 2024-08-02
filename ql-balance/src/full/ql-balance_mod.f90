!
  module grid_mod
    !implicit none
    !save
    integer :: npoib,npoic,npoi_der,nbaleqs,neqset,iboutype, npoimin, npoi
    integer :: mwind
    double precision :: rmin,rmax
    double precision :: gg_factor, gg_width, gg_r_res;
    double precision :: rb_cut_in,re_cut_in,rb_cut_out,re_cut_out
    integer,          dimension(:),   allocatable :: ipbeg,ipend
    double precision, dimension(:),   allocatable :: rb,rc,Sb,Sc
    double precision, dimension(:),   allocatable :: sqg_bthet_overc,Ercov
    double precision, dimension(:),   allocatable :: sqg_bthet_overcavg,Ercovavg
    double precision, dimension(:),   allocatable :: y,dery,dery_equisource
    double precision, dimension(:),   allocatable :: dae11,dae12,dae22
    double precision, dimension(:),   allocatable :: dai11,dai12,dai22
    double precision, dimension(:),   allocatable :: dni22,visca,polforce
    double precision, dimension(:),   allocatable :: dqle11,dqle12,dqle21,dqle22
    double precision, dimension(:),   allocatable :: dqli11,dqli12,dqli21,dqli22
    double precision, dimension(:),   allocatable :: de11,de12,de21,de22
    double precision, dimension(:),   allocatable :: di11,di12,di21,di22
    double precision, dimension(:),   allocatable :: d11_misalign !> diffusion due to misalignment of equipotential and flux surfaces
    double complex, dimension(:),   allocatable :: Es_pert_flux !> part of Es from perturbed flux surface
    double precision, dimension(:),   allocatable :: qlheat_e,qlheat_i
    double precision, dimension(:),   allocatable :: cneo,gpp_av,qsafb,qsaf
    double precision, dimension(:,:), allocatable :: deriv_coef,reint_coef
    double precision, dimension(:,:), allocatable :: params,dot_params
    double precision, dimension(:,:), allocatable :: ddr_params, ddr_params_nl
    double precision, dimension(:,:), allocatable :: fluxes_dif,fluxes_con,fluxes_con_nl
    double precision, dimension(:,:), allocatable :: params_b
    double precision, dimension(:,:), allocatable :: init_params
    double precision, dimension(:,:), allocatable :: params_lin,params_b_lin
    double precision, dimension(:),   allocatable :: alpha
    double precision, dimension(:),   allocatable :: source_term
    double precision, dimension(:),   allocatable :: Ercov_lin 
    double precision, dimension(:),   allocatable :: r_resonant
  end module grid_mod
!
  module baseparam_mod
    double precision, parameter :: pi=3.14159265358979d0
    double precision,parameter  :: c = 29979245800.0;
    double precision,parameter  :: e_charge=4.8032d-10
    double precision,parameter  :: e_mass=9.1094d-28
    double precision,parameter  :: p_mass=1.6726d-24
    double precision,parameter  :: ev=1.6022d-12
    double precision :: btor,rtor,dperp,Z_i,am,pertamp,omega,rsepar
  end module baseparam_mod
!
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

module diag_mod
logical :: write_diag,write_diag_b
integer :: iunit_diag,iunit_diag_b
integer :: i_mn_loop
logical :: toggle
double precision :: relax_ff
double precision, dimension(:), allocatable :: field_fac_old
end module diag_mod


module parallelTools
    use mpi
    implicit none
    integer :: irank
    integer :: ierror
end module parallelTools