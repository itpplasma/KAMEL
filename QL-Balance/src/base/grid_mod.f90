
module grid_mod

    use QLBalance_kinds, only: dp

    implicit none

    integer :: npoib,npoic,npoi_der,nbaleqs,neqset,iboutype, npoimin, npoi
    integer :: mwind
    real(dp) :: rmin,rmax
    real(dp) :: gg_factor, gg_width, gg_r_res;
    real(dp) :: rb_cut_in, re_cut_in, rb_cut_out, re_cut_out
    complex(dp) :: Ipar

    integer, dimension(:),   allocatable :: ipbeg,ipend
    real(dp), dimension(:),   allocatable :: rb,rc,Sb,Sc
    real(dp), dimension(:),   allocatable :: sqrt_g_times_B_theta_over_c,Ercov
    real(dp), dimension(:),   allocatable :: sqg_bthet_overcavg,Ercovavg
    real(dp), dimension(:),   allocatable :: y,dery,dery_equisource
    real(dp), dimension(:),   allocatable :: dae11,dae12,dae22
    real(dp), dimension(:),   allocatable :: dai11,dai12,dai22
    real(dp), dimension(:),   allocatable :: dni22,visca,polforce, polforce_ql
    real(dp), dimension(:), allocatable :: T_EM_phi_e, T_EM_phi_i
    real(dp), dimension(:), allocatable :: T_EM_phi_e_source, T_EM_phi_i_source
    ! ql RMP-induced diffusion
    real(dp), dimension(:),   allocatable :: dqle11,dqle12,dqle21,dqle22
    real(dp), dimension(:),   allocatable :: dqli11,dqli12,dqli21,dqli22
    ! 1/nu diffusion
    real(dp), dimension(:),   allocatable :: Donue11,Donue12,Donue21,Donue22
    real(dp), dimension(:),   allocatable :: Donui11,Donui12,Donui21,Donui22
    real(dp), dimension(:),   allocatable :: init_Donue11, init_Donue12, init_Donue21, init_Donue22
    real(dp), dimension(:),   allocatable :: init_Donui11, init_Donui12, init_Donui21, init_Donui22
    ! total diffusion
    real(dp), dimension(:),   allocatable :: de11,de12,de21,de22
    real(dp), dimension(:),   allocatable :: di11,di12,di21,di22
    real(dp), dimension(:),   allocatable :: d11_misalign !> diffusion due to misalignment of equipotential and flux surfaces
    double complex, dimension(:),   allocatable :: Es_pert_flux !> part of Es from perturbed flux surface
    real(dp), dimension(:),   allocatable :: qlheat_e,qlheat_i
    real(dp), dimension(:),   allocatable :: cneo,gpp_av
    real(dp), dimension(:,:), allocatable :: deriv_coef,reint_coef

    real(dp), dimension(:,:), allocatable :: fluxes_dif,fluxes_con,fluxes_con_nl
    real(dp), dimension(:),   allocatable :: source_term
    real(dp), dimension(:),   allocatable :: Ercov_lin 
    real(dp), dimension(:),   allocatable :: r_resonant

    real(dp), dimension(:), allocatable :: dummy


    contains

    subroutine set_boundary_condition

        implicit none

        ! boundary condition
        if (iboutype .eq. 1) then
            npoi = npoic - 1
        else
            npoi = npoic
        end if

    end subroutine

end module grid_mod
