
module grid_mod

    implicit none

    integer :: npoib,npoic,npoi_der,nbaleqs,neqset,iboutype, npoimin, npoi
    integer :: mwind
    double precision :: rmin,rmax
    double precision :: gg_factor, gg_width, gg_r_res;
    double precision :: rb_cut_in, re_cut_in, rb_cut_out, re_cut_out
    integer,          dimension(:),   allocatable :: ipbeg,ipend
    double precision, dimension(:),   allocatable :: rb,rc,Sb,Sc
    double precision, dimension(:),   allocatable :: sqg_bthet_overc,Ercov
    double precision, dimension(:),   allocatable :: sqg_bthet_overcavg,Ercovavg
    double precision, dimension(:),   allocatable :: y,dery,dery_equisource
    double precision, dimension(:),   allocatable :: dae11,dae12,dae22
    double precision, dimension(:),   allocatable :: dai11,dai12,dai22
    double precision, dimension(:),   allocatable :: dni22,visca,polforce
    ! ql RMP-induced diffusion
    double precision, dimension(:),   allocatable :: dqle11,dqle12,dqle21,dqle22
    double precision, dimension(:),   allocatable :: dqli11,dqli12,dqli21,dqli22
    ! 1/nu diffusion
    double precision, dimension(:),   allocatable :: Donue11,Donue12,Donue21,Donue22
    double precision, dimension(:),   allocatable :: Donui11,Donui12,Donui21,Donui22
    double precision, dimension(:),   allocatable :: init_Donue11, init_Donue12, init_Donue21, init_Donue22
    double precision, dimension(:),   allocatable :: init_Donui11, init_Donui12, init_Donui21, init_Donui22
    ! total diffusion
    double precision, dimension(:),   allocatable :: de11,de12,de21,de22
    double precision, dimension(:),   allocatable :: di11,di12,di21,di22
    double precision, dimension(:),   allocatable :: d11_misalign !> diffusion due to misalignment of equipotential and flux surfaces
    double complex, dimension(:),   allocatable :: Es_pert_flux !> part of Es from perturbed flux surface
    double precision, dimension(:),   allocatable :: qlheat_e,qlheat_i
    double precision, dimension(:),   allocatable :: cneo,gpp_av
    double precision, dimension(:,:), allocatable :: deriv_coef,reint_coef

    double precision, dimension(:,:), allocatable :: fluxes_dif,fluxes_con,fluxes_con_nl
    double precision, dimension(:),   allocatable :: source_term
    double precision, dimension(:),   allocatable :: Ercov_lin 
    double precision, dimension(:),   allocatable :: r_resonant

    double precision, dimension(:), allocatable :: dummy


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