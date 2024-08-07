!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module wave_code_data

    implicit none; 
    integer, parameter :: pp = 8; !pp = 4 for 32 bit and pp = 8 for 64 bit

    integer(pp), allocatable, dimension(:) :: flre_cd_ptr; !pointer to kilca's core data object for each mode
    integer(pp), allocatable, dimension(:) :: vac_cd_ptr; !pointer to kilca's core data object for each mode

    character(1024) :: flre_path; 
    character(1024) :: vac_path; 
    integer :: flre_call_ind = 0; 
    integer :: vac_call_ind = 0; 
    integer :: dim_mn; 
    integer, allocatable, dimension(:) :: m_vals, n_vals; 
!radial grid:
    integer :: dim_r; !radial grid dimension
    real(8), allocatable, dimension(:) :: r; !radial grid array

!wave fields (Gauss system units) for a mode:
    complex(8), allocatable, dimension(:) :: Er; !E field in rsp coordinates
    complex(8), allocatable, dimension(:) :: Es; !E field in rsp coordinates
    complex(8), allocatable, dimension(:) :: Ep; !E field in rsp coordinates
    complex(8), allocatable, dimension(:) :: Et; !Etheta field in cyl coordinates
    complex(8), allocatable, dimension(:) :: Ez; !Ez field in cyl coordinates

    complex(8), allocatable, dimension(:) :: Br; !B field in rsp coordinates
    complex(8), allocatable, dimension(:) :: Bs; !B field in rsp coordinates
    complex(8), allocatable, dimension(:) :: Bp; !B field in rsp coordinates
    complex(8), allocatable, dimension(:) :: Bt; !Btheta field in cyl coordinates
    complex(8), allocatable, dimension(:) :: Bz; !Bz field in cyl coordinates

!current densities (Gauss system units) for a mode:
    complex(8), allocatable, dimension(:) :: Jri; !Jr for ions
    complex(8), allocatable, dimension(:) :: Jsi; !Js for ions
    complex(8), allocatable, dimension(:) :: Jpi; !Jp for ions
    complex(8), allocatable, dimension(:) :: Jre; !Jr for electrons
    complex(8), allocatable, dimension(:) :: Jse; !Js for electrons
    complex(8), allocatable, dimension(:) :: Jpe; !Jp for electrons

!background profiles:
    real(8), allocatable, dimension(:) :: q; !safety factor
    real(8), allocatable, dimension(:) :: n; !1/cm^3
    real(8), allocatable, dimension(:) :: Ti; !eV
    real(8), allocatable, dimension(:) :: Te; !eV
    real(8), allocatable, dimension(:) :: Vth; !cm/c
    real(8), allocatable, dimension(:) :: Vz; !cm/c
    real(8), allocatable, dimension(:) :: dPhi0; !electric potential, Gauss units

    real(8), allocatable, dimension(:) :: nui; !ions collision frequency
    real(8), allocatable, dimension(:) :: nue; !electrons collision frequency

    real(8), allocatable, dimension(:) :: B0t; 
    real(8), allocatable, dimension(:) :: B0z; 
    real(8), allocatable, dimension(:) :: B0; 
!misc data for a mode:
    real(8), allocatable, dimension(:) :: kp; 
    real(8), allocatable, dimension(:) :: ks; 
    real(8), allocatable, dimension(:) :: om_E; 
!for a spectrum:
    real(8), allocatable, dimension(:) :: diss_pow_dens; 
! initial background profiles to be loaded from files:
    real(8), allocatable, dimension(:) :: rq, iq; !safety factor
    real(8), allocatable, dimension(:) :: rn, in; !1/cm^3
    real(8), allocatable, dimension(:) :: rTi, iTi; !eV
    real(8), allocatable, dimension(:) :: rTe, iTe; !eV
    real(8), allocatable, dimension(:) :: rVth, iVth; !cm/c
    real(8), allocatable, dimension(:) :: rVz, iVz; !cm/c
    real(8), allocatable, dimension(:) :: rep, idPhi0; !electric potential, Gauss units
    double precision :: antenna_factor

    character(1024) :: path2profs = './profiles/'; ! path to background profiles
end module