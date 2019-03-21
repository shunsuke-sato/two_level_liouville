include "parallel.f90"
include "communication.f90"
include "random_number/random_number_module.f90" 
module global_variables
  use parallel
  use communication
  use random_number
  implicit none


! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! model paramters
  real(8) :: Egap

! relaxation paramters
  real(8) :: T1_relax, T2_relax
  real(8) :: gamma1, gamma2

  complex(8) :: zH_mat(2,2)

  complex(8) :: zpsi(2)

! parameters for time-propagation
  integer :: nt, nt_cycle
  real(8) :: dt, Tprop
  real(8),allocatable :: tt(:)
  
  integer :: ntraj

! laser paraemter
  real(8) :: E0, omega0, T0
  real(8),allocatable :: Et(:),Et_dt2(:)
  

! Floquet
  logical,parameter :: if_floquet_state = .true.
  integer,parameter :: nmax_floquet = 32
  complex(8),allocatable :: zpsi_F(:,:)
  real(8) :: eps_F(2)

! Photoelectron spectroscopy
  logical,parameter :: if_calc_PES = .true.
  integer :: it_PES_ini, it_PES_fin
  integer :: NE_PES
  real(8),allocatable :: eps_PES(:)
  complex(8),allocatable :: zCt_PES(:)
  real(8),allocatable :: Et_env_PES(:)
  real(8) :: omega_PES, omega_range_PES
  real(8) :: Tpulse_PES

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call init_parallel
  call initialize_random_number_generator_mpi(comm_id_global, comm_nproc_global)

  call input
  call initialize

  if(if_calc_pes)call initialize_PES

  call time_propagation

  call fin_parallel

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  ntraj = 64 !1024

  Egap = 1d0
  T1_relax = 30d0
  T2_relax = 30d0

! L1 =  sigma_- =(0, 0), L2 = -sigma_z = (-1, 0)
!                (1, 0),                 (0, 1)
  gamma1 = 1d0/T1_relax
  gamma2 = 0.5d0*(1d0/T2_relax - 0.5d0*gamma1)

  omega0 = Egap
  E0 = 0.3d0
  T0 = 20d0*2d0*pi/omega0


  Tprop = 60d0*2d0*pi/omega0
  dt = 0.1d0


  write(*,"(A,2x,e26.16e3)")"input   dt=",dt
  nt_cycle = (2d0*pi/omega0)/dt + 1
  dt = (2d0*pi/omega0)/nt_cycle

  write(*,"(A,2x,e26.16e3)")"refined dt=",dt
  nt = Tprop/dt + 1



end subroutine input
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: it

  allocate(tt(0:nt+1))
  do it = 0, nt+1
    tt(it) = it*dt
  end do

  if(if_floquet_state) call calc_floquet_state
  call initialize_laser

end subroutine initialize
!-------------------------------------------------------------------------------
subroutine initialize_laser
  use global_variables
  implicit none
  integer :: it
  real(8) :: xx

  allocate(Et(0:nt+1),Et_dt2(0:nt+1))
  Et = 0d0; Et_dt2 = 0d0

  if(if_floquet_state)then
    do it = 0, nt+1
      xx = tt(it)
      Et(it) = E0*sin(omega0*xx)

      xx = tt(it) + dt*0.5d0
      Et_dt2(it) = E0*sin(omega0*xx)

    end do

  else

    do it = 0, nt+1
      xx = tt(it)
      if(xx < T0)then
        Et(it) = E0*sin(omega0*xx)*sin(0.5d0*pi*xx/T0)**2
      else
        Et(it) = E0*sin(omega0*xx)
      end if

      xx = tt(it) + dt*0.5d0
      if(xx < T0)then
        Et_dt2(it) = E0*sin(omega0*xx)*sin(0.5d0*pi*xx/T0)**2
      else
        Et_dt2(it) = E0*sin(omega0*xx)
      end if

    end do
  end if

end subroutine initialize_laser
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  real(8) :: prob
  integer :: it, ipes
  integer :: itraj
  real(8) :: ss
  real(8),allocatable :: pop(:,:), dip(:)
  real(8),allocatable :: pop_PES(:)

  allocate(pop(0:nt+1,2), dip(0:nt+1))
  if(if_calc_pes)then
    allocate(pop_PES(0:NE_PES))
    pop_PES = 0d0
  end if
  pop = 0d0; dip = 0d0

  do itraj = 1, ntraj
    if(mod(itraj, comm_nproc_global) /= comm_id_global)cycle

    zpsi = 0d0
    zpsi(2) = 1d0
    if(if_calc_pes)zCt_PES = 0d0

    call ranlux_double_mpi(prob)
    it = 0
    ss = sum(abs(zpsi)**2)
    pop(it,1) = pop(it,1) + abs(zpsi(1))**2/ss
    pop(it,2) = pop(it,2) + abs(zpsi(2))**2/ss
    dip(it)   = dip(it)   + 2d0*real(conjg(zpsi(1))*zpsi(2))/ss

    do it = 0,nt

      if(if_calc_pes .and. it>=it_PES_ini .and. it<=it_PES_fin)call dt_evolve_PES_1st_half(it)
      call dt_evolve(it, prob)
      if(if_calc_pes .and. it>=it_PES_ini .and. it<=it_PES_fin)call dt_evolve_PES_2nd_half(it)
      ss = sum(abs(zpsi)**2)
      pop(it+1,1) = pop(it+1,1) + abs(zpsi(1))**2/ss
      pop(it+1,2) = pop(it+1,2) + abs(zpsi(2))**2/ss
      dip(it+1)   = dip(it+1)   + 2d0*real(conjg(zpsi(1))*zpsi(2))/ss

    end do

    if(if_calc_pes)pop_PES = pop_PES + abs(zCt_PES(:))**2

  end do

  call comm_allreduce(pop); pop=pop/ntraj
  call comm_allreduce(dip); dip=dip/ntraj
  if(if_calc_pes)then
    call comm_allreduce(pop_PES); pop_PES=pop_PES/ntraj
  end if

  if(if_root_global)then
    open(21,file='quantities_t.out')
    do it = 0,nt+1
      write(21,"(999e26.16e3)")tt(it),Et(it),pop(it,1),pop(it,2),dip(it)
    end do
    close(21)

    open(20,file='pes_spectrum.out')
    do ipes = 0,NE_PES
      write(20,"(999e26.16e3)")eps_PES(ipes)-omega_PES,pop_PES(ipes)
    end do
    close(20)
  end if

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it, prob)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(inout) :: prob
  complex(8) :: zpsi_t(2), zpsi_rk(2, 4)
  complex(8) :: zHam_mat(2,2)
  complex(8) :: zLpsi_t(2,2)
  real(8) :: Et_tmp
  real(8) :: ss, p1, p2


! t
  Et_tmp = Et(it)

  zHam_mat(1,1) =  0.5d0*Egap -0.5d0*zI*(gamma1+gamma2)
  zHam_mat(2,1) =  Et_tmp
  zHam_mat(1,2) =  Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap -0.5d0*zI*gamma2

! RK1 
  zpsi_t = zpsi
  zpsi_rk(:,1) = matmul(zHam_mat, zpsi_t)


! t+dt/2
  Et_tmp = Et_dt2(it)
  zHam_mat(1,1) =  0.5d0*Egap -0.5d0*zI*(gamma1+gamma2)
  zHam_mat(2,1) =  Et_tmp
  zHam_mat(1,2) =  Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap -0.5d0*zI*gamma2

!RK2
  zpsi_t = zpsi -zI*dt*0.5d0*zpsi_rk(:,1)
  zpsi_rk(:,2) = matmul(zHam_mat, zpsi_t)

!RK3
  zpsi_t = zpsi -zI*dt*0.5d0*zpsi_rk(:,2)
  zpsi_rk(:,3) = matmul(zHam_mat, zpsi_t)

! t
  Et_tmp = Et(it+1)
  zHam_mat(1,1) =  0.5d0*Egap -0.5d0*zI*(gamma1+gamma2)
  zHam_mat(2,1) =  Et_tmp
  zHam_mat(1,2) =  Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap -0.5d0*zI*gamma2

!RK4
  zpsi_t = zpsi -zI*dt*zpsi_rk(:,3)
  zpsi_rk(:,4) = matmul(zHam_mat, zpsi_t)

  zpsi = zpsi -zI*dt/6d0*(zpsi_rk(:,1)&
                     +2d0*zpsi_rk(:,2)&
                     +2d0*zpsi_rk(:,3)&
                         +zpsi_rk(:,4))

  ss = abs(zpsi(1))**2 + abs(zpsi(2))**2
  if(ss > prob)return

  zLpsi_t(1,1) =  0d0    ; zLpsi_t(2,1) =  zpsi(1)
  zLpsi_t(1,2) = -zpsi(1); zLpsi_t(2,2) =  zpsi(2)

  p1 = gamma1*(abs(zLpsi_t(1,1))**2 + abs(zLpsi_t(2,1))**2)
  p2 = gamma2*(abs(zLpsi_t(1,2))**2 + abs(zLpsi_t(2,2))**2)
  ss = p1 + p2
  p1 = p1/ss; p2 = p2/ss
  
  call ranlux_double_mpi(ss)
  if(ss < p1)then
    zpsi = zLpsi_t(:,1)
  else
    zpsi = zLpsi_t(:,2)
  end if
  ss = abs(zpsi(1))**2 + abs(zpsi(2))**2
  zpsi = zpsi/sqrt(ss)

  call ranlux_double_mpi(prob)

end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine diag_2x2(zmat,zvec,lambda)
  implicit none
  complex(8),intent(in)  :: zmat(2,2)
  complex(8),intent(out) :: zvec(2,2)
  real(8),intent(out)     :: lambda(2)
  real(8) :: a,b,ss
  complex(8) :: zc
  complex(8) :: zx,zy
  
  a = real(zmat(1,1))
  b = real(zmat(2,2))
  zc = zmat(2,1)

  if(a>b)then
     lambda(1) = 0.5d0*( (a+b)+sqrt((a-b)**2+4d0*abs(zc)**2) )
     lambda(2) = 0.5d0*( (a+b)-sqrt((a-b)**2+4d0*abs(zc)**2) )
  else
     lambda(1) = 0.5d0*( (a+b)-sqrt((a-b)**2+4d0*abs(zc)**2) )
     lambda(2) = 0.5d0*( (a+b)+sqrt((a-b)**2+4d0*abs(zc)**2) )
  end if
  zy = zc/(lambda(1)-b)
  zx = conjg(zc)/(lambda(2)-a)

! vector 1
  ss = 1d0/sqrt((1d0+abs(zy)**2))
  zvec(1,1) = ss
  zvec(2,1) = ss*zy
  
! vector 2
  ss = 1d0/sqrt((1d0+abs(zx)**2))
  zvec(1,2) = ss*zx
  zvec(2,2) = ss

  
end subroutine diag_2x2
!-------------------------------------------------------------------------------
subroutine initialize_PES
  use global_variables
  implicit none
  real(8) :: wi, wf, dw
  integer :: iw, it
  real(8) :: xx

  if(.not. if_calc_PES) stop 'Error in initialize_PES.'

  NE_PES = 512
  omega_PES = Egap * 100d0
  omega_range_PES = 3d0*Egap

  Tpulse_PES = 20d0*2d0*pi/omega0

  allocate(eps_PES(0:NE_PES), zCt_PES(0:NE_PES))
  allocate(Et_env_PES(0:nt+1))
  zCt_PES = 0d0

  wi = omega_PES - omega_range_PES
  wf = omega_PES + omega_range_PES
  dw = (wf-wi)/NE_PES

  do iw = 0, NE_PES
    
    eps_PES(iw) = wi + dw*iw
    
  end do

  it_PES_ini = -1
  it_PES_fin = -1


  Et_env_PES = 0d0
  do it = 0,nt+1
    if(tt(it) >= Tprop - Tpulse_PES)then
      if(it_PES_ini == -1)it_PES_ini = it
      xx = tt(it) - (Tprop - Tpulse_PES)
      Et_env_PES(it) = sin(pi*xx/tpulse_PES)**2
    else if(tt(it) >= Tprop)then
      if(it_PES_fin == -1)it_PES_fin = it
    end if
  end do
  if(it_PES_fin == -1)it_PES_fin = nt+1
  write(*,"(A,2x,I7)")"it_pes_ini=",it_pes_ini
  write(*,"(A,2x,I7)")"it_pes_fin=",it_pes_fin

  if(if_root_global)then
    open(20,file='Et_env_PES.out')
    do it = 0,nt+1
      write(20,"(999e26.16e3)")tt(it),Et_env_PES(it)
    end do
    close(20)
  end if
  

end subroutine initialize_PES
!-------------------------------------------------------------------------------
subroutine dt_evolve_PES_1st_half(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ipes
  complex(8) :: zs


  do ipes = 0, NE_PES
    zs = 0.5d0*dt*0.5d0*Et_env_PES(it)*zpsi(2) ! |g>

    zCt_PES(ipes) = zCt_PES(ipes) + zs*exp(zI*(eps_PES(ipes) - omega_PES)*tt(it))
  end do
  

end subroutine dt_evolve_PES_1st_half
!-------------------------------------------------------------------------------
subroutine dt_evolve_PES_2nd_half(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ipes
  complex(8) :: zs


  do ipes = 0, NE_PES
    zs = 0.5d0*dt*0.5d0*Et_env_PES(it+1)*zpsi(2) ! |g>

    zCt_PES(ipes) = zCt_PES(ipes) + zs*exp(zI*(eps_PES(ipes) - omega_PES)*tt(it))
  end do
  
end subroutine dt_evolve_PES_2nd_half
!-------------------------------------------------------------------------------
subroutine calc_floquet_state
  use global_variables
  implicit none
  integer :: nmat
  complex(8),allocatable :: zHam_F(:,:)
  complex(8) :: zHam00(2,2),zHam01(2,2),zHam10(2,2)
  integer :: ifloquet, jfloquet
  integer :: i,j

! for LAPACK    =====
  integer :: lwork,info
  complex(8),allocatable :: work(:)
  real(8),allocatable    :: rwork(:),w(:)

  nmat = (2*nmax_floquet + 1)*2

  lwork = min(64, nmat)*max(1,2*nmat-1)
  if(lwork < 1024)lwork = 1024
  allocate(w(nmat))
  allocate(work(lwork), rwork(max(1, 3*nmat-2)))
! for LAPACK    =====


  nmat = (2*nmax_floquet + 1)*2
  allocate(zHam_F(nmat,nmat))

  zHam00 = 0d0
  zHam00(1,1) =  0.5d0*Egap
  zHam00(2,2) = -0.5d0*Egap

  zHam01 = 0d0
  zHam01(2,1) = -0.5d0*zI*E0
  zHam01(1,2) = -0.5d0*zI*E0

  zHam10 = 0d0
  zHam10(2,1) =  0.5d0*zI*E0
  zHam10(1,2) =  0.5d0*zI*E0


  zHam_F = 0d0
  do ifloquet = -nmax_floquet,nmax_floquet

! (m,m)
    jfloquet = ifloquet
    i = 2*(ifloquet+nmax_floquet)+1
    j = 2*(jfloquet+nmax_floquet)+1

    zHam_F(i:i+1,j:j+1) = zHam00(1:2,1:2)
    zHam_F(i,i) = zHam_F(i,i) - ifloquet*omega0
    zHam_F(i+1,i+1) = zHam_F(i+1,i+1) - ifloquet*omega0

! (m,m+1)
    jfloquet = ifloquet + 1
    if(jfloquet <= nmax_floquet)then
      i = 2*(ifloquet+nmax_floquet)+1
      j = 2*(jfloquet+nmax_floquet)+1

      zHam_F(i:i+1,j:j+1) = zHam01(1:2,1:2)
    end if

! (m,m+1)
    jfloquet = ifloquet - 1
    if(jfloquet >= -nmax_floquet)then
      i = 2*(ifloquet+nmax_floquet)+1
      j = 2*(jfloquet+nmax_floquet)+1

      zHam_F(i:i+1,j:j+1) = zHam10(1:2,1:2)
    end if

  end do

  call zheev('V', 'U', nmat, zHam_F, nmat, w, work, lwork, rwork, info)


  allocate(zpsi_F(nmat,2))
  i = 2*(0+nmax_floquet)+1
  zpsi_F(:,1:2) = zHam_F(:,i:i+1)
  eps_F(1:2) = w(i:i+1)

  zpsi = 0d0
  do ifloquet = -nmax_floquet, nmax_floquet
    i = 2*(ifloquet+nmax_floquet)+1

    zpsi(:) = zpsi(:) &
      + zpsi_F(i:i+1,2) ! 1 or 2
  end do


end subroutine calc_floquet_state
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
