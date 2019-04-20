include "random_number/luxury.f90"
module random_number
  use luxury
  implicit none
  integer :: nproc_mpi

  contains
!--------------------------------------------------------------------
    subroutine initialize_random_number_generator_mpi(irank, nproc_in)
      implicit none
      integer,intent(in) :: irank, nproc_in
      integer :: i
      real(8) :: rvec(1)

      nproc_mpi = nproc_in

      call initialize_random_number_generator

      do i = 1, irank
        call ranlux_double(rvec,1)
      end do


    end subroutine initialize_random_number_generator_mpi
!--------------------------------------------------------------------
    subroutine ranlux_double_mpi(rand_num)
      implicit none
      real(8),intent(out) :: rand_num
      integer :: len4 = 2
      real :: rvec4(2)
      integer(8) :: int1,i
      

      do i = 1,nproc_mpi
        
        CALL RANLUX (rvec4, len4)
        
      end do

      int1 = aint(rvec4(1)*1d6)*10000000 + aint(rvec4(2)*1d7)
      rand_num =dble(int1)*1d-13

    end subroutine ranlux_double_mpi
!--------------------------------------------------------------------
    subroutine initialize_random_number_generator
      implicit none
! parameters for random_number generator
      integer :: lux_ran = 3, K1_ran = 0, K2_ran = 0
      integer :: INT_ran

      INT_ran = 1234567
      CALL RLUXGO(lux_ran,int_ran,K1_ran,K2_ran)


    end subroutine initialize_random_number_generator
!--------------------------------------------------------------------
    subroutine ranlux_double(rvec,len)
      implicit none
      integer,intent(in) :: len
      real(8),intent(out) :: rvec(len)
      integer :: len4 = 2
      real :: rvec4(2)
      integer(8) :: int1,i
      

      do i = 1,len
        
        CALL RANLUX (rvec4, len4)

        int1 = aint(rvec4(1)*1d6)*10000000 + aint(rvec4(2)*1d7)
        rvec(i) =dble(int1)*1d-13
        
      end do

    end subroutine ranlux_double
!--------------------------------------------------------------------
end module random_number


