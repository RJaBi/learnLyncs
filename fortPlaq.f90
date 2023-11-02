
!! Calculates whole, spatial, temporal plaquettes



module types
        !use, intrinsic ::  ISO_C_BINDING, only: dc=>C_DOUBLE_COMPLEX, dp=>C_DOUBLE
        implicit none
        integer, parameter :: dc = kind((1.0D0,1.0D0)) !! Double precision complex scalars
        integer, parameter :: dp = kind(1.0D0) !! Double precision real

        real(dc) :: r_dc = (1.0_dp, 1.0_dp)
        real(dp) :: r_dp = 1.0_dp

end module

module plaq

contains

  subroutine MultiplyMatMat(MM, left, right)
          use types
    complex(kind=dc), dimension(3, 3), intent(in) :: left, right
    complex(kind=dc), dimension(3, 3), intent(out) :: MM
    !"""
    !Multiple left by right. Assumes 3x3 (colour) (complex) matrices
    !"""
    !# do the maths for the colour matrices
    MM(1, 1) = left(1, 1) * right(1, 1) + left(1, 2) * right(2, 1) + left(1, 3) * right(3, 1)
    MM(2, 1) = left(2, 1) * right(1, 1) + left(2, 2) * right(2, 1) + left(2, 3) * right(3, 1)
    MM(3, 1) = left(3, 1) * right(1, 1) + left(3, 2) * right(2, 1) + left(3, 3) * right(3, 1)
    !# second index
    MM(1, 2) = left(1, 1) * right(1, 2) + left(1, 2) * right(2, 2) + left(1, 3) * right(3, 2)
    MM(2, 2) = left(2, 1) * right(1, 2) + left(2, 2) * right(2, 2) + left(2, 3) * right(3, 2)
    MM(3, 2) = left(3, 1) * right(1, 2) + left(3, 2) * right(2, 2) + left(3, 3) * right(3, 2)
    !# third index
    MM(1, 3) = left(1, 1) * right(1, 3) + left(1, 2) * right(2, 3) + left(1, 3) * right(3, 3)
    MM(2, 3) = left(2, 1) * right(1, 3) + left(2, 2) * right(2, 3) + left(2, 3) * right(3, 3)
    MM(3, 3) = left(3, 1) * right(1, 3) + left(3, 2) * right(2, 3) + left(3, 3) * right(3, 3)
  end subroutine MultiplyMatMat


  subroutine MultiplyMatdagMatdag(MM, left, right)
    !!use, intrinsic ::  ISO_C_BINDING, only: dc=>C_DOUBLE_COMPLEX, dp=>C_DOUBLE
    use types
    complex(kind=dc), dimension(3, 3), intent(in) :: left, right
    complex(kind=dc), dimension(3, 3), intent(out) :: MM
    !"""
    !#Multiplies two (3,3) complex matrices together. Takes conjugate
    !Does (left*right)^dagger
    !"""
    !# take transpose manually
    MM(1, 1) = conjg(left(1, 1) * right(1, 1) + left(2, 1) * right(1, 2) + left(3, 1) * right(1, 3))
    MM(2, 1) = conjg(left(1, 2) * right(1, 1) + left(2, 2) * right(1, 2) + left(3, 2) * right(1, 3))
    MM(3, 1) = conjg(left(1, 3) * right(1, 1) + left(2, 3) * right(1, 2) + left(3, 3) * right(1, 3))
    !# but take conjugate using np
    MM(1, 2) = conjg(left(1, 1) * right(2, 1) + left(2, 1) * right(2, 2) + left(3, 1) * right(2, 3))
    MM(2, 2) = conjg(left(1, 2) * right(2, 1) + left(2, 2) * right(2, 2) + left(3, 2) * right(2, 3))
    MM(3, 2) = conjg(left(1, 3) * right(2, 1) + left(2, 3) * right(2, 2) + left(3, 3) * right(2, 3))
    !# last index
    MM(1, 3) = conjg(left(1, 1) * right(3, 1) + left(2, 1) * right(3, 2) + left(3, 1) * right(3, 3))
    MM(2, 3) = conjg(left(1, 2) * right(3, 1) + left(2, 2) * right(3, 2) + left(3, 2) * right(3, 3))
    MM(3, 3) = conjg(left(1, 3) * right(3, 1) + left(2, 3) * right(3, 2) + left(3, 3) * right(3, 3))
  end subroutine MultiplyMatdagMatdag


  subroutine RealTraceMultMatMat(TrMM, left, right)
    !use, intrinsic ::  ISO_C_BINDING, only: dc=>C_DOUBLE_COMPLEX, dp=>C_DOUBLE
    use types
    complex(kind=dc), dimension(3, 3), intent(in) :: left, right
    real(kind=dc), intent(out) :: TrMM
    !"""
    !# !Takes the real trace of (3,3) complex numbers left, right multiplied together
    !Tr(left*right)
    !"""
    TrMM = real(left(1, 1) * right(1, 1) + left(1, 2) * right(2, 1) + left(1, 3) * right(3, 1) + &
                left(2, 1) * right(1, 2) + left(2, 2) * right(2, 2) + left(2, 3) * right(3, 2) + &
                left(3, 1) * right(1, 3) + left(3, 2) * right(2, 3) + left(3, 3) * right(3, 3), kind=dp)
  end subroutine RealTraceMultMatMat



  subroutine plaquette(data, muStart, muEnd, nuEnd, sumTrP, nP, time)
          use types
    !use, intrinsic ::  ISO_C_BINDING, only: dc=>C_DOUBLE_COMPLEX, dp=>C_DOUBLE
    complex(kind=dc), dimension(:, :, :, :, :, :, :), intent(in)  :: data
    integer,                                          intent(in)  :: muStart, muEnd, nuEnd
    real(kind=dp),                                    intent(out) :: sumTrP, time
    integer,                                          intent(out) :: nP
    !"""
    !Calculates the plaquette over muStart to muEnd
    !data is [nt, nx, ny, nz, mu, colour, colour] complex
    !the plaquette over all lattice is muStart=1, muEnd=4, nuEnd=4
    !the spatial plaquette is muStart=2, muEnd=4, nuEnd=4
    !the temporal plaquette is muStart=1, muEnd=1, nuEnd=4
    !returns the sum of plaquettes, number of plaquettes measured,
    !the average plaquette and the time taken to calculate it
    !"""
    integer, dimension(7) :: dataShape
    integer, dimension(4) :: muCoord, nuCoord, coordBase, coord
    integer :: mu, nu, nx, ny, nz, nt, cc  ! Counters
    ! For intermediate calculating plaquette
    complex(kind=dc), dimension(3, 3) :: Umu_x, Unu_xmu, UmuUnu
    complex(kind=dc), dimension(3, 3) :: Umu_xnu, Unu_x, UmudagUnudag
    real(kind=dp) :: P
    ! Timers
    real :: start, end
    call cpu_time(start)
    dataShape = shape(data)
    ! write(*,*) 'fort shape', shape(data)
    !# hold the sum
    sumTrP = 0.0
    !# hold the number measured
    nP = 0
    do mu = muStart, muEnd
       muCoord(:) = 0
       !# This is the shift in mu
       muCoord(mu) = 1
       do nu = mu + 1, nuEnd
       ! write(*,*) 'mu, nu', mu, nu
          nuCoord(:) = 0
          !# This is the shift in nu
          nuCoord(nu) = 1
          !# loop over all sites
          do nx = 1, dataShape(2)
             do ny = 1, dataShape(3)
                do nz = 1, dataShape(4)
                   do nt = 1, dataShape(1)
                      !# U_mu(x)
                      coordBase = (/ nt, nx, ny, nz /)
                      coord = coordBase
                      Umu_x = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                      ! write(*,*) 'Umu_x', Umu_x
                      !# U_nu(x + amu)
                      coord = coordBase + muCoord
                      !# respect periodic boundary conditions
                      do cc = 1, size(coord)
                         if (coord(cc) > dataShape(cc)) then
                            coord(cc) = 1
                         end if
                      end do
                      Unu_xmu = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                      ! write(*,*) 'Unu_xmu', Unu_xmu
                      !# U_mu(x + anu)
                      coord = coordBase + nuCoord
                      do cc = 1, size(coord)
                         if (coord(cc) > dataShape(cc)) then
                            coord(cc) = 1
                         end if
                      end do
                      Umu_xnu = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                      ! write(*,*) 'Umu_xnu', Umu_xnu
                      !# U_nu(x)
                      coord = coordBase
                      Unu_x = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                      ! write(*,*) 'Unu_x', Unu_x
                      !# Multiply bottom, right together
                      call MultiplyMatMat(UmuUnu, Umu_x, Unu_xmu)
                      ! write(*,*) 'UmuUnu', UmuUnu
                      !# Multiply left, top together, take dagger
                      call MultiplyMatdagMatdag(UmudagUnudag, Umu_xnu, Unu_x)
                      ! write(*,*) 'UmudagUnudag', UmudagUnudag
                      !# multiply two halves together, take trace
                      call RealTraceMultMatMat(P, UmuUnu, UmudagUnudag)
                      ! write(*,*) 'plaq', P
                      ! stop
                      sumTrP = sumTrP + P
                      nP = nP + 1
                   end do
                end do
             end do
          end do
       end do
    end do
    call cpu_time(end)
    time = end - start
  end subroutine plaquette


  end module plaq
