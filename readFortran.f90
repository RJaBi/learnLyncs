
!! Reads a single gaugefield in. Calculates whole, spatial, temporal plaquette

program readFortran
  implicit none
! Kinds
  integer, parameter :: dc = kind((1.0D0,1.0D0)) !! Double precision complex scalars.
  integer, parameter :: dp = kind(1.0D0) !! Double precision real
  ! For reading
  integer, parameter :: NS=24, NT=8, NC=3
  integer, parameter :: unit=303
  character(len=128) :: fileName
  ! What holds the gaugefield
  complex(kind=dc), dimension(NT, NS, NS, NS, 4, NC, NC) :: U
  ! For plaquette
  real(kind=dp) :: sumTrP, time
  integer :: nP
  ! Read this file in
  fileName = 'Gen2_8x24_gfAr0.fort'
  open(unit, file=fileName, form='unformatted', action='read', access='stream')
  read(unit) U
  close(unit)
  write(*,*) 'type, sumTrp, nP, ave, time (seconds)'
  ! Whole plaquette
  call plaquette(U, 1, 4, 4, sumTrP, nP, time)
  write(*,*) 'whole', sumTrP, nP, sumTrP / real(nP, kind=dp), time
  ! Spatial Plaquette
  call plaquette(U, 2, 4, 4, sumTrP, nP, time)
  write(*,*) 'spatial', sumTrP, nP, sumTrP / real(nP, kind=dp), time
  ! Temporal Plaquette
  call plaquette(U, 1, 1, 4, sumTrP, nP, time)
  write(*,*) 'temporal', sumTrP, nP, sumTrP / real(nP, kind=dp), time

contains

  subroutine MultiplyMatMat(MM, left, right)
    complex(kind=dc), dimension(:, :), intent(in) :: left, right
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
    complex(kind=dc), dimension(:, :), intent(in) :: left, right
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
    complex(kind=dc), dimension(:, :), intent(in) :: left, right
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
    implicit none
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
    
    !# hold the sum
    sumTrP = 0.0
    !# hold the number measured
    nP = 0
    do mu = muStart, muEnd
       muCoord(:) = 0
       !# This is the shift in mu
       muCoord(mu) = 1
       do nu = mu + 1, nuEnd
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
                      !# U_nu(x + amu)
                      coord = coordBase + muCoord
                      !# respect periodic boundary conditions
                      do cc = 1, size(coord)
                         if (coord(cc) > dataShape(cc)) then
                            coord(cc) = 1
                         end if
                      end do
                      Unu_xmu = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                      !# U_mu(x + anu)
                      coord = coordBase + nuCoord
                      do cc = 1, size(coord)
                         if (coord(cc) > dataShape(cc)) then
                            coord(cc) = 1
                         end if
                      end do
                      Umu_xnu = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                      !# U_nu(x)
                      coord = coordBase
                      Unu_x = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                      !# Multiply bottom, right together
                      call MultiplyMatMat(UmuUnu, Umu_x, Unu_xmu)
                      !# Multiply left, top together, take dagger
                      call MultiplyMatdagMatdag(UmudagUnudag, Umu_xnu, Unu_x)
                      !# multiply two halves together, take trace
                      call RealTraceMultMatMat(P, UmuUnu, UmudagUnudag)
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


  
  
end program readFortran
