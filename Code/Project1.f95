program Project
    implicit none
    integer :: xmax, ymax
    parameter (xmax = 500)
    parameter (ymax = 500)
    real :: x(xmax)
    real :: y(ymax)
    real :: x2(xmax), y2(ymax)
    real :: M(xmax,ymax)
    real :: T(xmax,ymax)
    integer :: i,j,k
    complex :: a
    character (len = 10) :: fname1,fname2

    fname1 = "data1.dat"
    fname2 = "data2.dat"

    call linspace(-1.0,1.0,x)
    call linspace(-1.0,1.0,y)
    call linspace(-1.8,1.8,x2)
    call linspace(-0.7,0.7,y2)

    call iterate(x,y,phi,fname1,M)
    call iterate(x2,y2,phi2,fname2,T)

contains

    subroutine linspace(from, to, array)
        implicit none
        real, intent(in) :: from, to
        real, intent(out) :: array(:)
        real :: range
        integer :: n, i
        n = size(array)
        range = to - from

        if (n == 0) return

        if (n == 1) then
            array(1) = from
            return
        end if

        do i=1, n
            array(i) = from + range * (i - 1) / (n - 1)
        end do
    end subroutine linspace

    subroutine iterate(x,y,phi,filename,M)
        real, intent(in) :: x(:), y(:)
        complex :: phi
        character(len=10) :: filename
        real, intent(out) :: M(:,:)
        complex :: a
        do i = 1, xmax
            do j = 1, ymax
                M(i,j) = 1
                a = cmplx(x(i),y(j))
                do k = 1, 100
                    a = phi(a)
                    if (abs(a) > 2) then
                        M(i,j) = 2
                    end if
                end do
            end do
        end do
        open (unit = 8, file=filename)
        do i = 1,xmax
            write(8,*) (M(i,j),j=1,ymax)
        end do
        close(8)
    end subroutine iterate

    complex function phi(z)
        complex, intent(in) :: z
        phi = z*z
    end function phi

    complex function phi2(z)
        complex, intent(in) :: z
        phi2 = z**2 - 1.25
    end function phi2

end program Project
