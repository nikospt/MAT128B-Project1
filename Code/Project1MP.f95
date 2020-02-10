program Project

    integer :: xmax, ymax
    parameter (xmax = 500)
    parameter (ymax = 500)
    real :: start_time
    real :: end_time
    real :: xl, yl, xu, yu
    character (len = 10) :: fname1
    character (len = 350) :: sysin
    character (len = 350) :: sysin1, sysin2, sysin3
    character (len = 4) :: test
    fname1 = "data1.dat"

    ! The 3 is the number of digits of result and may need to be changed if you increase ymax sufficiently
    write(test, '(I3)') ymax/2

    sysin = 'gnuplot -e "reset; set terminal png; set output \"test.png\"; set xrange [-1:1]; set yrange [-1:1]; &
            set key off; unset colorbox; set size ratio -1; &
            set palette maxcolors 2; set palette defined (1 \"#FF0000\", 2 \"#FFFFFF\"); &
            plot \"data1.dat\" using ((\$1-'//trim(test)//')/'//trim(test)//'):((\$2-'//trim(test)// &
            ')/'//trim(test)//'):3 matrix with image"'
    write(*,*) sysin

    call cpu_time(start_time)
    xl = -1.0
    xu = 1.0
    yl = -1.0
    yu = 1.0
    call iterate(xl,xu,yl,yu,phi,fname1)

    call cpu_time(end_time)
    write(*,*) "Process took ", end_time-start_time, "seconds"
    call system(sysin)

    ! Probabley can get rid of this
    ! sysin = 'gnuplot -e "reset; set terminal png; set output \"test.png\"; set xrange [-1:1]; set yrange [-1:1]; &
    !         set key off; unset colorbox; set size ratio -1; &
    !         set palette maxcolors 2; set palette defined (1 \"#FF0000\", 2 \"#FFFFFF\"); &
    !         plot \"data1.dat\" using ((\$1-250)/250):((\$2-250)/250):3 matrix with image"'

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

    subroutine iterate(xl,xu,yl,yu,phi,filename)
        character(len=10) :: filename
        real, intent(in) :: xl, xu, yl, yu
        real :: x(xmax), y(ymax)
        real :: M(xmax,ymax)
        real :: recc(100,500)
        complex :: phi
        complex :: a
        integer :: i, j, k

        call linspace(xl,xu,x)
        call linspace(yl,yu,y)
        !$OMP DO
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
        !$OMP END DO

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
