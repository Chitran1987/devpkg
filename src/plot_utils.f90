module plot_utils
    use iso_fortran_env, only: real64
    implicit none
    
contains
    subroutine plot_local(X, Y)
        implicit none
        real(real64), dimension(:) :: X, Y
        integer :: i, tmp
        open(newunit = tmp, file = 'data.txt', status = 'replace')
        do i = 1, size(X)
            write(tmp, *) X(i), Y(i)
        end do
        close(tmp)
        call system('graph -T X data.txt')
    end subroutine plot_local

    subroutine plot_gnu(X, Y)
        implicit none
        real(real64), dimension(:) :: X, Y
        integer :: i, tmp
        open(newunit = tmp, file = 'data.txt', status = 'replace')
        do i = 1, size(X)
            write(tmp, *) X(i), Y(i)
        end do
        close(tmp)
        call execute_command_line( &
            'gnuplot -p -e "set grid; set xlabel ''X''; set ylabel ''Y''; ' // &
            'plot ''data.txt'' using 1:2 with linespoints title ''Y vs X''"' )
    end subroutine plot_gnu

    subroutine plot_oct(X, Y)
        implicit none
        real(real64), dimension(:) :: X, Y
        integer :: i, tmp
        open(newunit = tmp, file = 'data.txt', status = 'replace')
        do i = 1, size(X)
            write(tmp, *) X(i), Y(i)
        end do
        close(tmp)
        call execute_command_line( &
        'octave-cli --quiet --persist --eval "D=dlmread(''data.txt''); ' // &
        'plot(D(:,1),D(:,2),''-o''); grid on; xlabel(''X''); ylabel(''Y''); ' // &
        'title(''Y vs X''); drawnow;"' )
    end subroutine plot_oct
    
end module plot_utils