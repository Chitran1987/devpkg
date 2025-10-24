program build_save_matrix
    use vec_utils
    implicit none
    
    !------------The matrix build declarations-----------!
    real(real64) :: M(3, 2)  !The matrix
    real(real64) :: del !Assign periodicity to the matrix
    integer :: i, j !Matrix indices

    !------------The datawriting build declarations-----------!
    integer :: tmp, r

    !---------------Fill the matrix with data-----------------!
    del = 2*asin(1.0_real64)/100_real64
    do j = 1, size(M, 2)
        do i = 1, size(M, 1)
            M(i,j) = cos(3*i*del)*sin(3*j*del)
        end do
    end do

    r = print_mat(M)
    !---------------Fill the matrix with data-----------------!

    !---------------Save the matrix data to data.txt-----------------!
    open(newunit=tmp, file='matrix.txt', status='replace', action='write')
     ! Loop over rows
    do i = 1, size(M,1)
        ! Loop over columns
        do j = 1, size(M,2)
            if (j < size(M,2)) then
                write(tmp,'(G0,1X)', advance='no') M(i,j)
            else
                write(tmp,'(G0)') M(i,j)
            end if
        end do
    end do
    close(tmp)

    !--------------Call the matrix data for plotting-------------------!
     

end program build_save_matrix