module lattice_utils
  use, intrinsic :: iso_c_binding
  use iso_fortran_env, only: real64
  use vec_utils
  use calculus_utils
  use fft_utils1
  use plot_utils
  implicit none
  
  !private :: int_search, latt_R_convert


contains

!Private procedure to this module
!Creates a matrix of integers such that m,n,p,q = -int_range to +int_range
!Ensures that all combinations of m,n,p,q exist in a matrix where
!.....1st column are values of m, 2nd col->n, 3rd col->p, 4th col->q
!5th col -> mq-np
!Hence, returns a (2*int_range + 1) row, 5 column matrix
!Used in functions local to module ---> private 
function int_search(int_range) result(list_mat)
    integer, intent(in) :: int_range
    integer :: list_mat1((2*int_range+1)**4, 5)
    integer, allocatable :: list_mat(:,:)
    integer :: m, n, p, q
    integer :: dumm, cnt, len
    logical :: mask((2*int_range+1)**4)
    !code
    m = int_range
    cnt = 1
    do m = -int_range, int_range
        do n = -int_range, int_range
            do p = -int_range, int_range
                do q = -int_range, int_range
                    dumm = m*q - n*p
                    list_mat1(cnt,1) = m
                    list_mat1(cnt, 2) = n
                    list_mat1(cnt, 3) = p
                    list_mat1(cnt, 4) = q
                    list_mat1(cnt, 5) = dumm
                    cnt = cnt+1
                end do
            end do
        end do
    end do
    mask = abs(list_mat1(:,5)) == 1
    len = count(mask)
    allocate(list_mat(len,4))
    list_mat(:,1) = pack(list_mat1(:,1), mask)
    list_mat(:,2) = pack(list_mat1(:,2), mask)
    list_mat(:,3) = pack(list_mat1(:,3), mask)
    list_mat(:,4) = pack(list_mat1(:,4), mask)
    return
end function int_search

!Private procedure to this module
!Converts provided lattices R1 and R2 to 1st quadrant R1 and R2 
!Returns data in a matrix M
!R1 is the 1st column, R2 is the 2nd column
function latt_R_convert(R1, R2, search_rad) result(col_R1_R2)
        !Inputs
        real(real64), intent(in) :: R1(2), R2(2)
        integer, intent(in) :: search_rad
        !Outputs
        real(real64) :: col_R1_R2(2,2)
        !Internals
        real(real64) :: latt_mat(2,2), det_latt
        real(real64), allocatable :: M1(:,:) !M2(:,:) 
        real(real64) :: S1(2), S2(2)
        integer :: m, n, p, q
        !logical, allocatable :: mask(:)
        integer :: i !mask_true
        !Error check
        latt_mat(:,1) = R1
        latt_mat(:,2) = R2
        det_latt = latt_mat(1,1)*latt_mat(2,2) - latt_mat(1,2)*latt_mat(2,1)
        if ( abs(det_latt) <= 10**(-10) ) then
            error stop "R1 and R2 are linearly dependent"
        end if
        M1 = int_search(int_range = search_rad)
        do i = 1, size(M1(:,1))
            !define m, n, p, q
            m = M1(i,1)
            n = M1(i,2)
            p = M1(i,3)
            q = M1(i,4)
            S1 = m*R1 + n*R2
            S2 = p*R1 + q*R2
            if ( S1(1) > 0 .and. S1(2) > 0 .and. S2(1) > 0 .and. S2(2) > 0 ) then
                if ( .not.(S1(1) == 0 .and. S1(2) == 0) .and. .not.(S2(2) ==0 .and. S2(1) == 0)) then
                    col_R1_R2(:,1) = S1
                    col_R1_R2(:,2) = S2
                    return 
                end if
            end if
        end do
        return
end function latt_R_convert

!Define a square lattice
function square_latt_sb(X, Y, R_latt, A, sig) result(tens)
  real(real64) :: X(:), Y(:), R_latt !Input
  real(real64) :: A, sig !Gaussian amplitude and sdev
  real(real64) :: tens(size(Y),size(X),3) !Output
  real(real64), allocatable :: tens_XY_00(:,:,:) !The same as tens with the X-Y grid in tens(:,:,:) set to X=0, Y=0, at the leftmost corner 
  real(real64) :: dmp_x, dmp_y !dummies 
  integer :: m_x, m_y
  integer :: i, j, k, cnt, dmp
  real(real64) :: xmin, xmax, ymin, ymax, gauss_dummy(size(Y),size(X)) ! gauss_dummy(size(Y),size(X),3)
  real(real64), allocatable  :: X_sh(:), Y_sh(:) !The shifted X and Y vectors
  real(real64) , allocatable :: pos(:,:)

  !error check
  if ( R_latt > maxval(X) - minval(X)  ) then
    error stop "R_latt too large for array"
  end if
  if ( R_latt > maxval(Y) - minval(Y) ) then
    error stop "R_latt too large for array"
  end if

  !allocate size of tens_XY_00 and build its X and Y and then shift the leftmost corner to (0,0)
  allocate(tens_XY_00(size(Y), size(X), 3))
  tens(:,:,2:3) = grid_2(X, Y)
  tens_XY_00(:,:,2:3) = tens(:,:,2:3)
  tens_XY_00(:,:,2) = tens_XY_00(:,:,2) - minval(X) !shift X to 0
  tens_XY_00(:,:,3) = tens_XY_00(:,:,3) - minval(Y) !shift y to 0
  tens_XY_00(:,:,1) = 0.0_real64

  !figure out the size of the loop and then the vector
  m_x = floor((maxval(X) - minval(X))/R_latt)
  m_y = floor((maxval(Y) - minval(Y))/R_latt) 
  allocate(pos((m_x+1)*(m_y+1), 2)) !Debug vector. Remove after debug

  !Run the loop
  xmin = minval(tens_XY_00(:,:,2))
  xmax = maxval(tens_XY_00(:,:,2))
  ymin = minval(tens_XY_00(:,:,3))
  ymax = maxval(tens_XY_00(:,:,3))
  X_sh = X - minval(X)
  Y_sh = Y - minval(Y)
  cnt = 1 !New indexing integer
  do k = 0, m_x
    do j = 0, m_y
      gauss_dummy = gauss_2D_nocorr_core(X=tens_XY_00(:,:,2), Y=tens_XY_00(:,:,3), A = 1.0_real64, x0 = k*R_latt, y0 = j*R_latt, sig_x = sig, sig_y = sig)
      tens_XY_00(:,:,1) = tens_XY_00(:,:,1) + gauss_dummy
    end do
  end do

  tens(:,:,1) = tens_XY_00(:,:,1)
end function square_latt_sb

!Define a rectangular lattice
function rect_latt_sb(X, Y, R_latt_x, R_latt_y, A, sig) result(tens)
  real(real64) :: X(:), Y(:), R_latt_x, R_latt_y !Input
  real(real64) :: A, sig !Gaussian amplitude and sdev
  real(real64) :: tens(size(Y),size(X),3) !Output
  real(real64), allocatable :: tens_XY_00(:,:,:) !The same as tens with the X-Y grid in tens(:,:,:) set to X=0, Y=0, at the leftmost corner
  integer :: m_x, m_y
  integer :: i, j, k, cnt, dmp
  real(real64) :: xmin, xmax, ymin, ymax, gauss_dummy(size(Y),size(X)) ! gauss_dummy(size(Y),size(X),3)
  real(real64), allocatable  :: X_sh(:), Y_sh(:) !The shifted X and Y vectors
  real(real64) , allocatable :: pos(:,:)

  !error check
  if ( R_latt_x > maxval(X) - minval(X)  ) then
    error stop "R_latt too large for array"
  end if
  if ( R_latt_y > maxval(Y) - minval(Y) ) then
    error stop "R_latt too large for array"
  end if

  !allocate size of tens_XY_00 and build its X and Y and then shift the leftmost corner to (0,0)
  allocate(tens_XY_00(size(Y), size(X), 3))
  tens(:,:,2:3) = grid_2(X, Y)
  tens_XY_00(:,:,2:3) = tens(:,:,2:3)
  tens_XY_00(:,:,2) = tens_XY_00(:,:,2) - minval(X) !shift X to 0
  tens_XY_00(:,:,3) = tens_XY_00(:,:,3) - minval(Y) !shift y to 0
  tens_XY_00(:,:,1) = 0.0_real64

  !figure out the size of the loop and then the vector
  m_x = floor((maxval(X) - minval(X))/R_latt_x)
  m_y = floor((maxval(Y) - minval(Y))/R_latt_y)
  allocate(pos((m_x+1)*(m_y+1), 2)) !Debug vector. Remove after debug

  !Run the loop
  X_sh = X - minval(X)
  Y_sh = Y - minval(Y)
  cnt = 1 !New indexing integer
  do k = 0, m_x
    do j = 0, m_y
      gauss_dummy = gauss_2D_nocorr_core(X=tens_XY_00(:,:,2), Y=tens_XY_00(:,:,3), A = 1.0_real64, x0 = k*R_latt_x, y0 = j*R_latt_y, sig_x = sig, sig_y = sig)
      tens_XY_00(:,:,1) = tens_XY_00(:,:,1) + gauss_dummy
    end do
  end do

  tens(:,:,1) = tens_XY_00(:,:,1)
end function rect_latt_sb

!Define an arbitrary lattice

!Define a hexagonal lattice

!Define a honeycomb lattice


  
end module lattice_utils