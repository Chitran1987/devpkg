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
        integer, allocatable :: M1(:,:) !M2(:,:)
        real(real64) :: S1(2), S2(2)
        integer :: m, n, p, q
        !logical, allocatable :: mask(:)
        integer :: i !mask_true
        !Error check
        latt_mat(:,1) = R1
        latt_mat(:,2) = R2
        det_latt = latt_mat(1,1)*latt_mat(2,2) - latt_mat(1,2)*latt_mat(2,1)
        if ( abs(det_latt) <= 1.0e-10_real64 ) then
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
                col_R1_R2(:,1) = S1
                col_R1_R2(:,2) = S2
                return
            end if
        end do
        error stop "No first-quadrant primitive basis found. Increase search_rad or provide a manually reoriented unit cell."
end function latt_R_convert

!create a latt_pnt_lst() function
!The function should return a list of lattice points of m*R1+n*R2 type in the X-Y space
!shft = 0 means that once the lattice points are created
    !!There shouldat least be one lattice point present at the bottom left corner (X_min, Y_min)
    !!shft /= 0 means one latt point at (X_min+shft(1), Y_min+shft(2))
function latt_pnt_lst(X, Y, R1, R2, shft) result(lst_mat)
    real(real64) :: R1(2), R2(2), shft(2) !Inputs
    real(real64), intent(in) :: X(:), Y(:)
    real(real64), allocatable :: lst_mat(:,:) !Outputs
    real(real64) :: R1R2_mat(2,2), X1(size(X)), Y1(size(Y)) !Internals
    real(real64) :: R1_int(2), R2_int(2), lenX, lenY, t1, t2, p(2) !Internals
    real(real64) :: V_low(2), V_up(2) !Internals
    real(real64) :: distX, distY, dist_R1, dist_R2, vec(2), c1, c2, min_vec(2)  !Internals
    real(real64), allocatable :: M1(:,:), Mx(:,:), My(:,:)  !Internals
    real(real64), allocatable :: dis_small(:) !Internals
    integer :: r, s, cnt, i, j, k, shft_idx !Internals
    logical, allocatable :: mask_X(:), mask_Y(:)  !Internals

    !code----------------------------------------------------------------
    !!!swap R1 and R2 if necessary----------
    X1 = X - minval(X)
    Y1 = Y - minval(Y)
    !Get R1 and R2 into first quadrants
    R1R2_mat = latt_R_convert(R1, R2, search_rad=5)
    R1_int = R1R2_mat(:,1)
    R2_int = R1R2_mat(:,2)
    !swap R1_int and R2_int if necessary
    if ( atan(R1_int(2)/R1_int(1)) > atan(R2_int(2)/R2_int(1)) ) then
        R1R2_mat(:,1) = R2_int
        R1R2_mat(:,2) = R1_int
    end if
    R1 = R1R2_mat(:,1)
    R2 = R1R2_mat(:,2)
    !!!Calculate the skewed shape--------------
    lenX = maxval(X1) - minval(X1)
    lenY = maxval(Y1) - minval(Y1)
    t1 = R1(2)/R1(1)
    t2 = R2(2)/R2(1)
    c1 = -(1.5_real64*lenX*t1 + 0.5_real64*lenY)
    c2 =  0.5_real64*lenX*t2 + 1.5_real64*lenY
    p(1) = -(c2 - c1)/(t2 - t1)
    p(2) = (t2*c1 - t1*c2)/(t2 - t1)
    V_low = [1.5_real64*lenX, -0.5_real64*lenY]
    V_up = [-0.5_real64*lenX, 1.5_real64*lenY]
    distX = sqrt(sum((p - V_low)*(p - V_low)))
    distY = sqrt(sum((p - V_up)*(p - V_up)))
    dist_R1 = sqrt( sum(R1*R1) )
    dist_R2 = sqrt( sum(R2*R2) )
    r = ceiling(distX/dist_R1)
    s = ceiling(distY/dist_R2)
    allocate(M1((r+1)*(s+1),2))
    vec = p
    cnt = 1
    !!! Populate the skewed shape-------------------------
    do i = 0, r
        do j = 0, s
            vec = p + i*R1 + j*R2
            M1(cnt,:) = vec
            cnt = cnt + 1
        end do
    end do
    !!!Snip the skewed shape to a large rectangle -----------------------------
    mask_X = ( M1(:,1)>=-0.5_real64*lenX ) .and. ( M1(:,1) <= 1.5_real64*lenX )
    allocate(Mx(count(mask_X),2))
    Mx(:,1) = pack(M1(:,1), mask_X)
    Mx(:,2) = pack(M1(:,2), mask_X)
    !Now Mx is succesfully populated
    !redefine the mask
    mask_Y = (Mx(:,2)>=-0.5_real64*lenY) .and. (Mx(:,2) <= 1.5_real64*lenY)
    allocate(My(count(mask_Y),2))
    My(:,1) = pack(Mx(:,1), mask_Y)
    My(:,2) = pack(Mx(:,2), mask_Y)
    !Now that My is the final list first stage before shifting
    lst_mat = My
    !!!Shift the list of point accordingly 
    !calculate the distance of each point 
    !from the lower left corner of the image
    !Then take the pt corresponding to the least distance
    !shift the entire image by the shift of the corresponding point 
    !to the co-ordinates of the lower-left image corner
    dis_small = sqrt(lst_mat(:,1)**2 + lst_mat(:,2)**2)
    shft_idx = minloc(dis_small, dim=1)
    min_vec = lst_mat(shft_idx,:)
    do k = 1, size(lst_mat(:,1))
        lst_mat(k,:) = lst_mat(k,:) - min_vec
    end do
    !Then shift the entire bunch by shft
    !Then shift it by the actual upper and lower corners of X and Y
    do k = 1, size(lst_mat(:,1))
        lst_mat(k,:) = lst_mat(k,:) + shft + [minval(X), minval(Y)] !shift by shft
    end do
    !Then shift it by 
    !Then return
    return
end function latt_pnt_lst

!Populate a lattice given The X, Y and lattice point list
!Returns a 3 dimensional tensor for it
function populate_latt(X, Y, A, sig, list) result(tens)
    real(real64), intent(in) :: X(:), Y(:)
    real(real64), intent(in) :: list(:,:)
    real(real64), intent(in) :: A, sig

    real(real64) :: tens(size(Y), size(X), 3)

    real(real64) :: gauss_dumm(size(Y), size(X))
    real(real64) :: pop_sum(size(Y), size(X))
    real(real64) :: Y_local(size(Y))

    integer :: i

    ! grid_2 reverses its Y argument in place,
    ! so give it a local copy.
    Y_local = Y

    ! Generate coordinate grids before entering OpenMP.
    tens(:,:,2:3) = grid_2(X, Y_local)

    pop_sum = 0.0_real64

    !$omp parallel do default(none) schedule(static)             &
    !$omp& shared(tens, list, A, sig)                            &
    !$omp& private(gauss_dumm) reduction(+:pop_sum)
    do i = 1, size(list, 1)
        gauss_dumm = gauss_2D_nocorr_core(X=tens(:,:,2), Y=tens(:,:,3), A=A, x0=list(i,1), y0=list(i,2), sig_x=sig, sig_y=sig)
        pop_sum = pop_sum + gauss_dumm
    end do
    !$omp end parallel do
    tens(:,:,1) = pop_sum
    return
end function populate_latt


!Define a square lattice
!single atomic basis(sb)
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
!single basis(sb)
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
function latt_arb(X, Y, R1, R2, A, sig, shft) result(tens)
    real(real64) :: X(:), Y(:), R1(2), R2(2), A, sig, shft(2) !Inputs
    real(real64) :: tens(size(Y), size(X), 3) !Outputs
    !Internal declarations
    real(real64), allocatable :: lst(:,:)
    lst = latt_pnt_lst(X=X, Y=Y, R1=R1, R2=R2, shft=shft)
    tens = populate_latt(X=X, Y=Y, A=A, sig=sig, list=lst)
    return
end function latt_arb


!Define a hexagonal lattice
!single atomic basis(sb)
function hex_latt_sb(X, Y, R_latt, A, sig) result(tens)
    real(real64) :: X(:), Y(:), R_latt, A, sig !Inputs
    real(real64) :: tens(size(Y), size(X), 3) !Outputs
    real(real64) :: R1(2), R2(2)
    !code
    R1 = [R_latt, 0.0_real64]
    R2 = [R_latt*cos(pi/3), R_latt*sin(pi/3)]
    tens = latt_arb(X=X, Y=Y, R1=R1, R2=R2, A=A, sig=sig, shft=[0.0_real64, 0.0_real64])
    return
end function hex_latt_sb

!Define a honeycomb lattice


  
end module lattice_utils