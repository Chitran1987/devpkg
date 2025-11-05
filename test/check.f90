program check
use iso_fortran_env, only: real64
use vec_utils
use stat_utils
implicit none
real(real64), allocatable :: X(:), Y1(:), Y2(:), Y(:)
real(real64), allocatable :: M(:,:), res(:,:)
real(real64) :: wind(1,2)
!build the dataset
X = seq(0.0_real64, 10.0_real64, 1000)
allocate(Y(size(X)))
allocate(Y1(size(X)))
Y = 2
Y1 = 0
Y = merge(Y, Y1, X <= 2)
print *, 'The original dataset check is ', sum(Y)/size(Y)
!add background
Y = Y + 3*X +1
print *, 'The background added dataset check is', sum(Y)/size(Y) 

!delete background
allocate(M(size(X),2))
M(:,1) = X
M(:, 2) = Y
print *, 'bla bla'
wind(1,1) = 4.0_real64
wind(1,2) = 10.0_real64
print *, 'bla bla'
print *, X(1)
res = lin_bg_sub_1D(M, wind)
print *, 'The background substracted dataset check is', sum(res(:,2))/size(res(:,2))
end program check
