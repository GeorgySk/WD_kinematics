module math

    implicit none
    private

    public :: PI, &
              multiplyMatrixByVector

    real*8, parameter :: PI = 4.d0 * datan(1.0d0)

contains

    function multiplyMatrixByVector(A,x) result(y)
        real*8,dimension(3,3),intent(in) :: A
        real*8,dimension(3),intent(in) :: x
        real*8,dimension(3) :: y
        y(1) = A(1, 1)*x(1) + A(1, 2)*x(2) + A(1, 3)*x(3)
        y(2) = A(2, 1)*x(1) + A(2, 2)*x(2) + A(2, 3)*x(3)
        y(3) = A(3, 1)*x(1) + A(3, 2)*x(2) + A(3, 3)*x(3)
    end function

end module