      subroutine toeplitz_cholesky_lower ( n, a, l )

c*********************************************************************72
c
cc TOEPLITZ_CHOLESKY_LOWER: lower Cholesky factor of a Toeplitz matrix.
c
c  Discussion:
c
c    The Toeplitz matrix must be positive semi-definite.
c
c    After factorization, A = L * L'.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Michael Stewart,
c    Cholesky factorization of semi-definite Toeplitz matrices.
c    Linear Algebra and its Applications,
c    Volume 254, pages 497-525, 1997.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the Toeplitz matrix.
c
c    Output, double precision L(N,N), the lower Cholesky factor.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision div
      double precision g(2,n)
      double precision g1j
      double precision g2j
      integer i
      integer j
      double precision l(n,n)
      double precision rho

      do j = 1, n
        do i = 1, n
          l(i,j) = 0.0D+00
        end do
      end do

      do j = 1, n
        g(1,j) = a(1,j)
      end do
      g(2,1) = 0.0D+00
      do j = 2, n
        g(2,j) = a(j,1)
      end do 

      do i = 1, n
        l(i,1) = g(1,i)
      end do
      do j = n, 2, -1
        g(1,j) = g(1,j-1)
      end do
      g(1,1) = 0.0D+00

      do i = 2, n

        rho = - g(2,i) / g(1,i)
        div = sqrt ( ( 1.0D+00 - rho ) * ( 1.0D+00 + rho ) )

        do j = i, n
          g1j = g(1,j)
          g2j = g(2,j)
          g(1,j) = (       g1j + rho * g2j ) / div
          g(2,j) = ( rho * g1j +       g2j ) / div
        end do

        do j = i, n
          l(j,i) = g(1,j)
        end do
        do j = n, i + 1, -1
          g(1,j) = g(1,j-1)
        end do
        g(1,i) = 0.0D+00

      end do

      return
      end
  