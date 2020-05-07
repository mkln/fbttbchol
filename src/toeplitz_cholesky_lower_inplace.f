      subroutine toeplitz_cholesky_lower_i (a, lda, bsize)

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
c    19 Feb 2020
c
c  Author:
c
c    John Burkardt
c    Inplace & upper block mod by mkln
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
c    Input, integer bsize, the order of the matrix.
c
c    Input/Output, double precision A(LDA, bsize), the block Toeplitz matrix's first column block.
c                  then, the lower Cholesky factor is on the lower triangular part
c                  the upper triangular part of A will be left untouched (like DPOTRF ?)
c
c    Input, integer LDA, number of rows of A (LDA = bsize * num_blocks)
c
c    Output integer IERR -- not used, set 0
c

      implicit none

      integer bsize
      integer lda

      double precision a(lda,bsize)
      double precision div
      double precision g(2,bsize)
      double precision g1j
      double precision g2j
      integer i
      integer j
C     double precision l(bsize,bsize)
      double precision rho

      do j = 1, bsize
        g(1,j) = a(1,j)
      end do
      g(2,1) = 0.0D+00
      do j = 2, bsize
        g(2,j) = a(j,1)
      end do 

      do i = 1, bsize
        a(i,1) = g(1,i)
      end do
      
      do j = bsize, 2, -1
        g(1,j) = g(1,j-1)
      end do
      g(1,1) = 0.0D+00

      do i = 2, bsize
        rho = - g(2,i) / g(1,i)
        div = sqrt ( ( 1.0D+00 - rho ) * ( 1.0D+00 + rho ) )

        do j = i, bsize
          g1j = g(1,j)
          g2j = g(2,j)
          g(1,j) = (       g1j + rho * g2j ) / div
          g(2,j) = ( rho * g1j +       g2j ) / div
        end do

C       do j = 1, i-1
C         a(j,i) = 0.0D+00
C       end do

        do j = i, bsize
          a(j,i) = g(1,j)
        end do
        
        do j = bsize, i + 1, -1
          g(1,j) = g(1,j-1)
        end do
        
        g(1,i) = 0.0D+00
        
      end do
      
      return
      end
  