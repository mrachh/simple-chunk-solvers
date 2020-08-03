      subroutine helm2d_comb_dir_mat(k,nch,n,srcvals,srccoefs,zpars,
     1   xmat)
c
c--------------
c  This subroutine generates the matrix corresponding to the 
c  Dirichlet data for the combined field representation 
c  
c  Rep:
c   u = \alpha S_{k}  + \beta D_{k}
c
c  Boundary data:
c   u  = f
c
c  Note: the matrix only includes the principal value part
c   of the discretization and not the identity term
c
c  Input arguments:
c    k: integer
c      order of discretization
c    nch: integer
c      number of chunks
c    n: integer
c      number of points
c    srcvals: real *8 (8,k,nch)
c      source info at discretization nodes
c        srcvals(1:2,j,ich) - xy data 
c        srcvals(3:4,j,ich) - dxdt,dydt data
c        srcvals(5:6,j,ich) - d2xdt2,d2ydt2 data
c        srcvals(7:8,j,ich) - normal data
c    srccoefs: real *8 (6,k,nch)
c      legendre expansion coeffs of the chunks
c        srccoefs(1:2,j,ich) - xy coefs 
c        srccoefs(3:4,j,ich) - dxdt,dydt coefs
c        srccoefs(5:6,j,ich) - d2xdt2,d2ydt2 coefs
c    zpars: complex *16(3)
c      parameters for problem prescription
c      zpars(1): k
c      zpars(2): alpha
c      zpars(3): beta
c
c  Output arguments:
c    xmat: complex *16 (n,n)
c      matrix corresponding to discretization
c  
      implicit none
      integer k,nch,n
      real *8 srcvals(8,k,nch),srccoefs(6,k,nch)
      complex *16 zpars(3),xmat(n,n)
      real *8 dpars
      integer ipars
      external h2d_comb

      call zgetmat_bdry(k,nch,n,srcvals,srccoefs,h2d_comb,8,0,
     1  dpars,3,zpars,0,ipars,xmat)

      return
      end
