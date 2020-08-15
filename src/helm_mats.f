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
c
c
c
c
c
c
c
c
c
      subroutine helm2d_dielec_mat(k,nch,n,srcvals,srccoefs,zpars,
     1   iregionl,iregionr,xmat)
c
c     This subroutine generates the matrix corresponding to the 
c     dielectric interface problem.
c     Unknowns are ordered by blocks corresponding to sigma/mu densities
c     one chunk after another. (All sigmas for a chunk, then all mus).
c     While there may be multiple regions involved, we are using a 
c     somewhat inefficient but robust global representation that makes 
c     the integral equation setup simple. Assuming regions are defined
c     by index i = 0,1,...,M, and u_i the restriction of the solution to 
c     domain i, we let
c
c     u_i = (-1/bi) S_ki[Gamma] sigma + (1/bi) D_ki[Gamma] mu .
c
c     wnere Gamma is the *total* boundary.
c     In domain i, the Helmholtz parameter is zk(i), and the
c     jump coefficients are bi = b(i) and ai = a(i).
c
c     The system matrix is 2n x 2n with n = nch*k.
c
c     The 11 block enforces the conditions
c     [au] = f1  ([.] meaning right - left) using the mu unknowns for 
c                all chunks.
c
c     The 12 block enforces the conditions
c     [au] = f1  ([.] meaning right - left) using the sigma unknowns for 
c                all chunks.
c
c     The 22 block enforces the conditions
c     [b dudn] = f2  ([.] meaning right - left), with the normal pointing
c                    from left region to right region, for sigma unknowns.
c
c     The 21 block enforces the conditions
c     [b dudn] = f2  ([.] meaning right - left), with the normal pointing
c                    from left region to right region, for mu unknowns.
c     
c
c     In forming ith row of the linear system for point on chunk j, let
c     iregionl(j) denote domain on the left and 
c     iregionr(j) denote domain on the right as one traverses the chunk.
c
c     Let al = a(iregionl(j)), ar = a(iregionr(j)) 
c         bl = b(iregionl(j)), br = b(iregionr(j)).
c         kl = zk(iregionl(j)), kr = zk(iregionr(j)).
c
c     [au] = (1/2)[ar/br + al/bl]mu + [(ar/br) D^*_kr - (al/bl)D^*_kl] mu 
c             - [ ar/br S_kr -  al/bl S_kl] sigma = f1
c     
c     In forming row of the linear system for point on chunk j,
c
c     [b dudn] = [D'_kr - D'_kl] mu  + sigma - [S'_kr - S'_kl] = f2
c
c     or (with q = (1/2)[ar/br+al/bl] )
c
c     [au]/q = mu + [(ar/br) D^*_kr - (al/bl)D^*_kl]/q mu 
c             - [ ar/br S_kr -  al/bl S_kl]/q sigma = f1/q
c     [b dudn] = [D'_kr - D'_kl] mu  + sigma - [S'_kr - S'_kl] = f2
c
c--------------
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
c      zpars(1): zkr
c      zpars(2): ar
c      zpars(3): br
c      zpars(4): zkl
c      zpars(5): al
c      zpars(6): bl
c    iregionl: integer (nch)  - left domain for chunk
c    iregionr: integer (nch)  - right domain for chunk
c
c  Output arguments:
c    xmat: complex *16 (2*n,2*n)
c      matrix corresponding to discretization
c  
      implicit none
      integer k,nch,n,i,j
      real *8 srcvals(8,k,nch),srccoefs(6,k,nch)
      complex *16 zpars(6),zsend(6),zq,ar
      complex *16 xmat(2*n,2*n)
      real *8 dpars
      integer ipars
      integer iregionl(nch)
      integer iregionr(nch)
      complex *16, allocatable :: xmatij(:,:)
      external h2d_transmission_dir
      external h2d_transmission_neu
c
      allocate(xmatij(n,n))
      zsend(1) = zpars(1)
      zsend(2) = zpars(4)
      zq = 0.5*(zpars(2)/zpars(3)+zpars(5)/zpars(6))
c
c     11 block  (mu to [au])
c
      call prinf(' doing 11 block *',nch,0)
      call prin2(' zpars are *',zpars,2*6)
      zsend(3) = 0.0d0
      zsend(4) = 0.0d0
      zsend(5) = (zpars(2)/zpars(3))/zq
      zsend(6) = -(zpars(5)/zpars(6))/zq
      call prin2(' zsend are *',zsend,2*6)
      call zgetmat_bdry(k,nch,n,srcvals,srccoefs,h2d_transmission_dir,
     1  8,0,dpars,6,zsend,0,ipars,xmatij)
c
      call prin2(' xmatij = *',xmatij,2*n*n)
      do i = 1,n
      do j = 1,n
         xmat(i,j) = xmatij(i,j)
      enddo
      enddo
c
c     12 block  (sigma to [au])
c
      call prinf(' doing 12 block *',nch,0)
      zsend(3) = -(zpars(2)/zpars(3))/zq
      zsend(4) = (zpars(5)/zpars(6))/zq
      zsend(5) = 0.0d0
      zsend(6) = 0.0d0
      call zgetmat_bdry(k,nch,n,srcvals,srccoefs,h2d_transmission_dir,
     1  8,0,dpars,6,zsend,0,ipars,xmatij)
c
      call prin2(' xmatij = *',xmatij,2*n*n)
      do i = 1,n
      do j = 1,n
         xmat(i,j+n) = xmatij(i,j)
      enddo
      enddo
c
c     22 block  (sigma to [b dudn])
c
      call prinf(' doing 22 block *',nch,0)
      zsend(3) = -1.0d0
      zsend(4) = 1.0d0
      zsend(5) = 0.0d0
      zsend(6) = 0.0d0
      call zgetmat_bdry(k,nch,n,srcvals,srccoefs,h2d_transmission_neu,
     1  8,0,dpars,6,zsend,0,ipars,xmatij)
c
      call prin2(' xmatij = *',xmatij,2*n*n)
      do i = 1,n
      do j = 1,n
         xmat(i+n,j+n) = xmatij(i,j)
      enddo
      enddo
c
c     21 block  (mu to [b dudn])
c
      call prinf(' doing 21 block *',nch,0)
      zsend(3) = 0.0d0
      zsend(4) = 0.0d0
      zsend(5) = 1.0d0
      zsend(6) = -1.0d0
      call zgetmat_bdry(k,nch,n,srcvals,srccoefs,h2d_transmission_neu,
     1  8,0,dpars,6,zsend,0,ipars,xmatij)
c
      call prin2(' xmatij = *',xmatij,2*n*n)
      do i = 1,n
      do j = 1,n
         xmat(i+n,j) = xmatij(i,j)
      enddo
      enddo
c
      return
      end
