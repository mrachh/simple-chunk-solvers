c
c
c
c
c
c**********************************************************************
      subroutine projectonpoly(xs,ys,npts,sval,xcoeffs,ycoeffs,w1,w2)
c**********************************************************************
      implicit none
      integer *4 i,npts,norder
      real *8 xs(npts),ys(npts)
      real *8 sval(npts),xcoeffs(npts),ycoeffs(npts)
      real *8 w1(npts),w2(npts,npts)
c:::  local variables
      real *8 deltax, deltay, dlength
c
c     This subroutine takes a user-specified sequence of points
c     (xs(i),ys(i)), i = 1,npts  and computes a smooth curve
c     of order norder which interpolates them.
c
c     INPUT:
c
c     xs(i), ys(i) : real *8      ith user point   
c     npts         : integer *4   number of points
c     w1           : real *8      workspace of dimension npts
c     w2           : real *8      workspace of dimension npts*npts
c
c     OUTPUT:
c
c     sval(i)      : real *8   parameter value of projection of
c                              ith point onto segment
c                              [xs(1),ys(1)] -> [xs(npts),ys(npts)]
c     xcoeffs(j)   : real *8   array of expansion coefficients
c                              for x coordinates of curve
c        x(t) = sum_{i=1}^{npts} xcoeffs(i)*T_{i-1}(t).
c
c     ycoeffs(j)   : real *8   array of expansion coefficients
c                              for y coordinates of curve
c        y(t) = sum_{i=1}^{npts} ycoeffs(i)*T_{i-1}(t).
c
c     Version 1: Curve is the interpolating polynomial of order npts-1.
c
c---------------------------------------------------------------
      deltax = xs(npts) - xs(1)
      deltay = ys(npts) - ys(1)
      dlength = dsqrt(deltax*deltax+deltay*deltay)
c
c     compute projection of all given smooth curve pts onto straight 
c     segment from pt 1 to pt npts. sval(i) is distance along that 
c     segment of projection. If sval is not monotonic, then the curve 
c     cannot be viewed as a function of length along segment from pt 1 
c     to pt npts.
c
c     In that event, we will stop at the point where monotonicity fails
c     and start another segment.
c
c
      sval(1) = 0.0d0 
      do i = 2,npts-1
         sval(i) = (xs(i)-xs(1))*deltax + (ys(i)-ys(1))*deltay
         sval(i) = sval(i)/dlength
      enddo
      sval(npts) = dlength
ccc      call prin2( 'sval array is *',sval,npts)
c
      call getpinterp(npts,xs,sval,xcoeffs,w1,w2)
      call getpinterp(npts,ys,sval,ycoeffs,w1,w2)
ccc      call prin2( 'xcoeffs array is *',xcoeffs,npts)
ccc      call prin2( 'ycoeffs array is *',ycoeffs,npts)
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine getpinterp(npts,xs,sval,coeffs,xstemp,awork)
c**********************************************************************
      implicit none
c
c     the following subroutine takes npts values (sval(i),xs(i)) and
c     returns a simple interpolating polynomial defined by 
c     the coefficient array coeffs:
c
c        xs(i) = sum_{j=1}^{npts}  coeffs(j) * T_{j-1}(sval(i))
c
c     where T_j is the Chebyshev polynomial scaled to the interval
c     [sval(1),sval(npts)].
c
c     INPUT PARAMETERS:
c
c     npts       (integer *4):  number of points
c     xs(npts)      (real *8):  xs values to be interpolated
c     sval(npts)    (real *8):  s coordinates of pts to be interpolated
c                               on interval [sval(1),sval(npts)].
c     w1(npts)      (real *8):  workspace of dimension npts
c     w2(npts,npts) (real *8):  workspace of dimension npts*npts
c
c     OUTPUT PARAMETERS:
c
c     coeffs(npts) (real *8):  coefficients of interpolating polynomial
c
c----------------------------------------------------------------------
c
      integer  i, j, npts
c
      real *8  xs(npts), sval(npts), coeffs(npts), awork(npts,npts)
      real *8  dlength, xstemp(npts), rcond
c
      dlength = sval(npts) - sval(1)
c
c     write xs into work array xstemp
c
      do i = 1,npts
         xstemp(i) = xs(i)
      enddo
      call mkinterpmat(awork,npts,sval)
      call qrsolv(awork,npts,xstemp,rcond)
ccc      call prin2(' rcond = *',rcond,1)
      do i = 1,npts
         coeffs(i) = xstemp(i)
      enddo
      return
      end

c
c
c
c
c
c**********************************************************************
      subroutine mkinterpmat(a,npts,sval)
c**********************************************************************
      implicit real *8 (a-h,o-z)
      real *8 a(npts,npts)
      real *8 sval(npts)
c
c     create Vandermonde matrix for interpolation coefficients.
c
c
      pi = 4*datan(1.0d0)
      dlength = sval(npts)-sval(1)
      do i = 1,npts
         a(i,1) = 1.0d0
         u = 2/dlength
         v = 1 - u*sval(npts) 
         xscal = u*sval(i) + v
         tjm2 = 1
         tjm1 = xscal
c-------------------------------------------
ccc         if (xscal.gt.1.0d0) then
ccc            arg = 0.0d0
ccc         else if (xscal.lt.-1.0d0) then
ccc            arg = pi
ccc         else
ccc            arg = dacos(xscal)
ccc         endif
c-------------------------------------------
c
         do j = 2,npts
ccc            tj = dcos( (j-1)* arg )
            a(i,j) = tjm1
            tj = 2*xscal*tjm1-tjm2
            tjm2 = tjm1
            tjm1 = tj
         enddo
      enddo
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine cheval(x,val,texp,n,a,b)
c**********************************************************************
C--------------
C
C     Subroutine evaluates Chebyshev expansion with coefficients 
C     TEXP at point X in interval [A,B].
C
C     A = left endpoint of interval
C     B = right endpoint of interval
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion (0 ,..., N-1)
C     VAL = computed value
C
C     In this implementation, if X is outside of [A,B] it moves X to 
C     the nearest endpoint. This prevents nonsense from happening when
C     calling DACOS(X). 
C
      implicit real *8 (a-h,o-z)
      integer *4 n,j
      real *8 texp(1),a,b,x,xscal,u,v,val,tj
C
      pi = 4*datan(1.0d0)
      u = 2/(b-a)
      v = 1 - u*b
      xscal = u*x + v
      call chebexev(xscal,val,texp,n-1)
ccc      call prin2(' chebexev val is *',val,1)
ccc      if (2.ne.3) goto 111
c
c     hack - make sure dacos doesn't mess up near 1,-1
c
ccc      if (xscal.gt.1.0d0) then
ccc         arg = 0.0d0
ccc      else if (xscal.lt.-1.0d0) then
ccc         arg = pi
ccc      else
ccc         arg = dacos(xscal)
ccc      endif
ccc      val = 0
ccc      do 600 j = 0,n-1
cccCC	 tj = dcos( j* dacos(xscal) )
ccc         tj = dcos( j* arg)
ccc         val = val + tj*texp(j+1)
ccc      call prin2(' in dumb loop val is *',val,1)
ccc600   continue
ccc      call prin2(' dumb loop val is *',val,1)
ccc111   continue
      return
      end
C
c
c
      subroutine chbdif(texp,tdiff,n,a,b)
C
C     Subroutine computes Chebyshev coefficients of derivative of
C     Chebyshev expansion. In other words, the (N-1)st degree expansion
C     TEXP is mapped to the (N-2) degree expansion TDIFF which
C     represents the derivative of the original expansion.
C
C     N = number of Chebyshev nodes
C     TEXP = array of expansion coefficients
C     TDIFF = array of expansion coefficients of derivative
C
      implicit real *8 (a-h,o-z)
      integer *4 n
      real *8 tdiff(n),texp(n),a,b,u,sc
c
c-----local variables
c
      integer *4 k
C
      u = 2/(b-a)
      sc = 2*u
      if (n .eq. 1) then
	 tdiff(1) = 0
      else if (n .eq. 2) then
	 tdiff(2) = 0
	 tdiff(1) = texp(2)*sc/2
      else
	 tdiff(n) = 0
	 tdiff(n-1) = (n-1)*texp(n)*sc
	 do k = 2,n-1
	    tdiff(n-k) = tdiff(n-k+2) + (n-k)*texp(n-k+1)*sc
         enddo
	 tdiff(1) = tdiff(1)/2
      endif
CCC   call prin2(' differentiated expansion is *',tdiff,n)
      return
      end
C
