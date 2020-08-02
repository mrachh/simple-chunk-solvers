      implicit real *8 (a-h,o-z)

      parameter (maxc = 10000)
 
      complex *16 zk,eye
      complex *16, allocatable :: xmat(:,:),xmat3(:,:)
      complex *16, allocatable :: xmat2(:,:)
      real *8, allocatable :: wgeo(:)
      real *8 pars(1000)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:)
      real *8, allocatable :: whts(:,:)
      integer, allocatable :: adjs(:,:)
      real *8, allocatable :: srccoefs(:,:,:),srcvals(:,:,:)

      complex *16 a2,b2,pars1(10),pars2(10)

      real *8 xyin(2)
      complex *16 grad(100)

      complex *16, allocatable :: zrhs(:)
      complex *16, allocatable :: zsoln(:)
      real *8 errs(1000),targ(2),src_normal(2),targ_normal(2)
      complex *16 zz,pot,potex,zpars(3)

      external gk2d_new,zkernel_dlp,zkernel_cfie,multa,fcurve
      external h2d_comb
      data eye/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = 4*atan(done)

      zk = 1.2d0 + 0.1*eye

      eps = 1.0d-12
      ifclosed = 1
      ier = 0

      lused = 0
      k = 16
      ibell = 1

      nover = 2
      ta = 0
      tb = 2*pi
      pars(1) = 1.0
      pars(2) = 0.3
      ier = 0
      ifclosed = 1
      chsmall = 1000

      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))


      nch = 0
      call prin2('ta=*',ta,1)
      call prin2('tb=*',tb,1)


      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurve,pars,nover,k,
     1       nch,chunks,adjs,ders,ders2,hs)

       call prinf('nch=*',nch,1)

       lused = 100*k*nch+10
       allocate(wgeo(lused))

       lused = 0

      call chunkpack(k,nch,chunks,adjs,ders,ders2,hs,wgeo,lused)

      call chunkunpack1(wgeo,k,nch,ichunks,iadjs,iders,iders2,ihs)

      ntot = k*nch
      allocate(xmat(ntot,ntot),xmat2(ntot,ntot))
      allocate(whts(k,nch))
      
      call chunkwhts(k,nch,wgeo(ichunks),wgeo(iders),wgeo(ihs),whts)

      norder = k

cc    form matrix corresponding to the double layer

      call prin2('zk=*',zk,2)
      call prinf('norder=*',norder,1)
      
      t1 = second()

      a2 = eye*(0.8 + 1.2*zk) 
      b2 = 1

      call zbuildmat(norder,wgeo,zkernel_cfie,a2,b2,gk2d_new,zk,
     1      pars1,pars2,nout,xmat)
      print *, "after getting first matrix"

      t2 = second()
      call prin2('total time=*',t2-t1,1)

      allocate(srcvals(8,k,nch),srccoefs(6,k,nch))

      call chunks_to_srcinfo(k,nch,chunks,ders,ders2,hs,srcvals,
     1  srccoefs)

      n = k*nch

      zpars(1) = zk
      zpars(2) = a2
      zpars(3) = b2
      
      call zgetmat_bdry(k,nch,n,srcvals,srccoefs,h2d_comb,8,0,
     1  dpars,3,zpars,0,ipars,xmat2)

      ichtest = 5
      erra = 0
      ra = 0
      istart = (ichtest-1)*k
      do i=1,n
        do j=1,n
          ra = ra + abs(xmat(j,i))**2
          erra = erra + abs(xmat2(j,i)-xmat(j,i))**2
        enddo
      enddo

      ii = k+1
      j = 1
      print *, xmat2(17,1)
      print *, xmat(17,1)


      call prin2('erra=*',sqrt(erra),1)
      call prin2('ra=*',sqrt(ra),1)

cc      call prin2('xmat2=*',xmat2(istart+1:istart+k,istart+1:istart+k),
cc     1  2*k*k)
cc      call prin2('xmat=*',xmat(istart+1:istart+k,istart+1:istart+k),
cc     1  2*k*k)

      erra = sqrt(erra/ra)
      call prin2('error in matrix entries=*',erra,1)
 
      return
      end
c-----------------------------------------------------------


      subroutine gk2d_new(zk, src, targ, pars1, pars2, val,grad, hess)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: src(2), targ(2)
        complex *16, intent(in) :: zk, pars1(*), pars2(*)
        complex *16, intent(out) :: val, grad(2), hess(2,2)

        complex *16 :: z, h0, h1, h1p
        complex *16, parameter :: ima = (0,1)

        !
        ! return the value of h0, its gradient, and its hessian using
        ! a call to hank103.
        !
        ! input:
        !   zk - complex wavenumber
        !   src - source location
        !   targ - target location
        !   pars1 - dummy variables
        !   pars2 - dummy variables
        !
        ! output:
        !   val - the value of H_0(kr)
        !   grad - the gradient *at the target* of H_0(kr)
        !   hess - the hessian *at the target* of H_0(kr)
        !
        ! NOTE: the source gradient and hessian, and mixed hessian
        !   can be computed from the target derivatives. they are
        !   left out of the calling sequence to simplify things.
        !

        ifexpon = 1
        dx = targ(1)-src(1)
        dy = targ(2)-src(2)
        r = sqrt(dx**2 + dy**2)

        z = zk*r
        call hank103(z, h0, h1, ifexpon)
        val = ima*h0/4
        grad(1) = -dx*ima/4*h1*zk/r
        grad(2) = -dy*ima/4*h1*zk/r

        h1p = h0 - h1/z
        hess(1,1) = ima*zk/4/r*(-h1 + dx**2/r/r*h1 - 
     1       zk*dx**2*h1p/r)
        hess(1,2) = ima*zk/4/r*(dx*dy/r/r*h1 - zk*dx*dy*h1p/r)
        hess(2,1) = ima*zk/4/r*(dx*dy/r/r*h1 - zk*dx*dy*h1p/r)
        hess(2,2) = ima*zk/4/r*(-h1 + dy**2/r/r*h1 - 
     1       zk*dy**2*h1p/r)

        return
      end subroutine
 
      
c-----------------------------------------------------------


        subroutine gk2d(zk,xy,xy0,cval,ifgrad,cgrad)
        implicit real *8 (a-h,o-z)
        real *8 xy(2),xy0(2)
        complex *16 zk,ima,cval,cgrad(2),z,h0,h1
c
c       evaluates the 2d green's function for the helmholtz equation,
c       which we define to be i/4 * H_0(x-y) (note the absense of a minus)
c
c
        ima=(0,1)
c
        r=(xy(1)-xy0(1))**2+(xy(2)-xy0(2))**2
        r=sqrt(r)
        z=r*zk
        ifexpon=1
        call hank103(z,h0,h1,ifexpon)
        cval=ima/4*h0
c
        if (ifgrad .eq. 0) return
c
        cgrad(1)=-ima/4*h1*(xy(1)-xy0(1))*zk/r
        cgrad(2)=-ima/4*h1*(xy(2)-xy0(2))*zk/r
c
        if (ifgrad .eq. -1) then
        cgrad(1)=-cgrad(1)
        cgrad(2)=-cgrad(2)
        endif
c
        return
        end

c-----------------------------------------------------------------

       subroutine multa(a,p1,p2,p3,p4,x,y,n)
       implicit real *8 (a-h,o-z)

       complex *16 a(n,n),x(1),y(1),d

       do 1400 i=1,n

       d = 0

       do 1200 j=1,n

       d = d + a(i,j)*x(j)

 1200  continue

       y(i) = d

 1400  continue       

       return
       end
c------------------------------------------------------------
        subroutine fcurve(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)
c
c        radx=1
c        rady=1
c
c

        rad0 = pars(1)
        radin = pars(2)

        rr = rad0 + radin*cos(3*t)
        rp = -3*radin*sin(3*t)
        rpp = -9*radin*cos(3*t)
        
        x=rr*cos(t)
        y=rr*sin(t)
c
        dxdt=-rr*sin(t) + rp*cos(t)
        dydt=rr*cos(t) + rp*sin(t)
c
        dxdt2=(rpp-rr)*cos(t) - 2*rp*sin(t) 
        dydt2=(rpp-rr)*sin(t) + 2*rp*cos(t)
c
        return
        end
c
