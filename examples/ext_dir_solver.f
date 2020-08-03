      implicit real *8 (a-h,o-z)

      parameter (maxc = 10000)
 
      complex *16 zk,eye
      complex *16, allocatable :: xmat(:,:)
      real *8 pars(1000)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:)
      real *8, allocatable :: whts(:,:)
      integer, allocatable :: adjs(:,:)
      real *8, allocatable :: srccoefs(:,:,:),srcvals(:,:,:)

      complex *16 a2,b2,pars1(10),pars2(10)

      real *8 xyin(2)
      complex *16 grad(100)

      complex *16, allocatable :: zrhs(:,:)
      complex *16, allocatable :: zsoln(:,:)
      real *8 errs(1000),targ(2),src_normal(2),targ_normal(2)
      complex *16 zz,pot,potex,zpars(3),zid

      external fcurve,h2d_comb
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
      pars(1) = 1.0d0
      pars(2) = 0.3d0
      ier = 0
      ifclosed = 1
      chsmall = 1000

      xyin(1) = 0.01d0
      xyin(2) = -0.07d0

      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))


      nch = 0

      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurve,pars,nover,k,
     1       nch,chunks,adjs,ders,ders2,hs)

       call prinf('nch=*',nch,1)

      ntot = k*nch
      allocate(xmat(ntot,ntot))
      
      t1 = second()

      a2 = eye*(0.8 + 1.2*zk) 
      b2 = 1.0d0

      allocate(srcvals(8,k,nch),srccoefs(6,k,nch),whts(k,nch))

      call chunks_to_srcinfo(k,nch,chunks,ders,ders2,hs,srcvals,
     1  srccoefs,whts)

      ra = 0
      do ich=1,nch
        do j=1,k
          ra = ra + whts(j,ich)
        enddo
      enddo

      n = k*nch

      zpars(1) = zk
      zpars(2) = a2
      zpars(3) = b2
      
      call zgetmat_bdry(k,nch,n,srcvals,srccoefs,h2d_comb,8,0,
     1  dpars,3,zpars,0,ipars,xmat)
      
      allocate(zrhs(k,nch),zsoln(k,nch))

      do ich=1,nch
        do j=1,k
          call h2d_slp(xyin,8,srcvals(1,j,ich),ndd,dpars,1,zk,0,
     1       ipars,zrhs(j,ich))
        enddo
      enddo

      call prinf('n=*',n,1)

      call prin2('zrhs=*',zrhs,2*n)

      zid = b2*0.5d0
      numit = 200
      niter = 0

      eps = 1.0d-15

      call prin2('zid=*',zid,2)

      call zgmres_solver(n,xmat,zid,zrhs,numit,eps,niter,errs,rres,
     1   zsoln)
      call prin2('rres=*',rres,1)


c
c
c    test soln at an exterior point
c
      targ(1) = 5.1d0
      targ(2) = 3.1d0

      call h2d_slp(xyin,2,targ,0,dpars,1,zk,0,ipars,potex)
      call prin2('zpars=*',zpars,6)
      
      pot = 0
      do ich=1,nch
        do j=1,k
          call h2d_comb(targ,8,srcvals(1,j,ich),0,dpars,3,zpars,0,
     1      ipars,zz)
          pot = pot + zsoln(j,ich)*zz*whts(j,ich)
cc          call prin2('zz=*',zz,2)
        enddo
      enddo

      call prin2('pot=*',pot,2)
      call prin2('potex=*',potex,2)
      call prin2('abs error=*',abs(potex-pot),1)
      call prin2('rel error=*',abs(potex-pot)/abs(potex),1)
 
      return
      end
c-----------------------------------------------------------
c
c
c
c
c

        subroutine fcurve(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)

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
