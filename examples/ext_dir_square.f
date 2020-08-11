c
c     test driver that solves exterior Dirichlet problem
c     for a square. 
c
c
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

      external fcurve1,h2d_comb
      data eye/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = 4*atan(done)

      zk = 1.2d0 + 0.0d0*eye

      eps = 1.0d-7
      ifclosed = 1
      ier = 0

      lused = 0
      k = 16

      nover = 1
      ta = 0
      tb = 1.0d0
      ier = 0
      ifclosed = 0
      chsmall = 1.0d-7

      xyin(1) = 0.3d0
      xyin(2) = 0.2d0

      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))
c
c     chunk up each side of square (positively oriented).
c     routine fcurve1 defines a side based on 
c     x in [pars(1),pars(2), y in [pars(3)mparrs(4)] 
c
      nch = 0
      pars(1) = 0.0d0
      pars(2) = 1.0d0
      pars(3) = 0.0d0
      pars(4) = 0.0d0
      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurve1,pars,nover,k,
     1       nch1,chunks,adjs,ders,ders2,hs)
      nch = nch+nch1
c
      pars(1) = 1.0d0
      pars(2) = 1.0d0
      pars(3) = 0.0d0
      pars(4) = 1.0d0
      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurve1,pars,nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     2       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
      nch = nch+nch1
c
ccc      if (2.ne.3) goto 111
      pars(1) = 1.0d0
      pars(2) = 0.0d0
      pars(3) = 1.0d0
      pars(4) = 1.0d0
      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurve1,pars,nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     2       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
      nch = nch+nch1
c
      pars(1) = 0.0d0
      pars(2) = 0.0d0
      pars(3) = 1.0d0
      pars(4) = 0.0d0
      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurve1,pars,nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     2       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
      nch = nch+nch1
111   continue
       call prinf('nch=*',nch,1)
      open (unit=17,status='unknown',file='pts.m')
      write(17,*) ' xpts = ['
      do i = 1,nch
      do j = 1,k
         write(17,*) chunks(1,j,i), chunks(2,j,i)
      enddo
      enddo
      write(17,*) ' ];'
      write(17,*) ' plot(xpts(:,1),xpts(:,2),''.'')'
      ntot = k*nch
      allocate(xmat(ntot,ntot))
      
      t1 = second()

      a2 = eye*(0.8 + 1.2*zk)
      b2 = 1.0d0

      allocate(srcvals(8,k,nch),srccoefs(6,k,nch),whts(k,nch))

      call chunks_to_srcinfo(k,nch,chunks,ders,ders2,hs,srcvals,
     1  srccoefs,whts)

      n = k*nch

      zpars(1) = zk
      zpars(2) = a2
      zpars(3) = b2

c
c  get the matrix
c


      call cpu_time(t1)
C$       t1 = omp_get_wtime()

      call helm2d_comb_dir_mat(k,nch,n,srcvals,srccoefs,zpars,xmat)

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      
      call prin2('matrix generation time=*',t2-t1,1)

      allocate(zrhs(k,nch),zsoln(k,nch))

      do ich=1,nch
        do j=1,k
          call h2d_slp(xyin,8,srcvals(1,j,ich),0,dpars,1,zk,0,
     1       ipars,zrhs(j,ich))
        enddo
      enddo

      zid = b2/2
      numit = 200
      niter = 0

      eps = 1.0d-15

      call zgmres_solver(n,xmat,zid,zrhs,numit,eps,niter,errs,rres,
     1   zsoln)

c
c
c    test soln at an exterior point
c
      targ(1) = 5.31d0
      targ(2) = 3.33d0

      call h2d_slp(xyin,2,targ,0,dpars,1,zk,0,ipars,potex)
      
      pot = 0
      do ich=1,nch
        do j=1,k
          zz = 0
          call h2d_comb(srcvals(1,j,ich),2,targ,0,dpars,3,zpars,0,
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

        subroutine fcurve1(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)
c
c       polygonal segment defined by 
c       x in [pars(1),pars(2), y in [pars(3)mparrs(4)] 
c       assuming t in [0,1].
c
        x = pars(1) + t*(pars(2)-pars(1))
        y = pars(3) + t*(pars(4)-pars(3))

        dxdt=pars(2)-pars(1)
        dydt=pars(4)-pars(3)
c
        dxdt2=0.0d0
        dydt2=0.0d0
c
        return
        end
c
