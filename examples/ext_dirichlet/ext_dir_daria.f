      implicit real *8 (a-h,o-z)

      parameter (maxc = 10000)
      parameter (maxp = 100)
 
      complex *16 zk,eye,uin
      complex *16, allocatable :: xmat(:,:)
      real *8 pars(1000)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:)
      real *8, allocatable :: whts(:,:)
      integer, allocatable :: idt(:)
      integer, allocatable :: iregionl(:)
      integer, allocatable :: iregionr(:)
      integer, allocatable :: adjs(:,:)
      real *8, allocatable :: srccoefs(:,:,:),srcvals(:,:,:)
      real *8 parsall(maxp+3,4)

      complex *16 a2,b2,pars1(10),pars2(10)

      real *8 xyin(2)
      complex *16 grad(100)

      complex *16, allocatable :: zrhs(:,:)
      complex *16, allocatable :: zsoln(:,:)
      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: xtarg(:),ytarg(:)
      real *8 errs(1000),targ(2),src_normal(2),targ_normal(2)
      integer, allocatable :: isin(:)
      complex *16, allocatable :: pottarg(:),pottargex(:)
      complex *16 zz,pot,potex,zpars(3),zid

      external fcurve,h2d_comb,fcurve1
      data eye/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = 4*atan(done)

      call getspecs(iscat,itot,depth,zk,parsall,maxp,norder,thetain)

      write(6,*)'zk = ',zk
      write(6,*)'thetain = ',thetain
      eps = 1.0d-9
      ier = 0
c
c     k is order of accuracy of chunk discretization
c     chsmall is dyadic refinement parameter for corners
c
      k = 16
      chsmall = 1.0d-3
      ta = 0
      tb = 1.0d0
c
c     xyin is an interior ppoint
c
      xyin(1) = 0.1d0
      xyin(2) = -1.52d0

      nover = 1
      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))
c
c     npw = points per wavelength
c     ireflinel,irefiner = 1 means refine to corner
c
      npw = 10
      irefinel = 1
      irefiner = 1
      ifclosed = 0
      nch = 0
c
c     straight sides
c
      do ii = 1,3
         call chunkfunczk(zk,npw,eps,ifclosed,irefinel,
     1     irefiner,chsmall,ta,tb,fcurve1,parsall(1,ii),nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     1       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
         nch = nch+nch1
      call prinf('nch=*',nch,1)
      enddo
c
c     top side defined by sine series
c
      call chunkfunczk(zk,npw,eps,ifclosed,irefinel,
     1     irefiner,chsmall,ta,tb,fcurve,parsall(1,4),nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     1       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
      nch = nch+nch1
      call prinf('nch=*',nch,1)
c
      ntot = k*nch
      allocate(xmat(ntot,ntot))
      allocate(iregionl(nch))
      allocate(iregionr(nch))
c      
      t1 = second()
      a2 = eye*(0.8 + 1.2*zk)
      b2 = 1.0d0
c
      allocate(srcvals(8,k,nch),srccoefs(6,k,nch),whts(k,nch))
      call chunks_to_srcinfo(k,nch,chunks,ders,ders2,hs,srcvals,
     1  srccoefs,whts)
c
      n = k*nch
      zpars(1) = zk
      zpars(2) = a2
      zpars(3) = b2
c
c  get the matrix
c
      call cpu_time(t1)
C$       t1 = omp_get_wtime()

ccc      if (2.ne.3) goto 111

      call helm2d_comb_dir_mat(k,nch,n,srcvals,srccoefs,zpars,xmat)

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      
      call prin2('matrix generation time=*',t2-t1,1)

      allocate(zrhs(k,nch),zsoln(k,nch))

      do ich=1,nch
        do j=1,k
           if (iscat.eq.0) then
             call h2d_slp(xyin,8,srcvals(1,j,ich),0,dpars,1,zk,0,
     1       ipars,zrhs(j,ich))
           else 
             uin = cdexp(-eye*zk*srcvals(2,j,ich))
             zrhs(j,ich) = -uin
           endif
        enddo
      enddo
      zid = b2/2
      numit = 200
      niter = 0
      eps = 1.0d-12
      call cpu_time(t1)
      call zgmres_solver(n,xmat,zid,zrhs,numit,eps,niter,errs,rres,
     1   zsoln)
      call cpu_time(t2)
      
      call prin2('solve time=*',t2-t1,1)

c
111   continue
      open(unit=17,file='bdrypts.m',status='unknown')
      write(17,*) ' xx = ['
      do ich=1,nch
        do j=1,k
          write(17,*) srcvals(1,j,ich), srcvals(2,j,ich)
        enddo
      enddo
      write(17,*) ' ];'
      write(17,*) ' plot(xx(:,1),xx(:,2),''.'')'
      write(17,*) ' hold on'
ccc      stop
      
      open(unit=9,file='targ.dat',status='unknown')
      read(9,*) ntarg
      write(6,*) ntarg

      allocate(targs(2,ntarg),pottargex(ntarg),
     1   pottarg(ntarg))
      allocate(xtarg(ntarg),ytarg(ntarg))
      allocate(idt(ntarg))
c
      do i = 1,ntarg
         read(9,*) xtarg(i), ytarg(i)
         targs(1,i) = xtarg(i)
         targs(2,i) = ytarg(i)
         write(6,*) xtarg(i), ytarg(i)
      enddo
c
      do i = 1,nch
         iregionl(i) = 1
         iregionr(i) = 0
      enddo
c
      if (iscat.eq.0) then
         do i = 1,ntarg
            call h2d_slp(xyin,2,targs(1,i),0,dpars,1,zk,0,ipars,
     1       pottargex(i))
         enddo
      endif
c
      call helm2d_comb_dir_targ(k,nch,n,srccoefs,srcvals,zpars,zsoln,
     1  ntarg,targs,pottarg)

      if ((iscat.eq.1).and.(itot.eq.1)) then
         do i = 1,ntarg
            uin = cdexp(-eye*zk*targs(2,i))
            pottarg(i) = pottarg(i) + uin
         enddo
      endif
      erra = 0
      ra = 0
      do i=1,ntarg
          err = abs(pottargex(i)-pottarg(i))
          ra = ra + abs(pottargex(i))**2
          erra = erra + abs(pottargex(i)-pottarg(i))**2
      enddo
      open(unit = 20,file = 'pot.m',status='unknown')
      write(20,*) 'zpot = ['
      do i = 1,ntarg
         write(20,*) pottarg(i)
      enddo
      write(20,*) ']'

      if (iscat.eq.0) then
         erra = sqrt(erra/ra)
         call prin2('relative l2 error on grid of targets=*',erra,1)
      endif
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
c
c
c       seg is   p1  ------- p2 , sin series 
c
        pi = 4.0d0*datan(1.0d0)
c
        x = pars(1) + t*(pars(2)-pars(1))
        dxdt= pars(2)-pars(1)
        dxdt2= 0.0d0
c
        norder = nint(pars(3))
ccc        y = x*(pi-x)
ccc        dydt=-2*x + pi
ccc        dydt2=-2.0d0
ccc        y = x*(pi-x)
ccc        dydt=(-2*x+pi)*dxdt
ccc        dydt2=-2.0d0*dxdt*dxdt
        y = 0.d0
        dydt= 0.0d0
        dydt2= 0.0d0
        do i = 1,norder
           rinc = pars(3+i)*dsin(i*x)
           rincp = pars(3+i)*dcos(i*x)
           y= y+ rinc
           dydt=dydt+ i*rincp*dxdt
           dydt2=dydt2- rinc*(i*dxdt)**2
        enddo
c
        return
        end
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
