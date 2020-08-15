c
c     test driver that solves dielectric interface problem
c     on square for TE mode 
c
c
      implicit real *8 (a-h,o-z)

      parameter (maxc = 10000)
 
      complex *16 zk0,zk1,eye
      complex *16, allocatable :: xmat(:,:)
      real *8 pars(1000)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:)
      real *8, allocatable :: whts(:,:)
      integer, allocatable :: adjs(:,:)
      integer, allocatable :: iregionl(:)
      integer, allocatable :: iregionr(:)
      real *8, allocatable :: srccoefs(:,:,:),srcvals(:,:,:)

      complex *16 a2,b2,pars1(10),pars2(10)

      real *8 xyin(2),xyout(2)
      real *8 srcinfoin(8)
      real *8 srcinfoout(8)
      complex *16 grad(100)

      complex *16, allocatable :: zrhs(:,:,:)
      complex *16, allocatable :: zsoln(:,:,:)
      real *8 errs(1000),targ(2),src_normal(2),targ_normal(2)
      complex *16 zz,pot,potex,zpars(6),zid
      complex *16 pots,potd,gx,gy
      complex *16 zpar3(3),zq

      external fcurve,h2d_comb
      data eye/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = 4*atan(done)

      zk0 = 1.0d0 + 0.0d0*eye
      zk1 = 1.3d0 + 0.0d0*eye

      eps = 1.0d-7

      k = 16

      nover = 2
      ta = 0
      tb = 2*pi
      pars(1) = 1.0d0
      pars(2) = 0.3d0

      ier = 0
      ifclosed = 1
      chsmall = 1000

      xyin(1) = 0.3d0
      xyin(2) = 0.2d0
      xyout(1) = 6.3d0
      xyout(2) = 6.2d0

      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))
      allocate(iregionl(maxc),iregionr(maxc))
c
c     chunk up each side of square (positively oriented).
c     routine fcurve1 defines a side based on 
c     x in [pars(1),pars(2), y in [pars(3)mparrs(4)] 
c
      nch = 0

      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurve,pars,nover,k,
     1       nch,chunks,adjs,ders,ders2,hs)

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
      allocate(xmat(2*ntot,2*ntot))
      
      t1 = second()


      allocate(srcvals(8,k,nch),srccoefs(6,k,nch),whts(k,nch))

      call chunks_to_srcinfo(k,nch,chunks,ders,ders2,hs,srcvals,
     1  srccoefs,whts)

      n = k*nch
c
      a0 = 1.0d0
      b0 = 1.2d0
      a1 = 1.3d0
      b1 = 1.4d0
      zpars(1) = zk0
      zpars(2) = a0
      zpars(3) = b0
      zpars(4) = zk1
      zpars(5) = a1
      zpars(6) = b1
      zq = 0.5*(zpars(2)/zpars(3)+zpars(5)/zpars(6))
      call prin2(' zq is *',zq,2)
c
c  get the matrix
c


      call cpu_time(t1)
C$       t1 = omp_get_wtime()

      call helm2d_dielec_mat(k,nch,n,srcvals,srccoefs,zpars,
     1     iregionl,iregionr,xmat)

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      
      call prin2('matrix generation time=*',t2-t1,1)

      allocate(zrhs(k,nch,2),zsoln(k,nch,2))

      do ich=1,nch
        do j=1,k
           call h2d_slp(xyin,8,srcvals(1,j,ich),0,dpars,1,zk0,0,
     1       ipars,zz)
           zrhs(j,ich,1) = zz*zpars(2)/zq
           call h2d_slp(xyout,8,srcvals(1,j,ich),0,dpars,1,zk1,0,
     1       ipars,zz)
          zrhs(j,ich,1) = zrhs(j,ich,1) - zz*zpars(5)/zq
        enddo
      enddo

      do ich=1,nch
        do j=1,k
           call h2d_sprime(xyin,8,srcvals(1,j,ich),0,dpars,1,zk0,0,
     1       ipars,zz)
           zrhs(j,ich,2) = zz*zpars(3)
           call h2d_sprime(xyout,8,srcvals(1,j,ich),0,dpars,1,zk1,0,
     1       ipars,zz)
          zrhs(j,ich,2) = zrhs(j,ich,2) - zz*zpars(6)
        enddo
      enddo

ccc      zid = b2/2
      zid = 1.0d0
      numit = 200
      niter = 0

      eps = 1.0d-15

      call prin2(' zrhs1 is *',zrhs(1,1,1),k*nch*2)
      call prin2(' zrhs2 is *',zrhs(1,1,2),k*nch*2)
      call zgmres_solver(2*n,xmat,zid,zrhs,numit,eps,niter,errs,rres,
     1   zsoln)
      call prin2(' zsoln1 is *',zsoln(1,1,1),k*nch*2)
      call prin2(' zsoln2 is *',zsoln(1,1,2),k*nch*2)
      
c
c
c    test soln at an exterior point
c
      targ(1) = 5.31d0
      targ(2) = 3.33d0

      call h2d_slp(xyin,2,targ,0,dpars,1,zk0,0,ipars,potex)

      zpar3(1) = zk0
      zpar3(2) = -1.0d0/zpars(3)
      zpar3(3) = 1.0d0/zpars(3)
      
      pot = 0
      pots = 0
      potd = 0
      do ich=1,nch
        do j=1,k
          zz = 0
          call h2d_dlp(srcvals(1,j,ich),2,targ,0,dpars,1,zk0,0,
     1      ipars,zz)
          pot = pot + zsoln(j,ich,1)*zz*whts(j,ich)/zpars(3)
ccc             pot = pot + zrhs(j,ich,1)*zz*whts(j,ich)/zpars(3)
             potd = potd + zsoln(j,ich,1)*zz*whts(j,ich)/zpars(3)
          call h2d_slp(srcvals(1,j,ich),2,targ,0,dpars,1,zk0,0,
     1      ipars,zz)
          pot = pot - zsoln(j,ich,2)*zz*whts(j,ich)/zpars(3)
ccc          pot = pot - zrhs(j,ich,2)*zz*whts(j,ich)/zpars(3)
          pots = pots - zsoln(j,ich,2)*zz*whts(j,ich)/zpars(3)
cc          call prin2('zz=*',zz,2)
        enddo
      enddo

      call prin2('pots=*',pots,2)
      call prin2('potd=*',potd,2)
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
c
