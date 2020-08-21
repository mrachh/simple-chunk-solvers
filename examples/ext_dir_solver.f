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
      real *8, allocatable :: targs(:,:)
      real *8 errs(1000),targ(2),src_normal(2),targ_normal(2)
      integer, allocatable :: isin(:)
      complex *16, allocatable :: pottarg(:),pottargex(:)
      complex *16 zz,pot,potex,zpars(3),zid

      external fcurve,h2d_comb
      data eye/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = 4*atan(done)

      zk = 1.2d0 + 0.0d0*eye

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

      xyin(1) = 0.1d0
      xyin(2) = -0.02d0

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

      nlat = 300
      ntarg = nlat*nlat

      allocate(targs(2,ntarg),isin(ntarg),pottargex(ntarg),
     1   pottarg(ntarg))

      do ix=1,nlat
        do iy=1,nlat
          itt = (ix-1)*nlat + iy
          targs(1,itt) = -3  +6*(ix-1.0d0)/(nlat-1.0d0)
          targs(2,itt) = -3  +6*(iy-1.0d0)/(nlat-1.0d0)
          isin(itt) = 1
          rr = targs(1,itt)**2 + targs(2,itt)**2
          thet = atan2(targs(2,itt),targs(1,itt))
          rtest = (pars(1) + pars(2)*cos(3*thet))**2
          if(rr.gt.rtest) isin(itt) = 0
          call h2d_slp(xyin,2,targs(1,itt),0,dpars,1,zk,0,ipars,
     1       pottargex(itt))
        enddo
      enddo

      call helm2d_comb_dir_targ(k,nch,n,srccoefs,srcvals,zpars,zsoln,
     1  ntarg,targs,pottarg)

      erra = 0
      ra = 0
      do i=1,ntarg
        if(isin(i).eq.0) then
          err = abs(pottargex(i)-pottarg(i))
          ra = ra + abs(pottargex(i))**2
          erra = erra + abs(pottargex(i)-pottarg(i))**2

          write(33,*) i,err,isin(i),targs(1,i),targs(2,i) 
        endif
      enddo

      erra = sqrt(erra/ra)
      call prin2('relative l2 error on grid of targets=*',erra,1)
 
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
