      implicit real *8 (a-h,o-z)
      real *8, allocatable :: targs(:,:),parsall(:,:)
      complex *16 zk
      character *300 fsol,ftarg

cc      call prini(6,13)
      done = 1
      pi = atan(done)*4


      maxp = 100
      ntarg = 128
      allocate(parsall(maxp+3,4),targs(2,ntarg))
c
c   generate targets
c 
c
      do i=1,ntarg
        targs(1,i) = -pi + 3*pi*(i-1.0d0)/(ntarg-1.0d0)
        targs(2,i) = 3*pi/2
      enddo

      iscat = 1
      itot = 1
      depth = pi
      norder = 4

      nbatsize = 40
      nbat = 2
      nzk = 37
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibat,iconf_start,iconf_end)
C$OMP$PRIVATE(fsol,ftarg,iconf,parsall,rl1,rsc,izk,zk,i)
C$OMP$SCHEDULE(DYNAMIC)
      do ibat = 1,nbat
        print *, "starting ibat=",ibat
        iconf_start = (ibat-1)*nbatsize + 1
        iconf_end = ibat*nbatsize
        write(fsol,'(a,i4.4,a,i4.4,a,i2.2,a,i2.2,a,i1,a)') 
     1   'data/soln_iconf',iconf_start,'-',iconf_end,'_nzk',nzk,
     2   '_norder',norder,'_iscat',iscat,'.dat'
        write(fsol,'(a,i1,a)') 'sol_ibat',ibat,'.dat'
        write(ftarg,'(a,i4.4,a,i4.4,a,i2.2,a,i2.2,a,i1,a)') 
     1   'data/targ_iconf',iconf_start,'-',iconf_end,'_nzk',nzk,
     2   '_norder',norder,'_iscat',iscat,'.dat'
        write(ftarg,'(a,i1,a)') 'targ_ibat',ibat,'.dat'
        open(unit=33+ibat,file=trim(ftarg))
        open(unit=233+ibat,file=trim(fsol))
        write(33+ibat,*) ntarg,norder
        do i=1,ntarg
          write(33+ibat,*) targs(1,i),targs(2,i)
        enddo
        close(33+ibat)
        close(233+ibat)
        parsall  = 0
        parsall(1,1) =  0.0d0
        parsall(3,1) =  0.0d0
        parsall(2,1) =  0.0d0
        parsall(4,1) = -depth
c
        parsall(1,2) = 0.0d0
        parsall(3,2) = -depth
        parsall(2,2) = pi
        parsall(4,2) = -depth
c
        parsall(1,3) = pi
        parsall(3,3) = -depth
        parsall(2,3) = pi
        parsall(4,3) = 0.0d0
c

        do iconf = iconf_start,iconf_end
          if(mod(iconf,10).eq.1) print *, "starting iconf=",iconf
          parsall(1,4) = pi
          parsall(2,4) = 0.0d0
          parsall(3,4) = norder

          rl1 = 0
          do i=1,norder
            parsall(3+i,4) = hkrand(0) - 0.5d0
            rl1 = rl1 + abs(parsall(3+i,4))
          enddo
          rsc = hkrand(0)*0.5d0
          do i=1,norder
            parsall(3+i,4) = parsall(3+i,4)/rl1*rsc
          enddo

          do izk = 1,nzk
            if(mod(izk,9).eq.1) print *, iconf,izk
            zk = 1+(izk-1)*0.25d0
            open(unit=33+ibat,file=trim(ftarg),access='append')
            open(unit=233+ibat,file=trim(fsol),access='append')
            write(33+ibat,*) iconf,real(zk)
            write(233+ibat,*) iconf,real(zk)
            write(33+ibat,*) depth
            write(233+ibat,*) depth
            do i=1,norder
              write(33+ibat,*) parsall(3+i,4)
              write(233+ibat,*) parsall(3+i,4)
            enddo
            close(33+ibat)
            close(233+ibat)
            call daria_datgen(ibat,iscat,itot,depth,zk,parsall,maxp,
     1        norder,ntarg,targs,ftarg,fsol)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO      


      stop
      end



      subroutine daria_datgen(ibat,iscat,itot,depth,zk,parsall,maxp,
     1   norder,ntarg,targs,ftarg,fsol)
      implicit real *8 (a-h,o-z)

      parameter (maxc = 10000)
 
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
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      integer ipars

      complex *16 a2,b2,pars1(10),pars2(10)

      real *8 xyin(2)
      complex *16 grad(100)

      complex *16, allocatable :: zrhs(:,:)
      complex *16, allocatable :: zsoln(:,:)
      real *8 targs(2,ntarg)
      real *8, allocatable :: xtarg(:),ytarg(:)
      real *8 errs(1000),targ(2),src_normal(2),targ_normal(2)
      integer, allocatable :: ixys(:),norders(:),iptype(:)
      real *8, allocatable :: ab(:,:)
      integer, allocatable :: isin(:)
      complex *16, allocatable :: pottarg(:),pottargex(:)
      complex *16 zz,pot,potex,zpars(3),zid

      character (len=*) fsol,ftarg


      external fcurve,h2d_comb,fcurve1
      data eye/(0.0d0,1.0d0)/


      done = 1
      pi = 4*atan(done)

      eps = 1.0d-4
      ier = 0
c
c     k is order of accuracy of chunk discretization
c     chsmall is dyadic refinement parameter for corners
c
      k = 16
      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      chsmall = 1.0d-2*pi
      ta = 0
      tb = 1.0d0
c
c     xyin is an interior ppoint
c
      xyin(1) = 0.6d0
      xyin(2) = -0.77d0

      nover = 1
      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))
      allocate(srcvals(8,k,maxc),srccoefs(6,k,maxc),whts(k,maxc))
      allocate(ixys(maxc+1),norders(maxc),iptype(maxc),ab(2,maxc))
c
c     npw = points per wavelength
c     ireflinel,irefiner = 1 means refine to corner
c
      npw = 20
      irefinel = 1
      irefiner = 1
      ifclosed = 0
      nch = 0
c
c     straight sides
c
      ndd = 4
      ndz = 0
      ndi = 0
      nwav = real(zk)*rlmax/2/pi
      rlmax = (k+0.0d0)/real(zk)/npw*2*pi
      do ii = 1,3
         nch1 = 0
         ier = 0
         call chunkfunc_guru(eps,rlmax,ifclosed,irefinel,
     1    irefiner,chsmall,ta,tb,fcurve1,ndd,parsall(1,ii),
     2    ndz,zpars,ndi,ipars,nover,k,maxc,nch1,norders(nch+1),
     3    ixys,iptype(nch+1),n0,srcvals(1,1,nch+1),srccoefs(1,1,nch+1),
     4    ab(1,nch+1),adjs(1,nch+1),ier)

         nch = nch+nch1
      enddo
c
c     top side defined by sine series
c
      ndd0 = nint(parsall(3,4))
      ndd = ndd0 + 3
      ier = 0
      call chunkfunc_guru(eps,rlmax,ifclosed,irefinel,
     1  irefiner,chsmall,ta,tb,fcurve,ndd,parsall(1,4),
     2  ndz,zpars,ndi,ipars,nover,k,maxc,nch1,norders(nch+1),
     3  ixys,iptype(nch+1),n0,srcvals(1,1,nch+1),srccoefs(1,1,nch+1),
     4  ab(1,nch+1),adjs(1,nch+1),ier)

      nch = nch+nch1

cc      call prinf('nch=*',nch,1)
c
      ntot = k*nch
      allocate(xmat(ntot,ntot))
      allocate(iregionl(nch))
      allocate(iregionr(nch))
c      
      t1 = second()
      a2 = eye*(0.8 + 1.2*zk)
      b2 = 1.0d0
      do i=1,nch
        do j=1,k
          rr = sqrt(srcvals(3,j,i)**2 + srcvals(4,j,i)**2)
          whts(j,i) = rr*wts(j) 
        enddo
      enddo
c
      n = k*nch
      zpars(1) = zk
      zpars(2) = a2
      zpars(3) = b2
c
c  get the matrix
c

      call helm2d_comb_dir_mat(k,nch,n,srcvals,srccoefs,zpars,xmat)

      
cc      call prin2('matrix generation time=*',t2-t1,1)

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
      eps = 1.0d-8
      call zgmres_solver(n,xmat,zid,zrhs,numit,eps,niter,errs,rres,
     1   zsoln)
      
cc      call prin2('solve time=*',t2-t1,1)
      
      allocate(pottargex(ntarg),
     1   pottarg(ntarg))
c

      do i=1,ntarg
        pottargex(i) = 0
        pottarg(i) = 0

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

      if (iscat.eq.0) then
         erra = sqrt(erra/ra)
cc         call prin2('relative l2 error on grid of targets=*',erra,1)
      endif


      open(unit=33+ibat,file=trim(ftarg),access='append')
      open(unit=233+ibat,file=trim(fsol),access='append')
      do i=1,ntarg
        write(33+ibat,'(2(2x,e22.16))') real(pottarg(i)),
     1    imag(pottarg(i))
      enddo
      
      npts = nch*k
      write(233+ibat,*) nch,k
      do i=1,nch
        do j=1,k
          write(233+ibat,'(10(2x,e22.16))') srcvals(1,j,i),
     1     srcvals(2,j,i),
     1     srcvals(3,j,i),srcvals(4,j,i),srcvals(5,j,i),srcvals(6,j,i),
     2     real(zrhs(j,i)),imag(zrhs(j,i)),real(zsoln(j,i)),
     3     imag(zsoln(j,i))
        enddo
      enddo
      if(iscat.eq.0) write(233+ibat,*) erra
      close(233+ibat)
      close(33+ibat)


      return
      end
c-----------------------------------------------------------
c
c
c
c
c

        subroutine fcurve(t,ndd,dpars,ndz,zpars,ndi,ipars,
     1     x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 dpars(ndd)
        complex *16 zpars(ndz)
        integer ipars(ndi)
c
c
c       seg is   p1  ------- p2 , sin series 
c
        pi = 4.0d0*datan(1.0d0)
c
        x = dpars(1) + t*(dpars(2)-dpars(1))
        dxdt= dpars(2)-dpars(1)
        dxdt2= 0.0d0
c
        norder = nint(dpars(3))
        y = 0.d0
        dydt= 0.0d0
        dydt2= 0.0d0
        do i = 1,norder
           rinc = dpars(3+i)*dsin(i*x)
           rincp = dpars(3+i)*dcos(i*x)
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
        subroutine fcurve1(t,ndd,dpars,ndz,zpars,ndi,ipars,
     1     x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 dpars(ndd)
        complex *16 zpars(ndz)
        integer ipars(ndi)
c
c       polygonal segment defined by 
c       x in [pars(1),pars(2), y in [pars(3)mparrs(4)] 
c       assuming t in [0,1].
c
        x = dpars(1) + t*(dpars(2)-dpars(1))
        y = dpars(3) + t*(dpars(4)-dpars(3))

        dxdt=dpars(2)-dpars(1)
        dydt=dpars(4)-dpars(3)
c
        dxdt2=0.0d0
        dydt2=0.0d0
c
        return
        end
c
