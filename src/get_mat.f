      subroutine zgetmat_bdry(k,nch,n,srcinfo,srccoefs,
     1    fker,ndt,ndd,dpars,ndz,zpars,ndi,ipars,xmat)
c
c   this subroutine generates the matrix corresponding to the
c   discretization of a weakly singular integral operator
c   whose kernel is given by fker
c
c   The expected calling sequence of fker is
c      fker(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c
c   where 
c     - srcinfo: double precision(8) 
c         is the source info
c     - ndt: integer
c         <=8 is the expected
c       number of parameters to be used in the target info,
c     - targinfo: double precision(ndt) 
c         is the target info
c     - ndd: integer
c         number of double precision parameters
c     - dpars: double precision(ndd) 
c         the double precision parameters 
c     - ndz: integer
c         number of complex parameters
c     - zpars: double precision(ndz)
c         complex *16 parameters
c     - ndi: integer
c         number of integer parameters
c     - ipars: integer(ndi) 
c         integer parameters 
c
c   other input arguments:
c     - k: integer
c         order of discretization
c     - nch: integer
c         number of chunks
c     - n: integer
c        total number of points on the boundary = k*nch,
c        also leading dimension of output matrix
c     - srcinfo: real *8 (8,k,nch)
c         source information
c         * srcinfo(1:2,j,ich) - x,y info at jth node on ichth chunk
c         * srcinfo(3:4,j,ich) - dxdt,dydt info at jth node on ichth
c           chunk
c         * srcinfo(5:6,j,ich) - d2xdt2,d2ydt2 info at jth node on ichth
c           chunk
c         * srcinfo(7:8,j,ich) - normal information at jth node on ichth
c           chunk
c     - srccoefs: real *8 (6,k,nch)
c         Legendre expansion coeffs of the first 6 arrays in srcinfo
c         array
c   output parameters:
c     - xmat: complex *16(n,n)
c         discretized matrix corresponding to fker
c
c  This subroutine uses generalized gaussian quadratures
c  for self, and adaptive integration for the rest of the chunks
c
c-----

      implicit real *8 (a-h,o-z)
      integer, intent(in) :: k,nch,n
      real *8, intent(in) :: srcinfo(8,k,nch),srccoefs(6,k,nch)
      integer, intent(in) :: ndt,ndd,ndi,ndz
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: ipars(ndi)

      complex *16, intent(out) :: xmat(n,n)
      external fker

c
c    breemer \self quad variables
c
      integer nquad0,nquad1
      real *8, allocatable :: xs1(:),whts1(:),xs0(:,:),whts0(:,:)
      real *8, allocatable :: ainterp(:,:,:)
      real *8, allocatable :: ts0(:),w0(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: tsquad(:),wquad(:)
      real *8, allocatable :: tstmp(:)
      real *8, allocatable :: work(:)
      complex *16, allocatable :: xtmp(:,:),xtmp2(:)

c
c       temporary variables 
c
      integer i,j,l
      integer ich,jch,isch,itch
      integer ismat,itmat
      integer lwork

      call get_quads_info(k, nquad1, nquad0)
      allocate(xs1(nquad1), whts1(nquad1))
      allocate(xs0(k,k), whts0(k,k))
      call get_quads(k, nquad1, xs1, whts1, nquad0, xs0, whts0)

      itype = 2
      allocate(ts0(k),w0(k),umat(k,k),vmat(k,k))
      call legeexps(itype,k,ts0,umat,vmat,w0)

      allocate(tstmp(k+10))
      
      allocate(xtmp(k,k),xtmp2(k))
      allocate(ainterp(k,k,k))

      lwork = k*k*10+10
      allocate(work(lwork))
      do i=1,k
        call lematrin(k,k,xs0(1,i),ainterp(1,1,i),tstmp,work)
      enddo

      m = 20
      itype = 1
      allocate(tsquad(m),wquad(m))
      call legeexps(itype,m,tsquad,utmp,vtmp,wquad) 

      do isch=1,nch
        ismat = (isch-1)*k + 1
        do itch=1,nch
          itmat = (itch-1)*k + 1
          do i=1,k
            do j=1,k
              xtmp(j,i) = 0
            enddo
          enddo


          if(isch.eq.itch) then
            call zself_buildmat(k,srccoefs(1,1,isch),
     1         srcinfo(1,1,itch),fker,ndt,ndd,dpars,ndz,zpars,
     2         ndi,ipars,xs0,whts0,ainterp,xtmp)

          else
            do j=1,k
              xtmp2 = 0
              call zadapquad(k,srccoefs(1,1,isch),srcinfo(1,j,itch),
     1          fker,8,ndd,dpars,ndz,zpars,ndi,ipars,umat,m,tsquad,
     2          wquad,xtmp2)
              do i=1,k
                xtmp(j,i) = xtmp2(i)
              enddo
            enddo
          endif
          call zinsertmat(k,k,xtmp,itmat,ismat,n,n,xmat)
        enddo
      enddo

      return
      end
c
c
c
c
c
      subroutine zself_buildmat(k,srccoefs,srcvals,fker,ndt,ndd,
     1   dpars,ndz,zpars,ndi,ipars,xs0,whts0,ainterp,xmat)
      implicit real *8 (a-h,o-z)
      integer k
      real *8 srccoefs(6,k),srcvals(8,k)
      integer ndt,ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 xs0(k,k),whts0(k,k)
      complex *16 xmat(k,k)
      real *8 ainterp(k,k,k)
      complex *16 wtmp2(k),wtmp3(k)
      complex *16 val
      real *8 pols(k),srcinfo(8)

      real *8 tsrc
      integer i,j,l
      real *8 alpha,beta
      external fker

      alpha = 1.0d0
      beta = 0.0d0

      do i=1,k
        do j=1,k
          tsrc = xs0(j,i)
          call legepols(tsrc,k-1,pols)
          call dgemv('n',6,k,alpha,srccoefs,6,pols,1,beta,srcinfo,1)
          dst = sqrt(srcinfo(3)**2 + srcinfo(4)**2)
          srcinfo(7) = srcinfo(4)/dst
          srcinfo(8) = -srcinfo(3)/dst
          call fker(srcinfo,8,srcvals(1,i),ndd,dpars,ndz,zpars,ndi,
     1      ipars,wtmp2(j))
          wtmp2(j) = wtmp2(j)*whts0(j,i)*dst
        enddo

        call dgemm('n','n',2,k,k,alpha,wtmp2,2,ainterp(1,1,i),k,beta,
     1    wtmp3,2) 
        do j=1,k
          xmat(i,j) = wtmp3(j)
        enddo
      enddo
      
      return
      end
c
c
c
c
c
      subroutine zadapquad(k,srccoefs,targ,fker,ndt,ndd,dpars,ndz,
     1   zpars,ndi,ipars,umat,m,ts0,w0,xmat)
      implicit real *8 (a-h,o-z)
      real *8 srccoefs(6,k),targ(ndt)
      integer ndt,ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 umat(k,k),ts0(m),w0(m)
      complex *16 xmat(k)
      complex *16 rints(k)
      real *8, allocatable :: stack(:,:)
      complex *16, allocatable :: vals(:,:)
      complex *16 value2(k),value3(k)
      external fker

      a = -1.0d0
      b = 1.0d0

      maxrec = 0
      numit = 0
      
      maxdepth = 50
      allocate(stack(2,maxdepth),vals(k,maxdepth))

      nnmax = 100000
      eps = 1.0d-13
      do j=1,maxdepth
        stack(1,j) = 0
        stack(2,j) = 0
        vals(1:k,j) = 0
      enddo



      call adinrecm(ier,stack,a,b,k,srccoefs,targ,fker,
     1   ndt,ndd,dpars,ndz,zpars,ndi,ipars,m,ts0,w0,vals,nnmax,
     2   eps,rints,maxdepth,maxrec,numit,value2,value3)

c
c      convert rints to matrix entries
c
      call zrmatmatt(1,k,rints,k,umat,xmat)

      return
      end
c
c
c
c

      subroutine adinrecm(ier,stack,a,b,k,srccoefs,targ,fker,
     1   ndt,ndd,dpars,ndz,zpars,ndi,ipars,m,ts0,w0,vals,nnmax,
     2   eps,rints,maxdepth,maxrec,numit,value2,value3)
      implicit real *8 (a-h,o-z)
      real *8 stack(2,maxdepth),a,b
      integer k
      real *8 srccoefs(6,k),targ(ndt),dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 ts0(m),w0(m)
      complex *16 vals(k,maxdepth),rints(k),value2(k),value3(k)
      
      


      ier = 0
      stack(1,1) = a
      stack(2,1) = b

      call oneintmu(a,b,k,srccoefs,targ,fker,ndt,ndd,dpars,ndz,
     1   zpars,ndi,ipars,m,ts0,w0,vals(1,1))
      
c 
c       recursively integrate the thing
c 
      j=1
      do i=1,k
        rints(i)=0
      enddo
c 
      ier=0
      maxrec=0
      do 3000 i=1,nnmax
        numint=i
        if(j .gt. maxrec) maxrec=j
c 
c       subdivide the current subinterval
c 
         c=(stack(1,j)+stack(2,j))/2
         call oneintmu(stack(1,j),c,k,srccoefs,targ,fker,ndt,ndd,
     1      dpars,ndz,zpars,ndi,ipars,m,ts0,w0,value2)
         call oneintmu(c,stack(2,j),k,srccoefs,targ,fker,ndt,ndd,
     1      dpars,ndz,zpars,ndi,ipars,m,ts0,w0,value3)
c 
          dd=0
          do jj=1,k
            ddd=abs(value2(jj)+value3(jj)-vals(jj,j) )
            if(ddd .gt. dd) dd=ddd
          enddo
c 
          ifdone=0
          if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
          if(ifdone  .eq. 0) goto 2000
c 
          do jj=1,k
            rints(jj)=rints(jj)+value2(jj)+value3(jj)
          enddo
          j=j-1
c 
c        if the whole thing has been integrated - return
c 
          if(j .eq. 0) return
          goto 3000
 2000 continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
        vals(1:k,j+1) = value2
        vals(1:k,j) = value3
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
c 
        j=j+1
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if(j .lt. maxdepth) goto 3000
        ier=8
        return
 3000 continue
      ier=16

      return
      end
c
c
c
c
c
      subroutine oneintmu(a,b,k,srccoefs,targ,fker,ndt,ndd,dpars,
     1   ndz,zpars,ndi,ipars,m,t,w,vals)
      implicit real *8 (a-h,o-z)
      real *8 a,b
      real *8 srccoefs(6,k),targ(ndt)
      real *8 srctmp(8),pols(8)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 t(m),w(m)
      complex *16 vals(k),val
      external fker
      
      do i=1,k
        vals(i) = 0
      enddo

      u = (b-a)/2
      v = (b+a)/2
      alpha = 1.0d0
      beta = 0.0d0

      do i=1,m
        tt = u*t(i)+v
        call legepols(tt,k-1,pols) 
        call dgemv('n',6,k,alpha,srccoefs,6,pols,1,beta,srctmp,1)
        dst = sqrt(srctmp(3)**2 + srctmp(4)**2)
        srctmp(7) = srctmp(4)/dst
        srctmp(8) = -srctmp(3)/dst
        val = 0
        call fker(srctmp,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1      ipars,val)
        do j=1,k
          vals(j) = vals(j) + dst*w(i)*val*pols(j) 
        enddo
      enddo

      do i=1,k
        vals(i) = vals(i)*u
      enddo

      return
      end
c
c
c
c
c
      subroutine zrmatmatt(m, n, a, k, b, c)
      implicit double precision (a-h,o-z)
      complex *16 :: a(n,m),c(k,m)
      real *8 :: b(n,k)
      complex *16, allocatable :: bz(:,:)
      character *1 :: transa, transb
      complex *16 :: alpha, beta


      allocate(bz(n,k))
      bz = 0
      call dcopy(n*k,b,1,bz,2)

      transa = 'T'
      transb = 'N'
      alpha = 1
      beta = 0

      call zgemm(transa, transb, k, m, n, alpha, bz, n, a, n,
     1     beta, c, k)

      return
      end


