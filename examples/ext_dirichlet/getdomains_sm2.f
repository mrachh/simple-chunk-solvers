c
      subroutine domainflagsm(k,nch,n,srccoefs,srcinfo,
     1      iregionl,iregionr,nreg,ntarget,xtarg,ytarg,idt) 
c
c     evaluate SLP and/or DLP at off surface targets.
c   
c     INPUT:
c 
c     k:          nodes per chunk
c     nch:        number of chunks
c     n:          k*nch
c     srccoeffs:  chunk data struct
c     srcinfo:    chunk data struct
c     sigma:      slp density
c     mu :        dlp density
c     ntarget     number of target points 
c     xtarg,ytarg target points 
c
c     OUTPUT:
c
c     upotall  potential at target points
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real *8 (a-h,o-z)
      integer k,nch,n,nreg,ifin
      integer idt(ntarget) 
      integer iregionl(nch) 
      integer iregionr(nch) 
      real *8 xtarg(ntarget)
      real *8 ytarg(ntarget)
      real *8 srcinfo(8,k,nch),srccoefs(6,k,nch)
      real *8 thetain
      real *8, allocatable :: targ2(:,:)
      complex *16, allocatable :: sigma(:,:)
      complex *16, allocatable :: mu(:,:)
      complex *16, allocatable :: pot(:)
c
      allocate(pot(ntarget))
      allocate(targ2(2,ntarget))
      allocate(sigma(k,nch))
      allocate(mu(k,nch))
c
ccc      call prinf(' idt *',idt,ntarget)
ccc      call prinf(' nreg *',nreg,1)
ccc      call prinf(' nch *',nch,1)
      do i = 1,ntarget
         targ2(1,i) = xtarg(i)
         targ2(2,i) = ytarg(i)
         idt(i) = 0
      enddo
c
      do ich=1,nch
        do j=1,k
           sigma(j,ich) = 0.0d0
           mu(j,ich) = 0.0d0
        enddo
      enddo
ccc      goto 333
      do ich=1,nch
         do j=1,k
            mu(j,ich) = mu(j,ich) + 1.0d0*iregionl(ich)
         enddo
         do j=1,k
            mu(j,ich) = mu(j,ich) -1.0d0*iregionr(ich)
         enddo
      enddo
ccc         call prin2(' mu is *',mu,2*k*nch)
c
ccc         call prinf(' nch *',nch,1)
ccc         call prinf(' k *',k,1)
ccc         call prinf(' ntarget *',ntarget,1)
ccc         call prinf(' n *',n,1)
ccc         call prin2(' targ2 *',targ2,2*ntarget)
ccc         call helm2d_dielec_targ(k,nch,n,srccoefs,srcinfo,
ccc     1      zpars,sigma,mu,ntarget,targ2,pot) 
ccc         t1 = second()
      call lap2d_dlp_targ(k,nch,n,srccoefs,srcinfo,
     1      mu,ntarget,targ2,pot) 
ccc         t2 = second()
ccc         call prin2(' time for lap2d targ = *',t2-t1,1)
ccc         call prinf(' nts = *',nts,1)
ccc         call prin2(' pot = *',pot,2*ntarget)
      do i = 1,ntarget
         idt(i) = idt(i) + nint(abs(pot(i)))
      enddo
      return
      end
c
