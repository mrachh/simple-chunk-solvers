c
      subroutine evalpotvol(k,nch,n,srccoefs,srcinfo,zks,alphas,betas,
     1          sigma,mu,nreg,idt,ntarget,xtarg,ytarg,upotall)
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
c     zpars(3):   zk, weight of slp, weight of dlp
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
      integer k,nch,n,nreg
      integer idt(ntarget) 
      real *8 xtarg(ntarget)
      real *8 ytarg(ntarget)
      real *8 srcinfo(8,k,nch),srccoefs(6,k,nch)
      real *8, allocatable :: targ2(:,:)
      complex *16 zks(0:nreg),alphas(0:nreg),betas(0:nreg)
      complex *16 sigma(k,nch)
      complex *16 mu(k,nch)
      complex *16 upotall(ntarget)
      complex *16, allocatable :: pot(:)
      complex *16 zpars(3)
c
      allocate(pot(ntarget))
      allocate(targ2(2,ntarget))
      nts = 0
      zpars(1) = zks(0)
      zpars(2) = -1.0d0/betas(0)
      zpars(3) = 1.0d0/betas(0)
      call prinf(' idt *',idt,ntarget)
      call prinf(' nreg *',nreg,1)
      do i = 1,ntarget
         if (idt(i).eq.0) then
            write(33,*) i
            nts = nts+1
            targ2(1,nts) = xtarg(i)
            targ2(2,nts) = ytarg(i)
         endif
      enddo
      call prin2(' calling fmm 0 *',pot(1),0)

      if (nts.gt.0) call helm2d_dielec_targ(k,nch,n,srccoefs,srcinfo,
     1   zpars,sigma,mu,nts,targ2,pot) 
      call prin2(' after fmm 0 *',pot(1),0)

      nts = 0
      do i = 1,ntarget
         if (idt(i).eq.0) then
            nts = nts+1
            upotall(i) = pot(nts)
c            pot(i) = 0
ccc            upotall(i) = upotall(i) + cdexp(-eye*zk0*targs(2,i))
         endif
      enddo
ccc      call prin2(' pot(1) is *',pot(1),2)
c
ccc      goto 333
      do ii = 1,nreg
ccc         write(34,*) ' component', ii
         nts = 0
         zpars(1) = zks(ii)
         zpars(2) = -1.0d0/betas(ii)
         zpars(3) = 1.0d0/betas(ii)
ccc         call prinf(' ii reg  *',ii,1)
         do i = 1,ntarget
            if (idt(i).eq.ii) then
ccc               write(34,*) i
               nts = nts+1
               targ2(1,nts) = xtarg(i)
               targ2(2,nts) = ytarg(i)
            endif
         enddo
ccc         call prin2(' targ2 *',targ2,2*nts)
         if (nts.gt.0) call helm2d_dielec_targ(k,nch,n,srccoefs,srcinfo,
     1      zpars,sigma,mu,nts,targ2,pot) 
ccc         call prinf(' nts = *',nts,1)
ccc         call prin2(' pot = *',pot,2*nts)
         nts = 0
         do i = 1,ntarget
            if (idt(i).eq.ii) then
ccc               write(34,*) i
               nts = nts+1
               upotall(i) = pot(nts)
            endif
         enddo
      enddo
      return
      end
c
