cc Copyright (C) 2018-2019: Leslie Greengard, Zydrunas Gimbutas, 
cc and Manas Rachh
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$

      subroutine lfmm2dmain(nd,eps,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,
     $     ntarget,targetsort,nexpc,expcsort,
     $     iaddr,rmlexp,mptemp,lmptmp,
     $     itree,ipointer,ndiv,nlevels, 
     $     nboxes,boxsize,rscales,centers,laddr,mnlist1,mnlist2,
     $     mnlist3,mnlist4,nterms,ntj,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg,
     $     jsort,scjsort)
c   Helmholtz FMM in R^2: evaluate all pairwise particle
c   interactions (ignoring self-interaction) 
c   and interactions with targets.
c
c   We use log for the Green's function.
c   Self-interactions are not included
c
c   l2d: charge and dipstr are complex valued, x in \R^2
c
c   \phi(x_i) = \sum_{j\ne i} charge_j log(x_i-x_j)
c   + dipstr_j/(x_i - x_j)
c
c   This subroutine lfmm2dgqbx (helmholtz fmm implementation
c   of global qbx)
c   returns the local expansions of order ntj due to the
c   collection of sources given by the user at the expansion centers
c   and evaluates the potential and gradients
c   at the targets
c
c   Dependencies
c   helmrouts2d.f
c   hank103cc.f
c   l2dterms.f
c   cdjseval2d.f
c   levrtree2d.f (Only required for documentation purposes)
c
c   All the source/target/expansion center related quantities
c   are assumed to be tree-sorted
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:   number of charge densities
c
c   eps:  FMM precision requested
c
c   nsource:     integer:  number of sources
c   sourcesort: real *8 (2,ns):  source locations
c
c   ifcharge:  charge computation flag
c              ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c   chargesort: complex *16 (nsource): charge strengths
c
c   ifdipole:  dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c   dipstrsort: complex *16 (nsource): dipole strengths
c   ntarget: integer:  number of targets
c   targetsort: real *8 (2,ntarget):  target locations
c   nexpc: number of expansion centers
c   expcsort: real *8 (2,nexpc): expansion center locations
c   iaddr: (2,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local expansion of ibox
c  mptemp: (lmptmp): temporary multipole/local expansion
c                        (may not be needed in new setting)
c  lmptmp: length of temporary expansion
c   
c
c
c   itree    in: integer (ltree)
c             This array contains all the information
c             about the tree
c             Refer to levrtree2d.f
c
c   ltree    in: integer
c            length of tree
c
c    ipointer in: integer(30)
c             ipointer is a collection of pointers 
c             which points to where different elements 
c             of the tree are stored in the itree array
c
c             itree(ipointer(1):ipointer(2)-1) == laddr(2,0:nlevels)
c             laddr(1,i) is the first box on level i
c             laddr(2,i) is the last box on level i
c
c             itree(ipointer(2):ipointer(3)-1) == iparent(nboxes)
c             iparent(i) is the parent of box i
c
c             itree(ipointer(3):ipointer(4)-1) == nchild(nboxes)
c             nchild(i) is the number of children of box i
c
c             itree(ipointer(4):ipointer(5)-1) == child(4,nboxes)
c             child(i,j) is the ith child of box j
c
c             itree(ipointer(5):ipointer(6)-1) == isource(ns)
c             isource is the mapping to tree sort source chunks
c
c             itree(ipointer(6):ipointer(7)-1) == itarget(nt)
c             itarget is the mapping to tree sort targets
c
c             itree(ipointer(7):ipointer(8)-1) == iexpc(nexpc)
c             iexpc is the mapping to tree sort expansion centers
c
c             itree(ipointer(8):ipointer(9)-1) == ihsfirst(nboxes)
c             ihsfirst(i) is the location in isource of the first
c             hung chunk in box i
c 
c             itree(ipointer(9):ipointer(10)-1) == ihslast(nboxes)
c             ihslast(i) is the location in isource of the last hung
c             chunk in box i
c
c             itree(ipointer(10):ipointer(11)-1) == isfirst(nboxes)
c             isfirst(i) is the location in isource of the first
c             source in box i
c 
c             itree(ipointer(11):ipointer(12)-1) == islast(nboxes)
c             islast(i) is the location in isource of the last
c             source in box i
c
c             itree(ipointer(12):ipointer(13)-1) == itfirst(nboxes)
c             itfirst(i) is the location in itarget of the first
c             target in box i
c 
c             itree(ipointer(13):ipointer(14)-1) == itlast(nboxes)
c             itlast(i) is the location in itarget of the last
c             target in box i
c
c             itree(ipointer(14):ipointer(15)-1) == ihefirst(nboxes)
c             ihsfirst(i) is the location in iexpc of the first
c             hung expansion center in box i
c 
c             itree(ipointer(15):ipointer(16)-1) == ihelast(nboxes)
c             ihslast(i) is the location in iexpc of the last hung
c             expansion center in box i
c
c             itree(ipointer(16):ipointer(17)-1) == iefirst(nboxes)
c             iefirst(i) is the location in iexpc of the first
c             expansion center in box i
c 
c             itree(ipointer(17):ipointer(18)-1) == ielast(nboxes)
c             ielast(i) is the location in iexpc of the last
c             expansion center in box i
c
c             itree(ipointer(18):ipointer(19)-1) == nlist1(nboxes)
c             nlist1(i) is the number of boxes in list 1 of
c             a box. Max number of boxes in list 1 of a box is 13
c             (mnlist1 = 13)
c
c             itree(ipointer(19):ipointer(20)-1) == list1(mnlist1,
c                                                   nboxes)
c             list1(j,i) is the id of the jth box in list1 of
c             box i
c
c             itree(ipointer(20):ipointer(21)-1) == nlist2(nboxes)
c             nlist2(i) is the number of boxes in list 2 of
c             a box. Max number of boxes in list 2 of a box is 27
c             (mnlist2 = 27)
c
c             itree(ipointer(21):ipointer(22)-1) == list2(mnlist2,
c                                                   nboxes)
c             list2(j,i) is the id of the jth box in list 2 of
c             box i
c
c             itree(ipointer(22):ipointer(23)-1) == nlist3(nboxes)
c             nlist3(i) is the number of boxes in list 3 of
c             a box. Max number of boxes in list 3 of a box is 20
c             (mnlist3 = 20)
c
c             itree(ipointer(23):ipointer(24)-1) == list3(mnlist3,
c                                                   nboxes)
c             list3(j,i) is the id of the jth box in list3 of
c             box i
c
c             itree(ipointer(24):ipointer(25)-1) == nlist4(nboxes)
c             nlist4(i) is the number of boxes in list 4 of
c             a box. Max number of boxes in list 4 of a box is 5
c             (mnlist4 = 5)
c
c             itree(ipointer(25):ipointer(26)-1) == list4(mnlist4,
c                                                   nboxes)
c             list4(j,i) is the id of the jth box in list4 of
c             box i
c
c             itree(ipointer(26):ipointer(27)-1) == nhungsrc(nboxes)
c             nhungsrc(i) is the number of hung chunks in box i
c
c             itree(ipointer(27):ipointer(28)-1) == nhungexp(nboxes)
c             nhungexp(i) is the number of hung expansion centers 
c             in box i
c 
c             itree(ipointer(28):ipointer(29)-1) == nhunglistsrc(nboxes)
c             nhunglistsrc(i) is the number of hung sources relevant
c             to box i
c 
c             itree(ipointer(29):ipointer(30)-1) == 
c             ihunglistsrc(mhung,nboxes)
c             ihunglistsrc(i,j) = src id of ith hung source relevant
c             to box j
c
c     ndiv    in: integer
c             Max number of chunks per box
c
c     nlevels in: integer
c             number of levels in the tree
c
c     
c     nboxes  in: integer
c             number of boxes in the tree
c
c     boxsize in: real*8 (0:nlevels)
c             boxsize(i) is the size of the box from end to end
c             at level i
c
c     centers in: real *8(2,nboxes)
c                 array containing the centers of all the boxes
c
c     mnlist1 in: integer
c             max number of boxes in list1 of a given box
c
c     mnlist2 in: integer
c             max number of boxes in list2 of a given box
c
c     mnlist3 in: integer
c             max number of boxes in list3 of a given box
c
c     mnlist4 in: integer
c             max number of boxes in list4 of a given box
c
c     nterms: (0:nlevels) length of multipole and local expansions
c              at various levels
c     ntj     in: integer
c             order of the output expansions
c
c     ifpgh  in: integer
c             flag for evaluating potential/gradients/hessians 
c             at sources.
c             ifpgh = 1, only potentials will be evaluated
c             ifpgh = 2, potentials/gradients will be evaluated
c             ifpgh = 3, potentials/gradients/hessians will be evaluated
c
c     ifpghtarg  in: integer
c             flag for evaluating potential/gradients/hessians 
c             at targets.
c             ifpghtarg = 1, only potentials will be evaluated
c             ifpghtarg = 2, potentials/gradients will be evaluated
c             ifpghtarg = 3, potentials/gradients/hessians will be evaluated
c
c     flags      in: integer(nt)
c                potential/gradients/ at target i
c                will be computed using point fmm if flags(i).lt.0
c
c     nadd       in: integer
c                number of terms to be added to fmm p
c
c   OUTPUT
c
c   Expansions at the targets
c   jexps : coeffs for local expansion
c   scj: scaling parameter for the expansions
c
c   pot: potential at the source locations
c   grad: gradient at the source locations
c   hess: gradient at the source locations
c  
c   pottarg: potential at the target locations
c   gradtarg: gradient at the target locations
c   hesstarg: gradient at the target locations
c------------------------------------------------------------------

      implicit none

      integer nd

      integer nsource,ntarget,nexpc
      integer ndiv,nlevels,ntj

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 eps

      real *8 sourcesort(2,nsource)

      complex *16 chargesort(nd,*)
      complex *16 dipstrsort(nd,*)

      real *8 targetsort(2,ntarget)
      complex *16 jsort(nd,0:ntj,*)

      real *8 expcsort(2,*)

      complex *16 pot(nd,*)
      complex *16 grad(nd,*)
      complex *16 hess(nd,*)

      complex *16 pottarg(nd,*)
      complex *16 gradtarg(nd,*)
      complex *16 hesstarg(nd,*)

      integer iaddr(2,nboxes),lmptmp
      real *8 rmlexp(*)
      complex *16 mptemp(lmptmp)
       
      real *8 timeinfo(10)
      real *8 timelev(0:200)
      real *8 centers(2,*)

      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer ipointer(30)
      integer itree(*)
      integer nboxes
      real *8 rscales(0:nlevels),boxsize(0:nlevels)

      real *8 scjsort(*)

      real *8 thresh

      integer nterms_eval(4,0:200)

c     temp variables
      integer i,j,k,l,idim
      integer ibox,jbox,ilev,npts
      integer nchild,nlist1,nlist2,nlist3,nlist4
      integer mnlist1,mnlist2,mnlist3,mnlist4

      integer istart,iend,istarts,iends
      integer isstart,isend,jsstart,jsend
      integer jstart,jend
      integer istarte,iende,istartt,iendt

      integer ifprint

      integer ifhesstarg,nn
      real *8 d,time1,time2,omp_get_wtime
      real *8 tt1,tt2
      complex *16 pottmp,gradtmp,hesstmp
      
      real *8, allocatable :: carray(:,:)
      integer ldc


      
      double precision dlam, pi, boxlam
      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c      
        ifprint=0

        pi = 4*atan(1.0d0)
c

        do i=0,nlevels
          timelev(i) = 0
        enddo

        ldc = 100
        allocate(carray(0:ldc,0:ldc))

        call l2d_init_carray(carray,ldc)

cc        call prinf('iaddr=*',iaddr,2*nboxes)
cc        call prinf('nd = *',nd,1)
cc        call prinf('ifcharge=*',ifcharge,1)
cc        call prinf('ifdipole=*',ifdipole,1)
cc        call prin2('rscales=*',rscales,nlevels+1)
        
        



c
c
c     ... set the expansion coefficients to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim,i,j)
      do i=1,nexpc
         do j = 0,ntj
           do idim=1,nd
             jsort(idim,j,i)=0
           enddo
         enddo
      enddo
C$OMP END PARALLEL DO
C
c       
        do i=1,10
          timeinfo(i)=0
        enddo
c
c       ... set all multipole and local expansions to zero
c
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            call l2dmpzero_vec(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
            call l2dmpzero_vec(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
         enddo
C$OMP END PARALLEL DO         
       enddo

c     Set scjsort
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
         do ibox = laddr(1,ilev), laddr(2,ilev)
            nchild = itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then
                istart = itree(ipointer(16)+ibox-1)
                iend = itree(ipointer(17)+ibox-1)
                do i=istart,iend
                   scjsort(i) = rscales(ilev)
                enddo
            endif
         enddo
C$OMP END PARALLEL DO         
      enddo
       
c
c
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
        call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions

      do ilev = 2,nlevels
C
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             nchild = itree(ipointer(3)+ibox-1)
             istart = itree(ipointer(10)+ibox-1)
             iend = itree(ipointer(11)+ibox-1)
             npts = iend-istart+1
c              Check if current box is a leaf box            
             if(nchild.eq.0.and.npts.gt.0) then
                 call l2dformmpc_vec(nd,rscales(ilev),
     1             sourcesort(1,istart),npts,chargesort(1,istart),
     2             centers(1,ibox),nterms(ilev),
     3             rmlexp(iaddr(1,ibox)))
             endif
          enddo
C$OMP END PARALLEL DO 
        endif

        if(ifdipole.eq.1.and.ifcharge.eq.0) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             nchild = itree(ipointer(3)+ibox-1)
             istart = itree(ipointer(10)+ibox-1)
             iend = itree(ipointer(11)+ibox-1)
             npts = iend-istart+1
c              Check if current box is a leaf box            
             if(nchild.eq.0.and.npts.gt.0) then
                call l2dformmpd_vec(nd,rscales(ilev),
     1          sourcesort(1,istart),npts,dipstrsort(1,istart),
     2          centers(1,ibox),
     3          nterms(ilev),rmlexp(iaddr(1,ibox))) 
             endif
          enddo
C$OMP END PARALLEL DO 
        endif

        if(ifdipole.eq.1.and.ifcharge.eq.1) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             nchild = itree(ipointer(3)+ibox-1)
             istart = itree(ipointer(10)+ibox-1)
             iend = itree(ipointer(11)+ibox-1)
             npts = iend-istart+1
c             Check if current box is a leaf box            
             if(nchild.eq.0.and.npts.gt.0) then
                call l2dformmpcd_vec(nd,rscales(ilev),
     1             sourcesort(1,istart),npts,chargesort(1,istart),
     2             dipstrsort(1,istart),
     3             centers(1,ibox),
     4             nterms(ilev),rmlexp(iaddr(1,ibox))) 
             endif
          enddo
C$OMP END PARALLEL DO 
        endif
      enddo


      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1

      if(ifprint.ge.1)
     $      call prinf('=== STEP 2 (form lo) ====*',i,0)
      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nlist4 = itree(ipointer(24)+ibox-1)
            do i=1,nlist4
              jbox = itree(ipointer(25)+(ibox-1)*mnlist4+i-1)
              istart = itree(ipointer(10)+jbox-1)
              iend = itree(ipointer(11)+jbox-1)
              npts = iend-istart+1

              call l2dformtac_vec(nd,rscales(ilev),
     1            sourcesort(1,istart),npts,
     2            chargesort(1,istart),centers(1,ibox),
     3            nterms(ilev),rmlexp(iaddr(2,ibox)))
            enddo
          enddo
C$OMP END PARALLEL DO        
        endif
        if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nlist4 = itree(ipointer(24)+ibox-1)
            do i=1,nlist4
              jbox = itree(ipointer(25)+(ibox-1)*mnlist4+i-1)
              istart = itree(ipointer(10)+jbox-1)
              iend = itree(ipointer(11)+jbox-1)
              npts = iend-istart+1

              call l2dformtad_vec(nd,rscales(ilev),
     1          sourcesort(1,istart),npts,
     2          dipstrsort(1,istart),
     3          centers(1,ibox),nterms(ilev),rmlexp(iaddr(2,ibox)))
            enddo
          enddo
C$OMP END PARALLEL DO        
        endif
        if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nlist4 = itree(ipointer(24)+ibox-1)
            do i=1,nlist4
              jbox = itree(ipointer(25)+(ibox-1)*mnlist4+i-1)
              istart = itree(ipointer(10)+jbox-1)
              iend = itree(ipointer(11)+jbox-1)
              npts = iend-istart+1

              call l2dformtacd_vec(nd,rscales(ilev),
     1          sourcesort(1,istart),npts,
     2          chargesort(1,istart),dipstrsort(1,istart),
     3          centers(1,ibox),
     3          nterms(ilev),rmlexp(iaddr(2,ibox)))
            enddo
          enddo
C$OMP END PARALLEL DO        
        endif
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

cc      print *, ldc
cc      call prin2('carray=*',carray,(ldc+1)*(ldc+1))

      if(ifprint .ge. 1)
     $      call prinf('=== STEP 3 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do ilev=nlevels-1,1,-1

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,istart,iend,npts,mptemp)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(ipointer(3)+ibox-1)
          do i=1,nchild
            jbox = itree(ipointer(4)+4*(ibox-1)+i-1)
            istart = itree(ipointer(10)+jbox-1)
            iend = itree(ipointer(11)+jbox-1)
            npts = iend-istart+1
            if(npts.gt.0) then
              call l2dmpmp_vec(nd,rscales(ilev+1),
     1             centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2             nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3             rmlexp(iaddr(1,ibox)),nterms(ilev),carray,ldc)
            endif
          enddo
        enddo
C$OMP END PARALLEL DO    
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3)=time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (mp to loc) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 2,nlevels

       tt1 = second()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,mptemp,i,nlist2)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          npts = 0
          if(ifpghtarg.gt.0) then
            istart = itree(ipointer(12)+ibox-1)
            iend = itree(ipointer(13)+ibox-1)
            npts = npts + iend-istart+1
          endif

          istart = itree(ipointer(16)+ibox-1)
          iend = itree(ipointer(17)+ibox-1)
          npts = npts + iend-istart+1

          if(ifpgh.gt.0) then
            istart = itree(ipointer(10)+ibox-1)
            iend = itree(ipointer(11)+ibox-1)
            npts = npts + iend-istart+1
          endif

          if(npts.gt.0) then
            nlist2 = itree(ipointer(20)+ibox-1)
            do i=1,nlist2
              jbox = itree(ipointer(21)+mnlist2*(ibox-1)+i-1)
              call l2dmploc_vec(nd,rscales(ilev),
     $          centers(1,jbox),rmlexp(iaddr(1,jbox)),nterms(ilev),
     2          rscales(ilev),centers(1,ibox),rmlexp(iaddr(2,ibox)),
     3          nterms(ilev),carray,ldc)
            enddo
          endif
        enddo
C$OMP END PARALLEL DO        
       tt2 = second()
       timelev(ilev) = tt2-tt1
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(4) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 5 (split loc) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 1,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,istart,iend,npts,mptemp)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(ipointer(3)+ibox-1)
          istart = itree(ipointer(16)+ibox-1)
          iend = itree(ipointer(17)+ibox-1)
          npts = iend - istart + 1


          if(ifpghtarg.gt.0) then
            istart = itree(ipointer(12)+ibox-1)
            iend = itree(ipointer(13)+ibox-1)
            npts = npts + iend-istart+1
          endif

          if(ifpgh.gt.0) then
            istart = itree(ipointer(10)+ibox-1)
            iend = itree(ipointer(11)+ibox-1)
            npts = npts + iend-istart+1
          endif

          if(npts.gt.0) then
            do i=1,nchild
              jbox = itree(ipointer(4)+4*(ibox-1)+i-1)
              call l2dlocloc_vec(nd,rscales(ilev),centers(1,ibox),
     1          rmlexp(iaddr(2,ibox)),nterms(ilev),rscales(ilev+1),
     2          centers(1,jbox),rmlexp(iaddr(2,jbox)),nterms(ilev+1),
     3          carray,ldc)
            enddo
          endif
        enddo
C$OMP END PARALLEL DO        
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(5) = time2-time1

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      if(ifprint.ge.1)
     $    call prinf('=== Step 6 (mp eval) ===*',i,0)

cc      call prinf('ifpgh=*',ifpgh,1)
cc      call prinf('ifpghtarg=*',ifpghtarg,1)
cc      call prinf('laddr=*',laddr,2*(nlevels+1))
      do ilev=1,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,j,i,mptemp)
C$OMP$PRIVATE(jbox)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=laddr(1,ilev),laddr(2,ilev)
          nlist3 = itree(ipointer(22)+ibox-1)

          istart = itree(ipointer(16)+ibox-1)
          iend = itree(ipointer(17)+ibox-1)
          do j=istart,iend
            do i=1,nlist3
              jbox = itree(ipointer(23)+(ibox-1)*mnlist3+i-1)
c                 shift multipole expansion directly to box
c                 for all expansion centers
              call l2dmploc_vec(nd,rscales(ilev+1),
     $          centers(1,jbox),rmlexp(iaddr(1,jbox)),nterms(ilev+1),
     2          scjsort(j),expcsort(1,j),jsort(1,0,j),ntj,carray,ldc)
            enddo
          enddo

c              evalute multipole expansion at all targets
          istart = itree(ipointer(12)+ibox-1)
          iend = itree(ipointer(13)+ibox-1)
          npts = iend-istart+1

          if(ifpghtarg.eq.1) then
            do i=1,nlist3
              jbox = itree(ipointer(23)+(ibox-1)*mnlist3+i-1)
                  
              call l2dmpevalp_vec(nd,rscales(ilev+1),
     1         centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2         nterms(ilev+1),targetsort(1,istart),npts,
     3         pottarg(1,istart))
            enddo
          endif
          if(ifpghtarg.eq.2) then
            do i=1,nlist3
              jbox = itree(ipointer(23)+(ibox-1)*mnlist3+i-1)
              call l2dmpevalg_vec(nd,rscales(ilev+1),
     1          centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2          nterms(ilev+1),targetsort(1,istart),npts,
     3          pottarg(1,istart),gradtarg(1,istart))
            enddo
          endif
          if(ifpghtarg.eq.3) then
            do i=1,nlist3
              jbox = itree(ipointer(23)+(ibox-1)*mnlist3+i-1)

              call l2dmpevalh_vec(nd,rscales(ilev+1),
     1          centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2          nterms(ilev+1),targetsort(1,istart),npts,
     3          pottarg(1,istart),
     3          gradtarg(1,istart),hesstarg(1,istart))
            enddo
          endif


c              evalute multipole expansion at all sources
          istart = itree(ipointer(10)+ibox-1)
          iend = itree(ipointer(11)+ibox-1)
          npts = iend-istart+1
            

          if(ifpgh.eq.1) then
            do i=1,nlist3
              jbox = itree(ipointer(23)+(ibox-1)*mnlist3+i-1)
              call l2dmpevalp_vec(nd,rscales(ilev+1),
     1           centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2           nterms(ilev+1),sourcesort(1,istart),npts,
     3           pot(1,istart))
            enddo
          endif
          if(ifpgh.eq.2) then
            do i=1,nlist3
              jbox = itree(ipointer(23)+(ibox-1)*mnlist3+i-1)
              call l2dmpevalg_vec(nd,rscales(ilev+1),
     1           centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2           nterms(ilev+1),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,istart))
            enddo
          endif
          if(ifpgh.eq.3) then
            do i=1,nlist3
              jbox = itree(ipointer(23)+(ibox-1)*mnlist3+i-1)
              call l2dmpevalh_vec(nd,rscales(ilev+1),
     1           centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2           nterms(ilev+1),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,istart),hess(1,istart))
            enddo
          endif

        enddo
C$OMP END PARALLEL DO     
      enddo

 1000 continue    


      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(6) = time2-time1


      if(ifprint.ge.1)
     $    call prinf('=== step 7 (eval lo) ===*',i,0)

c     ... step 7, evaluate all local expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,mptemp,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(ipointer(3)+ibox-1)
          if(nchild.eq.0) then
            istart = itree(ipointer(16)+ibox-1)
            iend = itree(ipointer(17)+ibox-1)
            do i=istart,iend
              call l2dlocloc_vec(nd,rscales(ilev),
     $          centers(1,ibox),
     1          rmlexp(iaddr(2,ibox)),nterms(ilev),scjsort(i),
     2          expcsort(1,i),jsort(1,0,i),ntj,carray,ldc)
            enddo
c
cc               evaluate local expansion
c                at targets
            istart = itree(ipointer(12)+ibox-1)
            iend = itree(ipointer(13)+ibox-1)
            npts = iend-istart + 1
            if(ifpghtarg.eq.1) then
              call l2dtaevalp_vec(nd,rscales(ilev),
     1              centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2              nterms(ilev),targetsort(1,istart),npts,
     3              pottarg(1,istart))
            endif
            if(ifpghtarg.eq.2) then
              call l2dtaevalg_vec(nd,rscales(ilev),
     1          centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2          nterms(ilev),targetsort(1,istart),npts,
     3          pottarg(1,istart),gradtarg(1,istart))
            endif
            if(ifpghtarg.eq.3) then
              call l2dtaevalh_vec(nd,rscales(ilev),
     1          centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2          nterms(ilev),targetsort(1,istart),npts,
     3          pottarg(1,istart),gradtarg(1,istart),
     4          hesstarg(1,istart))
            endif

c
cc                evaluate local expansion at sources

            istart = itree(ipointer(10)+ibox-1)
            iend = itree(ipointer(11)+ibox-1)
            npts = iend-istart+1
            if(ifpgh.eq.1) then
              call l2dtaevalp_vec(nd,rscales(ilev),
     1           centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2           nterms(ilev),sourcesort(1,istart),npts,
     3           pot(1,istart))
            endif
            if(ifpgh.eq.2) then
              call l2dtaevalg_vec(nd,rscales(ilev),
     1           centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2           nterms(ilev),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,istart))
            endif
            if(ifpgh.eq.3) then
              call l2dtaevalh_vec(nd,rscales(ilev),
     1           centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2           nterms(ilev),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,istart),hess(1,istart))
            endif
          endif
        enddo
C$OMP END PARALLEL DO        
      enddo

      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(7) = time2 - time1

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (direct) =====*',i,0)

c
cc     set threshold for ignoring interactions with 
c      |r| < thresh
c
      thresh = boxsize(0)*1.0d-14

cc      call prin2('thresh=*',thresh,1)
c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime()  
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istartt,iendt,i,jstart,jend,istarte,iende)
C$OMP$PRIVATE(nlist1,istarts,iends)
C$OMP$SCHEDULE(DYNAMIC)  
         do ibox = laddr(1,ilev),laddr(2,ilev)

            istartt = itree(ipointer(12)+ibox-1)
            iendt = itree(ipointer(13)+ibox-1)


            istarte = itree(ipointer(16)+ibox-1)
            iende = itree(ipointer(17)+ibox-1)

            istarts = itree(ipointer(10)+ibox-1)
            iends = itree(ipointer(11)+ibox-1)


            nlist1 = itree(ipointer(18)+ibox-1)
            do i =1,nlist1
               jbox = itree(ipointer(19)+mnlist1*(ibox-1)+i-1)

               jstart = itree(ipointer(10)+jbox-1)
               jend = itree(ipointer(11)+jbox-1)

               call lfmm2dexpc_direct_vec(nd,jstart,jend,istarte,
     1         iende,rscales,nlevels, 
     2         sourcesort,ifcharge,chargesort,ifdipole,dipstrsort,
     3         expcsort,jsort,scjsort,ntj)

                
               call lfmm2dpart_direct_vec(nd,jstart,jend,istartt,
     1         iendt,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipstrsort,targetsort,ifpghtarg,pottarg,
     3         gradtarg,hesstarg,thresh)
         
               call lfmm2dpart_direct_vec(nd,jstart,jend,istarts,iends,
     1         sourcesort,ifcharge,chargesort,ifdipole,
     2         dipstrsort,sourcesort,ifpgh,pot,grad,hess,
     3         thresh)
            enddo   
         enddo
C$OMP END PARALLEL DO         
      enddo

      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(8) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,8)
      d = 0
      do i = 1,8
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)
      if(ifprint.ge.1) call prin2('timlev=*',timelev,nlevels+1)

      return
      end
c
      subroutine lfmm2dexpc_direct_vec(nd,istart,iend,jstart,jend,
     $     rscales,nlevels,source,ifcharge,charge,ifdipole,dipstr,
     $     targ,jexps,scj,ntj)
c--------------------------------------------------------------------
c     This subroutine adds the local expansions due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the existing local
c     expansions
c
c     INPUT arguments
c-------------------------------------------------------------------
c     nd           in: integer
c                   number of expansions
c
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in target array at which we
c                  wish to compute the expansions
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to compute the expansions
c 
c     rscales       in: real*8(0:nlevels)
c                  Scale of expansions formed at all levels
c
c     nlevels      in:Integer
c                  Number of levels in the tree structure
c
c     source       in: real *8(2,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: complex *16
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr        in: complex *16(ns)
c                   dip strengths at the source locations
c
c     targ        in: real *8(2,nexpc)
c                 Expansion center locations
c
c     scj         in: real *8(nexpc)
c                 scaling parameter for expansions
c
c     ntj         in: Integer
c                 Number of terms in expansion
c------------------------------------------------------------
c     OUTPUT
c
c   Updated expansions at the targets
c   jexps : coeffs for local expansions
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j
        integer ifcharge,ifdipole,ier,nd
        real *8 source(2,*)
        real *8 rscales(0:nlevels)
        complex *16 charge(nd,*),dipstr(nd,*)
        real *8 targ(2,*)
        real *8 scj(*)

        integer nlevels,ntj
c
        complex *16 jexps(nd,0:ntj,*)
        
c
        ns = iend - istart + 1
        do j=jstart,jend
           if(ifcharge.eq.1.and.ifdipole.eq.0) then
              call l2dformtac_vec(nd,scj(j),
     1        source(1,istart),charge(1,istart),ns,targ(1,j),
     2        ntj,jexps(1,0,j))
           endif

           if(ifdipole.eq.1.and.ifcharge.eq.0) then
               call l2dformtad_vec(nd,scj(j),
     1         source(1,istart),dipstr(1,istart),
     2         ns,targ(1,j),ntj,jexps(1,0,j))
           endif        
           if(ifdipole.eq.1.and.ifcharge.eq.1) then
               call l2dformtacd_vec(nd,scj(j),
     1         source(1,istart),charge(1,istart),dipstr(1,istart),
     2         ns,targ(1,j),ntj,jexps(1,0,j))
           endif        
        enddo
c
        return
        end
c------------------------------------------------------------------     
      subroutine lfmm2dpart_direct_vec(nd,istart,iend,jstart,jend,
     $     source,ifcharge,charge,ifdipole,dipstr,
     $     targ,ifpgh,pot,grad,hess,thresh)
c--------------------------------------------------------------------
c     This subroutine adds the contribuition due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the computed velocities
c     and gradients. Note that contributions for sources
c     within thresh of the targets are not added to the potential
c     
c
c     INPUT arguments
c-------------------------------------------------------------------
c     nd           in: integer
c                  number of charge densities
c
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in target array at which we
c                  wish to update the potential and gradients
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to update the potential and gradients
c
c     source       in: real *8(2,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: complex *16
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr        in: complex *16(ns)
c                 dipole strengths at the source locations
c
c     targ        in: real *8(2,nt)
c                 target locations
c
c     ifpgh        in: Integer
c                  Flag for computing the potential/gradient/hessian.
c                  ifpgh = 1, only potential is computed
c                  ifpgh = 2, potential/gradient are computed
c                  ifpgh = 3, potential/gradient/hessian are computed
c
c     thresh       in: real *8
c                  threshold for computing interactions
c                  if |r| < threshold, then interactions are
c                  not included
c
c
c------------------------------------------------------------
c     OUTPUT
c
c   Updated velocity and gradients at the targets
c   pot : potential at the targets
c   grad: gradient at the targets
c   hess: Hessian at the targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole

        integer nd



        real *8 source(2,*)
        complex *16 charge(nd,*),dipstr(nd,*)

        integer ifpgh
        real *8 targ(2,*),thresh
        
c
        complex *16 pot(nd,*)
        complex *16 grad(nd,*)
        complex *16 hess(nd,*)

c
        ns = iend - istart + 1
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          if(ifpgh.eq.1) then
             do j=jstart,jend
               call l2d_directcp_vec(nd,source(1,istart),ns,
     1            charge(1,istart),targ(1,j),pot(1,j),thresh)
             enddo
          endif

          if(ifpgh.eq.2) then
             do j=jstart,jend
               call l2d_directcg_vec(nd,source(1,istart),ns,
     1            charge(1,istart),targ(1,j),pot(1,j),grad(1,j),
     2            thresh)
             enddo
          endif
          if(ifpgh.eq.3) then
             do j=jstart,jend
               call l2d_directch_vec(nd,source(1,istart),ns,
     1            charge(1,istart),targ(1,j),pot(1,j),grad(1,j),
     2            hess(1,j),thresh)
             enddo
          endif
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             do j=jstart,jend
               call l2d_directdp_vec(nd,source(1,istart),ns,
     1            dipstr(1,istart),
     2            targ(1,j),pot(1,j),thresh)
             enddo
          endif

          if(ifpgh.eq.2) then
             do j=jstart,jend
               call l2d_directdg_vec(nd,source(1,istart),ns,
     1            dipstr(1,istart),
     2            targ(1,j),pot(1,j),grad(1,j),
     2            thresh)
             enddo
          endif
          if(ifpgh.eq.3) then
             do j=jstart,jend
               call l2d_directdh_vec(nd,source(1,istart),ns,
     1            dipstr(1,istart),targ(1,j),
     2            pot(1,j),grad(1,j),
     2            hess(1,j),thresh)
             enddo
          endif
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             do j=jstart,jend
               call l2d_directcdp_vec(nd,source(1,istart),ns,
     1            charge(1,istart),dipstr(1,istart),
     2            targ(1,j),pot(1,j),thresh)
             enddo
          endif

          if(ifpgh.eq.2) then
             do j=jstart,jend
               call l2d_directcdg_vec(nd,source(1,istart),ns,
     1            charge(1,istart),dipstr(1,istart),
     2            targ(1,j),pot(1,j),grad(1,j),
     2            thresh)
             enddo
          endif
          if(ifpgh.eq.3) then
             do j=jstart,jend
               call l2d_directcdh_vec(nd,source(1,istart),ns,
     1            charge(1,istart),dipstr(1,istart),
     2            targ(1,j),pot(1,j),grad(1,j),
     2            hess(1,j),thresh)
             enddo
          endif
        endif


c
        return
        end
c------------------------------------------------------------------    
      subroutine l2dmpalloc(nd,laddr,iaddr,nlevels,lmptot,
     1                          nterms)
c     This subroutine determines the size of the array
c     to be allocated for the multipole expansions
c     iaddr(1,i) points to the starting location of the multipole
c     expansion of box i and iaddr(2,i) points to the local
c     expansion of box i
c  
c     Input arguments
c     nd          in: integer
c                 number of expansions
c
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array providing access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     nterms      in: Integer(0:nlevels)
c                 Number of terms requried in expansions at each
c                 level
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr       out: Integer(2,nboxes)
c                 Points the multipole and local expansions in box i
c 
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,nterms(0:nlevels),nd,nsig,nt1,nt2,next235
      integer iaddr(2,*), lmptot, laddr(2,0:nlevels)
      integer ibox,i,iptr,istart,nn,itmp
      real *8 ddn
c
      istart = 1
      do i = 0,nlevels

         nn = (nterms(i)+1)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(1,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
       enddo
c
c            Allocate memory for the local expansion
c
       do i=0,nlevels
         nn = (nterms(i)+1)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
             itmp = ibox - laddr(1,i)
             iaddr(2,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
      lmptot = istart

      return
      end
c----------------------------------------------------------------     
