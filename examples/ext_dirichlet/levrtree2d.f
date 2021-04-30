c
c     d2hplratree.f - hung pruned level restricted adaptive
c                     tree
c     radsrc parameter (preventing sources with large radius 
c                         from being pushed to the finest level)
c     radexp parameter (preventing exp centers with large radius 
c                         from being pushed to the finest level)
c     ndiv parameter (usual refinement parameter extending tree
c                         depth until all leaf nodes hanve fewer than
c                         ndiv particles)
c     idivflag       (allows use of ndiv parameter for sources OR
c                         targets OR sources+targets OR sources+targets
c                         + expansion centers)
c
c     We first create a fully adaptive tree and then 
c     use the fix tree algorithm to make it a level restricted
c     tree. This algorithm is discussed in depth in Frank Ethridge's
c     thesis
c
c     This tree code MUST be used along with its memory code
c     maketreemem located at the end of this file. The memory
c     code precomputes the number of levels, the number of
c     boxes, the max number of hung sources and the length
c     of the tree.
c
      subroutine maketree2d(src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   mhung,nlevels,nboxes,
     $                   centers,boxsize,itree,ltree,
     $                   ipointer,mnlist1,mnlist2,mnlist3,mnlist4)
      implicit none
      integer ns,nt,nexpc,idivflag,ndiv,isep,mhung
      integer nlevels,lcenters,ltree
      integer nlmax,nbmax,nboxes, nlevtmp,nbtmp, mhungtmp
      integer itree(ltree)
      integer i,j
      integer, allocatable :: laddr(:,:)
      integer, allocatable :: ilevel(:)
      integer, allocatable :: iparenttemp(:)
      integer, allocatable :: nchild(:)
      integer, allocatable :: ichildtemp(:,:)
      integer, allocatable :: nnbors(:)
      integer, allocatable :: nbors(:,:)
      integer, allocatable :: isourcetemp(:)
      integer, allocatable :: itargettemp(:)
      integer, allocatable :: iexpctemp(:)
      integer, allocatable :: ihsfirsttemp(:)
      integer, allocatable :: ihslasttemp(:)
      integer, allocatable :: isfirsttemp(:)
      integer, allocatable :: islasttemp(:)
      integer, allocatable :: itfirsttemp(:)
      integer, allocatable :: itlasttemp(:)
      integer, allocatable :: ihefirsttemp(:)
      integer, allocatable :: ihelasttemp(:)
      integer, allocatable :: iefirsttemp(:)
      integer, allocatable :: ielasttemp(:)
      integer, allocatable :: nhungsrc(:)
      integer, allocatable :: nhungexp(:)
      integer, allocatable :: nhunglistsrc(:)
      integer, allocatable :: ihunglistsrc(:,:)

      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer, allocatable :: nlist1(:)
      integer, allocatable :: list1(:,:)
      integer, allocatable :: nlist2(:)
      integer, allocatable :: list2(:,:)
      integer, allocatable :: nlist3(:)
      integer, allocatable :: list3(:,:)
      integer, allocatable :: nlist4(:)
      integer, allocatable :: list4(:,:)

      integer ipointer(30)
      real *8 boxsize(0:nlevels)
      real *8 src(2,ns),radsrc(ns)
      real *8 trg(2,nt)
      real *8 centers(2,nboxes)
      real *8 expc(2,nexpc)
      real *8 radexp(nexpc)
c
      real *8 xmin,xmax,ymin,ymax,sizey
      integer ictr,ih,irefine,is,ie
      integer nss,nee

c
c     INPUT:
c
c     src           source locations        
c     ns            number of sources 
c     rads          source radii (determines deepest level that
c                   the source can reach) 
c     trg           target locations        
c     nt            number of targets
c
c     expc          expansion center locations
c     nexpc         number of expansion centers
c
c     idivflag      0 => divide on sources
c                   1 => divide on targets
c                   2 => divide on sources+targets
c                   3 => divide on sources+targets+expansion centers
c 
c     ndiv          refinement criterion - extend tree until all
c                   nodes at finest level have less than ndiv 
c                   source/targets/sources+targets depending on
c                   idivflag
c
c     ltree         length of ltree
c     nlevels       number of levels determined using the memory code
c     nboxes        number of boxes determined using mem code
c                   refer to mkplratreemem.
c
c
c     OUTPUT:
c     centers       array of box centers
c     boxsize       box dimensions at all levels
c
c     itree         tree array
c     
c     iladdr = 1
c     itree(iladdr) <-> laddr
c                   indexing array providing access to boxes at
c                   each level. 
c                   the first box on level i is laddr(1,i)
c                   the last box on level i is laddr(2,i)
c
c     iparent = iladdr + 2*(nlevels+1)
c     itree(iparent) <-> parent
c                   parent of box (set to -1 for root node)
c
c     inchild = iparent + nboxes
c     itree(inchild) <-> nchild
c                   number of children for each box 
c
c     ichild = inchild + nboxes
c     itree(ichild) <-> child
c                   only first nchild entries are defined 
c                   (others set to -1). The array should be viewed 
c                   as dimensioned (4,nboxes)
c
c     isource = ichild + 4*nboxes
c     itree(isource) <-> isource
c                   tree-ordered array of sources
c
c     itarget = isource + ns
c     itree(itarget) <-> itarget
c                   tree-ordered array of targets
c
c     iexpc = itarget + nt
c     itree(iexpc) <-> iexpc
c                      tree-ordered array of expansion centers
c
c     ihsfirst = itarget + nt
c     itree(ihsfirst) <-> ihsfirst
c                   ihsfirst(j) = location in isource of first hung 
c                                source for box j 
c
c     ihslast = ihsfirst + nboxes
c     itree(ihslast) <-> ihslast
c                   ihslast(j)  = location in isource of last hung
c                                source for box j 
c
c     isfirst = ihlast + nboxes
c     itree(isfirst) <-> isfirst
c                   isfirst(j) = location in isource of first 
c                                source for box j 
c
c     islast = isfirst + nboxes
c     itree(islast) <-> islast
c                   islast(j)  = location in isource of last 
c                                source for box j 
c
c     itfirst = islast + nboxes
c     itree(itfirst) <-> itfirst
c                   itfirst(j) = location in itarget of first 
c                                target for box j 
c
c     itlast = itfirst + nboxes
c     itree(itlast) <-> itlast
c                   itlast(j)  = location in itarget of last 
c                                target for box j
c
c     ihefirst = itlast + nboxes
c     itree(ihefirst) <-> ihefirst
c                   ihefirst(j) = location in iexpc of first hung 
c                                 expansion center for box j 
c
c     ihelast = ihefirst + nboxes
c     itree(ihelast) <-> ihelast
c                   ihelast(j)  = location in iexpc of last hung
c                                expansion center for box j 
c
c     iefirst = ihelast + nboxes
c     itree(iefirst) <-> iefirst
c                    iefirst(j) = location in iexpc of first target
c                                 in box j
c
c     ielast = iefirst + nboxes
c     itree(ielast) <-> ielast
c                       ielast(j) = location in iexpc of last target
c                                   in box j
c
c     inlist1 = ielast + nboxes
c     itree(inlist1) <-> nlist1
c                        nlist1(i) = number of boxes in list 1
c                        of box i. The list 1 of a box i, Ui, is
c                        the set of boxes which touch
c                        box i. In a level restricted tree,
c                        the maximum number of boxes in list
c                        1 of a box can be 13 (mnlist1 = 13)
c
c 
c     ilist1 = inlist1 + nboxes
c     itree(ilist1) <-> list1
c                       list1(j,i) is the id of the jth box
c                       in list 1 of box i
c
c     inlist2 = ilist1 + mnlist1*nboxes
c     itree(inlist2) <-> nlist2
c                       nlist2(i) = number of boxes in list 2 of
c                       box i. The list 2 of a box i, Vi, is
c                       the set of boxes which are descendants
c                       of the colleagues of the parent of box i
c                       which are well separated from box i
c                       at the scale of box i. In a level
c                       restricted tree, the maximum number of boxes
c                       in list 2 of a box can be 27 (mnlist2=27)
c
c     ilist2 = inlist2 + nboxes
c     itree(ilist2) <-> list2
c                       list2(j,i) is the id fo the jth box in 
c                       list 2 of box i
c
c     inlist3 = ilist2 + mnlist2*nboxes
c     itree(inlist3) <-> nlist3
c                        nlist3(i) = number of boxes in list3 of
c                        box i. The list 3 of a box i, Wi, is
c                        the set of boxes which are descendants
c                        of the colleagues of box i which
c                        are not in list 1 of box i. In a level
c                        restricted tree, the maximum number of
c                        boxes in list 3 of a box can be 
c                        20 (mnlist3 = 20). Note
c                        that list 3 of a box is non empty
c                        iff the box is childless
c
c     ilist3 = inlist3 + nboxes
c     itree(ilist3) <-> list3(j,i) is the id of the jth box in
c                       list 3 of box i
c
c     inlist4 = ilist3 + mnlist3*nboxes
c     itree(inlist4) <-> nlist4
c                        nlist4(i) = number of boxes in list4
c                        of box i. The list 4 of a box i, Xi,
c                        is dual to Wi. j in Xi if i in Wj, that
c                        is box j is in list 4 of box i if i is
c                        in list 3 of box j. In a level restricted
c                        tree the max number of boxes in list 4
c                        of a box is 5.
c       
c     ilist4 = inlist4 + nboxes
c     itree(ilist4) <-> list4
c                       list4(j,i) is the id of the jth box
c                       in list 4 of box i
c      
c     inhungsrc = ielast + nboxes
c     itree(inhungsrc) <-> inhungsrc
c                inhungsrc(j)  = number of hung chunks in box j
c
c     inhungexp = inhungsrc + nboxes
c     itree(inhungexp) <-> inhungexp
c                inhungexp(j) = number of hung expansion centers
c                               in box j
c 
c     nhunglistsrc = inhungexp + nboxes
c     itree(inhunglistsrc) <-> inhunglistsrc
c           inhunglistsrc(j) = Total number of hung sources
c                              relevant for box j
c
c     ihunglistsrc = inhungexp + nboxes
c     itree(ihunglistsrc) <-> ihunglistsrc
c                   ihunglistsrc(m,j) = src id of  mth hung src
c                   relevant to box j. A hung src is relevant to
c                   a box if the box is a descendant of the
c                   neighbor of the box in which the src is hung
c
c     ltree = ihunglistsrc + mhung*nboxes
c
c     ipointers is the collection of pointers
c     ipointer(1) = iladdr
c     ipointer(2) = iparent = ipointer(1) + 2*(nlevels+1)
c     ipointer(3) = inchild = ipointer(2) + nboxes
c     ipointer(4) = ichild = ipointer(3) + nboxes
c     ipointer(5) = isource = ipointer(4) + 4*nboxes
c     ipointer(6) = itarget = ipointer(5) + ns
c     ipointer(7) = iexpc = ipointer(6) + nt
c     ipointer(8) = ihsfirst = ipointer(7) + nexpc
c     ipointer(9) = ihslast = ipointer(8) + nboxes
c     ipointer(10) = isfirst = ipointer(9) + nboxes
c     ipointer(11) = islast = ipointer(10) + nboxes
c     ipointer(12) = itfirst = ipointer(11) + nboxes
c     ipointer(13) = itlast = ipointer(12) + nboxes
c     ipointer(14) = ihefirst = ipointer(13) + nboxes
c     ipointer(15) = ihelast = ipointer(14) + nboxes
c     ipointer(16) = iefirst = ipointer(15) + nboxes
c     ipointer(17) = ielast = ipointer(16) + nboxes
c     ipointer(18) = inlist1 = ipointer(17) + nboxes
c     ipointer(19) = ilist1 = ipointer(18) + nboxes
c     ipointer(20) = inlist2 = ipointer(19) + mnlist1*nboxes
c     ipointer(21) = ilist2 = ipointer(20) + nboxes
c     ipointer(22) = inlist3 = ipointer(21) + mnlist2*nboxes
c     ipointer(23) = ilist3 = ipointer(22) + nboxes
c     ipointer(24) = inlist4 = ipointer(23) + mnlist3*nboxes
c     ipointer(25) = ilist4 = ipointer(24) + nboxes
c     ipointer(26) = nhungsrc = ipointer(25) + mnlist4*nboxes
c     ipointer(27) = nhungexp = ipointer(26) + nboxes
c     ipointer(28) = nhunglistsrc = ipointer(27) + nboxes
c     ipointer(29) = hunglistsrc = ipointer(28) + nboxes
c     ipointer(30) = lentree = ipointer(29) + mhung*nboxes
c
c     mnlist1 : max number of box possible in list 1 of box
c     mnlist2 : max number of box possible in list 2 of box
c     mnlist3 : max number of box possible in list 3 of box
c     mnlist4 : max number of box possible in list 4 of box
c---------------------------------------------------------------------
c
c     Other notes:
c
c     refinement criterion w.r.t rads is:
c
c     hang if (rads .geq. boxsize)
c
c
      lcenters = nboxes

      allocate(laddr(2,0:nlevels))
      allocate(ilevel(nboxes))
      allocate(iparenttemp(nboxes))
      allocate(nchild(nboxes))
      allocate(ichildtemp(4,nboxes))
      allocate(isourcetemp(ns))
      allocate(itargettemp(nt))
      allocate(iexpctemp(nexpc))
      allocate(ihsfirsttemp(nboxes))
      allocate(ihslasttemp(nboxes))
      allocate(isfirsttemp(nboxes))
      allocate(islasttemp(nboxes))
      allocate(itfirsttemp(nboxes))
      allocate(itlasttemp(nboxes))
      allocate(ihefirsttemp(nboxes))
      allocate(ihelasttemp(nboxes))
      allocate(iefirsttemp(nboxes))
      allocate(ielasttemp(nboxes))
      allocate(nhungsrc(nboxes))
      allocate(nhungexp(nboxes))
      allocate(nbors(9,nboxes))
      allocate(nnbors(nboxes))

c     Step 1: Find enclosing box
c
      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
c
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
      enddo
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
      enddo

      do i=1,nexpc
        if(expc(1,i) .lt. xmin) xmin=expc(1,i)
        if(expc(1,i) .gt. xmax) xmax=expc(1,i)
        if(expc(2,i) .lt. ymin) ymin=expc(2,i)
        if(expc(2,i) .gt. ymax) ymax=expc(2,i)
      enddo
      boxsize(0)=xmax-xmin
      sizey=ymax-ymin
      if(sizey .gt. boxsize(0)) boxsize(0)=sizey
c
c     initialize arrays at level 0
c
      centers(1,1)=(xmin+xmax)/2
      centers(2,1)=(ymin+ymax)/2
      laddr(1,0) = 1
      laddr(2,0) = 1
      ilevel(1) = 0
      iparenttemp(1) = -1
      isfirsttemp(1) = 1
      nhungsrc(1) = 0
      nhungexp(1) = 0
c
c     count number of hung sources
c     and hang up "big" sources
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) nhungsrc(1)=nhungsrc(1)+1 
      enddo
      isfirsttemp(1) = nhungsrc(1)+1
      islasttemp(1) = ns
      if (nhungsrc(1) .gt. 0) then 
         ihsfirsttemp(1) = 1
         ihslasttemp(1) = nhungsrc(1)
      else
         ihsfirsttemp(1) = 0
         ihslasttemp(1) = -1
      endif

c     Count number of hung expansion centers      
c     and hang up "big" expansion centers
      do i=1,nexpc
         if (radexp(i).gt.boxsize(0)) nhungexp(1)=nhungexp(1)+1
      enddo
      iefirsttemp(1) = nhungexp(1)+1
      ielasttemp(1) = nexpc
      if (nhungexp(1).gt.0) then
          ihefirsttemp(1) = 1
          ihelasttemp(1) = nhungexp(1)
      else
         ihefirsttemp(1) = 0
         ihelasttemp(1) = -1
      endif

c
c     reorder isourcetemp to put hung sources in beginning
      ih = 0
      is = nhungsrc(1)
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) then
            ih = ih+1
            isourcetemp(ih) = i
         else
            is = is+1
            isourcetemp(is) = i
         endif
      enddo
c     reorder iexptemp to put hung expansion centers in beginning
      ih = 0
      ie = nhungexp(1)
      do i= 1,nexpc
         if(radexp(i).gt.boxsize(0)) then
            ih = ih+1
            iexpctemp(ih) = i
         else
            ie = ie+1
            iexpctemp(ie) = i
         endif
      enddo

c     initialize itargettemp 
      do i = 1,nt
         itargettemp(i) = i
      enddo
      itfirsttemp(1) = 1
      itlasttemp(1) = nt

      nchild(1) = 0
      ichildtemp(1,1) = -1
      ichildtemp(2,1) = -1
      ichildtemp(3,1) = -1
      ichildtemp(4,1) = -1

c
      irefine = 0
      nss = ns - nhungsrc(1)
      nee = nexpc - nhungexp(1)
      if ((idivflag .eq.0).and.(nss.gt.ndiv)) irefine=1
      if ((idivflag .eq.1).and.(nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.2).and.(nss+nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.3).and.(nss+nt+nee.gt.ndiv)) irefine=1

c     Reset nlevels, nboxes      
      nbmax = nboxes
      nlmax = nlevels
      nlevels = 0
      nboxes = 1
      do i = 1,nlmax
         if (irefine.eq.1) then
            call subdivide_adap(src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   nlevels,nboxes,
     $                   centers,lcenters,boxsize,nbmax,nlmax,
     $                   laddr,ilevel,iparenttemp,nchild,ichildtemp,
     $                   isourcetemp,itargettemp,iexpctemp,
     $                   ihsfirsttemp,ihslasttemp,
     $                   isfirsttemp,islasttemp,
     $                   itfirsttemp,itlasttemp,
     $                   ihefirsttemp,ihelasttemp,
     $                   iefirsttemp,ielasttemp,nhungsrc,
     $                   nhungexp,irefine)

         else
            exit
         endif
      enddo

      do i=1,nboxes
         nnbors(i) = 0
         do j=1,9
            nbors(j,i) = -1
         enddo
      enddo

c     compute colleagues
      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparenttemp,nchild,
     2                   ichildtemp,nnbors,nbors)
c      call prinf('nboxes=*',nboxes,1)


       if(nlevels.ge.2) 
     $  call d2hpfixtree(src,ns,radsrc,trg,nt,expc,nexpc,
     $              radexp,nlevels,nboxes,
     $              centers,boxsize,nbmax,
     $              laddr,ilevel,iparenttemp,nchild,ichildtemp,
     $              nnbors,nbors,  
     $              isourcetemp,itargettemp,iexpctemp,
     $              ihsfirsttemp,ihslasttemp,
     $              isfirsttemp,islasttemp,
     $              itfirsttemp,itlasttemp,
     $              ihefirsttemp,ihelasttemp,
     $              iefirsttemp,ielasttemp,nhungsrc,
     $              nhungexp)

c      call prinf('nboxes=*',nboxes,1)
       
      mnlist1 = 13
      mnlist2 = 27
      mnlist3 = 20
      mnlist4 = 5
      allocate(nlist1(nboxes))
      allocate(nlist2(nboxes))
      allocate(nlist3(nboxes))
      allocate(nlist4(nboxes))

      allocate(list1(mnlist1,nboxes))
      allocate(list2(mnlist2,nboxes))
      allocate(list3(mnlist3,nboxes))
      allocate(list4(mnlist4,nboxes))

      allocate(nhunglistsrc(nboxes))
      allocate(ihunglistsrc(mhung,nboxes))

      do i=1,nboxes
         nlist1(i) = 0
         nlist2(i) = 0
         nlist3(i) = 0
         nlist4(i) = 0
         nhunglistsrc(i) = 0
         do j=1,mnlist1
            list1(j,i) = -1
         enddo
         do j=1,mnlist2
            list2(j,i) = -1
         enddo
         do j=1,mnlist3
            list3(j,i) = -1
         enddo
         do j=1,mnlist4
            list4(j,i) = -1
         enddo
         do j=1,mhung
            ihunglistsrc(j,i) = -1
         enddo
      enddo
      call computelists(nlevels,nboxes,laddr,boxsize,centers,
     1                  iparenttemp,nchild,ichildtemp,nnbors,
     2                  nbors,mnlist1,nlist1,list1,mnlist2,nlist2,
     3                  list2,mnlist3,nlist3,list3,mnlist4,
     4                  nlist4,list4)

c     Compute hunglist  
      if(mhung.gt.0) then
         call computehunglist(mhung,nlevels,nboxes,laddr,
     1                        ns,isourcetemp,iparenttemp,nchild,
     2                        ihsfirsttemp,ihslasttemp,
     3                        nnbors,nbors,mnlist1,nlist1,
     4                        list1,nhungsrc,
     5                        nhunglistsrc,ihunglistsrc)

      endif


c     Move the output to itree
c     Store iladdr - first and last box at level ilev
      ictr = 1
      ipointer(1) = ictr
      do i=0,nlevels
         itree(ictr) = laddr(1,i)
         ictr = ictr + 1
         itree(ictr) = laddr(2,i)
         ictr = ictr + 1
      enddo

c     Store iparent - parent of box i
      ipointer(2)  = ictr
      do i = 1,nboxes
          itree(ictr) = iparenttemp(i)
          ictr = ictr + 1
      enddo

c     Store inchild - number of children of given box
      ipointer(3)  = ictr
      do i = 1,nboxes
         itree(ictr) = nchild(i)
         ictr = ictr + 1
      enddo 

c     Store ichild temp - children of given box i
      ipointer(4)  = ictr
      do i = 1,nboxes
         do j=1,4
             itree(ictr) = ichildtemp(j,i)
             ictr = ictr + 1
         enddo
      enddo

c     isource - tree sorted array of sources
      ipointer(5)  = ictr
      do i =1,ns
         itree(ictr) = isourcetemp(i)
         ictr = ictr + 1
      enddo

c     itarget - tree sorted array of targets
      ipointer(6)  = ictr
      do i = 1,nt
         itree(ictr) = itargettemp(i)
         ictr = ictr + 1
      enddo

      ipointer(7) = ictr
      do i= 1,nexpc
         itree(ictr) = iexpctemp(i)
         ictr = ictr + 1
      enddo

c     ihfirst - first hung chunk in box i
      ipointer(8)  = ictr
      do i=1,nboxes
         itree(ictr)  = ihsfirsttemp(i)
         ictr = ictr + 1
      enddo

c     ihlast - last hung chunk in box i
      ipointer(9)  = ictr
      do i=1,nboxes
         itree(ictr) = ihslasttemp(i)
         ictr = ictr + 1
      enddo

c     isfirst - first source in box i
      ipointer(10)  = ictr
      do i=1,nboxes
         itree(ictr)  = isfirsttemp(i)
         ictr = ictr + 1
      enddo

c     islast - last source in box i
      ipointer(11)  = ictr
      do i=1,nboxes
         itree(ictr) = islasttemp(i)
         ictr = ictr + 1
      enddo

c     itfirst - first target in box i
      ipointer(12)  = ictr
      do i=1,nboxes
         itree(ictr)  = itfirsttemp(i)
         ictr = ictr + 1
      enddo

c     itlast - last target in box i
      ipointer(13)  = ictr
      do i=1,nboxes
         itree(ictr) = itlasttemp(i)
         ictr = ictr + 1
      enddo

c     ihefirst - first hung exansion center in
c     box i
      ipointer(14) = ictr
      do i=1,nboxes
          itree(ictr) = ihefirsttemp(i)
          ictr = ictr + 1
      enddo

c     ihelast - last hung expansion center in
c     box i
      ipointer(15) = ictr
      do i=1,nboxes
         itree(ictr) = ihelasttemp(i)
         ictr = ictr + 1
      enddo

c     iefirst - first exansion center in
c     box i
      ipointer(16) = ictr
      do i=1,nboxes
          itree(ictr) = iefirsttemp(i)
          ictr = ictr + 1
      enddo

c     ielast - last hung expansion center in
c     box i
      ipointer(17) = ictr
      do i=1,nboxes
         itree(ictr) = ielasttemp(i)
         ictr = ictr + 1
      enddo

c     nlist1 - number of boxes in list 1 of a given box
      ipointer(18) = ictr
      do i=1,nboxes
         itree(ictr) = nlist1(i)
         ictr = ictr + 1
      enddo

c     list1(j,i) - id of jth box in list 1 of box i
      ipointer(19) = ictr
      do i=1,nboxes
         do j=1,mnlist1
            itree(ictr) = list1(j,i)
            ictr = ictr + 1
         enddo
      enddo
      
c     nlist2 - number of boxes in list 2 of a given box
      ipointer(20) = ictr
      do i=1,nboxes
         itree(ictr) = nlist2(i)
         ictr = ictr + 1
      enddo

c     list2(j,i) - id of jth box in list 1 of box i
      ipointer(21) = ictr
      do i=1,nboxes
         do j=1,mnlist2
            itree(ictr) = list2(j,i)
            ictr = ictr + 1
         enddo
      enddo
      
c     nlist3 - number of boxes in list 3 of a given box
      ipointer(22) = ictr
      do i=1,nboxes
         itree(ictr) = nlist3(i)
         ictr = ictr + 1
      enddo

c     list3(j,i) - id of jth box in list 3 of box i
      ipointer(23) = ictr
      do i=1,nboxes
         do j=1,mnlist3
            itree(ictr) = list3(j,i)
            ictr = ictr + 1
         enddo
      enddo
      
c     nlist4 - number of boxes in list 4 of a given box
      ipointer(24) = ictr
      do i=1,nboxes
         itree(ictr) = nlist4(i)
         ictr = ictr + 1
      enddo

c     list4(j,i) - id of jth box in list 4 of box i
      ipointer(25) = ictr
      do i=1,nboxes
         do j=1,mnlist4
            itree(ictr) = list4(j,i)
            ictr = ictr + 1
         enddo
      enddo

c     nhungsrc  - Number of hung chunks in a given box
      ipointer(26)  = ictr
      do i=1,nboxes
         itree(ictr) = nhungsrc(i)
         ictr = ictr + 1
      enddo
 
c     nhungexp - Number of hung expansion centers
c                in box i
      ipointer(27) = ictr
      do i=1,nboxes
         itree(ictr) = nhungexp(i)
         ictr = ictr + 1
      enddo

c     nhunglistsrc - Number of hung chunks for box i
c        (ordered by given ordering of targets, not tree sorted)
      ipointer(28)  = ictr
      do i=1,nboxes
         itree(ictr) = nhunglistsrc(i)
         ictr = ictr + 1
      enddo

c     hunglistsrc - chunk hunglist for each target
      ipointer(29)  = ictr
      do i=1,nboxes
         do j=1,mhung
            itree(ictr) = ihunglistsrc(j,i)
            ictr=ictr+1
         enddo
      enddo
      ipointer(30) = ictr

      return
      end
c-------------------------------------------------------------------
      subroutine subdivide_adap(src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   nlevels,nboxes,
     $                   centers,lcenters,boxsize,nbmax,nlmax,
     $                   laddr,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,irefine)
      implicit none
      integer ns,nt,nexpc,idivflag,ndiv,mhung
      integer nlevels,nboxes,lcenters,nbmax,nlmax
      integer irefine, irefinebox
      real *8 src(2,ns),radsrc(ns)
      real *8 trg(2,nt)
      real *8 expc(2,nexpc),radexp(nexpc)
      real *8 centers(2,lcenters)
      real *8 boxsize(0:nlmax)
      integer laddr(2,0:nlmax)
      integer ilevel(nbmax)
      integer iparent(nbmax)
      integer nchild(nbmax)
      integer ichild(4,nbmax)
      integer isource(ns)
      integer itarget(nt)
      integer iexpc(nexpc)
      integer ihsfirst(nbmax)
      integer ihslast(nbmax)
      integer isfirst(nbmax)
      integer islast(nbmax)
      integer itfirst(nbmax)
      integer itlast(nbmax)
      integer ihefirst(nbmax)
      integer ihelast(nbmax)
      integer iefirst(nbmax)
      integer ielast(nbmax)
      integer nhungsrc(nbmax)
      integer nhungexp(nbmax)
c     Temporary variables
      integer isrctmp(ns),itargtmp(nt),iexpctmp(nexpc)
      integer i,j,i12,i34,istart,jstart,kstart,ii,iii,nss,nee,ntt
      integer ibox,ifirstbox,ilastbox,nbfirst
      integer is,it,ie
      integer nsc(4),ntc(4),nh(4),nexpcc(4),nhc(4)
c
c     for every box at level nlevels,
c     sort into children, updating various arrays 
c     perhaps just build tree here paren/child/particle sorting...
c     lists in second call ???
c     
c     allocate temp array for isourcetemp2 itargtemp2
c     after all done, write back to isourcetemp, itargtemp
c     this is O(N) * nlevels work for rewriting.
c     can be fancier I suppose.
c    
      irefine = 0
      ifirstbox = laddr(1,nlevels)
      ilastbox =  laddr(2,nlevels)

c
      nbfirst = nboxes+1
      boxsize(nlevels+1) = boxsize(0)/2.0d0**(nlevels+1)

      do ibox = ifirstbox,ilastbox
c        Determine if current box needs to be subdivided
c        if current box needs to be subdivided then
c        set irefinebox = 1
         nss = islast(ibox) - isfirst(ibox) + 1
         ntt = itlast(ibox) - itfirst(ibox) + 1
         nee = ielast(ibox) - iefirst(ibox) + 1
         irefinebox = 0
         if ((idivflag .eq.0).and.(nss.gt.ndiv)) irefinebox = 1
         if ((idivflag .eq.1).and.(ntt.gt.ndiv)) irefinebox=1
         if ((idivflag .eq.2).and.(nss+ntt.gt.ndiv)) irefinebox=1
         if ((idivflag .eq.3).and.(nss+ntt+nee.gt.
     1                                    ndiv)) irefinebox=1  
         if(irefinebox.eq.1) then
c           Based on the subdivision criterion, the current
c           box needs to be divided. 
c
c           Allocate temporary array to figure out which child you
c           belong to
c           which child?  1,2,3,4? counter ns1,ns2,ns3,ns4
c           The box nomenclature is as follows
c           3   4
c           1   2
            i12 = isfirst(ibox)-1
            i34 = 0
            do is = isfirst(ibox),islast(ibox)
               if(src(2,isource(is)) - centers(2,ibox).lt.0) then
                  i12 = i12+1
                  isource(i12) = isource(is)
               else
                  i34 = i34 + 1
                  isrctmp(i34) = isource(is)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end

c           Reorder sources to include sources in 34 in the array
            do i=1,i34
               isource(i12+i) = isrctmp(i)
            enddo
            nsc(1) = 0
            nsc(2) = 0
            nsc(3) = 0
            nsc(4) = 0
c           Sort into boxes 1 and 2
            do is = isfirst(ibox),i12
               if(src(1,isource(is))-centers(1,ibox).lt.0) then
                  isource(isfirst(ibox)+nsc(1)) = isource(is)
                  nsc(1) = nsc(1) + 1
               else
                  nsc(2) = nsc(2) + 1
                  isrctmp(nsc(2)) = isource(is)
               endif
            enddo
c           Reorder sources so that sources in 2 are at the
c           end of the array
            do i=1,nsc(2)
               isource(isfirst(ibox)+nsc(1)+i-1) = isrctmp(i)
            enddo
c           Sort into boxes 3 and 4
            do is = i12+1, islast(ibox)
               if(src(1,isource(is))-centers(1,ibox).lt.0) then
                   isource(i12+1+nsc(3)) = isource(is)
                   nsc(3) = nsc(3) + 1
                else
                   nsc(4) = nsc(4)+1
                   isrctmp(nsc(4)) = isource(is)
                endif
            enddo
            do i=1,nsc(4)
               isource(i12+nsc(3)+i) = isrctmp(i)
            enddo
c           End of sorting sources

            istart = isfirst(ibox)-1
            do j=1,4
c           check hung -> counter nh1,nh2,nh3,nh4
               ii = 0
               nh(j) = 0
               do i=1,nsc(j)
                  if(radsrc(isource(istart+i)).gt.boxsize(nlevels+1))
     1            then     
                     nh(j) = nh(j) + 1
                     isource(istart+nh(j)) = isource(istart+i)
                   else
                      ii = ii+1
                      isrctmp(ii) = isource(istart+i)
                   endif
                enddo
c            Reorder sources to have hung chunks at the star
c            of the sorted sources in the box ibox
                do i=1,ii
                   isource(istart+nh(j)+i) = isrctmp(i)
                enddo
                istart = istart + nsc(j)
            enddo
         
c           which child?  1,2,3,4? counter nt(1),nt(2),ns(3),ns(4)
c           The box nomenclature is as follows
c           3   4
c           1   2
            i12 = itfirst(ibox)-1
            i34 = 0
            do it = itfirst(ibox),itlast(ibox)
               if(trg(2,itarget(it)) - centers(2,ibox).lt.0) then
                  i12 = i12+1
                  itarget(i12) = itarget(it)
               else
                  i34 = i34 + 1
                  itargtmp(i34) = itarget(it)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end

c           Reorder sources to include sources in 34 in the array
            do i=1,i34
               itarget(i12+i) = itargtmp(i)
            enddo
            ntc(1) = 0
            ntc(2) = 0
            ntc(3) = 0
            ntc(4) = 0
c           Sort into boxes 1 and 2
            do it = itfirst(ibox),i12
               if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
                  itarget(itfirst(ibox)+ntc(1)) = itarget(it)
                  ntc(1) = ntc(1) + 1
               else
                  ntc(2) = ntc(2) + 1
                  itargtmp(ntc(2)) = itarget(it)
               endif
            enddo
c           Reorder sources so that sources in 2 are at the
c           end of the array
            do i=1,ntc(2)
               itarget(itfirst(ibox)+ntc(1)+i-1) = itargtmp(i)
            enddo
c           Sort into boxes 3 and 4
            do it = i12+1, itlast(ibox)
               if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
                  itarget(i12+1+ntc(3)) = itarget(it)
                  ntc(3) = ntc(3) + 1
               else
                  ntc(4) = ntc(4)+1
                  itargtmp(ntc(4)) = itarget(it)
               endif
            enddo
            do i=1,ntc(4)
               itarget(i12+ntc(3)+i) = itargtmp(i)
            enddo
c           End of sorting targets

c           Sort expansion centers
c           which child?  1,2,3,4? counter nt(1),nt(2),ns(3),ns(4)
c           The box nomenclature is as follows
c           3   4
c           1   2
            i12 = iefirst(ibox)-1
            i34 = 0
            do ie = iefirst(ibox),ielast(ibox)
               if(expc(2,iexpc(ie)) - centers(2,ibox).lt.0) then
                  i12 = i12+1
                  iexpc(i12) = iexpc(ie)
               else
                  i34 = i34 + 1
                  iexpctmp(i34) = iexpc(ie)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end

c           Reorder sources to include sources in 34 in the array
            do i=1,i34
               iexpc(i12+i) = iexpctmp(i)
            enddo
            nexpcc(1) = 0
            nexpcc(2) = 0
            nexpcc(3) = 0
            nexpcc(4) = 0
c           Sort into boxes 1 and 2
            do ie = iefirst(ibox),i12
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(iefirst(ibox)+nexpcc(1)) = iexpc(ie)
                  nexpcc(1) = nexpcc(1) + 1
               else
                  nexpcc(2) = nexpcc(2) + 1
                  iexpctmp(nexpcc(2)) = iexpc(ie)
               endif
            enddo
c           Reorder sources so that sources in 2 are at the
c           end of the array
            do i=1,nexpcc(2)
               iexpc(iefirst(ibox)+nexpcc(1)+i-1) = iexpctmp(i)
            enddo
c           Sort into boxes 3 and 4
            do ie = i12+1, ielast(ibox)
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(i12+1+nexpcc(3)) = iexpc(ie)
                  nexpcc(3) = nexpcc(3) + 1
               else
                  nexpcc(4) = nexpcc(4)+1
                  iexpctmp(nexpcc(4)) = iexpc(ie)
               endif
            enddo
            do i=1,nexpcc(4)
               iexpc(i12+nexpcc(3)+i) = iexpctmp(i)
            enddo
c           End of sorting expanison centers

            istart = iefirst(ibox)-1
            do j=1,4
c           check hung -> counter nhc1,nhc2,nhc3,nhc4
               ii = 0
               nhc(j) = 0
               do i=1,nexpcc(j)
                  if(radexp(iexpc(istart+i)).gt.boxsize(nlevels+1))
     1            then     
                     nhc(j) = nhc(j) + 1
                     iexpc(istart+nhc(j)) = iexpc(istart+i)
                   else
                      ii = ii+1
                      iexpctmp(ii) = iexpc(istart+i)
                   endif
                enddo
c           Reorder sources to have hung chunks at the star
c           of the sorted sources in the box ibox
               do i=1,ii
                  iexpc(istart+nhc(j)+i) = iexpctmp(i)
               enddo
               istart = istart + nexpcc(j)
            enddo

            nchild(ibox) = 0
c           Create the required boxes
            istart = isfirst(ibox)
            jstart = itfirst(ibox)
            kstart = iefirst(ibox)
            do i=1,4
               ii = 2
               if(i.eq.1.or.i.eq.2) ii = 1
               if(nsc(i)+ntc(i)+nexpcc(i).ge.0) then
c                 Increment total number of boxes               
                  nboxes = nboxes + 1
c                 Increment number of children for the current box
                  nchild(ibox) = nchild(ibox)+1
c                 Update the array of children for the current box
                  ichild(nchild(ibox),ibox) = nboxes
c                 Update the array of levels for the child box
                  ilevel(nboxes) = nlevels+1
c                 Update the array of parents for the child box
                  iparent(nboxes) = ibox
c                 Compute center for the child box
                  centers(1,nboxes) = centers(1,ibox)+(-1)**i*
     1                                boxsize(nlevels+1)/2.0
                  centers(2,nboxes) = centers(2,ibox)+(-1)**ii*
     1                                boxsize(nlevels+1)/2.0
c                 Update children information for newly
c                 created child box
                  nchild(nboxes) = 0
                  ichild(1,nboxes) = -1
                  ichild(2,nboxes) = -1
                  ichild(3,nboxes) = -1
                  ichild(4,nboxes) = -1

c                 Update arrays ihsfirst,ihslast,isfirst,islast
                  ihsfirst(nboxes) = istart
                  ihslast(nboxes) = istart + nh(i) - 1
                  nhungsrc(nboxes) = nh(i)

                  isfirst(nboxes) = istart + nh(i)
                  islast(nboxes) = istart + nsc(i) - 1

c                 Update arrays itfirst, itlast
                  itfirst(nboxes) = jstart
                  itlast(nboxes) = jstart + ntc(i) - 1

c                 Update arrays ihefirst,ihelast,iefirst,ielast
                  ihefirst(nboxes) = kstart
                  ihelast(nboxes) = kstart + nhc(i)-1
                  nhungexp(nboxes) = nhc(i)

                  iefirst(nboxes) = kstart + nhc(i)
                  ielast(nboxes) = kstart + nexpcc(i) - 1

                  nss = islast(nboxes) - isfirst(nboxes)+1
                  nee = ielast(nboxes) - iefirst(nboxes)+1
c                 Check if further refinement required
                  if ((idivflag .eq.0).and.(nss.gt.ndiv)) irefine = 1
                  if ((idivflag .eq.1).and.(ntc(i).gt.ndiv)) irefine=1
                  if ((idivflag .eq.2).and.(nss+ntc(i).gt.ndiv)) 
     1                                          irefine=1
                  if ((idivflag .eq.3).and.(nss+ntc(i)+nee.gt.
     1                                    ndiv)) irefine=1  
                endif
                istart = istart + nsc(i)
                jstart = jstart + ntc(i)
                kstart = kstart + nexpcc(i)
            enddo
         endif
      enddo
      nlevels = nlevels+1
      laddr(1,nlevels) = nbfirst
      laddr(2,nlevels) = nboxes

      return
      end
c-------------------------------------------------------------      

      subroutine computecoll(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,
     2                       nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: real *8(0:;nlevels)
c                 Array of boxsizes
c 
c     centers     in: real *8(2,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(4,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes
      integer laddr(2,0:nlevels)
      real *8 boxsize(0:nlevels)
      real *8 centers(2,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(4,nboxes)
      integer nnbors(nboxes)
      integer nbors(9,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox


c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,nchild(jbox)
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
c               Check if kbox is a nearest neighbor or in list 2
                   if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                1.05*boxsize(ilev)).and.
     2                (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                1.05*boxsize(ilev))) then
                     
                     nnbors(ibox) = nnbors(ibox)+1
                     nbors(nnbors(ibox),ibox) = kbox
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
      enddo

      return
      end

      subroutine d2hpfixtree(src,ns,radsrc,trg,nt,expc,nexpc,
     $              radexp,nlevels,nboxes,
     $              centers,boxsize,nbmax,
     $              laddr,ilevel,iparent,nchild,ichild,
     $              nnbors,nbors, 
     $              isource,itarget,iexpc,
     $              ihsfirst,ihslast,
     $              isfirst,islast,
     $              itfirst,itlast,
     $              ihefirst,ihelast,
     $              iefirst,ielast,nhungsrc,
     $              nhungexp)
c     Give an adaptive tree, this subroutine fixes the tree
c     to make it level restricted. A level restricted tree
c     is one where no two boxes that contact each other are
c     more than one level apart. The subroutine corrects
c     a given adaptive tree by adding in new boxes. The
c     process involves flagging down bigger boxes and dividing
c     them and their children as necessary.
c
c     Input arguments
c     src         in: real *8(2,ns)
c                 x and y coordinates of the source locations
c    
c     ns          in: integer
c                 number of sources
c
c     radsrc      in: real *8(ns)
c                 radius of sources. A source i is hung
c                 at level ilev 
c                 if radsrc(i)>boxsize(ilev)
c
c     trg        in: real *8(2,nt)
c                 x and y coordinates of target locations
c
c     nt          in: integer
c                 number of targets
c
c     expc        in: real *8(2,nexpc)
c                 x and y coordinates of expansion centers
c
c     radexp      in: real *8(nexpc)
c                 radius of expansion centers. An expansion
c                 center i is hung at level ilev if
c                 radexp(i)>boxsize(ilev)
c
c     nlevels     in/out: integer
c                 current number of levels in the tree
c                 On output nlevels remains unchanged
c
c     nboxes      in/out: integer
c                 current number of boxes in the tree.
c                 On output it is the number of boxes
c                 in the level restricted tree
c
c     centers     in/out: real *8(2,nbmax)
c                 x and y coordinates of the box centers 
c                 in the tree
c
c     boxsize     in: real *8(0:nlevels)
c                 size of the box at any given level
c
c     nbmax       in: integer
c                 max number of boxes
c
c     laddr       in/out: integer(2,0:nlevels)
c                 laddr(1,i), laddr(2,i) are the first
c                 and the last box at level i
c
c     ilevel      in/out: integer(nbmax)
c                 ilevel(i) is the level of box i
c
c     iparent     in/out: integer(nbmax)
c                 iparent(i) is the parent of box i
c
c     nchild      in/out: integer(nbmax)
c                 nchild(i) is the number of children 
c                 of box i
c
c    ichild       in/out: integer(4,nbmax)
c                 ichild(j,i) is the jth child of box i
c
c    nnbors       in/out: integer(nbmax)
c                 nnbors(i) is the number of colleagues of box i
c
c    nbors        in/out: integer(9,nbmax)
c                 nbors(j,i) is the jth colleague of box i
c
c    isource      in/out: integer(ns)
c                 tree sorted array of sources
c
c    itarg        in/out: integer(nt)
c                 tree sorted array of targets
c
c    iexpc        in/out: integer(nexpc)
c                 tree sorted array of expansion centers
c
c    ihsfirst     in/out: integer(nbmax)
c                 ihsfirst(i) is the location in isource
c                 array for the first hung source in box i
c                 
c    ihslast      in/out: integer(nbmax)
c                 ihslast(i) is the location in isource
c                 array for the last hung source in box i
c                 
c    isfirst      in/out: integer(nbmax)
c                 isfirst(i) is the location in isource
c                 array for the first source in box i
c                 
c    islast       in/out: integer(nbmax)
c                 islast(i) is the location in isource
c                 array for the last source in box i
c                 
c    itfirst      in/out: integer(nbmax)
c                 itfirst(i) is the location in itarg
c                 array for the first target in box i
c                 
c    itlast       in/out: integer(nbmax)
c                 itlast(i) is the location in itarg
c                 array for the last target in box i
c                 
c    ihefirst     in/out: integer(nbmax)
c                 ihefirst(i) is the location in iexpc
c                 array for the first hung expansion center in box i
c                 
c    ihelast      in/out: integer(nbmax)
c                 ihelast(i) is the location in iexpc
c                 array for the last hung expansion center in box i
c                 
c    iefirst      in/out: integer(nbmax)
c                 iefirst(i) is the location in iexpc
c                 array for the first expansion center in box i
c                 
c    ielast       in/out: integer(nbmax)
c                 ielast(i) is the location in iexpc
c                 array for the last expansion center in box i
c                
c    nhungsrc     in/out: integer(nbmax)
c                 nhungsrc(i) is the number of hung sources in box i
c
c    nhungexp     in/out: integer(nbmax)
c                 nhungexp(i) is the number of hung 
c                 expansion centers in box i

      implicit none
c     Calling sequence variable declaration
      integer ns,nt,nexpc
      real *8 src(2,ns), trg(2,nt), expc(2,nexpc)
      real *8 radsrc(ns),radexp(nexpc)

      integer nlevels,nboxes,nbmax
      real *8 boxsize(0:nlevels), centers(2,nbmax)

      integer laddr(2,0:nlevels),ilevel(nbmax)
      integer iparent(nbmax)
      integer nchild(nbmax), ichild(4,nbmax)
      integer nnbors(nbmax), nbors(9,nbmax)
      integer isource(ns),itarget(nt),iexpc(nexpc)

      integer ihsfirst(nbmax), ihslast(nbmax)
      integer isfirst(nbmax), islast(nbmax)

      integer itfirst(nbmax), itlast(nbmax)

      integer ihefirst(nbmax), ihelast(nbmax)
      integer iefirst(nbmax), ielast(nbmax)
      
      integer nhungsrc(nbmax), nhungexp(nbmax)

      integer iflag(nbmax)
c     Temporary variables
      integer i,j,k,l
      integer ibox,jbox,kbox,ilev
      integer idad,igranddad
      real *8 xdis,ydis,distest

      integer laddrtail(2,0:nlevels)

c     Initialize flag array
      do i=1,nboxes
         iflag(i) = 0
      enddo

c     Flag boxes that violate level restriction by "1"
c     Violatioin refers to any box that is directly touching
c     a box that is more than one level finer
c
c     Method:
c     1) Carry out upward pass. For each box B, look at
c     the colleagues of B's grandparent
c     2) See if any of those colleagues are childless and in
c     contact with B.
c
c     Note that we only need to get up to level two, as
c     we will not find a violation at level 0 and level 1
c
c     For such boxes, we set iflag(i) = 1
c
      do ilev=nlevels,2,-1
c        This is the distance to test if two boxes separated
c        by two levels are touching
         distest = 1.05d0*(boxsize(ilev-1) + boxsize(ilev-2))/2.0d0
         do ibox = laddr(1,ilev),laddr(2,ilev) 
            idad = iparent(ibox)
            igranddad = iparent(idad)
            
c           Loop over colleagues of granddad            
            do i=1,nnbors(igranddad)
               jbox = nbors(i,igranddad)
c              Check if the colleague of grandad
c              is a leaf node. This automatically
c              eliminates the granddad
               if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                   xdis = centers(1,jbox) - centers(1,idad)
                   ydis = centers(2,jbox) - centers(2,idad)
                   if(abs(xdis).le.distest.and.abs(ydis).le.
     1                distest) then
                      iflag(jbox) = 1
                   endif
               endif
c              End of checking criteria for the colleague of
c              granddad
            enddo
c           End of looping over colleagues of
c           granddad
         enddo
c        End of looping over boxes at ilev         
      enddo
c     End of looping over levels and flagging boxes

c     Find all boxes that need to be given a flag+
c     A flag+ box will be denoted by setting iflag(box) = 2
c     This refers to any box that is not already flagged and
c     is bigger than and is contacting a flagged box
c     or another box that has already been given a flag +.
c     It is found by performing an upward pass and looking
c     at the flagged box's parents colleagues and a flag+
c     box's parents colleagues and seeing if they are
c     childless and present the case where a bigger box 
c     is contacting a flagged or flag+ box.

      do ilev = nlevels,2,-1
c        This is the distance to test if two boxes separated
c        by one level are touching
         distest = 1.05d0*(boxsize(ilev) + boxsize(ilev-1))/2.0d0
         do ibox = laddr(1,ilev),laddr(2,ilev)
            if(iflag(ibox).eq.1.or.iflag(ibox).eq.2) then
               idad = iparent(ibox)
c              Loop over dad's colleagues               
               do i=1,nnbors(idad)
                  jbox = nbors(i,idad)
c                 Check if the colleague of dad
c                 is a leaf node. This automatically
c                 eliminates the dad
                  if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                     xdis = centers(1,jbox) - centers(1,ibox)
                     ydis = centers(2,jbox) - centers(2,ibox)
                     if(abs(xdis).le.distest.and.abs(ydis).le.
     1                  distest) then
                        iflag(jbox) = 2
                      endif
                  endif
c                 End of checking criteria for the colleague of
c                dad
               enddo
c              End of looping over dad's colleagues               
            endif
c           End of checking if current box is relevant for
c           flagging flag+ boxes
         enddo
c        End of looping over boxes at ilev         
      enddo
c     End of looping over levels

c     Subdivide all flag and flag+ boxes. Flag all the children
c     of flagged boxes as flag++. Flag++ boxes are denoted
c     by setting iflag(box) = 3. The flag++ boxes need 
c     to be checked later to see which of them need further
c     refinement. While creating new boxes, we will
c     need to update all the tree structures as well.
c     Note that all the flagged boxes live between
c     levels 1 and nlevels - 2. We process the boxes via a
c     downward pass. We first determine the number of boxes
c     that are going to be subdivided at each level and 
c     everything else accordingly
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo
 
      do ilev = 1,nlevels-2
c        First subdivide all the flag and flag+
c        boxes with boxno nboxes+1, nboxes+ 2
c        and so on. In the second step, we reorganize
c        all the structures again to bring it back
c        in the standard format

         laddrtail(1,ilev+1) = nboxes+1

         call subdivide_flag(src,ns,radsrc,trg,nt,
     $                   expc,nexpc,radexp,
     $                   ilev,nboxes,
     $                   centers,boxsize,nbmax,nlevels,
     $                   laddr,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,iflag)
         laddrtail(2,ilev+1) = nboxes
      enddo
c     Reorganize the tree to get it back in the standard format

      call reorganizetree(nboxes,centers,nlevels,laddr,laddrtail,
     1                    ilevel,iparent,nchild,ichild,
     2                    ihsfirst,
     3                    ihslast,isfirst,islast,itfirst,itlast,
     4                    ihefirst,ihelast,iefirst,ielast,
     5                    nhungsrc,nhungexp,iflag)

c     Compute colleague information again      

      do i=1,nboxes
         nnbors(i) = 0
         do j=1,9
            nbors(j,i) = -1
         enddo
      enddo
      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,nnbors,nbors)

c     Processing of flag and flag+ boxes is done
c     Start processing flag++ boxes. We will use a similar
c     strategy as before. We keep checking the flag++
c     boxes that require subdivision if they still
c     violate the level restriction criterion, create
c     the new boxes, append them to the end of the list to begin
c     with and in the end reorganize the tree structure.
c     We shall accomplish this via a downward pass
c     as new boxes that get added in the downward pass
c     will also be processed simultaneously.
c     We shall additionally also need to keep on updating
c     the colleague information as we proceed in the 
c     downward pass

c     Reset the flags array to remove all the flag and flag+
c     cases. This is to ensure reusability of the subdivide
c     _flag routine to handle the flag++ case

      do ibox=1,nboxes
         if(iflag(ibox).ne.3) iflag(ibox) = 0
      enddo
 
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo

      do ilev = 3,nlevels-2

c     Step 1: Determine which of the flag++ boxes need
c     further division. In the even a flag++ box needs
c     further subdivision then flag the box with iflag(box) = 1
c     This will again ensure that the subdivide_flag routine
c     will take care of handling the flag++ case
         call updateflags(ilev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

         call updateflags(ilev,nboxes,nlevels,laddrtail,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)


         
c      Step 2: Subdivide all the boxes that need subdivision
c      in the laddr set and the laddrtail set as well
         laddrtail(1,ilev+1) = nboxes + 1
         call subdivide_flag(src,ns,radsrc,trg,nt,
     $                   expc,nexpc,radexp,
     $                   ilev,nboxes,
     $                   centers,boxsize,nbmax,nlevels,
     $                   laddr,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,iflag)

         call subdivide_flag(src,ns,radsrc,trg,nt,
     $                   expc,nexpc,radexp,
     $                   ilev,nboxes,
     $                   centers,boxsize,nbmax,nlevels,
     $                   laddrtail,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,iflag)
          laddrtail(2,ilev+1) = nboxes         
c      Step 3: Update the colleague information for the newly
c      created boxes
          do ibox = laddrtail(1,ilev+1),laddrtail(2,ilev+1)
            nnbors(ibox) = 0
c           Find the parent of the current box         
            idad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out colleagues
            do i=1,nnbors(idad)
                jbox = nbors(i,idad)
                do j=1,nchild(jbox)
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
c               Check if kbox is a nearest neighbor or in list 2
                   if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                1.05*boxsize(ilev+1)).and.
     2                (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                1.05*boxsize(ilev+1))) then
                     
                     nnbors(ibox) = nnbors(ibox)+1
                     nbors(nnbors(ibox),ibox) = kbox
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
      enddo

c     Reorganize tree once again and we are all done      
      call reorganizetree(nboxes,centers,nlevels,laddr,laddrtail,
     1                    ilevel,iparent,nchild,ichild,
     2                    ihsfirst,
     3                    ihslast,isfirst,islast,itfirst,itlast,
     4                    ihefirst,ihelast,iefirst,ielast,
     5                    nhungsrc,nhungexp,iflag)

c     Compute colleague information again      

      do i=1,nboxes
         nnbors(i) = 0
         do j=1,9
            nbors(j,i) = -1
         enddo
      enddo
      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,nnbors,nbors)
      
      return
      end

      subroutine subdivide_flag(src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,
     $                   curlev,nboxes,
     $                   centers,boxsize,nbmax,nlevels,
     $                   laddr,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,iflag)
      implicit none
      integer ns,nt,nexpc
      integer nlevels,nboxes,nbmax, curlev
      real *8 src(2,ns),radsrc(ns)
      real *8 trg(2,nt)
      real *8 expc(2,nexpc),radexp(nexpc)
      real *8 centers(2,nbmax)
      real *8 boxsize(0:nlevels)
      integer laddr(2,0:nlevels)
      integer ilevel(nbmax)
      integer iparent(nbmax)
      integer nchild(nbmax)
      integer ichild(4,nbmax)
      integer isource(ns)
      integer itarget(nt)
      integer iexpc(nexpc)
      integer ihsfirst(nbmax)
      integer ihslast(nbmax)
      integer isfirst(nbmax)
      integer islast(nbmax)
      integer itfirst(nbmax)
      integer itlast(nbmax)
      integer ihefirst(nbmax)
      integer ihelast(nbmax)
      integer iefirst(nbmax)
      integer ielast(nbmax)
      integer nhungsrc(nbmax)
      integer nhungexp(nbmax)
      integer iflag(nbmax)
c     Temporary variables
      integer isrctmp(ns),itargtmp(nt),iexpctmp(nexpc)
      integer i,j,i12,i34,istart,jstart,kstart,ii,iii,nss,nee,ntt
      integer ibox,ifirstbox,ilastbox,nbfirst
      integer is,it,ie
      integer nsc(4),ntc(4),nh(4),nexpcc(4),nhc(4)
c
c     for every box at level nlevels,
c     sort into children, updating various arrays 
c     perhaps just build tree here paren/child/particle sorting...
c     lists in second call ???
c     
c     allocate temp array for isourcetemp2 itargtemp2
c     after all done, write back to isourcetemp, itargtemp
c     this is O(N) * nlevels work for rewriting.
c     can be fancier I suppose.
c     
      ifirstbox = laddr(1,curlev)
      ilastbox =  laddr(2,curlev)

c
      do ibox = ifirstbox,ilastbox
c        The current box needs to be subdivided if iflag(ibox).gt.0       
         if(iflag(ibox).gt.0) then
c           Based on flagging criterion, the current
c           box needs to be divided. 
c
c           Allocate temporary array to figure out which child you
c     1      belong to
c           which child?  1,2,3,4? counter ns1,ns2,ns3,ns4
c           The box nomenclature is as follows
c           3   4
c           1   2
            i12 = isfirst(ibox)-1
            i34 = 0
            do is = isfirst(ibox),islast(ibox)
               if(src(2,isource(is)) - centers(2,ibox).lt.0) then
                  i12 = i12+1
                  isource(i12) = isource(is)
               else
                  i34 = i34 + 1
                  isrctmp(i34) = isource(is)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end

c           Reorder sources to include sources in 34 in the array
            do i=1,i34
               isource(i12+i) = isrctmp(i)
            enddo
            nsc(1) = 0
            nsc(2) = 0
            nsc(3) = 0
            nsc(4) = 0
c           Sort into boxes 1 and 2
            do is = isfirst(ibox),i12
               if(src(1,isource(is))-centers(1,ibox).lt.0) then
                  isource(isfirst(ibox)+nsc(1)) = isource(is)
                  nsc(1) = nsc(1) + 1
               else
                  nsc(2) = nsc(2) + 1
                  isrctmp(nsc(2)) = isource(is)
               endif
            enddo
c           Reorder sources so that sources in 2 are at the
c           end of the array
            do i=1,nsc(2)
               isource(isfirst(ibox)+nsc(1)+i-1) = isrctmp(i)
            enddo
c           Sort into boxes 3 and 4
            do is = i12+1, islast(ibox)
               if(src(1,isource(is))-centers(1,ibox).lt.0) then
                   isource(i12+1+nsc(3)) = isource(is)
                   nsc(3) = nsc(3) + 1
                else
                   nsc(4) = nsc(4)+1
                   isrctmp(nsc(4)) = isource(is)
                endif
            enddo
            do i=1,nsc(4)
               isource(i12+nsc(3)+i) = isrctmp(i)
            enddo
c           End of sorting sources

            istart = isfirst(ibox)-1
            do j=1,4
c           check hung -> counter nh1,nh2,nh3,nh4
               ii = 0
               nh(j) = 0
               do i=1,nsc(j)
                  if(radsrc(isource(istart+i)).gt.boxsize(curlev+1))
     1            then     
                     nh(j) = nh(j) + 1
                     isource(istart+nh(j)) = isource(istart+i)
                   else
                      ii = ii+1
                      isrctmp(ii) = isource(istart+i)
                   endif
                enddo
c            Reorder sources to have hung chunks at the star
c            of the sorted sources in the box ibox
                do i=1,ii
                   isource(istart+nh(j)+i) = isrctmp(i)
                enddo
                istart = istart + nsc(j)
            enddo
         
c           which child?  1,2,3,4? counter nt(1),nt(2),ns(3),ns(4)
c           The box nomenclature is as follows
c           3   4
c           1   2
            i12 = itfirst(ibox)-1
            i34 = 0
            do it = itfirst(ibox),itlast(ibox)
               if(trg(2,itarget(it)) - centers(2,ibox).lt.0) then
                  i12 = i12+1
                  itarget(i12) = itarget(it)
               else
                  i34 = i34 + 1
                  itargtmp(i34) = itarget(it)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end

c           Reorder sources to include sources in 34 in the array
            do i=1,i34
               itarget(i12+i) = itargtmp(i)
            enddo
            ntc(1) = 0
            ntc(2) = 0
            ntc(3) = 0
            ntc(4) = 0
c           Sort into boxes 1 and 2
            do it = itfirst(ibox),i12
               if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
                  itarget(itfirst(ibox)+ntc(1)) = itarget(it)
                  ntc(1) = ntc(1) + 1
               else
                  ntc(2) = ntc(2) + 1
                  itargtmp(ntc(2)) = itarget(it)
               endif
            enddo
c           Reorder sources so that sources in 2 are at the
c           end of the array
            do i=1,ntc(2)
               itarget(itfirst(ibox)+ntc(1)+i-1) = itargtmp(i)
            enddo
c           Sort into boxes 3 and 4
            do it = i12+1, itlast(ibox)
               if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
                  itarget(i12+1+ntc(3)) = itarget(it)
                  ntc(3) = ntc(3) + 1
               else
                  ntc(4) = ntc(4)+1
                  itargtmp(ntc(4)) = itarget(it)
               endif
            enddo
            do i=1,ntc(4)
               itarget(i12+ntc(3)+i) = itargtmp(i)
            enddo
c           End of sorting targets

c           Sort expansion centers
c           which child?  1,2,3,4? counter nt(1),nt(2),ns(3),ns(4)
c           The box nomenclature is as follows
c           3   4
c           1   2
            i12 = iefirst(ibox)-1
            i34 = 0
            do ie = iefirst(ibox),ielast(ibox)
               if(expc(2,iexpc(ie)) - centers(2,ibox).lt.0) then
                  i12 = i12+1
                  iexpc(i12) = iexpc(ie)
               else
                  i34 = i34 + 1
                  iexpctmp(i34) = iexpc(ie)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end

c           Reorder sources to include sources in 34 in the array
            do i=1,i34
               iexpc(i12+i) = iexpctmp(i)
            enddo
            nexpcc(1) = 0
            nexpcc(2) = 0
            nexpcc(3) = 0
            nexpcc(4) = 0
c           Sort into boxes 1 and 2
            do ie = iefirst(ibox),i12
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(iefirst(ibox)+nexpcc(1)) = iexpc(ie)
                  nexpcc(1) = nexpcc(1) + 1
               else
                  nexpcc(2) = nexpcc(2) + 1
                  iexpctmp(nexpcc(2)) = iexpc(ie)
               endif
            enddo
c           Reorder sources so that sources in 2 are at the
c           end of the array
            do i=1,nexpcc(2)
               iexpc(iefirst(ibox)+nexpcc(1)+i-1) = iexpctmp(i)
            enddo
c           Sort into boxes 3 and 4
            do ie = i12+1, ielast(ibox)
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(i12+1+nexpcc(3)) = iexpc(ie)
                  nexpcc(3) = nexpcc(3) + 1
               else
                  nexpcc(4) = nexpcc(4)+1
                  iexpctmp(nexpcc(4)) = iexpc(ie)
               endif
            enddo
            do i=1,nexpcc(4)
               iexpc(i12+nexpcc(3)+i) = iexpctmp(i)
            enddo
c           End of sorting expanison centers

            istart = iefirst(ibox)-1
            do j=1,4
c           check hung -> counter nhc1,nhc2,nhc3,nhc4
               ii = 0
               nhc(j) = 0
               do i=1,nexpcc(j)
                  if(radexp(iexpc(istart+i)).gt.boxsize(curlev+1))
     1            then     
                     nhc(j) = nhc(j) + 1
                     iexpc(istart+nhc(j)) = iexpc(istart+i)
                   else
                      ii = ii+1
                      iexpctmp(ii) = iexpc(istart+i)
                   endif
                enddo
c           Reorder sources to have hung chunks at the star
c           of the sorted sources in the box ibox
               do i=1,ii
                  iexpc(istart+nhc(j)+i) = iexpctmp(i)
               enddo
               istart = istart + nexpcc(j)
            enddo

            nchild(ibox) = 0
c           Create the required boxes
            istart = isfirst(ibox)
            jstart = itfirst(ibox)
            kstart = iefirst(ibox)
            do i=1,4
               ii = 2
               if(i.eq.1.or.i.eq.2) ii = 1
               if(nsc(i)+ntc(i)+nexpcc(i).ge.0) then
c                 Increment total number of boxes               
                  nboxes = nboxes + 1
c                 Increment number of children for the current box
                  nchild(ibox) = nchild(ibox)+1
c                 Update the array of children for the current box
                  ichild(nchild(ibox),ibox) = nboxes
c                 Update the array of levels for the child box
                  ilevel(nboxes) = curlev+1
c                 Update the array of parents for the child box
                  iparent(nboxes) = ibox
c                 Compute center for the child box
                  centers(1,nboxes) = centers(1,ibox)+(-1)**i*
     1                                boxsize(curlev+1)/2.0
                  centers(2,nboxes) = centers(2,ibox)+(-1)**ii*
     1                                boxsize(curlev+1)/2.0
c                 Update arrays ihsfirst,ihslast,isfirst,islast
                  ihsfirst(nboxes) = istart
                  ihslast(nboxes) = istart + nh(i) - 1
                  nhungsrc(nboxes) = nh(i)

                  isfirst(nboxes) = istart + nh(i)
                  islast(nboxes) = istart + nsc(i) - 1

c                 Update arrays itfirst, itlast
                  itfirst(nboxes) = jstart
                  itlast(nboxes) = jstart + ntc(i) - 1

c                 Update arrays ihefirst,ihelast,iefirst,ielast
                  ihefirst(nboxes) = kstart
                  ihelast(nboxes) = kstart + nhc(i)-1
                  nhungexp(nboxes) = nhc(i)

                  iefirst(nboxes) = kstart + nhc(i)
                  ielast(nboxes) = kstart + nexpcc(i) - 1

                  nchild(nboxes) = 0
                  ichild(1,nboxes) = -1
                  ichild(2,nboxes) = -1
                  ichild(3,nboxes) = -1
                  ichild(4,nboxes) = -1
                  if(iflag(ibox).eq.1) iflag(nboxes) = 3
                  if(iflag(ibox).eq.2) iflag(nboxes) = 0
                endif
                istart = istart + nsc(i)
                jstart = jstart + ntc(i)
                kstart = kstart + nexpcc(i)
            enddo
         endif
      enddo

      return
      end
c-------------------------------------------------------------      
      subroutine reorganizetree(nboxes,centers,nlevels,laddr,laddrtail,
     1                    ilevel,iparent,nchild,ichild,
     2                    ihsfirst,
     3                    ihslast,isfirst,islast,itfirst,itlast,
     4                    ihefirst,ihelast,iefirst,ielast,
     5                    nhungsrc,nhungexp,iflag)

c    This subroutine reorganizes the current data in all the tree
c    arrays to rearrange them in the standard format.
c    The boxes on input are assumed to be arranged in the following
c    format
c    boxes on level i are the boxes from laddr(1,i) to 
c    laddr(2,i) and also from laddrtail(1,i) to laddrtail(2,i)
c
c    At the end of the sorting, the boxes on level i
c    are arranged from laddr(1,i) to laddr(2,i)  
c
c    INPUT/OUTPUT arguments
c    nboxes         in: integer
c                   number of boxes
c
c    centers        in/out: real *8(2,nboxes)
c                   x and y coordinates of the center of boxes
c
c    nlevels        in: integer
c                   Number of levels in the tree
c
c    laddr          in/out: integer(2,0:nlevels)
c                   boxes at level i are numbered between
c                   laddr(1,i) to laddr(2,i)
c
c    laddrtail      in: integer(2,0:nlevels)
c                   new boxes to be added to the tree
c                   structure are numbered from
c                   laddrtail(1,i) to laddrtail(2,i)
c
c     ilevel      in/out: integer(nboxes)
c                 ilevel(i) is the level of box i
c
c     iparent     in/out: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c     nchild      in/out: integer(nboxes)
c                 nchild(i) is the number of children 
c                 of box i
c
c    ichild       in/out: integer(4,nboxes)
c                 ichild(j,i) is the jth child of box i
c
c    nnbors       in/out: integer(nboxes)
c                 nnbors(i) is the number of colleagues of box i
c
c    nbors        in/out: integer(9,nboxes)
c                 nbors(j,i) is the jth colleague of box i
c
c    ihsfirst     in/out: integer(nboxes)
c                 ihsfirst(i) is the location in isource
c                 array for the first hung source in box i
c                 
c    ihslast      in/out: integer(nboxes)
c                 ihslast(i) is the location in isource
c                 array for the last hung source in box i
c                 
c    isfirst      in/out: integer(nboxes)
c                 isfirst(i) is the location in isource
c                 array for the first source in box i
c                 
c    islast       in/out: integer(nboxes)
c                 islast(i) is the location in isource
c                 array for the last source in box i
c                 
c    itfirst      in/out: integer(nboxes)
c                 itfirst(i) is the location in itarg
c                 array for the first target in box i
c                 
c    itlast       in/out: integer(nboxes)
c                 itlast(i) is the location in itarg
c                 array for the last target in box i
c                 
c    ihefirst     in/out: integer(nboxes)
c                 ihefirst(i) is the location in iexpc
c                 array for the first hung expansion center in box i
c                 
c    ihelast      in/out: integer(nboxes)
c                 ihelast(i) is the location in iexpc
c                 array for the last hung expansion center in box i
c                 
c    iefirst      in/out: integer(nboxes)
c                 iefirst(i) is the location in iexpc
c                 array for the first expansion center in box i
c                 
c    ielast       in/out: integer(nboxes)
c                 ielast(i) is the location in iexpc
c                 array for the last expansion center in box i
c                
c    nhungsrc     in/out: integer(nboxes)
c                 nhungsrc(i) is the number of hung sources in box i
c
c    nhungexp     in/out: integer(nboxes)
c                 nhungexp(i) is the number of hung 
c                 expansion centers in box i
c
c    iflag        in/out: integer(nboxes)
c                 iflag(i) is a flag for box i required to generate
c                 level restricted tree from adaptive tree

      implicit none
c     Calling sequence variables and temporary variables
      integer nboxes,nlevels
      real *8 centers(2,nboxes), tcenters(2,nboxes)
      integer laddr(2,0:nlevels), tladdr(2,0:nlevels)
      integer laddrtail(2,0:nlevels)
      integer ilevel(nboxes), tilevel(nboxes)
      integer iparent(nboxes), tiparent(nboxes)
      integer nchild(nboxes), tnchild(nboxes)
      integer ichild(4,nboxes), tichild(4,nboxes)
      integer ihsfirst(nboxes), tihsfirst(nboxes)
      integer ihslast(nboxes), tihslast(nboxes)
      integer isfirst(nboxes), tisfirst(nboxes)
      integer islast(nboxes), tislast(nboxes)
      integer itfirst(nboxes), titfirst(nboxes)
      integer itlast(nboxes), titlast(nboxes)
      integer ihefirst(nboxes), tihefirst(nboxes)
      integer ihelast(nboxes), tihelast(nboxes)
      integer iefirst(nboxes), tiefirst(nboxes)
      integer ielast(nboxes), tielast(nboxes)
      integer nhungsrc(nboxes), tnhungsrc(nboxes)
      integer nhungexp(nboxes), tnhungexp(nboxes)
      integer iflag(nboxes),tiflag(nboxes)
      integer iboxtocurbox(nboxes)

c     Temporary variables
      integer i,j,k,l
      integer ibox,ilev, curbox

      do ilev = 0,nlevels
         tladdr(1,ilev) = laddr(1,ilev)
         tladdr(2,ilev) = laddr(2,ilev)
      enddo

      do ibox=1,nboxes
         tilevel(ibox) = ilevel(ibox)
         tcenters(1,ibox) = centers(1,ibox)
         tcenters(2,ibox) = centers(2,ibox)
         tiparent(ibox) = iparent(ibox)
         tnchild(ibox) = nchild(ibox)
         tichild(1,ibox) = ichild(1,ibox)
         tichild(2,ibox) = ichild(2,ibox)
         tichild(3,ibox) = ichild(3,ibox)
         tichild(4,ibox) = ichild(4,ibox)
         tihsfirst(ibox) = ihsfirst(ibox)
         tihslast(ibox) = ihslast(ibox)
         tisfirst(ibox) = isfirst(ibox)
         tislast(ibox) = islast(ibox)
         titfirst(ibox) = itfirst(ibox)
         titlast(ibox) = itlast(ibox)
         tihefirst(ibox) = ihefirst(ibox)
         tihelast(ibox) = ihelast(ibox)
         tiefirst(ibox) = iefirst(ibox)
         tielast(ibox) = ielast(ibox)
         tnhungsrc(ibox) = nhungsrc(ibox)
         tnhungexp(ibox) = nhungexp(ibox)
         tiflag(ibox) = iflag(ibox)
      enddo
     
c     Rearrange old arrays now

      do ilev = 0,1
         do ibox = laddr(1,ilev),laddr(2,ilev)
            iboxtocurbox(ibox) = ibox
         enddo
      enddo
      if(nlevels.ge.2) curbox = laddr(1,2)
      do ilev=2,nlevels
         laddr(1,ilev) = curbox
         do ibox = tladdr(1,ilev),tladdr(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            nchild(curbox) = tnchild(ibox)
            centers(1,curbox) = tcenters(1,ibox)
            centers(2,curbox) = tcenters(2,ibox)
            ihsfirst(curbox) = tihsfirst(ibox)
            ihslast(curbox) = tihslast(ibox)
            isfirst(curbox) = tisfirst(ibox)
            islast(curbox) = tislast(ibox)
            itfirst(curbox) = titfirst(ibox)
            itlast(curbox) = titlast(ibox)
            ihefirst(curbox) = tihefirst(ibox)
            ihelast(curbox) = tihelast(ibox)
            iefirst(curbox) = tiefirst(ibox)
            ielast(curbox) = tielast(ibox)
            nhungsrc(curbox) = tnhungsrc(ibox)
            nhungexp(curbox) = tnhungexp(ibox)
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         do ibox = laddrtail(1,ilev),laddrtail(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            centers(1,curbox) = tcenters(1,ibox)
            centers(2,curbox) = tcenters(2,ibox)
            nchild(curbox) = tnchild(ibox)
            ihsfirst(curbox) = tihsfirst(ibox)
            ihslast(curbox) = tihslast(ibox)
            isfirst(curbox) = tisfirst(ibox)
            islast(curbox) = tislast(ibox)
            itfirst(curbox) = titfirst(ibox)
            itlast(curbox) = titlast(ibox)
            ihefirst(curbox) = tihefirst(ibox)
            ihelast(curbox) = tihelast(ibox)
            iefirst(curbox) = tiefirst(ibox)
            ielast(curbox) = tielast(ibox)
            nhungsrc(curbox) = tnhungsrc(ibox)
            nhungexp(curbox) = tnhungexp(ibox)
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         laddr(2,ilev) = curbox-1
      enddo

c     Handle the parent children part of the tree 
c     using the mapping iboxtocurbox

      do ibox=1,nboxes
         if(tiparent(ibox).eq.-1) iparent(iboxtocurbox(ibox)) = -1
         if(tiparent(ibox).gt.0) 
     1    iparent(iboxtocurbox(ibox)) = iboxtocurbox(tiparent(ibox))
         if(tichild(1,ibox).eq.-1) ichild(1,iboxtocurbox(ibox)) = -1
         if(tichild(1,ibox).gt.0) 
     1    ichild(1,iboxtocurbox(ibox)) = iboxtocurbox(tichild(1,ibox))
         
         if(tichild(2,ibox).eq.-1) ichild(2,iboxtocurbox(ibox)) = -1
         if(tichild(2,ibox).gt.0) 
     1    ichild(2,iboxtocurbox(ibox)) = iboxtocurbox(tichild(2,ibox))
         if(tichild(3,ibox).eq.-1) ichild(3,iboxtocurbox(ibox)) = -1
         if(tichild(3,ibox).gt.0) 
     1    ichild(3,iboxtocurbox(ibox)) = iboxtocurbox(tichild(3,ibox))
         if(tichild(4,ibox).eq.-1) ichild(4,iboxtocurbox(ibox)) = -1
         if(tichild(4,ibox).gt.0) 
     1    ichild(4,iboxtocurbox(ibox)) = iboxtocurbox(tichild(4,ibox))
      enddo

      return
      end
c--------------------------------------------------------------------
      subroutine updateflags(curlev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

c      This subroutine is to check the boxes flagged as flag++
c      and determine which of the boxes need refinement. The flag
c      of the box which need refinement is updated to iflag(box)=1
c      and that of the boxes which do not need refinement is
c      updated to iflag(box) = 0
c
c      INPUT arguments
c      curlev         in: integer
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer
c                     total number of boxes
c
c      nlevels        in: integer
c                     total number of levels
c
c      laddr          in: integer(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer(4,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(9,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: real *8(2,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: real *8(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer curlev, nboxes, nlevels
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(4,nboxes)
      integer nnbors(nboxes), nbors(9,nboxes)
      integer iflag(nboxes)
      real *8 centers(2,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox
      real *8 distest,xdis,ydis

      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level      
      do ibox = laddr(1,curlev),laddr(2,curlev)
         if(iflag(ibox).eq.3) then
            iflag(ibox) = 0
c           Loop over colleagues of the current box      
            do i=1,nnbors(ibox)
c              Loop over colleagues of flag++ box        
               jbox = nbors(i,ibox)
c              Loop over the children of the colleague box
c              Note we do not need to exclude self from
c              the list of colleagues as a self box which
c              is flag++ does not have any children 
c              and will not enter the next loop
               do j=1,nchild(jbox)
                  kbox = ichild(j,jbox)
c                 Check if neighboring box has children
c                 Note that all flag++ boxes are childless,
c                 so this automatically eliminates self.
                  if(nchild(kbox).gt.0) then
c                    Check to see if the kbox touches ibox 
                     xdis = centers(1,ibox)-centers(1,kbox)
                     ydis = centers(2,ibox)-centers(2,kbox)
                     if(abs(xdis).le.distest.and.abs(ydis).le.
     1               distest) then
                         iflag(ibox) = 1
                         goto 1111
                      endif
                  endif
c                 End of looping over the children of the child
c                 of the colleague box
               enddo
c              End of looping over the children of the colleague box       
            enddo
c           End of looping over colleagues            
 1111       continue        
         endif
c        End of testing if the current box needs to checked for         
      enddo
c     End of looping over boxes at the current level      

      return
      end
c-----------------------------------------------------------------

      subroutine computelists(nlevels,nboxes,laddr,boxsize,
     1            centers,iparent,nchild,ichild,nnbors,nbors,
     2            mnlist1,nlist1,list1,
     3            mnlist2,nlist2,list2,
     4            mnlist3,nlist3,list3,
     5            mnlist4,nlist4,list4)
c
c     This subroutine computes the list1 and list2 for the tree
c     structure. For more details on the description
c     of the lists refer to the documentation of the main
c     routine maketree
c
c     
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: real *8(0:;nlevels)
c                 Array of boxsizes
c 
c     centers     in: real *8(2,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(4,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c     
c     nnbors      in: integer(nboxes)
c                 nnbors(i) is the number of colleagues of box i
c
c     nbors       in: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c
c     mnlist1     in: integer
c                 max number of boxes in list 1 of a box
c
c     mnlist2     in: integer
c                 max number of boxes in list 2 of a box
c
c     mnlist3     in: integer
c                 max number of boxes in list 3 of a box
c
c     mnlist4     in: integer
c                 max numbwe of boxes in list 4 of a box
c--------------------------------------------------------------
c     OUTPUT arguments:
c     nlist1      out: integer(nboxes)
c                 nlist1(i) is the number of boxes in list 1 
c                 of box i
c
c     list1       out: integer(mnlist1,nboxes)
c                 list1(j,i) is the box id of the jth box in 
c                 list1 of box i
c                      
c     nlist2      out: integer(nboxes)
c                 nlist2(i) is the number of boxes in the list 2
c                 of box i
c 
c     list2       out: integer(mnlist2,nboxes)
c                 list2(j,i) is the box id of the jth box in 
c                 list2 of box i
c
c     nlist3      out: integer(nboxes)
c                 nlist3(i) is the number of boxes in list 3
c                 of box i
c
c     list3       out: integer(mnlist3,nboxes)
c                 list3(j,i) is the box id of the jth box in 
c                 list3 of box i
c                      
c     nlist4      out: integer(nboxes)
c                 nlist4(i) is the number of boxes in the list 2
c                 of box i
c 
c     list4       out: integer(mnlist4,nboxes)
c                 list4(j,i) is the box id of the jth box in 
c                 list4 of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes
      integer laddr(2,0:nlevels)
      real *8 boxsize(0:nlevels)
      real *8 centers(2,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(4,nboxes)
      integer nnbors(nboxes)
      integer nbors(9,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer nlist1(nboxes), list1(mnlist1,nboxes)
      integer nlist2(nboxes), list2(mnlist2,nboxes)
      integer nlist3(nboxes), list3(mnlist3,nboxes)
      integer nlist4(nboxes), list4(mnlist4,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox
      real *8 distest,xdis,ydis


c     Setting parameters for level = 0
      if(nchild(1).eq.0) then
         nlist1(1) = 1 
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif
      nlist2(1) = 0
      nlist3(1) = 0
      nlist4(1) = 0
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,nchild(jbox)
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if((abs(centers(1,kbox)-centers(1,ibox)).ge.
     1                1.05*boxsize(ilev)).or.
     2                (abs(centers(2,kbox)-centers(2,ibox)).ge.
     3                1.05*boxsize(ilev))) then
                     
                     nlist2(ibox) = nlist2(ibox)+1
                     list2(nlist2(ibox),ibox) = kbox
                  endif
                enddo
            enddo

c           Compute list 1 and list 3 of ibox if ibox is childless
            if(nchild(ibox).eq.0) then
c              Loop over all colleagues of ibox              
               do i=1,nnbors(ibox)
                  jbox = nbors(i,ibox)
c                 If the colleague box is childless, then
c                 colleague box is in list 1
                  if(nchild(jbox).eq.0) then
                     nlist1(ibox) = nlist1(ibox)+1
                     list1(nlist1(ibox),ibox) = jbox

c                 If colleague box is not childless, then
c                 test to see if children of colleague
c                 box are in list1 or list 3. Since
c                 the tree is level restricted, we do not
c                 need to go any level deeper to determine
c                 list 1 and list 3
                  else
                     distest = 1.05d0*(boxsize(ilev) + 
     1                                 boxsize(ilev+1))/2.0d0
c                    Loop over children of colleague box              
                     do j=1,nchild(jbox)
                        kbox = ichild(j,jbox)
                        xdis = dabs(centers(1,kbox)-centers(1,ibox))
                        ydis = dabs(centers(2,kbox)-centers(2,ibox))
c                       Test to see if child of colleague box
c                       is in list1
                        if(xdis.lt.distest.and.ydis.lt.distest) then
                           nlist1(ibox) = nlist1(ibox)+1
                           list1(nlist1(ibox),ibox)=kbox

                           nlist1(kbox) = nlist1(kbox)+1
                           list1(nlist1(kbox),kbox) = ibox

c                       If it is not in list 1 of ibox then it 
c                       is in list3
                        else
                           nlist3(ibox) = nlist3(ibox)+1
                           list3(nlist3(ibox),ibox)=kbox

                           nlist4(kbox) = nlist4(kbox)+1
                           list4(nlist4(kbox),kbox)=ibox
                        endif
c                       End of figuring out whether child 
c                       of colleague box is in list 1 or list3
                     enddo
c                    End of looping over of children of colleague
c                    box
                  endif
c                 End of checking if colleague box has children
               enddo
c              End of looping over colleague boxes
            endif 
c           End of checking of current box is childless
         enddo
c        End of looping over boxes at level ilev         
      enddo
c     End of looping over levels      

      return
      end
c------------------------------------------------------------------

      subroutine computemhung(nlevels,nboxes,laddr,iparent,nchild,
     1                        nnbors,nbors,mnlist1,nlist1,list1,
     2                        nhungsrc,nhunglistsrc,mhung)    
c     This subroutine computes mhung for a given tree 
c     and the number of hung sources per box
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c    
c     iparent     in: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     nnbors      in: integer(nboxes)
c                 nnbors(i) is the number of boxes in list 1 of box i
c
c     nbors       in: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth box in
c                 list 1 of box i
c
c     mnlist1     in: integer
c                 max number of boxes in list 1 of a box
c
c     nlist1      in: integer(nboxes)
c                 nlist1(i) is the number of boxes in list1
c                 of box i
c
c     list1       in: integer(mnlist1,nboxes)
c                 list1(j,i) is the id of the jth box
c                 in list 1 of box i
c
c     nhungsrc   in: integer(nboxes)
c                 nhung(i) is the number of hung sources in box i
c 
c      OUTPUT
c      nhunglistsrc  out:integer(nboxes)
c                    nhunglistsrc(i) is the number of hung sources
c                    relevant to box i, nhunglistsrc(i) = 
c                    nhunglistsrc(idad) + \sum_{j=1}{nnbors(i)}
c                    nhungsrc(nbors(j,i))
c
c      mhung         out: integer
c                    max(nhunglistsrc)
      implicit none
      integer nlevels,nboxes
      integer laddr(2,0:nlevels)
      integer iparent(nboxes),nnbors(nboxes), nbors(9,nboxes)
      integer nchild(nboxes)
      integer mnlist1,mnlist4
      integer nlist1(nboxes),list1(mnlist1,nboxes)
      integer nhungsrc(nboxes), nhunglistsrc(nboxes)
      integer ilevel(nboxes)
      integer mhung

c     Temporary variables
      integer ilev,ibox,i,dad,jbox

c     initialize for level 0      
      
      do ilev = 0,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
            ilevel(ibox) = ilev
         enddo
      enddo
 
      nhunglistsrc(1) = nhungsrc(1)
      do ilev = 1,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
            dad = iparent(ibox)
            nhunglistsrc(ibox) = nhunglistsrc(dad)
            do i=1,nlist1(ibox)
               jbox = list1(i,ibox)
               if(ilevel(ibox).lt.ilevel(jbox)) then
                  nhunglistsrc(ibox) = nhunglistsrc(ibox) + 
     1                                nhungsrc(jbox)
               endif
            enddo

            do i=1,nnbors(ibox)
               jbox = nbors(i,ibox)
               nhunglistsrc(ibox) = nhunglistsrc(ibox) +
     1                                 nhungsrc(jbox)
            enddo
         enddo
      enddo

      mhung = 0
      do i=1,nboxes
         if(nhunglistsrc(i).gt.mhung) mhung = nhunglistsrc(i)
      enddo

      return
      end
c------------------------------------------------------------------
      subroutine computehunglist(mhung,nlevels,nboxes,laddr,ns,
     1           isource,iparent,nchild,ihsfirst,ihslast,nnbors,nbors,
     2           mnlist1,nlist1,list1,
     3           nhungsrc,nhunglistsrc,ihunglistsrc)
c     This subroutine computes the hung list for each box 
c     and stores it. The hung list of sources of a box is
c     the hunglist of the parent + the sources hung in colleagues
c     of the box + the sources hung in list 4 of the box
c 
c     INPUT arguments
c     mhung         in: Integer
c                   max number of hung chunks for a box
c
c     nlevels       in: Integer
c                   number of levels in the tree
c
c     nboxes        in: Integer
c                   number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     ns          in: integer
c                 number of sources
c
c     isource     in: integer(ns)
c                 Mapping from tree sorted array of sources
c                 to the user prescribed ordering
c  
c     iparent     in: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c     
c     ihsfirst    in: integer(nboxes)
c                 ihsfirst(i) points to the first hung source
c                 in box i in the sorted array
c
c      ihslast    in: integer(nboxes)
c                 ihslast(i) points to the last hung source
c                 in box i in the sorted array
c
c     nnbors      in: integer(nboxes)
c                 nnbors(i) is the number of boxes in list 1 of box i
c
c     nbors       in: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth box in
c                 list 1 of box i
c
c     mnlist1     in: integer
c                 max number of boxes in list 1 of a box
c
c     nlist1      in: integer(nboxes)
c                 nlist1(i) is the number of boxes in list1
c                 of box i
c
c     list1       in: integer(mnlist1,nboxes)
c                 list1(j,i) is the id of the jth box
c                 in list 1 of box i
c
c     nhungsrc   in: integer(nboxes)
c                 nhung(i) is the number of hung sources in box i
c 
c      OUTPUT
c      nhunglistsrc  out:integer(nboxes)
c                    nhunglistsrc(i) is the number of hung sources
c                    relevant to box i 
c
c      ihunglistsrc  out: integer(mhung,nboxes)
c                    ihunglistsrc(j,i) is the id of the jth
c                    hung source relevant to box i
c-----------------------------------------------------------------
       implicit none
       integer mhung,nlevels,nboxes,ns
       integer laddr(2,0:nlevels)
       integer iparent(nboxes), nchild(nboxes)
       integer isource(ns), ihsfirst(nboxes),ihslast(nboxes)
       integer nnbors(nboxes), nbors(9,nboxes)
       integer mnlist1, nlist1(nboxes), list1(mnlist1,nboxes)
       integer nhungsrc(nboxes),nhunglistsrc(nboxes)
       integer ilevel(nboxes)
       integer ihunglistsrc(mhung,nboxes)
c      Temp variables
       integer i,j,ibox,jbox,ilev,dad


       do ilev = 0,nlevels
          do ibox = laddr(1,ilev),laddr(2,ilev)
             ilevel(ibox) = ilev
          enddo
       enddo

c      Compute hung list for root node       
       nhunglistsrc(1) = nhungsrc(1)
       do i=1,nhungsrc(1)
          ihunglistsrc(i,1) = isource(ihsfirst(1)+i-1)
       enddo

       do ilev=1,nlevels
          do ibox = laddr(1,ilev),laddr(2,ilev)
             dad = iparent(ibox)
             do i=1,nhunglistsrc(dad)
                ihunglistsrc(i,ibox) = ihunglistsrc(i,dad)
             enddo
             nhunglistsrc(ibox) = nhunglistsrc(dad)
             do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                if(ilevel(ibox).lt.ilevel(jbox)) then
                   do j=1,nhungsrc(jbox)
                      ihunglistsrc(nhunglistsrc(ibox)+j,ibox)=
     1                isource(ihsfirst(jbox)+j-1)
                   enddo
                   nhunglistsrc(ibox) = nhunglistsrc(ibox)+
     1             nhungsrc(jbox)
                endif
             enddo

             do i=1,nnbors(ibox)
                jbox = nbors(i,ibox)
                do j=1,nhungsrc(jbox)
                   ihunglistsrc(nhunglistsrc(ibox)+j,ibox)=
     1                     isource(ihsfirst(jbox)+j-1)
                enddo
                nhunglistsrc(ibox) = nhunglistsrc(ibox)+
     1                  nhungsrc(jbox)
             enddo
          enddo
       enddo

       return
       end
c------------------------------------------------------------------
      subroutine maketree2dmem(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,nlmax,nbmax,nlevels,
     $                   nboxes,mhung,ltree)   
      implicit none
      integer ier,ns,nt,nexpc,idivflag,ndiv,isep,mhung
      integer nlevels,lcenters,ltree
      integer nlmax,nbmax,nboxes,nlevtmp,nbtmp, mhungtmp
      integer ntot
      integer i,j
      integer, allocatable :: laddr(:,:)
      integer, allocatable :: ilevel(:),ilevel2(:)
      integer, allocatable :: iparenttemp(:),iparenttemp2(:)
      integer, allocatable :: nchild(:),nchild2(:)
      integer, allocatable :: ichildtemp(:,:),ichildtemp2(:,:)
      integer, allocatable :: nnbors(:)
      integer, allocatable :: nbors(:,:)
      integer, allocatable :: isourcetemp(:)
      integer, allocatable :: itargettemp(:)
      integer, allocatable :: iexpctemp(:)
      integer, allocatable :: ihsfirsttemp(:),ihsfirsttemp2(:)
      integer, allocatable :: ihslasttemp(:),ihslasttemp2(:)
      integer, allocatable :: isfirsttemp(:),isfirsttemp2(:)
      integer, allocatable :: islasttemp(:),islasttemp2(:)
      integer, allocatable :: itfirsttemp(:),itfirsttemp2(:)
      integer, allocatable :: itlasttemp(:),itlasttemp2(:)
      integer, allocatable :: ihefirsttemp(:),ihefirsttemp2(:)
      integer, allocatable :: ihelasttemp(:),ihelasttemp2(:)
      integer, allocatable :: iefirsttemp(:),iefirsttemp2(:)
      integer, allocatable :: ielasttemp(:),ielasttemp2(:)
      integer, allocatable :: nlist1(:)
      integer, allocatable :: list1(:,:)
      integer, allocatable :: nlist2(:)
      integer, allocatable :: list2(:,:)
      integer, allocatable :: nlist3(:)
      integer, allocatable :: list3(:,:)
      integer, allocatable :: nlist4(:)
      integer, allocatable :: list4(:,:)
      integer, allocatable :: nhungsrc(:),nhungsrc2(:)
      integer, allocatable :: nhungexp(:),nhungexp2(:)
      integer, allocatable :: nhunglistsrc(:)

      real *8 boxsize(0:nlmax)
      real *8 src(2,ns),radsrc(ns)
      real *8 trg(2,nt)
      real *8, allocatable :: centers(:,:),centers2(:,:)
      real *8 expc(2,nexpc)
      real *8 radexp(nexpc)
c
      real *8 xmin,xmax,ymin,ymax,sizey
      integer ictr,ih,irefine,is,ie
      integer nss,nee,ntt
      integer ibox,ifirstbox,ilastbox,nnew,nbtot
      integer mnlist1, mnlist2, mnlist3, mnlist4


c
c     INPUT:
c
c     src           source locations        
c     ns            number of sources 
c     rads          source radii (determines deepest level that
c                   the source can reach) 
c     trg           target locations        
c     nt            number of targets
c
c     expc          expansion center locations
c     nexpc         number of expansion centers
c
c     idivflag      0 => divide on sources
c                   1 => divide on targets
c                   2 => divide on sources+targets
c                   3 => divide on sources+targets+expansion centers
c 
c     ndiv          refinement criterion - extend tree until all
c                   nodes at finest level have less than ndiv 
c                   source/targets/sources+targets depending on
c                   idivflag
c
c     nlmax         max number of levels expected in the tree
c     nbmax         max number of boxes expected in the tree
c                   (optional) The code by default
c                   assumes nbmax = 8 N log_{4} (N) boxes
c                   The main subdivision is called with
c                   nboxes = max(nbmax,8N log_{4} N)
c
c     OUTPUT:
c     ier           error code
c                   ier = 0, if memory determination successful
c                   ier = 1, if nbmax < nboxes
c
c     nlevels       number of levels
c     nboxes        number of boxes
c
c     mhung         max number of hung chunks
c     ltree         length of tree = (89 + mhung)*nboxes + 
c                   2*(nlevels+1) + ns + nt + nexpc
c
c
c     Other notes:
c
c     refinement criterion w.r.t rads is:
c
c     hang if (rads .geq. boxsize)
c
c

      ntot = ns + nt + nexpc
      nbtmp = 8*ntot*dlog(ntot+0.0d0)/dlog(4+0.0d0)
cc      nbmax = max(nbmax,nbtmp)
      nbmax = 10000
      nboxes = nbmax
      lcenters = nboxes
      nlevels = nlmax

      ier =0

      allocate(centers(2,nboxes))
      allocate(laddr(2,0:nlevels))
      allocate(ilevel(nboxes))
      allocate(iparenttemp(nboxes))
      allocate(nchild(nboxes))
      allocate(ichildtemp(4,nboxes))
      allocate(isourcetemp(ns))
      allocate(itargettemp(nt))
      allocate(iexpctemp(nexpc))
      allocate(ihsfirsttemp(nboxes))
      allocate(ihslasttemp(nboxes))
      allocate(isfirsttemp(nboxes))
      allocate(islasttemp(nboxes))
      allocate(itfirsttemp(nboxes))
      allocate(itlasttemp(nboxes))
      allocate(ihefirsttemp(nboxes))
      allocate(ihelasttemp(nboxes))
      allocate(iefirsttemp(nboxes))
      allocate(ielasttemp(nboxes))
      allocate(nhungsrc(nboxes))
      allocate(nhungexp(nboxes))

c     Step 1: Find enclosing box
c
      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
c
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
      enddo
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
      enddo

      do i=1,nexpc
        if(expc(1,i) .lt. xmin) xmin=expc(1,i)
        if(expc(1,i) .gt. xmax) xmax=expc(1,i)
        if(expc(2,i) .lt. ymin) ymin=expc(2,i)
        if(expc(2,i) .gt. ymax) ymax=expc(2,i)
      enddo
      boxsize(0)=xmax-xmin
      sizey=ymax-ymin
      if(sizey .gt. boxsize(0)) boxsize(0)=sizey
c
c     initialize arrays at level 0
c
      centers(1,1)=(xmin+xmax)/2
      centers(2,1)=(ymin+ymax)/2
      laddr(1,0) = 1
      laddr(2,0) = 1
      ilevel(1) = 0
      iparenttemp(1) = -1
      isfirsttemp(1) = 1
c
c     count number of hung sources
c     and hang up "big" sources
      nhungsrc(1) = 0
      nhungexp(1) = 0
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) nhungsrc(1)=nhungsrc(1)+1 
      enddo
      isfirsttemp(1) = nhungsrc(1)+1
      islasttemp(1) = ns
      if (nhungsrc(1) .gt. 0) then 
         ihsfirsttemp(1) = 1
         ihslasttemp(1) = nhungsrc(1)
      else
         ihsfirsttemp(1) = 0
         ihslasttemp(1) = -1
      endif
c     Count number of hung expansion centers      
c     and hang up "big" expansion centers
      do i=1,nexpc
         if (radexp(i).gt.boxsize(0)) nhungexp(1)=nhungexp(1)+1
      enddo
      iefirsttemp(1) = nhungexp(1)+1
      ielasttemp(1) = nexpc
      if (nhungexp(1).gt.0) then
          ihefirsttemp(1) = 1
          ihelasttemp(1) = nhungexp(1)
      else
         ihefirsttemp(1) = 0
         ihelasttemp(1) = -1
      endif

c
c     reorder isourcetemp to put hung sources in beginning
      ih = 0
      is = nhungsrc(1)
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) then
            ih = ih+1
            isourcetemp(ih) = i
         else
            is = is+1
            isourcetemp(is) = i
         endif
      enddo
c     reorder iexptemp to put hung expansion centers in beginning
      ih = 0
      ie = nhungexp(1)
      do i= 1,nexpc
         if(radexp(i).gt.boxsize(0)) then
            ih = ih+1
            iexpctemp(ih) = i
         else
            ie = ie+1
            iexpctemp(ie) = i
         endif
      enddo

c     initialize itargettemp 
      do i = 1,nt
         itargettemp(i) = i
      enddo
      itfirsttemp(1) = 1
      itlasttemp(1) = nt

      nchild(1) = 0
      ichildtemp(1,1) = -1
      ichildtemp(2,1) = -1
      ichildtemp(3,1) = -1
      ichildtemp(4,1) = -1

c
      irefine = 0
      nss = ns - nhungsrc(1)
      nee = nexpc - nhungexp(1)
      if ((idivflag .eq.0).and.(nss.gt.ndiv)) irefine=1
      if ((idivflag .eq.1).and.(nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.2).and.(nss+nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.3).and.(nss+nt+nee.gt.ndiv)) irefine=1

c     Reset nlevels, nboxes      
      nlevels = 0
      nboxes = 1
      do i = 1,nlmax
         if (irefine.eq.1) then

c
c
c      estimate number of new boxes to be created
c
            nnew = 0
            ifirstbox = laddr(1,nlevels)
            ilastbox = laddr(2,nlevels)

            do ibox = ifirstbox,ilastbox
               nss = islasttemp(ibox)-isfirsttemp(ibox)+1
               ntt = itlasttemp(ibox)-itfirsttemp(ibox)+1
               nee = ielasttemp(ibox)-iefirsttemp(ibox)+1
               if ((idivflag .eq.0).and.(nss.gt.ndiv)) nnew = nnew + 1
               if ((idivflag .eq.1).and.(nt.gt.ndiv)) nnew = nnew + 1
               if ((idivflag .eq.2).and.(nss+nt.gt.ndiv))  
     1            nnew = nnew + 1
               if ((idivflag .eq.3).and.(nss+nt+nee.gt.ndiv))
     1              nnew = nnew + 1          
            enddo

            nbtot = nboxes + 4*nnew

c
c            if current memory is not sufficient, 
c            delete previous arrays and allocate more memory
c            allocate more memory
c
            if(nbtot.gt.nbmax) then
              nbmax = nbtot
              lcenters = nbtot
              allocate(centers2(2,nboxes))
              allocate(ilevel2(nboxes),iparenttemp2(nboxes))
              allocate(nchild2(nboxes),ichildtemp2(4,nboxes))
              allocate(ihsfirsttemp2(nboxes),ihslasttemp2(nboxes))
              allocate(isfirsttemp2(nboxes),islasttemp2(nboxes))
              allocate(itfirsttemp2(nboxes),itlasttemp2(nboxes))
              allocate(ihefirsttemp2(nboxes),ihelasttemp2(nboxes))
              allocate(iefirsttemp2(nboxes),ielasttemp2(nboxes))
              allocate(nhungsrc2(nboxes),nhungexp2(nboxes))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)
              do ibox = 1,nboxes
                ilevel2(ibox) = ilevel(ibox)
                iparenttemp2(ibox) = iparenttemp(ibox)
                nchild2(ibox) = nchild(ibox)
                ihsfirsttemp2(ibox) = ihsfirsttemp(ibox)
                ihslasttemp2(ibox) = ihslasttemp(ibox)
                isfirsttemp2(ibox) = isfirsttemp(ibox)
                islasttemp2(ibox) = islasttemp(ibox)
                itfirsttemp2(ibox) = itfirsttemp(ibox)
                itlasttemp2(ibox) = itlasttemp(ibox)
                ihefirsttemp2(ibox) = ihefirsttemp(ibox)
                ihelasttemp2(ibox) = ihelasttemp(ibox)
                iefirsttemp2(ibox) = iefirsttemp(ibox)
                ielasttemp2(ibox) = ielasttemp(ibox)
                nhungsrc2(ibox) = nhungsrc(ibox)
                nhungexp2(ibox) = nhungexp(ibox)
                do j=1,2
                  centers2(j,ibox) = centers(j,ibox)
                enddo

                do j=1,4
                  ichildtemp2(j,ibox) = ichildtemp(j,ibox)
                enddo

              enddo
C$OMP END PARALLEL DO              
              deallocate(centers,ilevel,iparenttemp,nchild,ichildtemp)
              deallocate(ihsfirsttemp,ihslasttemp)
              deallocate(isfirsttemp,islasttemp)
              deallocate(itfirsttemp,itlasttemp)
              deallocate(ihefirsttemp,ihelasttemp)
              deallocate(iefirsttemp,ielasttemp)
              deallocate(nhungsrc,nhungexp)


              allocate(centers(2,nbmax))
              allocate(ilevel(nbmax),iparenttemp(nbmax))
              allocate(nchild(nbmax),ichildtemp(4,nbmax))
              allocate(ihsfirsttemp(nbmax),ihslasttemp(nbmax))
              allocate(isfirsttemp(nbmax),islasttemp(nbmax))
              allocate(itfirsttemp(nbmax),itlasttemp(nbmax))
              allocate(ihefirsttemp(nbmax),ihelasttemp(nbmax))
              allocate(iefirsttemp(nbmax),ielasttemp(nbmax))
              allocate(nhungsrc(nbmax),nhungexp(nbmax))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)
              do ibox = 1,nboxes
                ilevel(ibox) = ilevel2(ibox)
                iparenttemp(ibox) = iparenttemp2(ibox)
                nchild(ibox) = nchild2(ibox)
                ihsfirsttemp(ibox) = ihsfirsttemp2(ibox)
                ihslasttemp(ibox) = ihslasttemp2(ibox)
                isfirsttemp(ibox) = isfirsttemp2(ibox)
                islasttemp(ibox) = islasttemp2(ibox)
                itfirsttemp(ibox) = itfirsttemp2(ibox)
                itlasttemp(ibox) = itlasttemp2(ibox)
                ihefirsttemp(ibox) = ihefirsttemp2(ibox)
                ihelasttemp(ibox) = ihelasttemp2(ibox)
                iefirsttemp(ibox) = iefirsttemp2(ibox)
                ielasttemp(ibox) = ielasttemp2(ibox)
                nhungsrc(ibox) = nhungsrc2(ibox)
                nhungexp(ibox) = nhungexp2(ibox)
                do j=1,2
                  centers(j,ibox) = centers2(j,ibox)
                enddo
                do j=1,4
                  ichildtemp(j,ibox) = ichildtemp2(j,ibox)
                enddo
              enddo
C$OMP END PARALLEL DO              



              deallocate(centers2,ilevel2,iparenttemp2,nchild2)
              deallocate(ichildtemp2)
              deallocate(ihsfirsttemp2,ihslasttemp2)
              deallocate(isfirsttemp2,islasttemp2)
              deallocate(itfirsttemp2,itlasttemp2)
              deallocate(ihefirsttemp2,ihelasttemp2)
              deallocate(iefirsttemp2,ielasttemp2)
              deallocate(nhungsrc2,nhungexp2)
            endif


            call subdivide_adap(src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   nlevels,nboxes,
     $                   centers,lcenters,boxsize,nbmax,nlmax,
     $                   laddr,ilevel,iparenttemp,nchild,ichildtemp,
     $                   isourcetemp,itargettemp,iexpctemp,
     $                   ihsfirsttemp,ihslasttemp,
     $                   isfirsttemp,islasttemp,
     $                   itfirsttemp,itlasttemp,
     $                   ihefirsttemp,ihelasttemp,
     $                   iefirsttemp,ielasttemp,nhungsrc,
     $                   nhungexp,irefine)

         else
            exit
         endif
      enddo


c
c
c        check with Leslie/Dhairya/Alex/Zydrunas if 8 is enough
c

      nbtot = 8*nboxes
c
c            if current memory is not sufficient, 
c            delete previous arrays and allocate more memory
c            allocate more memory
c
      if(nbtot.gt.nbmax) then
        nbmax = nbtot
        lcenters = nbtot

        allocate(centers2(2,nboxes))
        allocate(ilevel2(nboxes),iparenttemp2(nboxes))
        allocate(nchild2(nboxes),ichildtemp2(4,nboxes))
        allocate(ihsfirsttemp2(nboxes),ihslasttemp2(nboxes))
        allocate(isfirsttemp2(nboxes),islasttemp2(nboxes))
        allocate(itfirsttemp2(nboxes),itlasttemp2(nboxes))
        allocate(ihefirsttemp2(nboxes),ihelasttemp2(nboxes))
        allocate(iefirsttemp2(nboxes),ielasttemp2(nboxes))
        allocate(nhungsrc2(nboxes),nhungexp2(nboxes))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)
        do ibox = 1,nboxes
          ilevel2(ibox) = ilevel(ibox)
          iparenttemp2(ibox) = iparenttemp(ibox)
          nchild2(ibox) = nchild(ibox)
          ihsfirsttemp2(ibox) = ihsfirsttemp(ibox)
          ihslasttemp2(ibox) = ihslasttemp(ibox)
          isfirsttemp2(ibox) = isfirsttemp(ibox)
          islasttemp2(ibox) = islasttemp(ibox)
          itfirsttemp2(ibox) = itfirsttemp(ibox)
          itlasttemp2(ibox) = itlasttemp(ibox)
          ihefirsttemp2(ibox) = ihefirsttemp(ibox)
          ihelasttemp2(ibox) = ihelasttemp(ibox)
          iefirsttemp2(ibox) = iefirsttemp(ibox)
          ielasttemp2(ibox) = ielasttemp(ibox)
          nhungsrc2(ibox) = nhungsrc(ibox)
          nhungexp2(ibox) = nhungexp(ibox)
          do j=1,2
            centers2(j,ibox) = centers(j,ibox)
          enddo
           do j=1,4
            ichildtemp2(j,ibox) = ichildtemp(j,ibox)
          enddo
         enddo
C$OMP END PARALLEL DO              
        deallocate(centers,ilevel,iparenttemp,nchild,ichildtemp)
        deallocate(ihsfirsttemp,ihslasttemp)
        deallocate(isfirsttemp,islasttemp)
        deallocate(itfirsttemp,itlasttemp)
        deallocate(ihefirsttemp,ihelasttemp)
        deallocate(iefirsttemp,ielasttemp)
        deallocate(nhungsrc,nhungexp)

        allocate(centers(2,nbmax))
        allocate(ilevel(nbmax),iparenttemp(nbmax))
        allocate(nchild(nbmax),ichildtemp(4,nbmax))
        allocate(ihsfirsttemp(nbmax),ihslasttemp(nbmax))
        allocate(isfirsttemp(nbmax),islasttemp(nbmax))
        allocate(itfirsttemp(nbmax),itlasttemp(nbmax))
        allocate(ihefirsttemp(nbmax),ihelasttemp(nbmax))
        allocate(iefirsttemp(nbmax),ielasttemp(nbmax))
        allocate(nhungsrc(nbmax),nhungexp(nbmax))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)
        do ibox = 1,nboxes
          ilevel(ibox) = ilevel2(ibox)
          iparenttemp(ibox) = iparenttemp2(ibox)
          nchild(ibox) = nchild2(ibox)
          ihsfirsttemp(ibox) = ihsfirsttemp2(ibox)
          ihslasttemp(ibox) = ihslasttemp2(ibox)
          isfirsttemp(ibox) = isfirsttemp2(ibox)
          islasttemp(ibox) = islasttemp2(ibox)
          itfirsttemp(ibox) = itfirsttemp2(ibox)
          itlasttemp(ibox) = itlasttemp2(ibox)
          ihefirsttemp(ibox) = ihefirsttemp2(ibox)
          ihelasttemp(ibox) = ihelasttemp2(ibox)
          iefirsttemp(ibox) = iefirsttemp2(ibox)
          ielasttemp(ibox) = ielasttemp2(ibox)
          nhungsrc(ibox) = nhungsrc2(ibox)
          nhungexp(ibox) = nhungexp2(ibox)
          do j=1,2
            centers(j,ibox) = centers2(j,ibox)
          enddo
          do j=1,4
            ichildtemp(j,ibox) = ichildtemp2(j,ibox)
          enddo
        enddo
C$OMP END PARALLEL DO              

        deallocate(centers2,ilevel2,iparenttemp2,nchild2)
        deallocate(ichildtemp2)
        deallocate(ihsfirsttemp2,ihslasttemp2)
        deallocate(isfirsttemp2,islasttemp2)
        deallocate(itfirsttemp2,itlasttemp2)
        deallocate(ihefirsttemp2,ihelasttemp2)
        deallocate(iefirsttemp2,ielasttemp2)
        deallocate(nhungsrc2,nhungexp2)
      endif


c     Set up computation of list1 and list2      
      allocate(nnbors(nbmax))
      allocate(nbors(9,nbmax))

      do i=1,nboxes
         nnbors(i) = 0
         do j=1,9
            nbors(j,i) = -1
         enddo
      enddo

c     compute colleagues
      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparenttemp,nchild,
     2                   ichildtemp,nnbors,nbors)

       if(nlevels.ge.2) 
     $  call d2hpfixtree(src,ns,radsrc,trg,nt,expc,nexpc,
     $              radexp,nlevels,nboxes,
     $              centers,boxsize,nbmax,
     $              laddr,ilevel,iparenttemp,nchild,ichildtemp,
     $              nnbors,nbors,  
     $              isourcetemp,itargettemp,iexpctemp,
     $              ihsfirsttemp,ihslasttemp,
     $              isfirsttemp,islasttemp,
     $              itfirsttemp,itlasttemp,
     $              ihefirsttemp,ihelasttemp,
     $              iefirsttemp,ielasttemp,nhungsrc,
     $              nhungexp)

       if(nboxes.gt.nbmax) then
           ier = 1
           call prinf('Not enough max boxes, set nbmax
     1                 to higher value. Default value
     2                 used is 8*n*log(n), where n=
     3                 ns+nt+nexpc',ier,1)
           return
       endif

       
      mnlist1 = 13
      mnlist2 = 27
      mnlist3 = 20
      mnlist4 = 5
      allocate(nlist1(nboxes))
      allocate(nlist2(nboxes))
      allocate(nlist3(nboxes))
      allocate(nlist4(nboxes))
      allocate(nhunglistsrc(nboxes))

      allocate(list1(mnlist1,nboxes))
      allocate(list2(mnlist2,nboxes))
      allocate(list3(mnlist3,nboxes))
      allocate(list4(mnlist4,nboxes))

      do i=1,nboxes
         nlist1(i) = 0
         nlist2(i) = 0
         nlist3(i) = 0
         nlist4(i) = 0
         nhunglistsrc(i) = 0
         do j=1,mnlist1
            list1(j,i) = -1
         enddo
         do j=1,mnlist2
            list2(j,i) = -1
         enddo
         do j=1,mnlist3
            list3(j,i) = -1
         enddo
         do j=1,mnlist4
            list4(j,i) = -1
         enddo
      enddo

      call computelists(nlevels,nboxes,laddr,boxsize,centers,
     1                  iparenttemp,nchild,ichildtemp,nnbors,
     2                  nbors,mnlist1,nlist1,list1,mnlist2,nlist2,
     3                  list2,mnlist3,nlist3,list3,mnlist4,
     4                  nlist4,list4)

       call computemhung(nlevels,nboxes,laddr,iparenttemp,nchild,
     1                   nnbors,nbors,mnlist1,nlist1,list1,
     2                   nhungsrc,nhunglistsrc,mhung)
      ltree = (88+mhung)*nboxes + 2*(nlevels+1) + ns+nexpc+nt

      return
      end
c-------------------------------------------------------------------
c
c
c
c
c
      subroutine print_tree(itree,ltree,nboxes,centers,boxsize,nlevels,
     1   iptr,fname)
c
c        this subroutine writes the tree info to a file
c
c        input arguments:
c          itree - integer (ltree)
c             packed array containing tree info
c          ltree - integer
c            length of itree
c          nboxes - integer
c             number of boxes
c          centers - real *8 (2,nboxes)
c             xy coordinates of box centers in tree hierarchy
c          boxsize - real *8 (0:nlevels)
c             size of box at various levels
c          nlevels - integer
c             number of levels
c          iptr - integer(30)
c            pointer to various arrays inside itree
c          fname - character *
c            file name to which tree info is to be written
c 
c          output
c            file with name fname, which contains the tree info
c            file can be plotted using the python script
c              tree_plot.py containted in src/common
c


      implicit real *8 (a-h,o-z)
      integer itree(ltree),ltree,nboxes,nlevels,iptr(30)
      integer, allocatable :: ilevel(:)
      real *8 centers(2,nboxes),boxsize(0:nlevels)
      character (len=*) fname

      allocate(ilevel(nboxes))
      
      do ilev=0,nlevels
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          ilevel(ibox) = ilev
        enddo
      enddo

      open(unit=33,file=trim(fname))
      

 1111 format(10(2x,e11.5))      

      do ibox=1,nboxes
         if(itree(iptr(3)+ibox-1).eq.0) then
           ilev = ilevel(ibox) 
           bs = boxsize(ilev)
           x1 = centers(1,ibox) - bs/2
           x2 = centers(1,ibox) + bs/2

           y1 = centers(2,ibox) - bs/2
           y2 = centers(2,ibox) + bs/2
           
           write(33,1111) x1,x2,x2,x1,x1,y1,y1,y2,y2,y1
         endif
      enddo

      close(33)

      return
      end
c
c
c
c
c
c
c
c
c
      subroutine print_tree_flag(itree,ltree,nboxes,centers,boxsize,
     1   nlevels,iptr,iflag,fname)
c
c        this subroutine writes the tree info to a file
c
c        input arguments:
c          itree - integer (ltree)
c             packed array containing tree info
c          ltree - integer
c            length of itree
c          nboxes - integer
c             number of boxes
c          centers - real *8 (2,nboxes)
c             xy coordinates of box centers in tree hierarchy
c          boxsize - real *8 (0:nlevels)
c             size of box at various levels
c          nlevels - integer
c             number of levels
c          iptr - integer(30)
c            pointer to various arrays inside itree
c          iflag - integer(nboxes)
c            print box if iflag(ibox)=1
c
c          fname - character *
c            file name to which tree info is to be written
c 
c          output
c            file with name fname, which contains the tree info
c            file can be plotted using the python script
c              tree_plot.py containted in src/common
c


      implicit real *8 (a-h,o-z)
      integer itree(ltree),ltree,nboxes,nlevels,iptr(30)
      integer iflag(nboxes)
      integer, allocatable :: ilevel(:)
      real *8 centers(2,nboxes),boxsize(0:nlevels)
      character (len=*) fname

      allocate(ilevel(nboxes))
      
      do ilev=0,nlevels
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          ilevel(ibox) = ilev
        enddo
      enddo

      open(unit=33,file=trim(fname))
      

 1111 format(10(2x,e11.5))      

      do ibox=1,nboxes
         if(iflag(ibox).eq.1) then
           ilev = ilevel(ibox) 
           bs = boxsize(ilev)
           x1 = centers(1,ibox) - bs/2
           x2 = centers(1,ibox) + bs/2

           y1 = centers(2,ibox) - bs/2
           y2 = centers(2,ibox) + bs/2
           
           write(33,1111) x1,x2,x2,x1,x1,y1,y1,y2,y2,y1
         endif
      enddo

      close(33)

      return
      end
c
c
c
c
