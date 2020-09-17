cc Copyright (C) 2010-2011: Leslie Greengard, Zydrunas Gimbustas 
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



c       
c   Helmholtz FMM in R^2: evaluate all pairwise particle
c   interactions (ignoring self-interaction) 
c   and interactions with targets.
c
c   We use log(r) for the Green's function.
c   Self-interactions are not included
c
c   h2d: charge and dipstr are complex valued, x in \R^2
c
c   \phi(x_i) = \sum_{j\ne i} charge_j log|x_i-x_j}
c   + dipstr_j/x_i - x_j
c
c

      subroutine lfmm2d_st_c_p(eps,ns,sources,
     1            charge,pot,nt,targ,pottarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   pottarg(nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns)

      complex *16 pot(ns)
      complex *16 pottarg(nt)

c
cc     temporary variables
c
      complex *16 dipstr
      complex *16 grad,gradtarg
      complex *16 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      nd = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------


      subroutine lfmm2d_st_c_g(eps,ns,sources,
     1            charge,pot,grad,nt,targ,pottarg,gradtarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(ns)      : gradients at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(nt)  : gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns)

      complex *16 pot(ns),grad(ns)
      complex *16 pottarg(nt),gradtarg(nt)

c
cc     temporary variables
c
      complex *16 dipstr
      complex *16 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      nd = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine lfmm2d_st_c_h(eps,ns,sources,
     1            charge,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(ns)      : gradients at the source locations
c   hess(ns)      : hessian at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(nt)  : gradient at the target locations
c   hesstarg(nt)  : hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns)

      complex *16 pot(ns),grad(ns),hess(ns)
      complex *16 pottarg(nt),gradtarg(nt),hesstarg(nt)

c
cc     temporary variables
c
      complex *16 dipstr
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 3

      nd = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end


c-------------------------------      

      subroutine lfmm2d_st_d_p(eps,ns,sources,
     1            dipstr,pot,nt,targ,pottarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   pottarg(nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(ns)

      complex *16 pot(ns)
      complex *16 pottarg(nt)

c
cc     temporary variables
c
      complex *16 charge
      complex *16 grad,gradtarg
      complex *16 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      nd = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------


      subroutine lfmm2d_st_d_g(eps,ns,sources,
     1            dipstr,pot,grad,nt,targ,pottarg,gradtarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(ns)      : gradients at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(nt)  : gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(ns)

      complex *16 pot(ns),grad(ns)
      complex *16 pottarg(nt),gradtarg(nt)

c
cc     temporary variables
c
      complex *16 charge
      complex *16 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      nd = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine lfmm2d_st_d_h(eps,ns,sources,
     1            dipstr,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(ns)      : gradients at the source locations
c   hess(ns)      : hessian at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(nt)  : gradient at the target locations
c   hesstarg(nt)  : hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(ns)

      complex *16 pot(ns),grad(ns),hess(ns)
      complex *16 pottarg(nt),gradtarg(nt),hesstarg(nt)

c
cc     temporary variables
c
      complex *16 charge
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      nd = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end


c-------------------------------      

      subroutine lfmm2d_st_cd_p(eps,ns,sources,charge,
     1            dipstr,pot,nt,targ,pottarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   pottarg(nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pot(ns)
      complex *16 pottarg(nt)

c
cc     temporary variables
c
      complex *16 grad,gradtarg
      complex *16 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      nd = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------


      subroutine lfmm2d_st_cd_g(eps,ns,sources,charge,
     1            dipstr,pot,grad,nt,targ,pottarg,gradtarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(ns)      : gradients at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(nt)  : gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pot(ns),grad(ns)
      complex *16 pottarg(nt),gradtarg(nt)

c
cc     temporary variables
c
      complex *16 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      nd = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine lfmm2d_st_cd_h(eps,ns,sources,charge,
     1            dipstr,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(ns)      : gradients at the source locations
c   hess(ns)      : hessian at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(nt)  : gradient at the target locations
c   hesstarg(nt)  : hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pot(ns),grad(ns),hess(ns)
      complex *16 pottarg(nt),gradtarg(nt),hesstarg(nt)

c
cc     temporary variables
c
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      nd = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end


      subroutine lfmm2d_st_c_p_vec(nd,eps,ns,sources,
     1            charge,pot,nt,targ,pottarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------


      subroutine lfmm2d_st_c_g_vec(nd,eps,ns,sources,
     1            charge,pot,grad,nt,targ,pottarg,gradtarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine lfmm2d_st_c_h_vec(nd,eps,ns,sources,
     1            charge,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   hess(nd,3,ns)    : hessian at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c   hesstarg(nd,3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 3

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end


c-------------------------------      

      subroutine lfmm2d_st_d_p_vec(nd,eps,ns,sources,
     1            dipstr,pot,nt,targ,pottarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------


      subroutine lfmm2d_st_d_g_vec(nd,eps,ns,sources,
     1            dipstr,pot,grad,nt,targ,pottarg,gradtarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine lfmm2d_st_d_h_vec(nd,eps,ns,sources,
     1            dipstr,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   hess(nd,3,ns)    : hessian at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c   hesstarg(nd,3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end


c-------------------------------      

      subroutine lfmm2d_st_cd_p_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,nt,targ,pottarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------


      subroutine lfmm2d_st_cd_g_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,grad,nt,targ,pottarg,gradtarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine lfmm2d_st_cd_h_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   hess(nd,3,ns)    : hessian at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c   hesstarg(nd,3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

c
cc     temporary variables
c
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      call lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end


c-------------------------------      
      subroutine lfmm2dpart(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns)    : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   dipstr(nd,ns)    : dipole strengths
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(2,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,2,*)    : gradients at the source locations
c   hess(nd,3,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,2,*): gradient at the target locations
c   hesstarg(nd,3,*): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 omp_get_wtime
      real *8 eps
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,*),dipstr(nd,*)

      complex *16 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
      complex *16 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)

c
cc      Tree variables
c
      integer, allocatable :: itree(:)
      integer ipointer(30)
      real *8, allocatable :: tcenters(:,:),boxsize(:)
      integer nexpc,ntj
      real *8 expc(2)
      real *8 scj
      complex *16 jexps(100)
      integer ier,idivflag,nlevels,nboxes,ndiv,nlmax,nbmax
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer mhung,ltree

      real *8, allocatable :: radsrc(:)
      real *8 radexp

c
cc     sorted arrays
c
      real *8, allocatable :: sourcesort(:,:)
      real *8, allocatable :: targsort(:,:)
      complex *16, allocatable :: chargesort(:,:),dipstrsort(:,:)
      complex *16, allocatable :: potsort(:,:),gradsort(:,:,:),
     1                             hesssort(:,:,:)
      complex *16, allocatable :: pottargsort(:,:),gradtargsort(:,:,:),
     1                              hesstargsort(:,:,:)

c
cc     additional fmm variables

      integer lmptot
      real *8, allocatable :: rscales(:)
      integer, allocatable :: nterms(:),iaddr(:,:)
      real *8, allocatable :: rmlexp(:)
      complex *16, allocatable :: mptemp(:)

c
cc      temporary variables
c
      integer i,ilev,lmptmp,nmax,idim
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg,ifprint
      real *8 time1,time2,pi,done

      done = 1
      pi = atan(done)*4.0d0


      nexpc = 0
      ier = 0

c
c    Need to fix ndiv in Helmmholtz FMM
c   
c

      nlevels = 0
      nboxes = 0
      idivflag =0
      ndiv = 20
      nlmax = 200
      nbmax = 0
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0

      mhung = 0
      ltree = 0

      ifprint = 0

      allocate(radsrc(ns))

      do i=1,ns
        radsrc(i) = 0
      enddo
      radexp = 0

c
cc      call the tree memory management
c       code to determine number of boxes,
c       number of levels and length of tree
c

      call maketree2dmem(ier,sources,ns,radsrc,targ,nt,expc,nexpc,
     1     radexp,idivflag,ndiv,nlmax,nbmax,nlevels,nboxes,mhung,
     2     ltree)

      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(tcenters(2,nboxes))

c
cc      call the tree code
c
      call maketree2d(sources,ns,radsrc,targ,nt,expc,nexpc,radexp,
     1      idivflag,ndiv,mhung,nlevels,nboxes,tcenters,boxsize,itree,
     2      ltree,ipointer,mnlist1,mnlist2,mnlist3,mnlist4)

      allocate(sourcesort(2,ns))
      allocate(targsort(2,nt))


      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        allocate(chargesort(nd,ns),dipstrsort(nd,1))
      endif
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        allocate(chargesort(nd,1),dipstrsort(nd,ns))
      endif
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        allocate(chargesort(nd,ns),dipstrsort(nd,ns))
      endif

      if(ifpgh.eq.1) then
        allocate(potsort(nd,ns),gradsort(nd,2,1),hesssort(nd,3,1))
      else if(ifpgh.eq.2) then
        allocate(potsort(nd,ns),gradsort(nd,2,ns),hesssort(nd,3,1))
      else if(ifpgh.eq.3) then
        allocate(potsort(nd,ns),gradsort(nd,2,ns),hesssort(nd,3,ns))
      else
        allocate(potsort(nd,1),gradsort(nd,2,1),hesssort(nd,3,1))
      endif

      
      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,2,1),
     1     hesstargsort(nd,3,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,2,nt),
     1      hesstargsort(nd,3,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,2,nt),
     1     hesstargsort(nd,3,nt))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,2,1),
     1     hesstargsort(nd,3,1))
      endif
      
c
cc      initialize potentials,hessians,gradients
c


      if(ifpgh.eq.1) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
          enddo
        enddo
      endif

      if(ifpgh.eq.2) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
          enddo
        enddo
      endif

      if(ifpgh.eq.3) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            hesssort(idim,1,i) = 0
            hesssort(idim,2,i) = 0
            hesssort(idim,3,i) = 0
          enddo
        enddo
      endif


      if(ifpghtarg.eq.1) then
        do i=1,nt
          do idim=1,nd
            pottarg(idim,i) = 0
            pottargsort(idim,i) = 0
          enddo
        enddo
      endif

      if(ifpghtarg.eq.2) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
          enddo
        enddo
      endif

      if(ifpghtarg.eq.3) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            hesstargsort(idim,1,i) = 0
            hesstargsort(idim,2,i) = 0
            hesstargsort(idim,3,i) = 0
          enddo
        enddo
      endif


c
cc      compute scaling factor for multipole/local expansions
c       and lengths of multipole and local expansions
c
      allocate(rscales(0:nlevels),nterms(0:nlevels))

      nmax = 0
      do i=0,nlevels
        rscales(i) = boxsize(i)
        call l2dterms(eps,nterms(i),ier)
        nterms(i) = nterms(i) 
        if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo

      if(ifprint.eq.1) call prinf('nmax=*',nmax,1)
      if(ifprint.eq.1) call prinf('nterms=*',nterms,nlevels+1)

c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(3,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(2,nboxes))

      lmptmp = (nmax+1)*nd
      allocate(mptemp(lmptmp))

c     reorder sources
c
      call dreorderf(2,ns,sources,sourcesort,itree(ipointer(5)))
      if(ifcharge.eq.1) 
     1    call dreorderf(2*nd,ns,charge,chargesort,itree(ipointer(5)))
      if(ifdipole.eq.1) then
         call dreorderf(2*nd,ns,dipstr,dipstrsort,itree(ipointer(5)))
      endif

c
cc     reorder targets
c
      call dreorderf(2,nt,targ,targsort,itree(ipointer(6)))



c
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
      call l2dmpalloc(nd,itree(ipointer(1)),iaddr,nlevels,lmptot,
     1    nterms)
      if(ifprint .eq. 1) call prinf(' lmptot is *',lmptot,1)

      allocate(rmlexp(lmptot),stat=ier)


c
cc     call the main fmm routine
c

c     Memory allocation is complete. 
c     Call main fmm routine
c
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call lfmm2dmain(nd,eps,
     $   ns,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipstrsort,
     $   nt,targsort,nexpc,expc,
     $   iaddr,rmlexp,mptemp,lmptmp,
     $   itree,ipointer,ndiv,nlevels,
     $   nboxes,boxsize,rscales,tcenters,itree(ipointer(1)),
     $   mnlist1,mnlist2,
     $   mnlist3,mnlist4,nterms,ntj,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,
     $   hesstargsort,jexps,scj)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)


c
cc      resort the output arrays in input order
c

      if(ifpgh.eq.1) then
        call dreorderi(2*nd,ns,potsort,pot,itree(ipointer(5)))
      endif

      if(ifpgh.eq.2) then
        call dreorderi(2*nd,ns,potsort,pot,itree(ipointer(5)))
        call dreorderi(2*nd,ns,gradsort,grad,itree(ipointer(5)))
      endif

      if(ifpgh.eq.3) then
        call dreorderi(2*nd,ns,potsort,pot,itree(ipointer(5)))
        call dreorderi(2*nd,ns,gradsort,grad,itree(ipointer(5)))
        call dreorderi(2*nd,ns,hesssort,hess,itree(ipointer(5)))
      endif

cc      call prini(6,13)
cc      call prin2('eps = *', eps, 1)
cc      call prin2('after lfmm2dmain, pottargsort = *', pottargsort, 30)
cc      stop
      
      if(ifpghtarg.eq.1) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itree(ipointer(6)))
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itree(ipointer(6)))
        call dreorderi(2*nd,nt,gradtargsort,gradtarg,itree(ipointer(6)))
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itree(ipointer(6)))
        call dreorderi(2*nd,nt,gradtargsort,gradtarg,itree(ipointer(6)))
        call dreorderi(2*nd,nt,hesstargsort,hesstarg,itree(ipointer(6)))
      endif


      return
      end
