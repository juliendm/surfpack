c _______________________________________________________________________
c
c Surfpack: A Software Library of Multidimensional Surface Fitting Methods
c Copyright (c) 2006, Sandia National Laboratories.
c This software is distributed under the GNU Lesser General Public License.
c For more information, see the README file in the top Surfpack directory.
c _______________________________________________________________________
c
c Sample user coded driver program for running MARS 3.6
c
c This is a simple program to run MARS on a small data set for testing
c purposes. The input data (provided below) is assumed to be on a file
c named 'mars.data' in the local directory. Printed output (also provided
c below) is sent to the standard output (Fortran file unit 6). Note that
c this output might be slightly different than that obtained by running
c this program on a computer with different floating point arithmetic.
c The resulting model ought to be (nearly) equivalent.
c
c This program provides only the minimal input necessary to run MARS. It
c thereby makes maximal use of internally supplied defaults, and therefore
c does not take advantage of the many user options available for guiding
c the analysis (see MARS documentation).
c
c set up problem dependent parameters for this (small) example:
c  50 observations (cases).
c   5 predictor variables (inputs).
c   2 variable interactions maximum.
c  15 basis functions maximum.
c
      parameter (n=50, np=5, mi=2, nk=15)
c
c set up (generic) graphics parameters:
c 100 raster points for curves.
c  40 raster points (on each axis) for surfaces.
c  plot piecewise-cubic (continuous derivative) model.
c  show surfaces inside respective bivariate convex hulls.
c 
      parameter (ngc=100, ngs=40, m=2, icx=1)
c
c  set up working storage:
c  constant dimensioned arrays may need to be set to larger values
c  for much larger problems (see MARS documentation).
c

      real x(n,np),y(n),w(n),fm(5000),sp(20000)
      real crv(ngc,2,nk),srf(ngs,ngs,nk)
      integer lx(np),im(5000),mm(5000)
      double precision dp(20000)

ccc 
ccc added by AAG:
ccc
      parameter (neval=4)
      real f(4), xeval(neval,np)
     
c
c set predictor variable flags to indicate all are unrestricted
c ordinal variables, and set observation weights:
c
      data lx,w /np*1,n*1.0/
c
c write output header:
c
      write(6,'(/,''  simple driver to test MARS 3.6. '')')
c
c read in data:
c
      open(10,file='mars.data',status='unknown')
      do 1 i=1,n
      read(10,*) (x(i,j),j=1,np),y(i)
    1 continue
c
c invoke MARS:
c
      call mars (n,np,x,y,w,nk,mi,lx,fm,im,sp,dp,mm)
c
c construct plots for interpreting resulting model:
c

      call plot (m,x,fm,im,ngc,ngs,icx,nc,crv,ns,srf,sp,mm)

ccc
ccc added by AAG: "fmod" subroutine call to get a 2-d slice of the 5-d space
ccc
c
c call fmod (m,n,x,fm,im,f,sp)
c
c calculates mars model response estimates for sets of covariate vectors.
c
c input:
c m = model flag:
c   = 1 => piecewise-linear mars model.
c   = 2 => piecewise-cubic mars model. (Ref[2] Sec. 3.7)
c n = number of covariate vectors.
c x(n,p) = covariate vectors.
c fm,im = same as in mars (see above).
c
c output:
c f(n) = value of the mars model estimate for each covariate vector.
c
c workspace:
c sp(n,2) : real.
c

      xeval(1,1) = 0.500
      xeval(1,2) = 0.333
      xeval(1,3) = 0.20
      xeval(1,4) = 0.143
      xeval(1,5) = 0.091

      xeval(2,1) = 0.500
      xeval(2,2) = 0.667
      xeval(2,3) = 0.20
      xeval(2,4) = 0.143
      xeval(2,5) = 0.091

      xeval(3,1) = 0.250
      xeval(3,2) = 0.333
      xeval(3,3) = 0.20
      xeval(3,4) = 0.143
      xeval(3,5) = 0.091

      xeval(4,1) = 0.250
      xeval(4,2) = 0.667
      xeval(4,3) = 0.20
      xeval(4,4) = 0.143
      xeval(4,5) = 0.091
      
c      call fmod (m,neval,xeval,fm,im,f,sp)                              SLB05
      call fmodm(m,neval,xeval,fm,im,f,sp)
      open(11,file='mars_2d.surf',status='unknown')
      write(11,*) (f(i),i=1,neval)


c
c write plots to output files for plotting with local graphics package:
c
      if(nc.le.0) go to 2
      open(11,file='mars.curves',status='unknown')
      write(11,*) (((crv(i,j,k),i=1,ngc),j=1,2),k=1,nc)
    2 continue
      if(ns.le.0) go to 3
      open(12,file='mars.surfs',status='unknown')
      write(12,*) (((srf(i,j,k),i=1,ngs),j=1,ngs),k=1,ns)
    3 continue
      stop
      end
c

