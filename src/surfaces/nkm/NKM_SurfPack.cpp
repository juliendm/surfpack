#include "NKM_SurfPack.hpp"
#include <cmath>

//the purpose of this file is to contain generic functions usable by everything

//#define __SURFPACK_ERR_CHECK__
//#define __SURFPACK_UNIT_TEST__

namespace nkm {

//not sure which of these usings are necessary, achieve compilation by random guessing
using std::fabs;
using std::pow;
using std::cerr;
using std::endl;
using std::ifstream;
using std::istream;
using std::ios;
//using std::numeric_limits;
using std::ofstream;
using std::ostream;
using std::setw;
using std::string;
//using std::vector;

int if_close_enough(double a, double b)
{
  //  if(std::fabs(a-b)>1.0e-5){
  //    std::printf("a=%20.14f b=%20.14f\n",a,b);  std::fflush(stdout);
  //    assert((std::fabs(a-b)<=1.0e-5));
  //  }
  return (std::fabs(a-b)<=1.0e-5);
};

///this function should be moved to surfpack.cpp
int nchoosek(int n, int k){
  int nck=1;

  int nloop=(k<n-k)?k:n-k;
  if(nloop>0) {
    nck=n;
    for(int iloop=1; iloop<nloop; iloop++)
      nck=(nck*(n-iloop))/(iloop+1);
  }
    
  return nck;
}

///this function should be moved to surfpack.cpp
MtxDbl& gen_rot_mat(MtxDbl& Rot, const MtxDbl& EulAng, int nvarsr){
#ifdef __SURFPACK_ERR_CHECK__
  assert((EulAng.getNRows()==(nvarsr*(nvarsr-1)/2))&&(EulAng.getNCols()==1));
#endif
  MtxDbl I(nvarsr,nvarsr), R(nvarsr,nvarsr), RotTemp(nvarsr,nvarsr);
  I.zero();
  int ivarr;
  for(ivarr=0; ivarr<nvarsr; ivarr++) 
    I(ivarr,ivarr)=1.0;
  Rot=I;
  int nang=nvarsr, iang, Iang=0;
  double c, s;
  for(ivarr=0;ivarr<nvarsr-1;ivarr++) {
    nang--;
    for(iang=0; iang<nang; iang++) {
      c=std::cos(EulAng(Iang,0));
      s=std::sin(EulAng(Iang,0));
      R=I;
      R(iang  ,iang  )= c;
      R(iang  ,iang+1)=-s;
      R(iang+1,iang  )= s;
      R(iang+1,iang+1)= c;
      matrix_mult(RotTemp,Rot,R,0,1.0);
      Rot=RotTemp;
      Iang++;
    }
  }
  return Rot;
}

//all coordinates are between -1 and 1
MtxDbl& gen_rand_rot_mat(MtxDbl& rot,int nvarsr) 
{
  int n_eul_ang=nchoosek(nvarsr, 2);
  //printf("n_eul_ang=%d\n",n_eul_ang);
  MtxDbl eul_ang(n_eul_ang,1);
  double pi=2.0*std::acos(0.0);
  int mymod = 1048576; //2^20 instead of 10^6 to be kind to the computer
  for(int i=0; i<n_eul_ang; ++i)
    eul_ang(i,0)=(std::rand() % mymod)*pi/mymod;
  rot.newSize(nvarsr,nvarsr);
  gen_rot_mat(rot, eul_ang, nvarsr);
  return rot;
}


///generates 2*nvarsr random samples between 0 and 1, the sample design is stored in a matrix with nvarsr rows and 2*nvars columns (one point per column), the sample design is binning optimal with the bins being chosen as the end points of a randomly rotated set of axes (so the BINS, but not the points, are maximin spaced) the design is NOT a latin hypercube and it is not symmetric, opposite octants are stored in sequential columns of the xr matrix.
MtxDbl& gen_rand_axis_bin_opt_samples_0to1(MtxDbl& xr, int nvarsr) 
{
  gen_rand_rot_mat(xr,nvarsr);
  xr.resize(nvarsr,2*nvarsr);
  int mymod = 1048576; //2^20 instead of 10^6 to be kind to the computer
  for(int j=nvarsr-1; j>=0; --j) {
    //printf("surfpack.cpp: i=%d",i);
    for(int i=0; i<nvarsr; ++i) {
      //printf(" j=%d",j); fflush(stdout);
      xr(i,2*j  )=2.0*std::floor(1.0+xr(i,j))-1.0;
      xr(i,2*j+1)=0.5*((-xr(i,2*j)*(std::rand() % mymod))/mymod+1.0);
      xr(i,2*j  )=0.5*(( xr(i,2*j)*(std::rand() % mymod))/mymod+1.0);
    }
    //printf("\n");
  }
  return xr;
}

/***********************************************************************/
/// num_multi_dim_poly_coef(int Nvarsr, int Ndeg) is for use with multi_dim_poly_power(), says how many coefficients are needed for a Nvarsr-dimensional polynomial "of degree abs(Ndeg)" if Ndeg>0 the polynomial contains all multidimensional monomials UPTO (inclusive) degree Ndeg, if Ndeg<0 it contains only those multidimensional monomials OF EXACT DEGREE abs(Ndeg)
int num_multi_dim_poly_coef(int Nvarsr, int Ndeg){
  int Npoly;
  if(Ndeg<0)
    Npoly=nchoosek(-Ndeg-1+Nvarsr,-Ndeg);
  else
    Npoly=nchoosek(Ndeg+Nvarsr,Ndeg);
  return Npoly;
}

/***********************************************************************/
/**** multi_dim_poly_power() is a rather long function so I gave it ****/
/**** this special "boundary" comment so you can easily tell where  ****/
/**** it begins and ends                                            ****/
/***********************************************************************/
///poly=multi_dim_poly_power(poly,Nvarsr,Ndeg,istart,jstart) returns a matrix of size={Npoly,Nvarsr}; if Ndeg>=0 this function returns the mixed partial powers of all Nvarsr-dimensional polynomials of degree less than or equal to zero (there are Npoly=nchoosek(Ndeg+Nvarsr,Ndeg) of these), if Ndeg<=0 this function returns the mixed partial powers of all Nvarsr-dimensional polynomials of exactly degree abs(Ndeg) (there are Npoly=nchoosek(abs(Ndeg)-1+Nvarsr,abs(Ndeg)) of these); istart and jstart are offsets from the beginning of the matrix that say where to start the next "group" of powers (istart and jstart are there to avoid needing to RECURSIVELY allocate, fill, copy contents of this "poly" to the next larger "poly" and then deallocate this poly (it's a performance/speed thing), when any function other than multi_dim_poly_power() calls multi_dim_poly_power(), istart and jstart should be left unspecified, which makes them default to zero; if the "user" specifies non-zero istart and jstart then it "voids the warranty" on multi_dim_poly_power(), specifically you could end up going beyond the bounds of poly (a memory error) 
MtxInt& multi_dim_poly_power(MtxInt& poly, int Nvarsr, int Ndeg, int istart, int jstart, int iffirst){
  //printf("istart=%d jstart=%d iffirst=%d Nvarsr=%d Ndeg=%d poly.NRows()=%d\n",istart,jstart,iffirst,Nvarsr,Ndeg,poly.getNRows());
  int Npoly, npoly;
#ifdef __SURFPACK_ERR_CHECK__
  int jstartorig=jstart;
#endif
  //determine the total number of polynomials that this call of multi_dim_poly_power() is supposed to find mixed partial powers for
  if(Ndeg<0)
    Npoly=nchoosek(-Ndeg-1+Nvarsr,-Ndeg);
  else
    Npoly=nchoosek(Ndeg+Nvarsr,Ndeg);

  //istart=jstart=0 (the default when istart and jstart are not specified) is the flag for this function being called by the "user" (as opposed to a bigger multi_dim_poly_power()) and we don't want to require the user to know how big poly should be so we will right size it for him/her, if the user specified nonzero istart and/or jstart then it is his/her responsibility for making sure that poly is already big enough but since I'm a nice guy I put in an assert for when the user runs this in debug mode
  if((istart==0)&&(jstart==0)&&(iffirst==1)) {
    //printf("newSizing: Npoly=%d\n",Npoly);
    poly.newSize(Nvarsr,Npoly);
  }
  else{
#ifdef __SURFPACK_ERR_CHECK__
    if(!((jstart+Npoly<=poly.getNCols())&&(istart+Nvarsr<=poly.getNRows()))) {
      printf("Error in multi_dim_poly_power(): you asked me to fill in mixed partial polynomial powers beyond the bounds of the MtxInt& poly that you gave me.  If you want me to resize poly don't specify a nonzero istart or jstart\n, jstart=%d Npoly=%d poly.NCols=%d istart=%d Nvarsr=%d poly.NRows=%d\n",jstart,Npoly,poly.getNCols(),istart,Nvarsr,poly.getNRows());
      assert((jstart+Npoly<=poly.getNCols())&&(istart+Nvarsr<=poly.getNRows()));
    }
#endif
  }

  int ivar, j;
  int ndeg; //this is for recursive calls to smaller multi_dim_poly_power()'s that's why it starts with little "n" instead of big "N"
  if(Ndeg==0) {
    // a Nvarsr-dimensional polynomial of total degree 0 has all mixed partial powers equal to zero, istart and jstart should be 0 but the user could have done something strange
    for(ivar=0;ivar<Nvarsr; ivar++)
      poly(istart+ivar,jstart)=0;
  }
  else if(Nvarsr==1) {
    //this is a failsafe for direct user input, it isn't necessary for recursive calls
    if(Ndeg>0){
      //all powers less than or equal to Ndeg for 1 variable are
      // [   0; ...
      //     1; ...
      //     2; ...
      //     3; ...
      //     .
      //     : 
      //  Ndeg];    
      // make sure to offset by istart and jstart (it shouldn't be necessary but who knows the user might have actually specified istart and jstart)
      for(j=0;j<=Ndeg;j++)
	poly(istart,jstart+j)=j;      
    }
    else{
      //Ndeg<0 says you only want powers of exactly abs(Ndeg) for 1 variable this is [-Ndeg]; istart and jstart "should be" zero but who knows what the user gave us as input
      poly(istart,jstart)=-Ndeg;
    }
  }
  else if(Ndeg==-1) {
    //Ndeg = -1 says you want all polynomials of EXACTLY degree 1, for which the "mixed" partial powers is just an identity matrix, but make sure to offset it by istart and jstart  
    for(j=0; j<Nvarsr; j++) {
      for(ivar=0; ivar<Nvarsr; ivar++)
	poly(istart+ivar,jstart+j)=0;    
      poly(istart+j,jstart+j)=1;
    }
  }
  else if(Ndeg>0) {
    //this logic should only be executed when this multi_dim_poly_power() is called by the "user" (it should not be executed when this multi_dim_poly_power() is called by a bigger multi_dim_poly_power()); Ndeg > 0 says "I want all polynomials of degree less than or equal to Ndeg" so we will do this in batches of exactly degree 0, exactly degree 1, exactly degree 2, ..., exactly degree Ndeg
    for(ndeg=0; ndeg<=Ndeg; ndeg++) {
      npoly=nchoosek(ndeg-1+Nvarsr,ndeg);
      multi_dim_poly_power(poly, Nvarsr, -ndeg, istart, jstart, 0);
      jstart+=npoly;
    }
#ifdef __SURFPACK_ERR_CHECK__
    assert(jstart-jstartorig==Npoly); //jstartorig should be zero but who knows what the user gave us
#endif      
  }
  else if(Nvarsr==2) { 
    //we know that Ndeg < -1, so we want all 2-dimensional polynomials of exaclty degree abs(Ndeg) and the answer for 2 dimensions is easy
    // [ abs(Ndeg)        0      ;...
    //   abs(Ndeg)-1      1      ;...
    //   abs(Ndeg)-2      2      ;...
    //      .             .      ;...
    //      :             :      ;...
    //      2        abs(Ndeg)-2 ;...
    //      1        abs(Ndeg)-1 ;...
    //      0        abs(Ndeg)   ]^T
    for(j=0; j<=-Ndeg; j++) {
      poly(istart  ,jstart+j)=-Ndeg-j;
      poly(istart+1,jstart+j)=j;
    }
  }
  else{ //we know that Ndeg < -1 and Nvarsr > 2
    //we want all Nvarsr-dimensional polynomials of exaclty degree abs(Ndeg) so we are going to start with the first row being ideg=abs(Ndeg) decrementing that and pairing it with ALL (Nvarsr-1)-dimensional polynomials of exactly degree ndeg=abs(Ndeg)-ideg
    int nvarsr=Nvarsr-1; //this is for recursive calls to smaller multi_dim_poly_power()'s that's why it starts with little "n" instead of big "N"
    for(int ideg=-Ndeg; ideg>=0; ideg--){
      ndeg=-Ndeg-ideg;
      npoly=nchoosek(ndeg-1+nvarsr,ndeg);
      for(j=0; j<npoly; j++)
	poly(istart,jstart+j)=ideg;
      multi_dim_poly_power(poly,nvarsr,-ndeg,istart+1,jstart,0);
      jstart+=npoly;
    }
#ifdef __SURFPACK_ERR_CHECK__
    assert(jstart-jstartorig==Npoly); //don't expect jstartorig to be zero 
#endif
  }
  
  return poly;  
}

/***********************************************************************/
/**** multi_dim_poly_power() ends here                              ****/
/***********************************************************************/

/// if ndeg>=0 generates a matrix "poly" of nvarsr-dimensional polynomial powers that for all polynomials WITHOUT MIXED POWERS up to (and including) degree ndeg.  If ndeg<0 it generates a nvarsr by nvarsr diagonal matrix "poly" of non-mixed polynomials of exact degree abs(ndeg), this diagonal matrix has abs(ndeg) for all its diagonal elements.
MtxInt& main_effects_poly_power(MtxInt& poly, int nvarsr, int ndeg) {

#ifdef __SURFPACK_ERR_CHECK__
  assert(nvarsr>0);
#endif

  if(ndeg<0) {
    int abs_ndeg=std::abs(ndeg);
    poly.newSize(nvarsr,nvarsr);
    poly.zero();
    for(int ivarsr=0; ivarsr<nvarsr; ++ivarsr)
      poly(ivarsr,ivarsr)=abs_ndeg;
    return poly;
  }
  else if(ndeg==0) {
    poly.newSize(nvarsr,1);
    poly.zero();
    return poly;
  }
  
  poly.newSize(nvarsr,1+nvarsr*ndeg);
  poly.zero();
  int ipoly=0;
  for(int ideg=1; ideg<=ndeg; ++ideg)
    for(int ivarsr=0; ivarsr<nvarsr; ++ivarsr)
      poly(ivarsr,++ipoly)=ideg;

  return poly;
}


MtxInt& poly_to_flypoly(MtxInt& flypoly, const MtxInt& poly, int maxorder) {
  int npoly=poly.getNCols();
  int nvars=poly.getNRows();
  flypoly.newSize(maxorder+1,npoly);
  for(int ipoly=0; ipoly<npoly; ++ipoly) {
    int nmult=0;
    for(int ivar=0; ivar<nvars; ++ivar)
      for(int i=0; i<poly(ivar,ipoly); ++i) 
	flypoly(++nmult,ipoly)=ivar;
    flypoly(0,ipoly)=nmult;
#ifdef __SURFPACK_ERR_CHECK__
    assert(nmult<=maxorder);
#endif
  }
  return flypoly;
}

//this function modifies flycoef as needed so make sure you pass in a COPY of 
//the original coefficients instead of the original coefficients themselves
void poly_der_to_flypoly(MtxInt& flypoly, MtxDbl& flycoef, const MtxInt& poly, 
			 const MtxInt& der, int ider, int maxorder) {
  int npoly=poly.getNCols();
  int nvars=poly.getNRows();
#ifdef __SURFPACK_ERR_CHECK__
  assert((flycoef.getNRows()==npoly)&&(flycoef.getNCols()==1)&&
	 (der.getNRows()==nvars)&&(0<=ider)&&(ider<der.getNCols()));
#endif
  flypoly.newSize(maxorder+1,npoly); //you really should have sized flypoly
  //appropriately outside of this function for efficiency, this is just a 
  //fail safe, you shouldn't be using it

  for(int ipoly=0; ipoly<npoly; ++ipoly) {
    double tempcoef=flycoef(ipoly,0);
    if(tempcoef==0.0) 
      flypoly(0,ipoly)=0;
    else {
      int nmult=0;
      for(int ivar=0; ivar<nvars; ++ivar) {
	//determine the derivative of this multidimensional monomial and
	//store it as a flypoly representation
	int thisder=der(ivar,ider);
	int powafterder=poly(ivar,ipoly)-thisder;
	if(powafterder<0) {
	  tempcoef=0.0;
	  nmult=0;
	  break; //break out of ivar loop
	}
	//since the derivative polynomial's coefficient isn't zero (so far)
	//we need to know what it is
	for(int jder=0; jder<thisder; ++jder) 
	  tempcoef*=(poly(ivar,ipoly)-jder);	  

	for(int i=0; i<powafterder; ++i)
	  flypoly(++nmult,ipoly)=ivar;
      }
      flycoef(ipoly,0)=tempcoef;
      flypoly(0,ipoly)=nmult;
    }
  }
  return;
}


MtxDbl& evaluate_flypoly(MtxDbl& y, const MtxInt& flypoly, const MtxDbl& coef, const MtxDbl& xr) {
  int npts=xr.getNCols();
  int npoly=flypoly.getNCols();
#ifdef __SURFPACK_ERR_CHECK__
  assert((0<npts)&&(0<xr.getNRows())&&(0<npoly)&&(0<flypoly.getNRows())&&
	 (npoly==coef.getNRows())&&(1==coef.getNCols()));
  {
    int maxorder=flypoly(0,0);
    for(int i=1; i<npoly; ++i)
      maxorder=(flypoly(0,i)<=maxorder)?maxorder:flypoly(0,i);
    assert(maxorder<flypoly.getNRows());
  }
#endif
  y.newSize(1,npts);
  for(int ipt=0; ipt<npts; ++ipt) {
    double tempy=0.0;
    for(int ipoly=0; ipoly<npoly; ++ipoly) {
      int nmult=flypoly(0,ipoly);
      double term=coef(ipoly,0);
      for(int imult=1; imult<=nmult; ++imult)
	term*=xr(flypoly(imult,ipoly),ipt);
      tempy+=term;
    }
    y(0,ipt)=tempy;
  }
  return y;
}

MtxDbl& evaluate_poly(MtxDbl& y, MtxInt& flypoly, const MtxInt& poly, const MtxDbl& coef, const MtxDbl& xr) {
  int nvars=poly.getNRows();
  int npoly=poly.getNCols();
#ifdef __SURFPACK_ERR_CHECK__
  assert((0<npoly)&&(0<nvars)&&(nvars==xr.getNRows()));
#endif  
  int maxorder=0;
  for(int ipoly=0; ipoly<npoly; ++ipoly) {
    int thisorder=poly(0,ipoly);
    for(int ivar=1; ivar<nvars; ++ivar)
      thisorder+=poly(ivar,ipoly);
    maxorder=(thisorder<=maxorder)?maxorder:thisorder;
  }
  poly_to_flypoly(flypoly, poly, maxorder);
  return evaluate_flypoly(y, flypoly, coef, xr);
}


MtxDbl& evaluate_poly_der(MtxDbl& dy, MtxInt& flypoly, MtxDbl& flycoef, 
			  const MtxInt& poly, const MtxInt& der, 
			  const MtxDbl& coef, const MtxDbl& xr) {
  int nvars=poly.getNRows();
  int npoly=poly.getNCols();
  int npts =xr.getNCols();
  int nder =der.getNCols();
#ifdef __SURFPACK_ERR_CHECK__
  assert((0<nvars)&&(nvars==xr.getNRows())&&(nvars==der.getNRows())&&
	 (0<npoly)&&(npoly==coef.getNRows())&&(1==coef.getNCols())&&
	 (0<npts)&&(0<nder));
#endif

  //determine the maximum total order of any multidimensional monomial in poly
  //so that we know how many rows flypoly will need.
  int maxorder=0;
  for(int ipoly=0; ipoly<npoly; ++ipoly) {
    int thisorder=poly(0,ipoly);
    for(int ivar=1; ivar<nvars; ++ivar)
      thisorder+=poly(ivar,ipoly);
    maxorder=(thisorder<=maxorder)?maxorder:thisorder;
  }
  

  dy.newSize(nder,npts); //we want all the information (derivatives) 
  //associated with a point to be contiguous in memory (this means in the
  //same column of the dy); The nkm::SurfMat  templated matrix class 
  //uses Column Major Ordering (the same convention as MATLAB and FORTRAN) 
  //to facilitate the use of BLAS and LAPACK (which is written in FORTRAN)
  //C++ normally uses Row Major Ordering (each column is contiguous, but 
  //rows are not);

  for(int ider=0; ider<nder; ++ider) { 
    //we are deliberately not going to assign the result in contiguous 
    //order BECAUSE we only want to determine each derivative of the 
    //polynomial once, instead of once for each point, and we don't want 
    //to have to allocate enough memory to store all the derivatives of 
    //all the polynomials at the same time.

    //determine the flypoly representation of the polynomial that is the 
    //der(:,ider) mixed multidimensional derivative of the polynomial in poly
    flycoef.copy(coef); //need to copy coef so we don't modify the original
    poly_der_to_flypoly(flypoly, flycoef, poly, der, ider, maxorder);
    
   
    for(int ipt=0; ipt<npts; ++ipt) { //loop over points
      double sumofterms=0.0; //this is workspace to compute/hold the 
      //current derivative at the current point so we don't need to 
      //matrix access "dy" so many times
      for(int ipoly=0; ipoly<npoly; ++ipoly) { //loop over multidimensional 
	//monomials in the derivative polynomial
	double term=flycoef(ipoly,0); //this is workspace to evaluate
	//the current multidimensional monomial
	int nmult=flypoly(0,ipoly); //this is the total order of the current
	//multidimensional monomial
	for(int imult=1; imult<=nmult; ++imult)
	  term*=xr(flypoly(imult,ipoly),ipt);
	sumofterms+=term;
      }
      dy(ider,ipt)=sumofterms;
    }
  }
  return dy;
}



//this function assumes that flypoly (the fast and "on the fly" polynomial representation) is already known and uses it to evaluate g at xr.
MtxDbl& evaluate_flypoly_basis(MtxDbl& g, const MtxInt& flypoly, const MtxDbl& xr) {
  int npts=xr.getNCols();
  int npoly=flypoly.getNCols();
#ifdef __SURFPACK_ERR_CHECK__
  assert((0<npts)&&(0<xr.getNRows())&&(0<npoly)&&(0<flypoly.getNRows()));
  {
    int maxorder=flypoly(0,0);
    for(int i=1; i<npoly; ++i)
      maxorder=(flypoly(0,i)<=maxorder)?maxorder:flypoly(0,i);
    assert(maxorder<flypoly.getNRows());
  }
#endif
  g.newSize(npoly,npts);
  for(int ipt=0; ipt<npts; ++ipt)
    for(int ipoly=0; ipoly<npoly; ++ipoly) {
      int nmult=flypoly(0,ipoly);
      double tempg=1.0;
      for(int imult=1; imult<=nmult; ++imult)
	tempg*=xr(flypoly(imult,ipoly),ipt);
      g(ipoly,ipt)=tempg;
    }
  return g;
}

//this function determines flypoly (the "on the fly" polynomial representation) from poly and then evaluates g at xr using flypoly
MtxDbl& evaluate_poly_basis(MtxDbl& g, MtxInt& flypoly, const MtxInt& poly, const MtxDbl& xr) {
  int nvars=poly.getNRows();
  int npoly=poly.getNCols();
#ifdef __SURFPACK_ERR_CHECK__
  assert(nvars==xr.getNRows());
#endif  
  int maxorder=0;
  for(int ipoly=0; ipoly<npoly; ++ipoly) {
    int thisorder=poly(0,ipoly);
    for(int ivar=1; ivar<nvars; ++ivar)
      thisorder+=poly(ivar,ipoly);
    maxorder=(thisorder<=maxorder)?maxorder:thisorder;
  }
  poly_to_flypoly(flypoly, poly, maxorder);
  return evaluate_flypoly_basis(g, flypoly, xr);
}



MtxDbl& evaluate_poly_der_basis(MtxDbl& dg, MtxInt& flypoly, MtxDbl& flycoef, 
				const MtxInt& poly, const MtxInt& der, 
				const MtxDbl& xr) {
  int nvars=poly.getNRows();
  int npoly=poly.getNCols();
  int nder=der.getNCols();
  int npts=xr.getNCols();
#ifdef __SURFPACK_ERR_CHECK__
  assert((0<nvars)&&(nvars==der.getNRows())&&(nvars==xr.getNRows())&&
	 (0<npoly)&&(0<nder)&&(0<npts));
#endif

  int maxorder=0;
  for(int ipoly=0; ipoly<npoly; ++ipoly) {
    int thisorder=poly(0,ipoly);
    for(int ivar=1; ivar<nvars; ++ivar)
      thisorder+=poly(ivar,ipoly);
    maxorder=(thisorder<=maxorder)?maxorder:thisorder;
  }

  flycoef.newSize(npoly,1);
  int ncol=nder*npts;
  dg.newSize(npoly,ncol); //layout in memory: polynomial is contigous 
  //(i.e. a single column), derivatives are outside that, points are outside 
  //that.  If der was a gradient, then using matlab notation...
  //dy=reshape(coef'*dg,nder,npts) would be a matrix with a gradient in each
  //column and different columns being different points.

  for(int ider=0; ider<nder; ++ider) {
    //we are deliberately not going to assign the result in contiguous 
    //order BECAUSE we only want to determine each derivative of the 
    //polynomial once, instead of once for each point, and we don't want 
    //to have to allocate enough memory to store all the derivatives of 
    //all the polynomials at the same time.

    //determine the flypoly representation of the polynomial that is the 
    //der(:,ider) mixed multidimensional derivative of the polynomial in poly
    for(int ipoly=0; ipoly<npoly; ++ipoly)
      flycoef(ipoly,0)=1.0;
    poly_der_to_flypoly(flypoly, flycoef, poly, der, ider, maxorder);
   
    int ipt, icol;
    for(ipt=0, icol=ider; ipt<npts; ++ipt, icol+=nder) 
      for(int ipoly=0; ipoly<npoly; ++ipoly) {
	double tempdg=flycoef(ipoly,0);
	int nmult=flypoly(0,ipoly);
	for(int imult=1; imult<=nmult; ++imult)
	  tempdg*=xr(flypoly(imult,ipoly),ipt);	  
	dg(ipoly,icol)=tempdg;
      }
  }
  return dg;
}






/// Throw an exception if end-of-file has been reached 
void surfpack::checkForEOF(istream& is)
{
  if (is.eof()) {
    throw surfpack::io_exception("End of file reached unexpectedly.");
  }
}

/// Return true if the file specified by parameter file name has the extension specified by parameter extension
bool surfpack::hasExtension(const string& filename, const string extension)
{
  return (filename.find(extension) == filename.size() - extension.size());
}

} // end namespace nkm

#ifdef __SURFPACK_UNIT_TEST__
//needs surf77_config.h  surfpack_LAPACK_wrappers.h  surfpack_system_headers.h everything else is in the nkm directory
//g++ -o NKM_SurfPack_UnitMain NKM_SurfPack.cpp NKM_SurfMat.o /usr/lib64/libblas.so /usr/lib64/liblapack.so

using namespace nkm;
int main(){
  for(int n=0; n<5; n++)
    for(int k=0; k<=n; k++)
      printf("nchoosek(%d,%d)=%d\n",n,k,nchoosek(n,k));

  MtxInt poly; 
  multi_dim_poly_power(poly,4,4);
  int npoly=poly.getNCols();
  int nvarsr=poly.getNRows();
  printf("poly.getNCols()=%d\npoly=...\n",poly.getNCols());
  
  for(int ipoly=0; ipoly<npoly; ipoly++) {
    printf("%8d [",ipoly);
    for(int ixr=0; ixr<nvarsr; ixr++)
      printf(" %d",poly(ixr,ipoly));
    printf(" ]^T\n");
  }
  
  double pi=std::acos(0.0)*2.0;
  MtxDbl EulAng(6,1), Rot(4,4), RotShouldBe(4,4); 
  EulAng(0,0)=pi/6.0;
  EulAng(1,0)=pi/4.0;
  EulAng(2,0)=pi/3.0;
  EulAng(3,0)=pi/6.0;
  EulAng(4,0)=pi/4.0;
  EulAng(5,0)=pi/3.0;
  //This RotShouldBe was evaluated for these EulAng by matlab code that I know works
  RotShouldBe(0,0)= -0.057800215120219;
  RotShouldBe(1,0)=  0.353765877365274;
  RotShouldBe(2,0)=  0.768283046242747;
  RotShouldBe(3,0)=  0.530330085889911;
  RotShouldBe(0,1)= -0.695272228311384;
  RotShouldBe(1,1)= -0.649306566066328;
  RotShouldBe(2,1)=  0.035320133098213;
  RotShouldBe(3,1)=  0.306186217847897;
  RotShouldBe(0,2)=  0.647692568794007;
  RotShouldBe(1,2)= -0.414729655649473;
  RotShouldBe(2,2)= -0.183012701892219;
  RotShouldBe(3,2)=  0.612372435695795;
  RotShouldBe(0,3)= -0.306186217847897;
  RotShouldBe(1,3)=  0.530330085889911;
  RotShouldBe(2,3)= -0.612372435695795;
  RotShouldBe(3,3)=  0.500000000000000;
  gen_rot_mat(Rot,EulAng,4);
  printf("Rot=...\n");
  for(int ixr=0; ixr<nvarsr; ixr++){
    printf("    [ ");
    for(int jxr=0; jxr<nvarsr; jxr++) {
      printf("%9f ",Rot(ixr,jxr));
      if(std::fabs(Rot(ixr,jxr)-RotShouldBe(ixr,jxr))>1.0e-15) {
	printf("\n(%d,%d):Rot=%.16f RotShouldBe=%.16f\n",ixr,jxr,Rot(ixr,jxr),RotShouldBe(ixr,jxr));
	fflush(stdout);
	assert(Rot(ixr,jxr)==RotShouldBe(ixr,jxr));
      }
    }
    printf("]\n");
  }
  printf("We know Rot is ok because it didn't assert(false)\n");

  MtxDbl should_be_identity(4,4);
  matrix_mult(should_be_identity,Rot,Rot,0.0,1.0,'N','T');
  printf("should_be_identity=Rot*Rot^T=...\n");
  for(int ixr=0; ixr<nvarsr; ixr++){
    printf("    [ ");
    for(int jxr=0; jxr<nvarsr; jxr++) {
      printf("%9f ",should_be_identity(ixr,jxr));
    }
    printf("]\n");
  }

  matrix_mult(should_be_identity,Rot,Rot,0.0,1.0,'T','N');
  printf("should_be_identity=Rot^T*Rot=...\n");
  for(int ixr=0; ixr<nvarsr; ixr++){
    printf("    [ ");
    for(int jxr=0; jxr<nvarsr; jxr++) {
      printf("%9f ",should_be_identity(ixr,jxr));
    }
    printf("]\n");
  }


  int rotnvarsr=5;
  MtxDbl BaseAxis(rotnvarsr,2*rotnvarsr); BaseAxis.zero();
  for(int ixr=0; ixr<rotnvarsr; ixr++) {
    BaseAxis(ixr,ixr)=1.0;
    BaseAxis(ixr,ixr+rotnvarsr)=-1.0;
  }
  MtxDbl Axis;

  int noct=static_cast<int>(std::pow(2,rotnvarsr));
  MtxInt InOct(noct,1); InOct.zero();
  int NDV=(rotnvarsr*(rotnvarsr-1))/2; //nchoosek(nvarsr,2);
  MtxDbl lowerBounds(NDV,1); lowerBounds.zero();
  MtxDbl upperBounds(NDV,1);
  int mymod=1048576; //2^20 instead of 10^6 to be kind to the computer
  for(int jxr=0; jxr<NDV; jxr++) upperBounds(jxr,0)=pi;

  int nguess=2*rotnvarsr*noct*100;

  EulAng.newSize(NDV,1);
  for(int iguess=0; iguess<nguess; iguess++) {
    for(int jxr=0; jxr<NDV; jxr++)
      EulAng(jxr,0)=(rand()%mymod)*upperBounds(jxr,0)/mymod;
    gen_rot_mat(Rot,EulAng,rotnvarsr);
    matrix_mult(Axis,Rot,BaseAxis,0.0,1.0,'T');
    for(int k=0; k<2*rotnvarsr; k++){
      int ioct=(Axis(0,k)>=0.0);
      for(int jxr=1; jxr<rotnvarsr; jxr++)
	ioct+=(Axis(jxr,k)>=0)*(static_cast<int>(std::pow(2.0,jxr)));
      InOct(ioct,0)++;
    }
  }
  printf("relative # of axis endpoints per octant\n");
  for(int ioct=0; ioct<noct; ioct++)
    printf("  InOct(%d/%d)=%g\n",ioct,noct,InOct(ioct,0)/(2.0*(nguess*rotnvarsr/noct)));

  NDV=3;
  int npts=2;
  int ndeg=3;
  int nder=6;
  multi_dim_poly_power(poly,NDV,ndeg);
  npoly=poly.getNCols();
  assert(npoly==num_multi_dim_poly_coef(NDV,ndeg));

  
  MtxDbl xr(NDV,npts);
  MtxDbl y(1,npts), y2(1,npts), y3(1,npts), dy(nder,npts), dy2(1,nder*npts);
  y3.zero();
  for(int ixr=0; ixr<NDV; ++ixr) {
    xr(ixr,0)=1.0;
    xr(ixr,1)=ixr+1;
  }
  
  MtxDbl g(npoly,npts), g2(npoly,npts), dg(npoly,nder*npts);
  MtxDbl coef1(npoly,1), coefs(npoly,1), flycoef;
  MtxInt flypoly, flypoly2, der(NDV,nder);
  der.zero();
  der(0,1)=1;
  der(1,2)=1;
  der(2,3)=1;
  der(1,4)=2;
  der(0,5)=1;
  der(1,5)=2;


  poly_to_flypoly(flypoly,poly,ndeg);  
  //for(int ipoly=0; ipoly<npoly; ++ipoly) {    
  //for(int ifly=0; ifly<=ndeg; ++ifly)
  //  printf("%d ",flypoly(ifly,ipoly));
  //printf("\n");
  //}
  //printf("\n");
  evaluate_flypoly_basis(g, flypoly, xr);
  evaluate_poly_basis(g2, flypoly2, poly, xr);
  printf("Testing Polynomials\n");
  double tempdouble;
  for(int ipoly=0; ipoly<npoly; ++ipoly) {
    coef1(ipoly,0)=1.0;
    coefs(ipoly,0)=ipoly+1;
    int ixr=0;
    printf("g%d(x)=x%d^%d",ipoly,ixr,poly(ixr,ipoly));
    for(ixr=1; ixr<NDV; ++ixr)
      printf("*x%d^%d",ixr,poly(ixr,ipoly));
    printf("=1");
    int nmult=flypoly(0,ipoly);
    assert(flypoly(0,ipoly)==flypoly2(0,ipoly));
    for(int imult=1; imult<=nmult; ++imult) {
      printf("*x%d",flypoly(imult,ipoly));
      assert(flypoly(imult,ipoly)==flypoly2(imult,ipoly));
    }
    printf(";");
    for(int ipt=0; ipt<npts; ++ipt) {
      ixr=0;
      printf(" g%d(%g",ipoly,xr(ixr,ipt));
      tempdouble=coefs(ipoly,0)*std::pow(xr(ixr,ipt),poly(ixr,ipoly));
      for(ixr=1; ixr<NDV; ++ixr) {
	printf(",%g",xr(ixr,ipt));
	tempdouble*=std::pow(xr(ixr,ipt),poly(ixr,ipoly));
      }
      printf(")=%g;",g(ipoly,ipt));
      assert(g(ipoly,ipt)==g2(ipoly,ipt));
      y3(0,ipt)+=tempdouble;
    }
    printf("\n");
  }
  matrix_mult(y,coefs,g,0.0,1.0,'T','N');
  evaluate_poly(y2, flypoly, poly, coefs, xr);
  assert((y.getNRows()==1)&&(y.getNCols()==npts)&&
	 (y2.getNRows()==1)&&(y2.getNCols()==npts)&&
	 (y3.getNRows()==1)&&(y3.getNCols()==npts));
  for(int ipt=0; ipt<npts; ++ipt) {
    printf("ipt=%d y=%g y2=%g y3=%g\n",ipt,y(0,ipt),y2(0,ipt),y3(0,ipt));
    assert((y(0,ipt)==y3(0,ipt))&&(y2(0,ipt)==y3(0,ipt)));
  }
  printf("we know that polynomials evaluated correcty because it didn't fail the assert assert.\n");
  
  printf("\n Begin Testing Polynomial Derivatives\n");
  evaluate_poly_der(dy, flypoly, flycoef, poly, der, coef1, xr);
  evaluate_poly_der_basis(dg, flypoly, flycoef, poly, der, xr);
  matrix_mult(dy2,coef1,dg,0.0,1.0,'T','N');
  dy2.reshape(nder,npts);

  for(int ider=0; ider<nder; ++ider) {
    printf("ider=%d\n",ider);
    flycoef.copy(coef1);
    poly_der_to_flypoly(flypoly, flycoef, poly, der, ider, ndeg);    
    evaluate_flypoly_basis(g, flypoly, xr); //this g is really a dg but 
    //only for derivative ider
    matrix_mult(y,flycoef,g,0.0,1.0,'T','N');  //this y really holds a dy    
    for(int ipt=0; ipt<npts; ++ipt)
      assert((dy(ider,ipt)==y(0,ipt))&&(dy2(ider,ipt)==y(0,ipt)));
    int totalderorder=0;    
    for(int ixr=0; ixr<NDV; ++ixr)
      totalderorder+=der(ixr,ider);
    for(int ipoly=0; ipoly<npoly; ++ipoly) {
      int ixr=0;
      printf("d^%d/dx^[%d",totalderorder,der(ixr,ider));
      for(ixr=1; ixr<NDV; ++ixr)
	printf(",%d",der(ixr,ider));
      ixr=0;
      printf("](x%d^%d",ixr,poly(ixr,ipoly));
      for(ixr=1; ixr<NDV; ++ixr)
	printf("*x%d^%d",ixr,poly(ixr,ipoly));
      printf(")=%g",flycoef(ipoly,0));
      int nmult=flypoly(0,ipoly);
      for(int imult=1; imult<=nmult; ++imult)
	printf("*x%d",flypoly(imult,ipoly));
      printf("\n");
    }
    printf("\n");
  }
  return 0;
}
#endif

