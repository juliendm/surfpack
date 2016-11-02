#include "NKM_SurfData.hpp"
#include "NKM_SurfPack.hpp"
#include <iostream>
#include <iomanip>

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT(nkm::SurfData) 
BOOST_CLASS_EXPORT(nkm::SurfDataScaler) 
#endif

namespace nkm {

using namespace std;
using std::ios;
using std::istringstream;
using std::istream;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;


/***********************************************************************/
/***********************************************************************/
/**** Unit Test functions for SurfData start here                   ****/
/***********************************************************************/
/***********************************************************************/

#ifdef __SURFDATA_SCALING_UNIT_TEST__
  //this will only work when void SurfData_scaling_unit_test() is a friend of the nkm::SurfData class.  void SurfData_scaling_unit_test() is only a friend of the nkm::SurfData class if __SURFDATA_SCALING_UNIT_TEST__ is #define'd within the NKM_SurfData.hpp file (otherwise only nkm::SurfDataScaler will be the friend of the nkm::SurfData class and only nkm::SurfPackModel 's will be friends of the nkm::SurfDataScaler (and therefore be able to access the scaling functionality)) 

void SurfData_scaling_unit_test(){
  int npts=100;
  int nvarsr=7; //this is hard coded below
  int nout=3; //this is hard coded below
  MtxInt lockxr(nvarsr,1);
  MtxInt sortedlockxr(nvarsr,2); //what it should be after being sorted
  //original     sorted into order    original position
  lockxr(0,0)=1; sortedlockxr(0,0)=0; sortedlockxr(0,1)=1; 
  lockxr(1,0)=0; sortedlockxr(1,0)=1; sortedlockxr(1,1)=0;
  lockxr(2,0)=2; sortedlockxr(2,0)=1; sortedlockxr(2,1)=3;
  lockxr(3,0)=1; sortedlockxr(3,0)=1; sortedlockxr(3,1)=5; 
  lockxr(4,0)=2; sortedlockxr(4,0)=2; sortedlockxr(4,1)=2;
  lockxr(5,0)=1; sortedlockxr(5,0)=2; sortedlockxr(5,1)=4;
  lockxr(6,0)=3; sortedlockxr(6,0)=3; sortedlockxr(6,1)=6; 
  double SingDimValue=-5.3; //real input dimension 3 is singular
  MtxDbl minmaxd(nvarsr,2);
  MtxDbl minmaxD(nvarsr,2);
  MtxDbl XR(nvarsr,npts);  
  MtxDbl Y(nout,npts);
  MtxDbl domain_new(nvarsr,2);
  domain_new(0,0) = 0.3; domain_new(0,1)= 0.7;
  domain_new(1,0) =-0.3; domain_new(1,1)= 5.0;
  domain_new(2,0) = 5.0; domain_new(2,1)= 8.3;
  domain_new(3,0) =-7.0; domain_new(3,1)=-4.9;
  domain_new(4,0) = 7.0; domain_new(4,1)= 9.4;
  domain_new(5,0) = 4.9; domain_new(5,1)= 7.0;
  domain_new(6,0) =-2.0; domain_new(6,1)= 2.0;

  int maxi, mini, modi;
  double tempd, mind, maxd;
  for(int ixr=0; ixr<nvarsr; ixr++) {
    mini=rand()%10+5;
    maxi=rand()%20+20;
    modi=5*(maxi-mini);
    for(int ipt=0; ipt<npts; ipt++)
      XR(ixr,ipt)=(rand()%modi)*0.2+mini;
  }
  //make one of the real input dimensions (dimension 3) singular
  for(int ipt=0; ipt<npts; ipt++)
    XR(3,ipt)=SingDimValue;

  for(int iy=0; iy<nout; iy++) {
    mini=rand()%10+5;
    maxi=rand()%20+20;
    modi=5*(maxi-mini);
    for(int ipt=0; ipt<npts; ipt++)
      Y(iy,ipt)=(rand()%modi)*0.2+mini;
  }


  //test1 is for default scaling WITHOUT locking aspect ratios of dimensions 
  SurfData test1a(XR,Y); test1a.scaleToDefault();
  SurfData test1b=test1a; test1b.unScale();
  SurfData test1c=test1b; test1c.scaleToDefault();
  SurfData test1d(XR,Y); //leave it unscaled
  SurfData test1e=test1a; test1e.scaleXrToDomain(domain_new);
  SurfData test1f=test1d; test1f.scaleXrToDomain(domain_new);
  SurfData test1g=test1d; test1g.scaleToFactors(test1e.unscalexr,test1e.unscaley);
  SurfData test1h=test1a; test1h.scaleXrToFactor(test1e.unscalexr);

  //test the xr scaling
  for(int ixr=0; ixr<nvarsr; ixr++) {
    if(ixr==3){
      //check to make sure that dimension 3 is flagged as being singular
      assert((test1a.unscalexr(3,0)==-1.0)&&
	     (test1a.unscalexr(3,1)==SingDimValue));
      for(int ipt=0; ipt<npts; ipt++)
	assert(test1a.xr(3,ipt)==0.0);
    }
    else{
      assert(test1a.unscalexr(ixr,0)>0.0);
      mind=maxd=test1a.xr(ixr,0);
      for(int ipt=1; ipt<npts; ipt++){
	if(test1a.xr(ixr,ipt)<mind) mind=test1a.xr(ixr,ipt);
	if(test1a.xr(ixr,ipt)>maxd) maxd=test1a.xr(ixr,ipt);
      }
      //assert(if_close_enough(mind,-0.5)&&if_close_enough(maxd,0.5)); //commented because the tolerance on if_close_enough was changed to be a lot smaller
      mind=maxd=test1e.xr(ixr,0);
      for(int ipt=1; ipt<npts; ipt++){
	if(test1e.xr(ixr,ipt)<mind) mind=test1e.xr(ixr,ipt);
	if(test1e.xr(ixr,ipt)>maxd) maxd=test1e.xr(ixr,ipt);
      }
      //assert(if_close_enough(mind,domain_new(ixr,0))&&
      //     if_close_enough(maxd,domain_new(ixr,1))); //commented because the tolerance on if_close_enough was changed to be a lot smaller
    }
    assert(if_close_enough(test1a.unscalexr(ixr,0),test1c.unscalexr(ixr,0))&&
	   if_close_enough(test1a.unscalexr(ixr,1),test1c.unscalexr(ixr,1))&&
	   if_close_enough(test1b.unscalexr(ixr,0),test1d.unscalexr(ixr,0))&&
	   if_close_enough(test1b.unscalexr(ixr,1),test1d.unscalexr(ixr,1))&&
	   (test1d.unscalexr(ixr,0)==1.0)&&(test1d.unscalexr(ixr,1)==0.0)&&
	   if_close_enough(test1e.unscalexr(ixr,0),test1f.unscalexr(ixr,0))&&
	   if_close_enough(test1e.unscalexr(ixr,1),test1f.unscalexr(ixr,1))&&
	   if_close_enough(test1e.unscalexr(ixr,0),test1g.unscalexr(ixr,0))&&
	   if_close_enough(test1e.unscalexr(ixr,1),test1g.unscalexr(ixr,1)));
	   
    for(int ipt=0; ipt<npts; ipt++)
      assert(if_close_enough(test1a.xr(ixr,ipt),test1c.xr(ixr,ipt))&&
	     if_close_enough(XR(ixr,ipt),test1b.xr(ixr,ipt))&&
	     (XR(ixr,ipt)==test1d.xr(ixr,ipt))&&
	     if_close_enough(test1e.xr(ixr,ipt),test1f.xr(ixr,ipt))&&
	     if_close_enough(test1e.xr(ixr,ipt),test1g.xr(ixr,ipt))&&
	     if_close_enough(test1e.xr(ixr,ipt),test1h.xr(ixr,ipt)));
  }
  
  //test the y scaling
  for(int iy=0; iy<nout; iy++) {
    mind=maxd=test1a.y(iy,0);
    for(int ipt=1; ipt<npts; ipt++){
      if(test1a.y(iy,ipt)<mind) mind=test1a.y(iy,ipt);
      if(test1a.y(iy,ipt)>maxd) maxd=test1a.y(iy,ipt);
    }
    //assert(if_close_enough(mind,-0.5)&&if_close_enough(maxd,0.5)); //commented because the tolerance on if_close_enough was changed to be a lot smaller
    assert(if_close_enough(test1a.unscaley(iy,0),test1c.unscaley(iy,0))&&
	   if_close_enough(test1a.unscaley(iy,1),test1c.unscaley(iy,1))&&
	   if_close_enough(test1b.unscaley(iy,0),test1d.unscaley(iy,0))&&
	   if_close_enough(test1b.unscaley(iy,1),test1d.unscaley(iy,1))&&
	   (test1d.unscaley(iy,0)==1.0)&&(test1d.unscaley(iy,1)==0.0));

    for(int ipt=0; ipt<npts; ipt++){
      assert(if_close_enough(test1a.y(iy,ipt),test1c.y(iy,ipt)));
      assert(if_close_enough(Y(iy,ipt),test1b.y(iy,ipt)));
      assert((Y(iy,ipt)==test1d.y(iy,ipt)));
      assert(if_close_enough(test1a.y(iy,ipt),test1g.y(iy,ipt)));
    }

  }

  //test2 is for default scaling WITH locking of aspect ratios of groups of dimensions
  SurfData test2a(lockxr,XR,Y); test2a.scaleToDefault();
  SurfData test2b=test2a; test2b.unScale();
  SurfData test2c=test2b; test2c.scaleToDefault();
  SurfData test2d(lockxr,XR,Y); //leave it unscaled

  for(int j=0; j<2; j++)
    for(int ixr=0; ixr<nvarsr; ixr++)
      assert((sortedlockxr(ixr,j)==test2a.lockxr(ixr,j))&&
	     (sortedlockxr(ixr,j)==test2b.lockxr(ixr,j))&&
	     (sortedlockxr(ixr,j)==test2c.lockxr(ixr,j))&&
	     (sortedlockxr(ixr,j)==test2d.lockxr(ixr,j)));

  //test the xr scaling
  for(int ixr=0; ixr<nvarsr; ixr++) {
    if(ixr==3){
      //check to make sure that dimension 3 is flagged as being singular
      assert((test2a.unscalexr(3,0)==-1.0)&&
	     (test2a.unscalexr(3,1)==SingDimValue));
      for(int ipt=0; ipt<npts; ipt++)
	assert(test2a.xr(3,ipt)==0.0);
    }
    else assert(test2a.unscalexr(ixr,0)>0.0);

    minmaxd(ixr,0)=minmaxd(ixr,1)=test2a.xr(ixr,0);
    minmaxD(ixr,0)=minmaxD(ixr,1)=XR(ixr,0);
    for(int ipt=1; ipt<npts; ipt++){
      if(test2a.xr(ixr,ipt)<minmaxd(ixr,0)) minmaxd(ixr,0)=test2a.xr(ixr,ipt);
      if(test2a.xr(ixr,ipt)>minmaxd(ixr,1)) minmaxd(ixr,1)=test2a.xr(ixr,ipt);
      if(XR(ixr,ipt)<minmaxD(ixr,0)) minmaxD(ixr,0)=XR(ixr,ipt);
      if(XR(ixr,ipt)>minmaxD(ixr,1)) minmaxD(ixr,1)=XR(ixr,ipt);
    }
  
    assert(if_close_enough(test2a.unscalexr(ixr,0),test2c.unscalexr(ixr,0))&&
	   if_close_enough(test2a.unscalexr(ixr,1),test2c.unscalexr(ixr,1)));
    for(int ipt=0; ipt<npts; ipt++)
      assert(if_close_enough(test2a.xr(ixr,ipt),test2c.xr(ixr,ipt))&&
	     if_close_enough(XR(ixr,ipt),test2b.xr(ixr,ipt))&&
	     (XR(ixr,ipt)==test2d.xr(ixr,ipt)));

  }
  assert(if_close_enough(1.0,(minmaxd(0,1)-minmaxd(0,0))*
			 (minmaxd(5,1)-minmaxd(5,0)))&&
	 if_close_enough((minmaxD(0,1)-minmaxD(0,0))/(minmaxd(0,1)-minmaxd(0,0)),
			 (minmaxD(5,1)-minmaxD(5,0))/(minmaxd(5,1)-minmaxd(5,0)))&&
	 if_close_enough(1.0,(minmaxd(2,1)-minmaxd(2,0))*
			 (minmaxd(4,1)-minmaxd(4,0)))&&
	 if_close_enough((minmaxD(2,1)-minmaxD(2,0))/(minmaxd(2,1)-minmaxd(2,0)),
			 (minmaxD(4,1)-minmaxD(4,0))/(minmaxd(4,1)-minmaxd(4,0)))&&
	 if_close_enough(1.0,(minmaxd(1,1)-minmaxd(1,0)))&&
	 if_close_enough(1.0,(minmaxd(6,1)-minmaxd(6,0)))&&
	 if_close_enough(0.0,(minmaxd(3,1)-minmaxd(3,0))));
	 

  //test the y scaling
  for(int iy=0; iy<nout; iy++) {
    mind=maxd=test2a.y(iy,0);
    for(int ipt=1; ipt<npts; ipt++){
      if(test2a.y(iy,ipt)<mind) mind=test2a.y(iy,ipt);
      if(test2a.y(iy,ipt)>maxd) maxd=test2a.y(iy,ipt);
    }
    //assert(if_close_enough(mind,-0.5)&&if_close_enough(maxd,0.5)); //commented because the tolerance on if_close_enough was changed to be a lot smaller
    assert(if_close_enough(test2a.unscaley(iy,0),test2c.unscaley(iy,0))&&
	   if_close_enough(test2a.unscaley(iy,1),test2c.unscaley(iy,1)));
    for(int ipt=0; ipt<npts; ipt++)
      assert(if_close_enough(test2a.y(iy,ipt),test2c.y(iy,ipt))&&
	     if_close_enough(Y(iy,ipt),test2b.y(iy,ipt)));
  }

  printf("class SurfData; passed the scaling unit tests\n");
  return;
}

#endif

/***********************************************************************/
/***********************************************************************/
/**** Unit Test functions for SurfData end here                     ****/
/***********************************************************************/
/***********************************************************************/



/***********************************************************************/
/***********************************************************************/
/**** SurfData member functions start here                          ****/
/***********************************************************************/
/***********************************************************************/
  
///assign all mixed partial derivatives of exact total order "der_order" for output jy, der_order can be up to 1 order higher than that already stored for this variable (e.g. you can't assign hessians until after you have assigned gradient but you can re-assign gradients after you have assigned hessians) you can use this function to assign function zeroth order derivatives (i.e. the value itself not a derivative) but it will be stored in the MtxDbl variable y rather than derY
void SurfData::putDerY(const MtxDbl& dny, int der_order, int iy) {
  if(iy==-99999) //default value
    iy=iout; //means use the user set default output (iout)
  assert((0<=iy)&&(iy<nout)&&(0<=der_order)&&(npts==dny.getNCols()));
  int nrows_dny_should_have=num_multi_dim_poly_coef(nvarsr,-der_order);
  assert((der_order<=derOrder(iy,0)+1)&&
	 (nrows_dny_should_have==dny.getNRows()));
  if(der_order>derOrder(iy,0)) {
    derY[iy].resize(der_order+1);
    derOrder(iy,0)=der_order;
  }
  if(der_order==0)
    y.putRows(dny,iy);
  else
    derY[iy][der_order].copy(dny);
  return;
}

///retrieve a copy of all mixed partial derivatives of exact total order "der_order" for output jy, you can use this function to retrieve zeroth order derivatives (i.e. the value itself not a derivative) 
MtxDbl& SurfData::getDerY(MtxDbl& dny, int der_order, int iy) const {
  if(iy==-99999) //default value
    iy=iout; //means use the user set default output (iout)
  assert((0<=iy)&&(iy<nout)&&(0<=der_order));
  assert(der_order<=derOrder(iy,0));
  if(der_order==0)
    y.getRows(dny,iy);
  else
    dny.copy(derY[iy][der_order]);
  return dny;
}

///assign all mixed partial derivatives up to total order "der_order" for output jy, this includes the zeroth order derivative (the value itself not a derivative) but the zeroth order derivative will be stored in the separate variable "y"
void SurfData::putUpToDerY(const MtxDbl& dny, int der_order, int iy) {
  if(iy==-99999) //default value
    iy=iout; //means use the user set default output (iout)
  assert((0<=iy)&&(iy<nout)&&(0<=der_order)&&(npts==dny.getNCols()));
  int nrows_dny_should_have=num_multi_dim_poly_coef(nvarsr,der_order);
  assert(nrows_dny_should_have==dny.getNRows());
  MtxDbl temp_vector(1,npts);
  dny.getRows(temp_vector,0);
  y.putRows(temp_vector,iy);
  int nrows_so_far=1;
  if(der_order>derOrder(iy,0)) {
    derY[iy].resize(der_order+1);
    derOrder(iy,0)=der_order;
  }
  MtxInt irows;
  for(int ider=1; ider<=der_order; ++ider) {
    int nder_this_exact_total_order=num_multi_dim_poly_coef(nvarsr,-ider);
    irows.newSize(nder_this_exact_total_order,1);
    for(int irow=0; irow<nder_this_exact_total_order; ++irow, ++nrows_so_far)
      irows(irow,0)=nrows_so_far;
    dny.getRows(derY[iy][ider],irows);
  }
  return;
}
 
///retrieve a copy of all mixed partial derivatives up to total order "der_order" for output iy, this includes function values themselves (zeroth order derivatives) 
MtxDbl& SurfData::getUpToDerY(MtxDbl& dny, int der_order, int iy) const {
  if(iy==-99999) //default value
    iy=iout; //means use the user set default output (iout)
  assert((0<=iy)&&(iy<nout)&&(0<=der_order));
  assert(der_order<=derOrder(iy,0));
  int nrows_dny_should_have=num_multi_dim_poly_coef(nvarsr,der_order);
  dny.reshape(nrows_dny_should_have,npts);
  if(nout==1)
    dny.putRows(y,0);
  else{
    MtxDbl temp_vector(1,npts);
    y.getRows(temp_vector,iy);
    dny.putRows(temp_vector,0);
  }
  int nrows_so_far=1;
  MtxInt irows;
  for(int ider=1; ider<=der_order; ++ider) {
    int nder_this_exact_total_order=derY[iy][ider].getNRows();
    irows.newSize(nder_this_exact_total_order,1);
    for(int irow=0; irow<nder_this_exact_total_order; ++irow, ++nrows_so_far)
      irows(irow,0)=nrows_so_far;
    dny.putRows(derY[iy][ider],irows);
  }
  assert(nrows_so_far==nrows_dny_should_have);
  return dny;
}
 
///a constructor for when there are real input variables that we don't want the model to group scale and there are no integer input variables, this will result in the model scaling its real input variables and output variable(s) to a hypercube of volume 1
  SurfData::SurfData(const MtxDbl& XR, const MtxDbl& Y, int iout_set) : npts(XR.getNCols()), nvarsr(XR.getNRows()), nvarsi(0), nout(Y.getNRows()), iout(iout_set), derOrder(nout,1), derY(nout), ifHaveMinMaxXr(false)
{
  assert(Y.getNCols()==npts);
  
  if(0<npts) {
    assert((0<=iout)&&(iout<nout));
    xr=XR;
    y=Y;
    dontScale();
    derOrder.zero();

  }
  else{
    iout=0;
    std::cerr << "Warning: SurfData() constructor was passed empty data matrices!!!" << std::endl;
  }
  defaultLabels();
  return;
}

  
///a constructor for when you pass in real inputs, and (arbitrarily high order) derivatives of the output with respect them, as well as the output, but don't pass in integer input variables.
  SurfData::SurfData(const MtxDbl& XR, const MtxDbl& Y, const MtxInt& der_order_in, const std::vector<std::vector<MtxDbl> > & derY_in, int iout_set) : npts(XR.getNCols()), nvarsr(XR.getNRows()), nvarsi(0), nout(Y.getNRows()), iout(iout_set), derOrder(der_order_in), derY(derY_in), ifHaveMinMaxXr(false), xr(XR), y(Y)
{
  assert(Y.getNCols()==npts);
  
  if(0<npts) {
    assert((0<=iout)&&(iout<nout)&&(1==derOrder.getNCols())&&
	   (nout==derOrder.getNRows())&&(nout==static_cast<int>(derY.size())));
    for(int iy=0; iy<nout; ++iy) {
      assert(derOrder(iy,0)>=0);
      if(1<=derOrder(iy,0)) {
	assert(static_cast<int>(derY[iy].size())==derOrder(iy,0)+1);
	for(int ider=1; ider<=derOrder(iy,0); ++ider)
	  assert((derY[iy][ider].getNCols()==npts)&&
		 (derY[iy][ider].getNRows()==
		  num_multi_dim_poly_coef(nvarsr,-ider)));
      }
    }

    //xr=XR;
    //y=Y;
    dontScale();

  }
  else{
    iout=0;
    std::cerr << "Warning: SurfData() constructor was passed empty data matrices!!!" << std::endl;
  }
  defaultLabels();
  return;
}




///a constructor for when there are real input variables that we want the model to group scale (if it is appropriate to the model) and there are no integer input variables.  If it is appropriate to the model to do so, the model will automatically scale the real input variables to a hyper-rectangle of volume 1 and the output variable(s) to a hypercube of volume 1
  SurfData::SurfData(const MtxInt& LOCKXR, const MtxDbl& XR, const MtxDbl& Y, int iout_set) : npts(XR.getNCols()), nvarsr(XR.getNRows()), nvarsi(0), nout(Y.getNRows()), iout(iout_set), derOrder(nout,1), derY(nout), ifHaveMinMaxXr(false)
{
  assert((LOCKXR.getNCols()==1)&&(LOCKXR.getNRows()==nvarsr)&&
	 (Y.getNCols()==npts));
  
  if(0<npts) {
    assert((0<=iout)&&(iout<nout));
    xr=XR;
    y=Y;
    lockxr.newSize(nvarsr,2);
    for(int ixr=0; ixr <nvarsr; ixr++) {
      lockxr(ixr,0)=LOCKXR(ixr,0);
      lockxr(ixr,1)=ixr; //need to retain original order of dimensions, because the rows of lockxr are going to be sorted by elements in column 0 of lockxr
    }
    //you want us to scale by groups, and although there are group flags for each dimension (lockxr), those dimensions aren't arranged so that member dimensions of a group are listed sequentionally so we need to sort lockxr by groups.
    lockxr.sortRows();
    dontScale();

    derOrder.zero();
  }
  else{
    iout=0;
    std::cerr << "Warning: SurfData() constructor was passed empty data matrices!!!" << std::endl;
  }
  defaultLabels();
  return;
}


///a constructor for when there is real input variables (that we don't want to group scale) and integer input variables.  Models that use the default scaling will automatically scale the real input variables and output variable(s) to a hypercube of volume 1, they don't scale the integer input variables 
SurfData::SurfData(const MtxDbl& XR, const MtxInt& XI, const MtxDbl& Y, int iout_set) : npts(XR.getNCols()), nvarsr(XR.getNRows()), nvarsi(XI.getNRows()), nout(Y.getNRows()), iout(iout_set), derOrder(nout,1), derY(nout), ifHaveMinMaxXr(false)
{
  assert((XI.getNCols()==npts)&&(Y.getNCols()==npts));

  if(0<npts) {
    assert((0<=iout)&&(iout<nout));
    xr=XR;
    y=Y;
    dontScale();
    derOrder.zero();

    xi=XI; //note that integer input variables will never be scaled
  }
  else{
    iout=0;
    std::cerr << "Warning: SurfData() constructor was passed empty data matrices!!!" << std::endl;
  }
  defaultLabels();
  return;
}


///a constructor for when there are real input variables (that we DO want the model to group scale if it is appropriate for the model) and integer input variables.  Models that use the default scaling will automatically scale the real input variables to a hyper-rectangle of volume 1 and the output variable(s) to a hypercube of volume 1, the integer input variables won't be scaled
  SurfData::SurfData(const MtxInt& LOCKXR, const MtxDbl& XR, const MtxInt& XI, const MtxDbl& Y, int iout_set) : npts(XR.getNCols()), nvarsr(XR.getNRows()), nvarsi(XI.getNRows()), nout(Y.getNRows()), iout(iout_set), derOrder(nout,1), derY(nout), ifHaveMinMaxXr(false)
{
  assert((LOCKXR.getNCols()==1)&&(LOCKXR.getNRows()==nvarsr)&&
	 (XI.getNCols()==npts)&&(Y.getNCols()==npts));

  if(0<npts) {
    assert((0<=iout)&&(iout<nout));
    xr=XR;
    y=Y;
    lockxr.newSize(nvarsr,2);
    for(int ixr=0; ixr<nvarsr; ixr++) {
      lockxr(ixr,0)=LOCKXR(ixr,0);
      lockxr(ixr,1)=ixr; //need to retain original order of dimensions, because the rows of lockxr are going to be sorted by elements in column 0 of lockxr
    }
    //you want us to scale by groups, and although there are group flags for each dimension (lockxr), those dimensions aren't arranged so that member dimensions of a group are listed sequentionally so we need to sort lockxr by groups.
    lockxr.sortRows();
    dontScale();

    derOrder.zero();

    xi=XI; //note that integer input variables will never be scaled
  }
  else{
    iout=0;
    std::cerr << "Warning: SurfData() constructor was passed empty data matrices!!!" << std::endl;
  }
  defaultLabels();
  return;
}

///a constructor that reads data from a file for when you don't want the model to group scale the real input variables
SurfData::SurfData(const string& filename, int nvarsr_in, int nvarsi_in, int nout_in, int iout_in, int der_order_in, int skip_columns) : nvarsr(nvarsr_in), nvarsi(nvarsi_in), nout(nout_in), iout(iout_in), derOrder(nout,1), derY(nout), ifHaveMinMaxXr(false)
{
  //std::cout << "filename=[" << filename <<"]\nnvarsr=" << nvarsr 
  //   <<"\nnvarsi=" << nvarsi <<"\nnout="<< nout << "\niout=" << iout 
  //   <<"\nskip_columns=" << skip_columns << "\nifscale=[" << ifscale <<"]\n";

  assert((0<nvarsr)&&(0<=nvarsi)&&(0<nout)&&(0<=iout));

  for(int iy=0; iy<nout; ++iy) {
    derOrder(iy,0)=der_order_in;
    derY[iy].resize(derOrder(iy,0)+1);
  }

  //lockxr is an empy matrix
  dontScale(); //set the initial scaling to identity
  read(filename,skip_columns);

  return;
}

SurfData::SurfData(const string& filename, int nvarsr_in, int nvarsi_in, int nout_in, int iout_in, const MtxInt& der_order_in, int skip_columns) : nvarsr(nvarsr_in), nvarsi(nvarsi_in), nout(nout_in), iout(iout_in), derOrder(der_order_in), derY(nout), ifHaveMinMaxXr(false)
{
  //std::cout << "filename=[" << filename <<"]\nnvarsr=" << nvarsr 
  //   <<"\nnvarsi=" << nvarsi <<"\nnout="<< nout << "\niout=" << iout 
  //   <<"\nskip_columns=" << skip_columns << "\nifscale=[" << ifscale <<"]\n";

  assert((0<nvarsr)&&(0<=nvarsi)&&(0<nout)&&(0<=iout));

  for(int iy=0; iy<nout; ++iy) 
    derY[iy].resize(derOrder(iy,0)+1);

  //lockxr is an empy matrix
  dontScale(); //set the initial scaling to identity
  read(filename,skip_columns);

  return;
} 


///a constructor that reads data from a file for when you DO want the model to group scale the real input variables if it is appropriate for the model
SurfData::SurfData(const string& filename, int nvarsr_in, int nvarsi_in, int nout_in, int iout_in, int der_order_in, int skip_columns, const MtxInt& LOCKXR) : nvarsr(nvarsr_in), nvarsi(nvarsi_in), nout(nout_in), iout(iout_in), derOrder(nout,1), derY(nout), ifHaveMinMaxXr(false)
{
  assert((0<nvarsr)&&(0<=nvarsi)&&(0<nout)&&(0<=iout)&&
	 (LOCKXR.getNCols()==1)&&(LOCKXR.getNRows()==nvarsr));

  lockxr.newSize(nvarsr,2);

  for(int iy=0; iy<nout; ++iy) {
    derOrder(iy,0)=der_order_in; 
    derY[iy].resize(derOrder(iy,0)+1);
  }

  for(int ixr=0; ixr<nvarsr; ++ixr) {
    lockxr(ixr,0)=LOCKXR(ixr,0);
    lockxr(ixr,1)=ixr; //need to retain original order of dimensions, because the rows of lockxr are going to be sorted by elements in column 0 of lockxr
  }
  //you want us to scale by groups, and although there are group flags for each dimension (lockxr), those dimensions aren't arranged so that member dimensions of a group are listed sequentionally so we need to sort lockxr by groups.
  lockxr.sortRows();

  dontScale(); //set the initial scaling to identity  
  read(filename,skip_columns);

  return;
}

///a constructor that reads data from a file for when you DO want the model to group scale the real input variables if it is appropriate for the model
  SurfData::SurfData(const string& filename, int nvarsr_in, int nvarsi_in, int nout_in, int iout_in, const MtxInt& der_order_in, int skip_columns, const MtxInt& LOCKXR) : nvarsr(nvarsr_in), nvarsi(nvarsi_in), nout(nout_in), iout(iout_in), derOrder(der_order_in), derY(nout), ifHaveMinMaxXr(false)
{
  assert((0<nvarsr)&&(0<=nvarsi)&&(0<nout)&&(0<=iout)&&
	 (LOCKXR.getNCols()==1)&&(LOCKXR.getNRows()==nvarsr));

  for(int iy=0; iy<nout; ++iy) 
    derY[iy].resize(derOrder(iy,0)+1);

  lockxr.newSize(nvarsr,2);

  for(int ixr=0; ixr<nvarsr; ++ixr) {
    lockxr(ixr,0)=LOCKXR(ixr,0);
    lockxr(ixr,1)=ixr; //need to retain original order of dimensions, because the rows of lockxr are going to be sorted by elements in column 0 of lockxr
  }
  //you want us to scale by groups, and although there are group flags for each dimension (lockxr), those dimensions aren't arranged so that member dimensions of a group are listed sequentionally so we need to sort lockxr by groups.
  lockxr.sortRows();

  dontScale(); //set the initial scaling to identity  
  read(filename,skip_columns);

  return;
}

///copy constructor performs a deep copy
SurfData::SurfData(const SurfData& other) : npts(other.npts), nvarsr(other.nvarsr), nvarsi(other.nvarsi), nout(other.nout), iout(other.iout), derOrder(other.derOrder), derY(other.derY), ifHaveMinMaxXr(false), lockxr(other.lockxr), unscalexr(other.unscalexr), unscaley(other.unscaley), xr(other.xr), xi(other.xi), y(other.y)
 //effective c++ says to initialize rather than assign 
{
  //I don't know if vectors have copy constructors so use the assignment operator which I know does the right thing.
  xrLabels=other.xrLabels;
  xiLabels=other.xiLabels;
  yLabels =other.yLabels;
  //copy(other);
  return;
}

///deep copy constructor that keeps only one column of output
SurfData::SurfData(const SurfData& other, int iout_keep) : npts(other.npts), nvarsr(other.nvarsr), nvarsi(other.nvarsi), nout(1), iout(0), ifHaveMinMaxXr(false), lockxr(other.lockxr), unscalexr(other.unscalexr), xr(other.xr), xi(other.xi) //effective c++ says to initialize rather than assign 
{
  //printf("inside SurfData::SurfData(const SurfData& other, int iout_keep)\n");  fflush(stdout);

  assert((0<=iout_keep)&&(iout_keep<other.nout)); //we shouldn't "throw" exceptions during construction but there is no prohibition of asserts during construction
  if(iout_keep==-1)
    iout_keep=other.iout;
  y.newSize(nout,npts);
  //printf("other.getNOut()=%d other.unscaley.getNRows()=%d\n", other.getNOut(),other.unscaley.getNRows());  fflush(stdout);
  other.y.getRows(y,iout_keep);
  other.unscaley.getRows(unscaley,iout_keep);

  derOrder.newSize(nout,1);
  derOrder(0,0)=other.derOrder(iout_keep,0);

  derY.resize(nout);
  derY[0]=other.derY[iout_keep];

  //printf("about to get labels\n"); fflush(stdout);
  xrLabels=other.xrLabels;
  xiLabels=other.xiLabels;
  //printf("got xrLabels and xiLabels, about to get yLabels\n"); fflush(stdout);
  yLabels.resize(1);
  yLabels[0]=other.yLabels[iout_keep];

  //printf("leaving SurfData::SurfData(const SurfData& other, int iout_keep)\n"); fflush(stdout);
  return;
}

void SurfData::clear(){
  npts=nvarsr=nvarsi=nout=iout=0;
  xr.clear();
  xi.clear();
  y.clear();
  unscalexr.clear();
  unscaley.clear();
  lockxr.clear();
  xrLabels.clear();
  xiLabels.clear();
  yLabels.clear();
  derOrder.clear();
  derY.clear();
}


///make a deep copy
SurfData& SurfData::copy(const SurfData& other) {
  npts      =other.npts;
  nvarsr    =other.nvarsr;
  nvarsi    =other.nvarsi;
  nout      =other.nout;
  iout      =other.iout;
  derOrder  =other.derOrder;
  xr        =other.xr;
  xi        =other.xi;
  y         =other.y;
  derY      =other.derY;
  unscalexr =other.unscalexr;
  unscaley  =other.unscaley;
  lockxr    =other.lockxr;
  xrLabels  =other.xrLabels;
  xiLabels  =other.xiLabels;
  yLabels   =other.yLabels;
  return *this;
}



/// Returns true if file has .bspd extension, false if it has .spd extension. Otherwise, an exception is thrown.
bool SurfData::hasBinaryFileExtension(const string& filename) const
{
  if (nkm::surfpack::hasExtension(filename,".bspd")) {
    return true;
  } else if (nkm::surfpack::hasExtension(filename,".spd")) {
    return false;
  } else if (nkm::surfpack::hasExtension(filename,".dat")) {
    return false;
  } else {
    throw nkm::surfpack::io_exception(
      "Unrecognized filename extension.  Use .bspd or .spd"
    );
  }
}

/// Read a set of points into SurfData from either a text or binary file.  Opens file and calls either readText() or readBinary()
void SurfData::read(const string& filename, int skip_columns)
{
  // Open file in binary or text mode based on filename extension (.bspd or .spd)
  bool binary = hasBinaryFileExtension(filename);
  //std::cout << "binary=" << binary << std::endl;
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  if (!infile) {
    //std::cout << "couldn't open file" << std::endl;
    throw nkm::surfpack::file_open_failure(filename);
  } else if (binary) {
    std::cout << "attempting to open a binary file" << std::endl;
    assert(0);
    readBinary(infile,skip_columns);
  } else {
    //std::cout << "attempting to open a text file" << std::endl;
    readText(infile,skip_columns);
  }
  //std::cout << "done reading from file\n";
  // Object may have already been created
  infile.close();
}

void SurfData::write(const string& filename) const
{
  //if (mapping.empty()) {
  //ostringstream errormsg;
  //errormsg << "Cannot write SurfData object to stream."
  //     << "  No active data points." << std::endl;
  //throw bad_surf_data(errormsg.str());
  //}
  bool binary = hasBinaryFileExtension(filename);
  ofstream outfile(filename.c_str(), 
    (binary ? ios::out|ios::binary : ios::out));
  if (!outfile) {
    throw nkm::surfpack::file_open_failure(filename);
  } else if (binary) {
    std::cout << "attempting to write a binary file" << std::endl;
    assert(0);
    //writeBinary(outfile);
  } else {
    // Write the header and label info for .spd, not for .dat
    bool write_labels = nkm::surfpack::hasExtension(filename,".spd");
    writeText(outfile, write_labels);
  }
  outfile.close();
}


/// Set real input variable, xr, labels to xr0 xr1, etc.; integer input variable, xi, labels to xi0, xi1, etc.; output variable, y, labels to y0 y1, etc.
void SurfData::defaultLabels()
{
  //real input variable labels
  xrLabels.resize(nvarsr);
  for(int ixr=0; ixr<nvarsr; ++ixr) {
    ostringstream os;
    os << "xr" << ixr ;
    xrLabels[ixr] = os.str();
  }

  //integer input variable labels
  xiLabels.resize(nvarsi);
  for(int ixi=0; ixi<nvarsi; ++ixi) {
    ostringstream os;
    os << "xi" << ixi;
    xiLabels[ixi] = os.str();
  }

  //output variable labels
  yLabels.resize(nout);
  for(int iy=0; iy<nout; ++iy) {
    ostringstream os;
    os << "y" << iy ;
    yLabels[iy] = os.str();
  }
}

bool SurfData::readLabelsIfPresent(string single_line, int skip_columns)
{
  if(! ((single_line[0] == '%')||(single_line[0] == '#'))) {
    defaultLabels();
    return false;
  } else { // use custom labels

    single_line[0] = ' ';
    string dummy;
    xrLabels.resize(nvarsr); //this is a vector of strings
    xiLabels.resize(nvarsi); //this is a vector of strings
    yLabels.resize(nout); //this is a vector of strings
    istringstream is(single_line);

    //labels for columns the user told us to skip
    for(int i=0; i<skip_columns; ++i) {
      is >> dummy;
      if(dummy == "") {
	// Not enough heading names in the line of column headings
	// Use the default headings and return
	defaultLabels();
	return false;
      }
    }

    //real input variable labels
    for(int ixr=0; ixr<nvarsr; ++ixr) {
      is >> xrLabels[ixr];
      if(xrLabels[ixr] == "") { 
	// Not enough heading names in the line of column headings
	// Use the default headings and return
        defaultLabels();
        return false;
      }
    } 

    //integer input variable labels
    for(int ixi=0; ixi<nvarsi; ++ixi) {
      is >> xiLabels[ixi];
      if(xiLabels[ixi] =="") {
	// Not enough heading names in the line of column headings
	// Use the default headings and return
	defaultLabels();
	return false;
      }
    }

    //output variable labels
    for(int iy=0; iy<nout; ++iy) {
      is >> yLabels[iy];
      if (yLabels[iy] == "") { 
	// Not enough heading names in the line of column headings
	// Use the default headings and return
        defaultLabels();
        return false;
      }
      //if derivatives are in the file they have to have column headings and they have to be in the "right" order, but we don't need to read them because we can regenerate "standard" versions of them when it's time to right the file
      int nder_up_to_derOrder=num_multi_dim_poly_coef(nvarsr,derOrder(iy,0));
      for(int ii=1; ii<nder_up_to_derOrder; ++ii) {
	//start at 1 instead of zero because we already read in the function value
	is >> dummy;
	if (dummy == "") { 
	  // Not enough heading names in the line of column headings
	  // Use the default headings and return
	  defaultLabels();
	  return false;
	}
      }
	  
    } 

  } // use custom labels
  return true;
}

///read one point into SurfData from a single_line taken from the text .spd file, skip_columns defaults to zero
void SurfData::readPointText(int ipt, const string& single_line, 
			     int skip_columns)
{
  int nvarsr_read=0, nvarsi_read=0, nout_read=0, nskip_read=0, ider=0, jder=0;
  string dummy;
  try {
    // read the point as text
    istringstream streamline(single_line);

    //cout << "skip_columns=" << skip_columns << std::endl;

    //skip leading columns
    for(nskip_read=0; nskip_read<skip_columns; ++nskip_read) {
      // Throw an exception if there are fewer values on this line that
      // expected.
      nkm::surfpack::checkForEOF(streamline);
      streamline >> dummy;
      //cout << "dummy=" << dummy << std::endl;
    }

    //read in real input variables "xr"
    for(nvarsr_read=0; nvarsr_read<nvarsr; ++nvarsr_read) {
      // Throw an exception if there are fewer values on this line that
      // expected.
      nkm::surfpack::checkForEOF(streamline);
      streamline >> xr(nvarsr_read,ipt);
      //cout << "xr(" << nvarsr_read << "," << ipt << ")=" << xr(nvarsr_read,ipt) << std::endl;
    }

    //read in integer input variables "xi"
    for(nvarsi_read=0; nvarsi_read<nvarsi; ++nvarsi_read) {
      // Throw an exception if there are fewer values on this line that
      // expected.
      nkm::surfpack::checkForEOF(streamline);
      streamline >> xi(nvarsi_read,ipt);
    } 

    //read in output variables "y"
    for(nout_read=0; nout_read<nout; ++nout_read) {
      // Throw an exception if there are fewer values on this line that
      // expected.
      ider=jder=0;
      nkm::surfpack::checkForEOF(streamline);
      streamline >> y(nout_read,ipt);
      for(ider=1; ider<=derOrder(nout_read,0); ++ider) {
	int nder_this_exact_total_order=derY[nout_read][ider].getNRows();
	for(jder=0; jder<nder_this_exact_total_order; ++jder) {
	    nkm::surfpack::checkForEOF(streamline);
	    streamline >> derY[nout_read][ider](jder,ipt);
	}
      }
    } 
  } catch(nkm::surfpack::io_exception&) {
    std::cerr << "Bad SurfPoint: " << single_line 
	      << "\nExpected on this line: " 
	      << "\n  " << skip_columns << " leading columns to skip and"
	      << "\n  " << nvarsr << " real input variable(s) and"
	      << "\n  " << nvarsi << " integer input variable(s) and"
	      << "\n  " << nout << " output variables(s) and"
	      << "\n  Order(" << derOrder(0,0);
    for(int iy=1; iy<nout; ++iy)
      std::cerr << "," << derOrder(iy,0);
    std::cerr << ") derivatives" 
	      << "so there should be (" << num_multi_dim_poly_coef(nvarsr,derOrder(0,0));
    for(int iy=1; iy<nout; ++iy)
      std::cerr << "," << num_multi_dim_poly_coef(nvarsr,derOrder(iy,0));
    std::cerr << ") columns for the respective output variables." << std::endl
	      << "Found: " 
	      << "\n  " << nskip_read << " leading columns to skip and"
	      << "\n  " << nvarsr_read << " real input variable(s) and"
	      << "\n  " << nvarsi_read << " integer input variable(s) and"
	      << "\n  " << nout_read << " output variables(s) and" 
	      << "\n did not finish reading in the " << jder 
	      << "-th mixed partial derivative of total order " << ider
	      << " of output variable " << nout_read << std::endl;
    throw;
  } catch (...) {
    std::cerr << "Exception caught and rethrown in SurfData::readPointText(...)" << std::endl;
    throw;
  }
  return;
}
    

///read one point into SurfData from a binary .bspd file, skip_columns defaults to zero
void SurfData::readPointBinary(int ipt, istream& is, int skip_columns)
{
  std::cout << "SurfData: reading from a binary file has not yet been implemented NEEDS MUCH WORK must deal with cross platform endian-ness variation\n";
  assert(0);

  int nvarsr_read=0, nvarsi_read=0, nout_read=0, nskip_read=0, ider=0, jder=0;

  try {
    // read the point in binary format
    for (nvarsr_read=0; nvarsr_read<nvarsr; ++nvarsr_read) {
       // Throw an exception if there are fewer values on this line that
       // expected.
       nkm::surfpack::checkForEOF(is);
       is.read(reinterpret_cast<char*>(xr.ptr(nvarsr_read,ipt)),sizeof(xr(nvarsr_read,ipt)));
    }
    for (nout_read=0; nout_read<nout; ++nout_read) {
       // Throw an exception if there are fewer values on this line that
       // expected.
      ider=jder=0;
      nkm::surfpack::checkForEOF(is);
      is.read(reinterpret_cast<char*>(y.ptr(nout_read,ipt)),sizeof(y(nout_read,ipt)));
      for(ider=1; ider<=derOrder(nout_read,0); ++ider) {
	int nder_this_exact_total_order=derY[nout_read][ider].getNRows();
	for(jder=0; jder<nder_this_exact_total_order; ++jder) {
	  nkm::surfpack::checkForEOF(is);
	  is.read(reinterpret_cast<char*>(derY[nout_read][ider].ptr(jder,ipt)),sizeof(derY[nout_read][ider].ptr(jder,ipt)));
	}
      }
    }
  } catch (nkm::surfpack::io_exception&) {
    std::cerr << "Bad SurfPoint: binary file"  
	      << "\nExpected on this line: " 
	      << "\n  " << skip_columns << " leading columns to skip and"
	      << "\n  " << nvarsr << " real input variable(s) and"
	      << "\n  " << nvarsi << " integer input variable(s) and"
	      << "\n  " << nout << " output variables(s) and"
	      << "\n  Order(" << derOrder(0,0);
    for(int iy=1; iy<nout; ++iy)
      std::cerr << "," << derOrder(iy,0);
    std::cerr << ") derivatives" 
	      << "so there should be (" << num_multi_dim_poly_coef(nvarsr,derOrder(0,0));
    for(int iy=1; iy<nout; ++iy)
      std::cerr << "," << num_multi_dim_poly_coef(nvarsr,derOrder(iy,0));
    std::cerr << ") columns for the respective output variables." << std::endl
	      << "Found: " 
	      << "\n  " << nskip_read << " leading columns to skip and"
	      << "\n  " << nvarsr_read << " real input variable(s) and"
	      << "\n  " << nvarsi_read << " integer input variable(s) and"
	      << "\n  " << nout_read << " output variables(s) and" 
	      << "\n did not finish reading in the " << jder 
	      << "-th mixed partial derivative of total order " << ider
	      << " of output variable " << nout_read << std::endl;
    throw;
  } catch (...) {
    std::cerr << "Exception rethrown in SurfData::readPointBinary(...)" 
	      << std::endl;
    throw;
  }
}



/** Read a set of points into SurfData from an input stream of a text file 
    "*.spd". At this point nvarsr, nvarsi, nout, and optionally lockxr and 
    skip_columns must have already been set by the user. skip_columns 
    defaults to zero
*/
void SurfData::readText(istream& is, int skip_columns) 
{
  //std::cout << "readText skip_columns=" << skip_columns << std::endl;

  string single_line;
  int nlines=0;
  npts=0;
  try {
    //determine how many lines are in the file, so we can allocate matrices if any lines are blank or commented (such as a variable label line) then we will need to resize() matrices (not reshape() not newSize()) at the end
    while(!(is.eof())) {
      nlines++;
      getline(is,single_line);
    }
    //cout << "nlines=" << nlines << std::endl;
    assert(nlines&&nvarsr&&nout); //replace with a
    //if(nlines==0) throw;
    xr.newSize(nvarsr,nlines);
    xi.newSize(nvarsi,nlines);
    y.newSize(nout,nlines);
    derY.resize(nout);
    for(int iy=0; iy<nout; ++iy) {
      derY[iy].resize(derOrder(iy,0)+1);
      for(int ider=1; ider<=derOrder(iy,0); ++ider) 
	derY[iy][ider].newSize(num_multi_dim_poly_coef(nvarsr,-ider),nlines);
    }
    //cout << "just newsized arrays\n";

    std::cout << "TODO in SurfData.cpp: void SurfData::readText(istream&is, int skip_columns)  need to check for \"failbit\" and \"badbit\" before doing \"is.clear()\"\n";
    is.clear();
    is.seekg(0, ios::beg);
    
    //std::cout << "I think I just reset file position to beginning of file\n";
    getline(is,single_line);
    //std::cout << "firstline is=[" << single_line << "]\n";
    istringstream streamline(single_line);
    if (!readLabelsIfPresent(single_line, skip_columns)) {
      if ((single_line != "") && 
	  (single_line != "\n") && 
	  (single_line[0] != '%')&&
	  (single_line[0] != '#')) {
	readPointText(npts, single_line, skip_columns);
        npts = 1;
      }
      //std::cout << "I didn't find labels; nout=" << nout << std::endl;
    }

    //we don't know whether npts=0 or 1 will occur at this line during 
    //execution (Could be either depending on input file) so we have to 
    //program generically => separate ilines and npts
    for(int iline=1; iline<nlines; ++iline) {
      //cout << "iline=" << iline << "  npts=" << npts << std::endl;
      getline(is,single_line);
      if ((single_line != "") && 
	  (single_line != "\n") && 
	  (single_line[0] != '%')&&
	  (single_line[0] != '#')) {
	readPointText(npts, single_line, skip_columns);
	++npts; //only increase npts when we add a new point
      }
    }
  } catch(nkm::surfpack::io_exception& exception) {
    //std::cout << "encountered an exception in readText()\n";
    std::cerr << exception.what() << std::endl;
    throw;
  } 
  if(npts < nlines){
    //std::cout << "downsizing arrays because npts=" <<npts << " < nlines=" << nlines << std::endl; 
    //there were one or more blank or commented lines so we will remove the extra lines from the input/output variable matrices with the resize() command
    xr.resize(nvarsr,npts);
    xi.resize(nvarsi,npts);
    y.resize(nout,npts);
    for(int iy=0; iy<nout; ++iy)
      for(int ider=1; ider<=derOrder(iy,0); ++ider) 
	derY[iy][ider].resize(num_multi_dim_poly_coef(nvarsr,-ider),npts);
  }else if(npts > nlines) {
    assert(0);  //replace with a throw, the only way to get here should be that whatever called this function passed it an istream that was not at the beginning of the .spd file, in which case we've read values into memory locations beyond the bounds of the matrices, i.e. it's a segfault waiting to happen
  }
}

/// Write a set of SurfPoints to an output stream
void SurfData::writeText(ostream& os, bool write_labels) const
{
  assert((0<=nvarsr)&&(0<=nvarsi)&&(0<nvarsr+nvarsi)&&(1<=nout));
  stringstream ss; //output to file is slow because of overhead of each 
  //output call, so "output" each line to a stringstream, and when done 
  //with the line output the stringstream to the file so we have fewer 
  //write to file calls, also since I declared stringstream ss; in this 
  //function it disappears when the function exits so I can format ss 
  //instead of os so I don't have to worry about restoring the previous
  //format settings for os
    if (write_labels) {
      ss.setf(ios_base::left, ios_base::adjustfield);

      ss << '%';
      int correction = 1;      
      for(int ixr=0; ixr<nvarsr; ++ixr) {
        ss << setw(nkm::surfpack::field_width - correction) << xrLabels[ixr] 
	   << " ";
        correction=0;
      }
      for(int ixi=0; ixi<nvarsi; ++ixi) {
        ss << setw(nkm::surfpack::field_width - correction) << xiLabels[ixi]
	   << " ";
        correction=0;
      }
      for(int iy=0; iy<nout; ++iy) {
	ss << setw(nkm::surfpack::field_width) << yLabels[iy] << " ";
	if(derOrder(iy,0)>0) {
	  MtxInt der;
	  multi_dim_poly_power(der,nvarsr,derOrder(iy,0));
	  int nder_up_to_derOrder=der.getNCols();
	  for(int jder=1; jder<nder_up_to_derOrder; ++jder) {
	    //start at jder=1 rather than 0 because we already handled the function/response value 
	    int this_der_order=0;
	    for(int ixr=0; ixr<nvarsr; ++ixr)
	      this_der_order+=der(ixr,jder);
	    ostringstream der_label_os;
	    der_label_os << "d^" << this_der_order << "/dxr^(" << der(0,jder);
	    for(int ixr=1; ixr<nvarsr; ++ixr)
	      der_label_os << "," << der(ixr,jder);
	    der_label_os <<")";
	    ss << setw(nkm::surfpack::field_width) << der_label_os.str() << " ";
	  }
	    

	}
      }
      //ss << yLabels[nout-1]; //left one extra space at the end of the labels line
      os << ss.str() << std::endl;
    }

  // ios_base::flags returns ios::fmtflags object, but OSF compiler doesn't 
  // like that.
  // Save the stream flags.  The output precision may be modified, but it 
  // will be restored to its old value before the method exits. 
    ss.precision(nkm::surfpack::output_precision);
    ss.setf(ios::scientific);
    //ss.width(nkm::surfpack::field_width);
    //ss.setw(nkm::surfpack::field_width);
    for(int ipt=0; ipt<npts; ++ipt) {
      ss.str("");
      if(nvarsr>0) {
	ss << setw(nkm::surfpack::field_width) << xr(0,ipt);
	for(int ixr=1; ixr<nvarsr; ++ixr) 
	  ss << " " << setw(nkm::surfpack::field_width) << xr(ixr,ipt);
	for(int ixi=0; ixi<nvarsi; ++ixi)
	  ss << " " << setw(nkm::surfpack::field_width) << xi(ixi,ipt); 
      }
      else{ 
	assert(nvarsi>0);
	ss << setw(nkm::surfpack::field_width) << xi(0,ipt); 
	for(int ixi=1; ixi<nvarsi; ++ixi)
	  ss << " " << setw(nkm::surfpack::field_width) << xi(ixi,ipt); 
      }
      for(int iy=0; iy<nout; ++iy) {
	ss << " " << setw(nkm::surfpack::field_width)<< y(iy,ipt);
	for(int ider=1; ider<=derOrder(iy,0); ++ider) {
	  int nder_this_exact_total_order=derY[iy][ider].getNRows();
	  for(int jder=0; jder<nder_this_exact_total_order; ++jder)
	    ss << " " << setw(nkm::surfpack::field_width)<< derY[iy][ider](jder,ipt);
	}
      }
      os << ss.str() << std::endl;
    }
}

/** Read a set of points into SurfData from an input stream of a binary file 
    "*.bspd". At this point nvarsr, nvarsi, nout, and optionally lockxr and 
    skip_columns must have already been set by the user. skip_columns 
    defaults to zero, need to do something better for variable labels (there should be a way to read them from the binary file)
*/
void SurfData::readBinary(istream& is, int skip_columns) 
{
  int npts_read = 0;
  try {
    is.read((char*)&npts  ,sizeof(npts));
    is.read((char*)&nvarsr,sizeof(nvarsr));
    is.read((char*)&nvarsi,sizeof(nvarsi));
    is.read((char*)&nout  ,sizeof(nout));
 
    derOrder.newSize(nout,1);
    for(int iy=0; iy<nout; ++iy)
      is.read((char*)derOrder.ptr(iy,0),sizeof(derOrder(iy,0)));

    xr.newSize(nvarsr,npts);
    xi.newSize(nvarsi,npts);
    y.newSize(nout,npts);
    derY.resize(nout);
    for(int iy=0; iy<nout; ++iy) {
      derY[iy].resize(derOrder(iy,0)+1);
      for(int ider=1; ider<=derOrder(iy,0); ++ider) 
	derY[iy][ider].newSize(num_multi_dim_poly_coef(nvarsr,-ider),npts);
    }
    defaultLabels(); //need to replace this with something better (there should be a way to store labels in the binary file)

    for (npts_read=0; npts_read < npts; ++npts_read) {
      // Throw an exception if we hit the end-of-file before we've
      // read the number of points that were supposed to be there.
      nkm::surfpack::checkForEOF(is);
      // True for fourth argument signals a binary read
      readPointBinary(npts_read, is, skip_columns);
    }
    
  } catch(nkm::surfpack::io_exception&) {
    std::cerr << "Expected: " << npts << " points.  "
	      << "Read: " << npts_read << " points." << std::endl;
    throw;
  } 
}

  ///this must be a private function that not even SurfDataScaler calls, only other scaling functions in SurfData, it scales the derivatives of Y (of total order greater than or equal to 1) according to the values in unscalexr and unscaley. There is an inline SurfData member function: inline void unScaleDerY(){scaleDerY(-1);return;};
void SurfData::scaleDerY(int scalepower) {
  assert((scalepower==1)||(scalepower==-1));
  MtxInt der;
  for(int iy=0; iy<nout; ++iy) 
    for(int ider=1; ider<=derOrder(iy,0); ++ider) {
      multi_dim_poly_power(der,nvarsr,-ider);
      int nder_this_exact_total_order=der.getNCols();
      for(int jder=0; jder<nder_this_exact_total_order; ++jder) {
	double temp_double=scaleFactorDerY(der,iy,jder);
	if(scalepower==-1)
	  temp_double=1.0/temp_double;
	for(int ipt=0; ipt<npts; ++ipt)
	  derY[iy][ider](jder,ipt)*=temp_double;
      }
    }
  return;
}


///scale the data (output variable and real input variables) Note that this function should only be called by the Model if it is appropriate for the model to do so, it should not be called by anything other than a model. Also note that copies of a SurfData object will use the same default scaling of the original SurfData object, i.e if during construction you told the original SurfData object to use group scaling so will the copy
void SurfData::scaleToDefault()
{

  //WE WANT TO CHECK IF IT IS OK TO SCALE BEFORE SCALING
  int i, ifalreadyscaled=0;
  MtxDbl minmaxgroup(0,0);
  if(unscalexr.getNCols()!=0) {
    //if this is a non empty matrix
    assert((unscalexr.getNCols()==2)&&(unscalexr.getNRows()==nvarsr)&&
	   (unscaley.getNCols()==2)&&(unscaley.getNRows()==nout ));

    for(int ixr=0; ixr<nvarsr; ixr++)
      if((std::fabs(unscalexr(ixr,0))!=1.0)||(unscalexr(ixr,1)!=0.0)) {
	ifalreadyscaled=1;
	break;
      }

    if(!ifalreadyscaled)
      for(int iy=0; iy<nout; iy++)
	if((std::fabs(unscaley(iy,0))!=1.0)||(unscaley(iy,1)!=0.0)) {
	  ifalreadyscaled=1;
	  break;
	}
  }
  else{
    printf("Warning: You just tried to scale() an empty surfdata object; ignoring request to scale!!!\n");
    assert(0);
    return;
  }

  if(ifalreadyscaled) unScale();
  
  //scale each output dimensions/variables individually
  //printf("scaleToDefault: before indivScale: y(0)=%g y(1)=%g",y(0),y(1));
  indivScale(y,unscaley,minmaxgroup,false);
  //printf(" after indivScale y(0)=%g y(1)=%g\n",y(0),y(1));
  if(lockxr.getNElems()) {
    //at this point we know that the user wants us to scale by groups and that the dimensions (in lockxr) have been sorted by groups, note that some groups may have only one dimension in them so individually scale the only-one-dimension-groups and only group-scale the groups with more than one dimension
    int igroupstart=0, igroupstop=0, ngroup,k;
    MtxInt igroup;    
    MtxDbl group;
    MtxDbl unscalegroup;

    for(int ixr=1; ixr<nvarsr; ixr++) {
      if(lockxr(ixr,0)==lockxr(ixr-1,0))
	igroupstop=ixr;
      
      if((igroupstop==nvarsr-1)&&(igroupstart==0)) {
	groupScale(xr,unscalexr,minMaxXr,ifHaveMinMaxXr);
	break;
      }

      if(((lockxr(ixr,0)!=lockxr(ixr-1,0))||
	  (igroupstop==nvarsr-1))&&
	 (igroupstop!=igroupstart)) {
	//group scale if statement, we've found the end of a group is lock(ixr,0)!==lock(ixr-1,0) or we've reached the last dimension and the start of the group is not the same element as the end of the group (i.e. don't "group scale" a single dimension)
	ngroup=igroupstop-igroupstart+1;
	igroup.newSize(ngroup,1);
	for(k=0,i=igroupstart; i<=igroupstop; i++,k++)
	  igroup(k,0)=lockxr(i,1);	  
	xr.getRows(group,igroup);
	if(ifHaveMinMaxXr==true)
	  minMaxXr.getRows(minmaxgroup,igroup);
	groupScale(group,unscalegroup,minmaxgroup,ifHaveMinMaxXr);
	xr.putRows(group,igroup);
	unscalexr.putRows(unscalegroup,igroup);
	igroupstop=igroupstart=igroupstop+1;
      }

      if((igroupstop<ixr)||(igroupstop==nvarsr-1)) {
	//if igroup<ixr then we're at the end of a group and I didn't enter the "group scale" if statement so there's only a single dimension in this group so I want to individual scale; however, if the final dimension is the only one in it's group, and there was a "real" (more than 1 dimension) group that ended on the dimension before the final dimension, i.e. if igroupstop was nvarsr-2 and was just increased (i.e. igroupstop=igroupstart=igroupstop+1) then igroupstart=igroupstop=nvarsr-1 and in that case we want to individual scale the final (i.e. nvarsr-1) dimension even though the preceding (group scale) if statement was just executed.
	i=lockxr(igroupstart,1);
	xr.getRows(group,i);
	if(ifHaveMinMaxXr==true)
	  minMaxXr.getRows(minmaxgroup,i);
	indivScale(group,unscalegroup,minmaxgroup,ifHaveMinMaxXr);
	//indivScale(group,unscalegroup);
	xr.putRows(group,i);
	unscalexr.putRows(unscalegroup,i);
	igroupstop=igroupstart=igroupstop+1;
      }
	  
    } //yay!!! we're done scaling all the groups.
  }
  else{
    //printf("scaleToDefault::calling indivScale(xr,unscalexr)\n");
    //the user wanted us to individually scale (i.e not group scale) all the real input variables/dimensions
    indivScale(xr,unscalexr,minMaxXr,ifHaveMinMaxXr);
  }
  scaleDerY();
  return;
}

  ///individually scale each dimension to length 1 centered at zero, (a "singular" dimension, i.e. a dimension in which all points share the same value, receives unscalea(0)=-1.0 and unscalea(1) is set to the single constant value, a negative unscalea(0) is a flag that says to exclude this dimension from the analysis) at the time of this writing this function is intended to only be called by the scaleToDefault() function
  void SurfData::indivScale(MtxDbl& a, MtxDbl& unscalea, const MtxDbl& minmaxa, bool if_have_minmaxa)
{
  int nptsa=a.getNCols();
  int ndimsa=a.getNRows();
  unscalea.newSize(ndimsa,2);
  
  //printf("indivScale:: nptsa=%d ndimsa=%d unscalea.getNCols()=%d unscalea.getNRows()=%d\n",nptsa,ndimsa,unscalea.getNCols(),unscalea.getNRows()); fflush(stdout);
  

  int ipta, idima;
  double mina, maxa, scalea;

  for(idima=0; idima<ndimsa; idima++) {
    if(if_have_minmaxa==true) {
      mina=minmaxa(idima,0);
      maxa=minmaxa(idima,1);
    }
    else 
      mina=maxa=a(idima,0);
    
    for(ipta=0; ipta<nptsa; ipta++) 
      if(a(idima,ipta)<mina) mina=a(idima,ipta);
      else if(a(idima,ipta)>maxa) maxa=a(idima,ipta);
  
    
    unscalea(idima,0)=maxa-mina;
    unscalea(idima,1)=0.5*(maxa+mina);

    if(unscalea(idima,0)==0.0){
      //I refer to this as a "singular" dimension, all of points in this 
      //dimension are the same so exclude this dimension from anaylsis  
      unscalea(idima,0)=-1.0; //a negative length scale is a flag to exclude this
      //dimension from the analysis; but we retain unscalea(idima,1) to tell us 
      //what that constant value is
      for(ipta=0;ipta<nptsa;ipta++)
	a(idima,ipta)=0.0;
    }
    else{
      scalea=1.0/unscalea(idima,0);
      for(ipta=0; ipta<npts; ipta++)
	a(idima,ipta)=(a(idima,ipta)-unscalea(idima,1))*scalea;
    }
  }
  return;
};

  ///scale _this_ group of dimensions to a hyperrectangle of volume 1 centered at zero while preserving/locking the aspect ratio of the group of dimensions (note that this should only be called for __a__ __group__ (as in 1 group at a time) of real input variables, (a "singular" dimension idima, i.e. a dimension in which all points share the same value, is not counted as part of the group, instead it is mapped to 0.0, has unscalea(idima,0)=-1.0 and unscalea(idima,1) is set to the single constant value) at the time of this writing this function is intended to only be called by the scaleToDefault() function
void SurfData::groupScale(MtxDbl& a, MtxDbl& unscalea, const MtxDbl& minmaxa, bool if_have_minmaxa)
{  
  int nptsa=a.getNCols();
  int ndimsa=a.getNRows();
  unscalea.newSize(ndimsa,2);

  int ipta, idima;
  double mina, maxa, scalea, unscale_all;
  int numdiffa=0;
  double proddiffa=1.0;
  
  for(idima=0; idima<ndimsa; idima++) {
    if(if_have_minmaxa==true) {
      mina=minmaxa(idima,0);
      maxa=minmaxa(idima,1);
    }
    else 
      mina=maxa=a(idima,0);
    
    for(ipta=0; ipta<nptsa; ipta++) 
      if(a(idima,ipta)<mina) mina=a(idima,ipta);
      else if(a(idima,ipta)>maxa) maxa=a(idima,ipta);
    

    unscalea(idima,0)=maxa-mina;
    unscalea(idima,1)=0.5*(maxa+mina);

    if(unscalea(idima,0)==0.0){
      //I refer to this as a "singular" dimension, all of points in this 
      //dimension are the same so exclude this dimension from anaylsis  
      unscalea(idima,0)=-1.0; //a negative length scale is a flag to exclude this
      //dimension from the analysis; but we retain unscalea(idima,1) to tell us 
      //what that constant value is
      for(ipta=0;ipta<nptsa;ipta++)
	a(idima,ipta)=0.0;
    }
    else{
      numdiffa++;
      proddiffa*=unscalea(idima,0);
    }
  }

  unscale_all=std::pow(proddiffa,1.0/numdiffa);
  scalea=1.0/unscale_all;

  for(idima=0; idima<ndimsa; idima++) 
    if(unscalea(idima,0)!=-1.0) {
      //exclude dimensions in which all points are the same from analysis
      unscalea(idima,0)=unscale_all;
      for(ipta=0; ipta<nptsa; ipta++)
	a(idima,ipta)=(a(idima,ipta)-unscalea(idima,1))*scalea;
    }
  return;
}

///scale real input variables to domain_new(ixr,0)<=xr(ixr,:)<=domain_new(ixr,1) for ixr=0,1,...,nvarsr;  For "singular" dimension ixr xr(ixr,:)=0.5*(domain_new(ixr,0)+domain_new(ixr,1)),unscalexr(ixr,0)=-1.0 and unscalexr(ixr,1) = the single value - xr(ixr,:)
void SurfData::scaleXrToDomain(MtxDbl& domain_new){
  MtxDbl unscale_xr(nvarsr,2);
  for(int ixr=0; ixr<nvarsr; ++ixr) {
    unscale_xr(ixr,1)=0.5*(domain_new(ixr,0)+domain_new(ixr,1));
    unscale_xr(ixr,0)=domain_new(ixr,1)-unscale_xr(ixr,1);
  }
  scaleXrToFactor(unscale_xr);
  return;
}

void SurfData::scaleXrToFactor(MtxDbl& unscale_xr) {
  assert((unscale_xr.getNCols()==2)&&(unscale_xr.getNRows()==nvarsr));
  unScaleDerY();
  double scaleby, shiftby;
  for(int ixr=0; ixr<nvarsr; ixr++) {
    assert(unscale_xr(ixr,0));
    scaleby=std::fabs(unscalexr(ixr,0))/std::fabs(unscale_xr(ixr,0));
    shiftby=(unscalexr(ixr,1)-unscale_xr(ixr,1))/std::fabs(unscale_xr(ixr,0));
    unscalexr(ixr,0)=unscale_xr(ixr,0);
    unscalexr(ixr,1)=unscale_xr(ixr,1);
    for(int ipt=0; ipt<npts; ipt++)
      xr(ixr,ipt)=xr(ixr,ipt)*scaleby+shiftby;
  }
  scaleDerY();
  return;
}

void SurfData::scaleYToFactor(MtxDbl& unscale_y) {
  assert((unscale_y.getNCols()==2)&&(unscale_y.getNRows()==nout));
  double scaleby, shiftby;
  for(int iy=0; iy<nout; ++iy) {
    assert(unscale_y(iy,0));
    scaleby=std::fabs(unscaley(iy,0))/std::fabs(unscale_y(iy,0));
    shiftby=(unscaley(iy,1)-unscale_y(iy,1))/std::fabs(unscale_y(iy,0));
    unscaley(iy,0)=unscale_y(iy,0);
    unscaley(iy,1)=unscale_y(iy,1);
    for(int ipt=0; ipt<npts; ++ipt)
      y(iy,ipt)=y(iy,ipt)*scaleby+shiftby;
    for(int ider=1; ider<=derOrder(iy,0); ++ider) {
      int nder_this_exact_total_order=derY[iy][ider].getNRows();
      for(int jder=0; jder<nder_this_exact_total_order; ++jder)
	for(int ipt=0; ipt<npts; ++ipt)
	  derY[iy][ider](jder,ipt)*=scaleby;
    }
  }
  return;
}

///unscale this surfdata object 
SurfData& SurfData::unScale()
{
  unScaleDerY();

  //unscale the real input variables
  double scaleby, shiftby;
  for(int ixr=0; ixr<nvarsr; ixr++) {
    scaleby=std::fabs(unscalexr(ixr,0)); //in case this dimension was singular
    shiftby=unscalexr(ixr,1);
    unscalexr(ixr,0)=1.0;
    unscalexr(ixr,1)=0.0;
    for(int ipt=0; ipt<npts; ipt++)
      xr(ixr,ipt)=xr(ixr,ipt)*scaleby+shiftby;
  }
  
  //unscale the output variables  
  for(int iy=0; iy<nout; iy++) {
    scaleby=std::fabs(unscaley(iy,0)); //in case this dimension was singular
    shiftby=unscaley(iy,1);
    unscaley(iy,0)=1.0;
    unscaley(iy,1)=0.0;
    for(int ipt=0; ipt<npts; ipt++)
      y(iy,ipt)=y(iy,ipt)*scaleby+shiftby;
  }
  
  return *this;
}

MtxDbl& SurfData::scaleYOther(MtxDbl& y_other, int iy)
{
  if(iy==-99999) iy=iout;
#ifdef __SURFDATA_ERR_CHECK__
  assert((0<=iy)&&(iy<nout));
#endif
  int npts_other=y_other.getNCols();
  int nout_other=y_other.getNRows();
  double scaleby, shiftby;
  if(nout_other==1) {
    scaleby=1.0/std::fabs(unscaley(iy,0));
    shiftby=unscaley(iy,1);
    for(int ipt=0; ipt<npts_other; ++ipt)
      y_other(0,ipt)=(y_other(0,ipt)-shiftby)*scaleby;
  }
  else if(nout_other==nout) {
    for(iy=0; iy<nout; ++iy) {
      scaleby=1.0/std::fabs(unscaley(iy,0));
      shiftby=unscaley(iy,1);
      for(int ipt=0; ipt<npts_other; ++ipt)
	y_other(iy,ipt)=(y_other(iy,ipt)-shiftby)*scaleby;
    }
  }
  else {
    printf("MtxDbl& SurfData::scaleYOther(MtxDbl& y_other, int iy=iout)... nout=%d & nout_other=%d but must equal 1 or nout\n",nout,nout_other);
#ifdef __SURFDATA_ERR_CHECK__
    fflush(stdin);
    assert((nout_other==1)||(nout_other==nout));
#endif
  }
  return y_other;
}
  
MtxDbl& SurfData::unScaleYOther(MtxDbl& y_other, int iy)
{
  if(iy==-99999) iy=iout;
#ifdef __SURFDATA_ERR_CHECK__
  assert((0<=iy)&&(iy<nout));
#endif
  int npts_other=y_other.getNCols();
  int nout_other=y_other.getNRows();
  double scaleby, shiftby;
  if(nout_other==1) {
    scaleby=std::fabs(unscaley(iy,0));
    shiftby=unscaley(iy,1);
    for(int ipt=0; ipt<npts_other; ++ipt)
      y_other(0,ipt)=y_other(0,ipt)*scaleby+shiftby;
  }
  else if(nout_other==nout) {
    for(iy=0; iy<nout; ++iy) {
      scaleby=std::fabs(unscaley(iy,0));
      shiftby=unscaley(iy,1);
      for(int ipt=0; ipt<npts_other; ++ipt)
	y_other(iy,ipt)=y_other(iy,ipt)*scaleby+shiftby;
    }
  }
  else{
    printf("MtxDbl& SurfData::unScaleYOther(MtxDbl& y_other, int iy=iout)... nout=%d & nout_other=%d but must equal 1 or nout\n",nout,nout_other);      
#ifdef __SURFDATA_ERR_CHECK__
    fflush(stdin);
    assert((nout_other==1)||(nout_other==nout));
#endif 
  }
  return y_other;
}




///this function should only be called by putPoints() it decides whether to recommend that this SurfData object be rescaled after newpoints2 are added to the current point and returns 0 or 1 to indicate the recommendation; newpoints2 must be non-empty and already have the same scale as this SurfData object.  * 0 means we do not recommend rescaling: 0 is returned if this SurfData object is unscaled, or if the all new points fall within this SurfData object's range of real inputs and outputs. * 1 means that we recommend rescaling: we recommend rescaling if this SurfData object is scaled and at least one of the points in newpoints2 is outside this SurfData object's range of real inputs and outputs.
int SurfData::ifRecommendRescale(SurfData& newpoints2) {
  int nnew=newpoints2.npts;
  //we need to determine if we should recommend rescaling to the 
  //user/programmer
  int ifrecommendrescale; //think about switching this from an int to a bool
  //First thing to check: if the current data is unscaled we will 
  //recommend to NOT rescale the data after adding the new points
  
  int ifunscaled=1;
  for(int ixr=0; ixr<nvarsr; ixr++)
    if((unscalexr(ixr,0)!=1.0)||(unscalexr(ixr,1)!=0.0)){
      ifunscaled=0;
      break;
    }
  if(ifunscaled)
    for(int iy=0; iy<nout; iy++)
      if((unscaley(iy,0)!=1.0)||(unscaley(iy,1)!=0.0)){
	ifunscaled=0;
	break;
      }
  if(ifunscaled) {
    //the current data is unscaled so we recommend to NOT rescale the data
    ifrecommendrescale=0;
  }
  else{
    ifrecommendrescale=0; //by default we recommend not to rescale
    //we need to check if any of the new points are outside the range of the 
    //current points, if so we will recommend to rescale if not we will 
    //recommend to NOT rescale
    MtxDbl minmaxold(1,2), minmaxnew(1,2);
    //first we check xr
    for(int ixr=0; ixr<nvarsr; ixr++) {
      minmaxold(0,0)=minmaxold(0,1)=xr(ixr,0);
      for(int ipt=1; ipt<npts; ipt++) {
	if(xr(ixr,ipt)<minmaxold(0,0)) minmaxold(0,0)=xr(ixr,ipt);
	if(xr(ixr,ipt)>minmaxold(0,1)) minmaxold(0,1)=xr(ixr,ipt);	       
      }
      minmaxnew(0,0)=minmaxnew(0,1)=newpoints2.xr(ixr,0);
      for(int ipt=1; ipt<nnew; ipt++){
	if(newpoints2.xr(ixr,ipt)<minmaxnew(0,0)) 
	  minmaxnew(0,0)=newpoints2.xr(ixr,ipt);
	if(newpoints2.xr(ixr,ipt)>minmaxnew(0,1)) 
	  minmaxnew(0,1)=newpoints2.xr(ixr,ipt);	       
      }
      if((minmaxnew(0,0)<minmaxold(0,0))||(minmaxold(0,1)<minmaxnew(0,1))){
	ifrecommendrescale=1;
	break;
      }
    }
    if(!ifrecommendrescale) {
      //now we check y
      for(int iy=0; iy<nout; iy++) {
	minmaxold(0,0)=minmaxold(0,1)=y(iy,0);
	for(int ipt=1; ipt<npts; ipt++) {
	  if(y(iy,ipt)<minmaxold(0,0)) minmaxold(0,0)=y(iy,ipt);
	  if(y(iy,ipt)>minmaxold(0,1)) minmaxold(0,1)=y(iy,ipt);	       
	}
	minmaxnew(0,0)=minmaxnew(0,1)=newpoints2.y(iy,0);
	for(int ipt=1; ipt<nnew; ipt++){
	  if(newpoints2.y(iy,ipt)<minmaxnew(0,0)) 
	    minmaxnew(0,0)=newpoints2.y(iy,ipt);
	  if(newpoints2.y(iy,ipt)>minmaxnew(0,1)) 
	    minmaxnew(0,1)=newpoints2.y(iy,ipt);	       
	}
	if((minmaxnew(0,0)<minmaxold(0,0))||(minmaxold(0,1)<minmaxnew(0,1))){
	  ifrecommendrescale=1;
	  break;
	}
      } 
    }
  }
  return ifrecommendrescale;
}

///int putPoints(SurfData& newpoints, int ipt) puts the single point in newpoints into this SurfData objects ipt-th point, ipt must be greater than or equal to zero, if ipt is less than the current number of points the previous ipt-th point is replaced with the new one, if ipt is greater than the current number of points this function will enlarge the matrices to hold ipt+1 points before inserting this point, if ipt is unspecified it appends the newpoint to the end of this SurfData object's current list of points (after enlarging the matrices), if newpoints contains multiple points and ipt is unspecified it appends all the newpoints to the end of this SurfData obeject's current list of points (if ipt is specified newpoints must contain exactly one point, otherwise you should call int putPoints(SurfData& newpoints, MtxInt& ipts) instead).  If this SurfData object contained any points before this function was called then the newpoints are automatically rescaled to match the current points before they are inserted.  putPoints returns an integer to recommend whether or not the user/programmer should scale this SurfData object after calling this function: * 0 means we do not recommend rescaling: 0 is returned if newpoints was empty, this SurfData object was unscaled before this function was called, or if the all new points fell within the range of real inputs and outputs that this SurfData object had before this function was called. * -1 means that this SurfData object was empty before this function was called so we kept the scaling in newpoints, * 1 means that we recommend rescaling: we recommend rescaling if this SurfData object was scaled before this function was called and at least one of the newpoints real inputs or outputs was outside the range of range of the old  points.  This function does not actually do the rescaling of this SurfData object itself in case the user/programmer specified had scaled to a factor or scaled to a domain, or if for the sake of consistency/comparison-to-the-values-before-the-new-datapoints-were-added he/she doesn't want to rescale.
int SurfData::putPoints(SurfData& newpoints, int ipt) {
  int nnew=newpoints.npts;
  if(nnew==0) {
    std::cerr << "Warning!!! in: 'int SurfData::putPoints(SurfData& newpoints, int ipt)' newpoints is empty." << std::endl;
    return 0; //recommend to NOT rescale this SurfData object because we didn't change it
  }

  if(npts>0) {
    //we already have some points so we need to make sure the newpoints
    //have the same number of each type of variable as the current points
    assert((nvarsr==newpoints.nvarsr)&&
	   (nvarsi==newpoints.nvarsi)&&
	   (nout  ==newpoints.nout  ));
    for(int iy=0; iy<nout; ++iy)
      assert(derOrder(iy,0)==newpoints.derOrder(iy,0));
  }
  else{
    if(ipt==-99999){
      *this=newpoints;
      return -1; //we kept the scaling in newpoints because we didn't have 
      //any points in this SurfData object previously
    }
    assert((nnew==1)&&(0<=ipt));
    npts=ipt+1;
    nvarsr   =newpoints.nvarsr;
    nvarsi   =newpoints.nvarsi;
    nout     =newpoints.nout;
    iout     =newpoints.iout;
    unscalexr=newpoints.unscalexr;
    unscaley =newpoints.unscaley;
    lockxr   =newpoints.lockxr;
    xrLabels =newpoints.xrLabels;
    xiLabels =newpoints.xiLabels;
    yLabels  =newpoints.yLabels;
    xr.newSize(nvarsr,npts); xr.putCols(newpoints.xr,ipt);
    xi.newSize(nvarsi,npts); xi.putCols(newpoints.xi,ipt);
    y.newSize( nout  ,npts); y.putCols( newpoints.y ,ipt);
    derOrder =newpoints.derOrder;
    derY.resize(nout);
    for(int iy=0; iy<nout; ++iy) {
      derY[iy].resize(derOrder(iy,0)+1);
      for(int ider=1; ider<=derOrder(iy,0); ++ider) {
	int nder_this_exact_total_order=newpoints.derY[iy][ider].getNRows();
	derY[iy][ider].newSize(nder_this_exact_total_order,npts);
	derY[iy][ider].putCols(newpoints.derY[iy][ider],ipt);
      }
    }
    return -1;  //we kept the scaling in newpoints because we didn't have 
    //any points in this SurfData object previously
  }



  //convert the newpoints to the same scale as the current points
  SurfData newpoints2(newpoints); //don't want to change this in case the
  //user is still interested in using it for something else
  newpoints2.scaleToFactors(unscalexr,unscaley);
  
  //decide whether or not to recommend that this SurfData object after 
  //the newpoints have been added to the current points
  int ifrecommendrescale=ifRecommendRescale(newpoints2);

  if(ipt==-99999) {
    //the flag to append at the end so we need to enlarge our matrices
    //this is triggered in by NOT passing in a valued for ipt, -99999
    //is a default value
    ipt=npts; //if it's a single point, set ipt to put it at the end
    npts=npts+nnew;
    xr.resize(nvarsr,npts);
    xi.resize(nvarsi,npts);
    y.resize( nout  ,npts);
    //derY does not need to be resized because nout hasn't changed
    //dery[:] don't need to be resized because derOrder hasn't changed


    if(nnew==1){
      //there's just one point so don't do anything special
      xr.putCols(newpoints2.xr,ipt);
      xi.putCols(newpoints2.xi,ipt);
      y.putCols( newpoints2.y ,ipt);
      for(int iy=0; iy<nout; ++iy) 
	for(int ider=1; ider<=derOrder(iy,0); ++ider) {
	  int nder_this_exact_total_order=
	    newpoints2.derY[iy][ider].getNRows();
	  derY[iy][ider].resize(nder_this_exact_total_order,npts);
	  derY[iy][ider].putCols(newpoints2.derY[iy][ider],ipt);
	}
    }
    else{
      //there's multiple new points, so we need to create an integer 
      //vector of indices to specify the rows of the matrices we want
      //to fill in
      MtxInt ipts(nnew,1);  
      for(int inew=0; inew<nnew; inew++)
	ipts(inew,0)=ipt+inew;
      xr.putCols(newpoints2.xr,ipts);
      xi.putCols(newpoints2.xi,ipts);
      y.putCols( newpoints2.y ,ipts);

      for(int iy=0; iy<nout; ++iy) 
	for(int ider=1; ider<=derOrder(iy,0); ++ider) {
	  int nder_this_exact_total_order=
	    newpoints2.derY[iy][ider].getNRows();
	  derY[iy][ider].resize(nder_this_exact_total_order,npts);
	  derY[iy][ider].putCols(newpoints2.derY[iy][ider],ipts);
	}
    }
  }
  else{
    //the user is specifying a particular col to fill in
    
    //make sure there is only one point, otherwise he would have needed 
    //to pass in a MtxInt& ipts containing the indices, and that would have
    //called the other (overloaded) function with the same name, i.e. 
    //void putPoints(SurfData& newpoints, MtxInt& ipts)
    assert((nnew==1)&&(0<=ipt));
    if(npts<=ipt) {
      //the user is specifying a point beyond the edge of our matrices so 
      //we need to enlarge them before we fill in the new point.
      //if the specified point to "replace" is "greater than" (or equal to) 
      //the number of points (i.e. the size of the matrices) we currently 
      //have, enlarge the matrices
      npts=npts+nnew;
      xr.resize(nvarsr,npts);
      xi.resize(nvarsi,npts);
      y.resize( nout  ,npts);

      for(int iy=0; iy<nout; ++iy) 
	for(int ider=1; ider<=derOrder(iy,0); ++ider) {
	  int nder_this_exact_total_order=
	    newpoints2.derY[iy][ider].getNRows();
	  derY[iy][ider].resize(nder_this_exact_total_order,npts);
	  derY[iy][ider].putCols(newpoints2.derY[iy][ider],ipt);
	}
    }
    else{
      for(int iy=0; iy<nout; ++iy) 
	for(int ider=1; ider<=derOrder(iy,0); ++ider) 
	  derY[iy][ider].putCols(newpoints2.derY[iy][ider],ipt);	
    }
    //fill in the new point
    xr.putCols(newpoints2.xr,ipt);
    xi.putCols(newpoints2.xi,ipt);
    y.putCols( newpoints2.y ,ipt);
  }

  return ifrecommendrescale; //make the recommendation about whether or
  //not to rescale
}


///int putPoints(SurfData& newpoints, MtxInt& ipts) puts the points in newpoints into the points in this SurfData object listed in ipts, all values contained in ipts must be greater than or equal to zero, if any point listed in ipts is less than the current number of points the previous point is replaced with the new one, if any point listed in ipts is greater than the current number of points this function will enlarge the matrices to hold ipts.max()+1 points before inserting these points, if ipts is unspecified then int putPoints(SurfData& newpoints, int ipt) is called instead (with ipt defaulting to -99999 which is the flag to append the new points to the current list of points).  If this SurfData object contained any points before this function was called then the newpoints are automatically rescaled to match the current points before they are inserted. putPoints returns an integer to recommend whether or not the user/programmer should scale this SurfData object after calling this function: * 0 means we do not recommend rescaling: 0 is returned if newpoints was empty, this SurfData object was unscaled before this function was called, or if the all new points fell within the range of real inputs and outputs that this SurfData object had before this function was called. * -1 means that this SurfData object was empty before this function was called so we kept the scaling in newpoints, * 1 means that we recommend rescaling: we recommend rescaling if this SurfData object was scaled before this function was called and at least one of the newpoints real inputs or outputs was outside the range of range of the old  points.  This function does not actually do the rescaling of this SurfData object itself in case the user/programmer specified had scaled to a factor or scaled to a domain, or if for the sake of consistency/comparison-to-the-values-before-the-new-datapoints-were-added he/she doesn't want to rescale.
int SurfData::putPoints(SurfData& newpoints, MtxInt& ipts) {
  int nnew=newpoints.npts;
  if(nnew==0) {
    std::cerr << "Warning!!! in: 'int SurfData::putPoints(SurfData& newpoints, MtxInt& ipts)' newpoints is empty." << std::endl;
    return 0; //recommend to NOT rescale this SurfData object because we didn't change it
  }
  assert((0<=ipts.minElem())&&
	 (newpoints.npts==ipts.getNElems()));
  int ifrecommendrescale;
  SurfData *newpts;
  if(npts==0) {
    nvarsr   =newpoints.nvarsr;
    nvarsi   =newpoints.nvarsi;
    nout     =newpoints.nout;
    iout     =newpoints.iout;
    unscalexr=newpoints.unscalexr;
    unscaley =newpoints.unscaley;
    lockxr   =newpoints.lockxr;
    xrLabels =newpoints.xrLabels;
    xiLabels =newpoints.xiLabels;
    yLabels  =newpoints.yLabels;
    derOrder =newpoints.derOrder;
    ifrecommendrescale=-1;
    newpts=&newpoints;
    derY.resize(nout);
    for(int iy=0; iy<nout; ++iy)
      derY[iy].resize(derOrder(iy,0)+1);
  }
  else{
    //check to make sure that we have valid input
    assert((nvarsr==newpoints.nvarsr)&&
	   (nvarsi==newpoints.nvarsi)&&
	   (nout  ==newpoints.nout  ));
    for(int iy=0; iy<nout; ++iy)
      assert(derOrder(iy,0)==newpoints.derOrder(iy,0));

    SurfData newpoints2=newpoints; 
    newpts=&newpoints2;

    newpoints2.scaleToFactors(unscalexr,unscaley);

    //decide whether or not to recommend that this SurfData object after 
    //the newpoints have been added to the current points
    ifrecommendrescale=ifRecommendRescale(newpoints2);
  }

  int iptsmax=ipts.maxElem();
  if(npts<=iptsmax){
    npts=iptsmax+1;
    xr.resize(nvarsr,npts);
    xi.resize(nvarsi,npts);
    y.resize( nout  ,npts);
    //derY does not need to be resized because nout hasn't changed
    //derY[:] don't need to be resized because derOrder hasn't changed    
  }
  xr.putCols(newpts->xr,ipts);
  xi.putCols(newpts->xi,ipts);
  y.putCols( newpts->y ,ipts);
  for(int iy=0; iy<nout; ++iy) 
    for(int ider=1; ider<=derOrder(iy,0); ++ider) {
      int nder_this_exact_total_order=
	newpts->derY[iy][ider].getNRows();
      derY[iy][ider].resize(nder_this_exact_total_order,npts);
      derY[iy][ider].putCols(newpts->derY[iy][ider],ipts);
    }


  return ifrecommendrescale;
}

///retrieve one point (the one with index ipt), return it as a reference to SurfData
SurfData& SurfData::getPoints(SurfData& result, int ipt) {
  assert((0<=ipt)&&(ipt<npts));
  result.npts     =1;
  result.nvarsr   =nvarsr;
  result.nvarsi   =nvarsi;
  result.nout     =nout;
  result.iout     =iout;
  result.unscalexr=unscalexr;
  result.unscaley =unscaley;
  result.lockxr   =lockxr;
  result.xrLabels =xrLabels;
  result.xiLabels =xiLabels;
  result.yLabels  =yLabels;
  result.derOrder =derOrder;
  result.derY.resize(nout);
  for(int iy=0; iy<nout; ++iy) {
    result.derY[iy].resize(derOrder(iy,0)+1);
    for(int ider=1; ider<=derOrder(iy,0); ++ider) 
      derY[iy][ider].getCols(result.derY[iy][ider],ipt);
  }

  xr.getCols(result.xr,ipt);
  xi.getCols(result.xi,ipt);
  y.getCols( result.y ,ipt);

  return result;
}

///retrieve multiple points (the ones whose indices are stored in ipts), pass them back as a reference to a SurfData object
SurfData& SurfData::getPoints(SurfData& result, MtxInt& ipts) {
  int nget=ipts.getNRows();
  //printf("ipts=[%d %d %d %d]\n",ipts(0,0),ipts(1,0),ipts(2,0),ipts(3,0));
  assert(ipts.getNCols()==1);
  ipts.uniqueElems();
  //printf("nget=%d ipts.getNRows=%d\n",nget,ipts.getNRows());
  //printf("ipts=[%d %d %d]\n",ipts(0,0),ipts(1,0),ipts(2,0));
  //fflush(stdout);
  assert(nget==ipts.getNRows());
  //assert(false);
  assert((0<=ipts(0,0))&&(ipts(nget-1,0)<npts));
  result.npts     =nget;
  result.nvarsr   =nvarsr;
  result.nvarsi   =nvarsi;
  result.nout     =nout;
  result.iout     =iout;
  result.unscalexr=unscalexr;
  result.unscaley =unscaley;
  result.lockxr   =lockxr;
  result.xrLabels =xrLabels;
  result.xiLabels =xiLabels;
  result.yLabels  =yLabels;
  result.derOrder =derOrder;
  result.derY.resize(nout);
  for(int iy=0; iy<nout; ++iy) {
    result.derY[iy].resize(derOrder(iy,0)+1);
    for(int ider=1; ider<=derOrder(iy,0); ++ider) 
      derY[iy][ider].getCols(result.derY[iy][ider],ipts);
  }

  xr.getCols(result.xr,ipts);
  xi.getCols(result.xi,ipts);
  y.getCols( result.y ,ipts);

  return result;
}

///retrieve all points except one (the one with index ipt), return them as a reference to a SurfData object
SurfData& SurfData::excludePoints(SurfData& result, int ipt) {
  assert((0<=ipt)&&(ipt<npts));
  result.npts     =npts-1;
  result.nvarsr   =nvarsr;
  result.nvarsi   =nvarsi;
  result.nout     =nout;
  result.iout     =iout;
  result.unscalexr=unscalexr;
  result.unscaley =unscaley;
  result.lockxr   =lockxr;
  result.xrLabels =xrLabels;
  result.xiLabels =xiLabels;
  result.yLabels  =yLabels;
  result.derOrder =derOrder;
  result.derY.resize(nout);
  for(int iy=0; iy<nout; ++iy) {
    result.derY[iy].resize(derOrder(iy,0)+1);
    for(int ider=1; ider<=derOrder(iy,0); ++ider)
      derY[iy][ider].excludeCols(result.derY[iy][ider],ipt);
  }

  xr.excludeCols(result.xr,ipt);
  xi.excludeCols(result.xi,ipt);
  y.excludeCols( result.y ,ipt);

  return result; 
}


///retrieve all points except the ones whose indices are listed in ipts, return them as a reference to a SurfData object
SurfData& SurfData::excludePoints(SurfData& result, MtxInt& ipts) {
  int nexclude=ipts.getNRows();
  assert(ipts.getNCols()==1);
  ipts.uniqueElems();
  assert(nexclude==ipts.getNRows());
  assert((0<=ipts(0,0))&&(ipts(nexclude-1,0)<npts));
  result.npts     =npts-nexclude;
  result.nvarsr   =nvarsr;
  result.nvarsi   =nvarsi;
  result.nout     =nout;
  result.iout     =iout;
  result.unscalexr=unscalexr;
  result.unscaley =unscaley;
  result.lockxr   =lockxr;
  result.xrLabels =xrLabels;
  result.xiLabels =xiLabels;
  result.yLabels  =yLabels;
  result.derOrder =derOrder;
  result.derY.resize(nout);
  for(int iy=0; iy<nout; ++iy) {
    result.derY[iy].resize(derOrder(iy,0)+1);
    for(int ider=1; ider<=derOrder(iy,0); ++ider)
      derY[iy][ider].excludeCols(result.derY[iy][ider],ipts);
  }

  xr.excludeCols(result.xr,ipts);
  xi.excludeCols(result.xi,ipts);
  y.excludeCols( result.y ,ipts);

  return result;
}

///create a split copy of the points stored in this surfdata, one point, the one whose index is listed in ipt, will be placed in "extracted", the rest will be placed in "rest", this is very useful for the PRESS metric
void SurfData::extractPoints(SurfData& rest, SurfData& extracted, int ipt) {
  assert((0<=ipt)&&(ipt<npts));
  /*
  assert((npts  ==xr.getNCols())&&
	 (npts  ==y.getNCols())&&
	 (nvarsr==xr.getNRows())&&
	 (nout  ==y.getNRows()));  
  printf("extractPoints(): npts=%d nvarsr=%d nout=%d\n",npts,nvarsr,nout); fflush(stdout);
  */

  extracted.npts     =1;
  extracted.nvarsr   =nvarsr;
  extracted.nvarsi   =nvarsi;
  extracted.nout     =nout;
  extracted.iout     =iout;
  extracted.unscalexr=unscalexr;
  //for(int ixr=0; ixr<nvarsr; ++ixr)
  //printf("unscalexr(:,%d)=[%g;%g] ",ixr,unscalexr(ixr,0), unscalexr(ixr,1));
  //printf("\n");
  extracted.unscaley =unscaley;
  //for(int iy=0; iy<nout; ++iy)
  //printf("unscaley(:,%d)=[%g;%g] ",iy,unscaley(iy,0), unscaley(iy,1));
  //printf("\n");
  extracted.lockxr   =lockxr;
  extracted.xrLabels =xrLabels;
  extracted.xiLabels =xiLabels;
  extracted.yLabels  =yLabels;
  extracted.derOrder =derOrder;
  extracted.derY.resize(nout);
  for(int iy=0; iy<nout; ++iy) {
    extracted.derY[iy].resize(derOrder(iy,0)+1);
    for(int ider=1; ider<=derOrder(iy,0); ++ider)
      derY[iy][ider].getCols(extracted.derY[iy][ider],ipt);
  }

  /*
  assert((npts ==xr.getNCols())&&
	 (npts ==y.getNCols())&&
	 (nvarsr==xr.getNRows())&&
	 (nout ==y.getNRows()));  
  printf("extractPoints(): npts=%d nvarsr=%d nout=%d\n",npts,nvarsr,nout); fflush(stdout);
  */

  xr.getCols(extracted.xr,ipt);
  if(nvarsi) xi.getCols(extracted.xi,ipt);
  y.getCols(extracted.y,ipt);

  rest.npts     =npts-1;
  rest.nvarsr   =nvarsr;
  rest.nvarsi   =nvarsi;
  rest.nout     =nout;
  rest.iout     =iout;
  rest.unscalexr=unscalexr;
  rest.unscaley =unscaley;
  rest.lockxr   =lockxr;
  rest.xrLabels =xrLabels;
  rest.xiLabels =xiLabels;
  rest.yLabels  =yLabels;
  rest.derOrder =derOrder;
  rest.derY.resize(nout);
  for(int iy=0; iy<nout; ++iy) {
    rest.derY[iy].resize(derOrder(iy,0)+1);
    for(int ider=1; ider<=derOrder(iy,0); ++ider)
      derY[iy][ider].excludeCols(rest.derY[iy][ider],ipt);
  }

  xr.excludeCols(rest.xr,ipt);
  if(nvarsi) xi.excludeCols(rest.xi,ipt);
  y.excludeCols(rest.y,ipt);

  return;
}

///create a split copy of the points stored in this surfdata, the ones whose indices are listed in ipts will be placed in "extracted", the rest will be placed in "rest", this is very useful for the Leave N Out cross Validation metric
void SurfData::extractPoints(SurfData& rest, SurfData& extracted, MtxInt& ipts) {
  int nextract=ipts.getNRows();
  assert(ipts.getNCols()==1);
  ipts.uniqueElems();
  assert(nextract==ipts.getNRows());
  assert((0<=ipts(0,0))&&(ipts(nextract-1,0)<npts));
  extracted.npts     =nextract;
  extracted.nvarsr   =nvarsr;
  extracted.nvarsi   =nvarsi;
  extracted.nout     =nout;
  extracted.iout     =iout;
  extracted.unscalexr=unscalexr;
  extracted.unscaley =unscaley;
  extracted.lockxr   =lockxr;
  extracted.xrLabels =xrLabels;
  extracted.xiLabels =xiLabels;
  extracted.yLabels  =yLabels;
  extracted.derOrder =derOrder;
  extracted.derY.resize(nout);
  for(int iy=0; iy<nout; ++iy) {
    extracted.derY[iy].resize(derOrder(iy,0)+1);
    for(int ider=1; ider<=derOrder(iy,0); ++ider)
      derY[iy][ider].getCols(extracted.derY[iy][ider],ipts);
  }

  xr.getCols(extracted.xr,ipts);
  if(nvarsi) xi.getCols(extracted.xi,ipts);
  y.getCols( extracted.y ,ipts);

  rest.npts     =npts-nextract;
  rest.nvarsr   =nvarsr;
  rest.nvarsi   =nvarsi;
  rest.nout     =nout;
  rest.iout     =iout;
  rest.unscalexr=unscalexr;
  rest.unscaley =unscaley;
  rest.lockxr   =lockxr;
  rest.xrLabels =xrLabels;
  rest.xiLabels =xiLabels;
  rest.yLabels  =yLabels;
  rest.derOrder =derOrder;
  rest.derY.resize(nout);
  for(int iy=0; iy<nout; ++iy) {
    rest.derY[iy].resize(derOrder(iy,0)+1);
    for(int ider=1; ider<=derOrder(iy,0); ++ider)
      derY[iy][ider].excludeCols(rest.derY[iy][ider],ipts);
  }

  xr.excludeCols(rest.xr,ipts);
  if(nvarsi) xi.excludeCols(rest.xi,ipts);
  y.excludeCols( rest.y ,ipts);

  return;
}


/***********************************************************************/
/***********************************************************************/
/**** SurfData member functions end here                            ****/
/***********************************************************************/
/***********************************************************************/

} // end namespace nkm


#ifdef __SURFDATA_SCALING_UNIT_TEST__
//needs surf77_config.h  surfpack_LAPACK_wrappers.h  surfpack_system_headers.h everything else is in the nkm directory
//g++ -g -o NKM_SurfData_UnitMain NKM_SurfData.cpp NKM_SurfPack.cpp NKM_SurfMat.cpp /usr/lib64/libblas.so /usr/lib64/liblapack.so

using namespace nkm;
int main(){
  SurfData_scaling_unit_test();
  return 0;
};
#endif

