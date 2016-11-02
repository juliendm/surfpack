#include "NKM_SurfData.hpp"
#include "NKM_KrigingModel.hpp"
#include "NKM_LinearRegressionModel.hpp"

#include <iostream>
#include <map>
#include <string>

using std::cout;
using std::endl;
using std::string;

void get_input_file_params(string& filename, int& npts, int& nvarsr, 
			   int& nvarsi, int& nout, int& jout, 
			   int& skip_columns, int test_number);
//void read_surf_data(nkm::SurfData& sd, int test_number);
void Press(nkm::SurfData& sd);
void validate();

// will change this to a "unit testing" main function later
int main(int argc, char* argv[])
{
  int test_number = (argc>1) ? atoi(argv[1]) : 1;
  //  validate();
  //  return 0;
  printf("Running test number %d.\n", test_number);

  string filename;
  int npts, nvarsr, nvarsi, nout, jout, skip_columns;
  get_input_file_params(filename, npts, nvarsr, nvarsi, nout, jout, 
			skip_columns, test_number);
  //filename="testout.spd"; //only used this to make sure input could read output

  //cout << "filename=[" << filename << "]\n";
  nkm::SurfData sd(filename, nvarsr, nvarsi, nout, jout, skip_columns);
  //nkm::SurfData sd;

  //read_surf_data(sd, test_number);
    
  //nkm::KrigingModel km(sd);

  Press(sd);
  return 0;
}

// this reads XR and Y, and scales only XR so that the result of Press will have the same units as Surfpack
//void read_surf_data(nkm::SurfData& sd, int test_number)
void get_input_file_params(string& filename, int& npts, int& nvarsr, 
			   int& nvarsi, int& nout, int& jout, 
			   int& skip_columns, int test_number)
{
  nvarsi=0;
  nout=1;
  jout=0;
  skip_columns=0;
  
  //int i, j;

  //printf("**********************************************************\n");
  //printf("**********************************************************\n");


  //string filename;

  //FILE* fp; 
  //int npts, nvarsr;
  //char yada[2048];
  switch(test_number) {
  case 1:
    filename="sat_scores.spd";
    //fp=fopen("sat_scores.spd","r");
    npts=50;
    nvarsr=6;
    //fgets(yada,2048,fp);
    printf("sat_scores.spd: Direct Performance with rcond condition number bounding\n  *Kriging Press with const trend function should be about press_Score=689.052;\n  *Kriging Press with linear trend should be about press_score=831.985;\n  *Kriging Press with (unroated) quadratic trend function should be about press_score=???5993.9???... possible explanations: A) the correlation matrix for this data set and polynomial degree has much larger than typical condition numbers (known to be true) and that is BECAUSE the system is very poorly behaved B)  there are 50 pieces of data, 28 polynomial coefficients and 6 correlation parameters, and the variance the general rule of thumb is you want twice as many pieces of data as variables you need to find the value of, probably both are true\n If it's higher than the above values, then something is wrong.\n");
    //\n  *Kriging Press with ROTATED quadratic trend function without condition number bounding should be about press_score=???3272.6???\n  *Kriging Press with (unrotated) quadratic trend function without condition number bounding should be about press_score=???3171.34???;
    break;
  case 2:
    filename="matlab_peaks.spd";
    //fp=fopen("matlab_peaks.spd","r");
    npts =20;
    nvarsr=2;
    //fgets(yada,2048,fp);
    printf("matlab_peaks.spd: After being unscaled, Kriging Press should be about press_score=5.89358; if it's higher, then something is wrong.\n");
    break;
  case 3: //don't use this one... there is far too much data, produces infinite negative log likelihoods, may be an issue of a singular correlation matrix but if so it didn't get far enough to tell
    filename="n26e171.spd";
    //fp=fopen("n26e171.spd","r");
    npts=1849;
    nvarsr=2;
    printf("n26e171.spd: has never worked... too big\n");
    break;
  case 4:
    filename="Grid_81_March2010_build_AREA.spd";
    //fp=fopen("Grid_81_March2010_build_AREA.spd","r");
    npts=81;
    nvarsr=2;
    //fgets(yada,2048,fp);
    printf("GoodYearFarmLugArea9x9grid: After being unscaled\n  *with self choosing Jitter (and correlation parameter optimization) to bound the condition number, max ChooseNug=0.1, and 2nd order polynomial trend functions Kriging Press should be about press_score=4.89787e-05\n  *with same except for linear trend functions Kriging Press should be about press_score=0.00127138\n  *with same except for a constant trend function Kriging Press should be about press_score=0.0742172\nif it's higher, then something is wrong.\n");
    break;
  case 5:
    filename="Grid_81_March2010_build_SED.spd";
    //fp=fopen("Grid_81_March2010_build_SED.spd","r");
    npts=81;
    nvarsr=2;
    //fgets(yada,2048,fp);
    printf("GoodYearFarmLugSED9x9grid: After being unscaled\n  *with self choosing Jitter (and correlation parameter optimization) to bound the condition number, max ChooseNug=0.1, and 2nd order polynomial trend functions Kriging Press should be about press_score=0.0103378\n  *with same except for linear trend functions Kriging Press should be about press_score=0.0336298\n  *with same except for a constant trend function Kriging Press should be about press_score=0.173744\nif it's higher, then something is wrong.\n");
    break;
  case 6:
    filename="MontserratFlowEastNorthH.spd";
    //fp=fopen("MontserratFlowEastNorthH.spd","r");
    npts=1172;
    nvarsr=2;
    //fgets(yada,2048,fp);
    printf("MontserratFlowEastNorthH.spd: has always resulted in a singular correlation matrix\n");
    break;
  case 7:
    filename="Paviani10D_200pts.spd";
    //fp=fopen("Paviani10D_200pts.spd","r");
    npts=200;
    nvarsr=10;
    //fgets(yada,2048,fp);
    printf("Paviani10D_200pts.spd: for direct optimizer with default settings and rcond condition number bounding, and non-rotated coordinates for the trend function, After being unscaled, \n *Kriging Press for constant trend should be about press_score=156.445\n *Kriging Press for linear trend should be about press_score= 153.575, and\n *Kriging Press for quadratic trend should be about press_score = 158.491.\n If it's higher then something is wrong.\n");
    break;
  default:
    assert(0);
  }

  

  return;

#ifdef nononono
  nkm::MtxDbl XR(npts,nvarsr),Y(npts);
  
  for(i=0; i<npts; i++) {
    for(j=0; j<nvarsr; j++)
      fscanf(fp,"%lf",XR.ptr(i,j));
    fscanf(fp,"%lf",Y.ptr(i));
    
    /*
      printf("XR(%d,:)={",i);
      for(j=0; j<nvarsr; j++) printf(" %g",XR(i,j));
      printf(" }; Y(%d)=%g\n",i,Y(i));
    */
  }
  fclose(fp);
  
  sd=nkm::SurfData(XR,Y);
  
  /*
    printf("minmax(XR)={%g,%g}; minmax(xr)={%g, %g}\n",XR.minElem(),XR.maxElem(),
    sd.xr.minElem(),sd.xr.maxElem());
    printf("minmax(Y)={%g,%g}; minmax(y)={%g, %g}\n",Y.minElem(),Y.maxElem(),
    sd.y.minElem(),sd.y.maxElem());
    for(j=0; j<nvarsr; j++)
    printf("unscalexr(:,%d)={%g, %g}\n",j,sd.unscalexr(0,j),sd.unscalexr(1,j));
    for(j=0; j<1; j++)
    printf("unscaley(:,%d)={%g, %g}\n",j,sd.unscaley(0,j),sd.unscaley(1,j));
  */
  
  return;
#endif
}

void Press(nkm::SurfData& sd)
{
  /*
    assert((sd.npts ==sd.xr.getNRows())&&
    (sd.npts ==sd.y.getNRows())&&
    (sd.nvarsr==sd.xr.getNCols())&&
    (sd.nout ==sd.y.getNCols()));  
  */
  int npts =sd.getNPts();
  //int nvarsr=sd.getNVarsr();
  //int nout =sd.getNOut();
  //printf("Press: npts=%d nvarsr=%d nout=%d\n",npts,nvarsr,nout); fflush(stdout);
  double temp_double;

  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  km_params["order"] = "linear";

  nkm::MtxDbl yeval(npts,1), y(npts,1);
  nkm::KrigingModel km(sd, km_params);
  nkm::SurfData sdeval(sd);
  km.create();
  nkm::MtxDbl d1y, d2y;
  //km.evaluate_d1y(d1y,sd.xr);
  //km.evaluate_d2y(d2y,sd.xr);
  km.evaluate(yeval,sd.xr);
  sdeval.y.putCols(yeval,sdeval.getJOut());
  string filename="testout.spd";
  sdeval.write(filename);
  sd.getY(y);
  double kmrms(0.0);
  
  for(int i=0; i<npts; i++){
    temp_double=y(i,0)-yeval(i,0);
    kmrms+=temp_double*temp_double;
  }
  kmrms=std::sqrt(kmrms/npts);
  
  nkm::LinearRegressionModel lrm(sd);
  
  double press_score_km =0.0;
  double press_score_lrm=0.0;

  
  nkm::SurfData rest, extracted;
  
  
  for(int ipt=0; ipt<npts; ipt++) {
  //for(int ipt=0; ipt<1; ipt++) {
    //printf("press ipt=%d/%d:\n",ipt,npts);
    /*
      assert((sd.npts ==sd.xr.getNRows())&&
      (sd.npts ==sd.y.getNRows())&&
      (sd.nvarsr==sd.xr.getNCols())&&
	   (sd.nout ==sd.y.getNCols()));  
	   printf("Press: S: ipt=%d npts=%d nvarsr=%d nout=%d\n",ipt,npts,nvarsr,nout); fflush(stdout);
    */
    
    sd.extractPoints(rest,extracted,ipt);
    //for(int jpt=0; jpt<npts-1; ++jpt) {
    //printf("yrest(%d,:)={%g",jpt,rest.y(jpt,0));
    //for(int jout=1; jout<nout; ++jout)
    //printf(", %g",rest.y(jpt,jout));
    //printf("}\n");
    //}

    //printf("ipt=%d: rest.minmax(xr)={%g, %g} rest.minmax(y)={%g, %g}, ex.minmax(xr)={%g, %g} ex.minmax(y)={%g, %g}\n",ipt,rest.xr.minElem(),rest.xr.maxElem(),rest.y.minElem(),rest.y.maxElem(),extracted.xr.minElem(),extracted.xr.maxElem(),extracted.y.minElem(),extracted.y.maxElem());
    /*
    //check that extraction is being done right
    for(int j=0; j<nvarsr; j++){
      int ir, is;
      for(is=0; is<ipt; is++)
	if(rest.xr(is,j)!=sd.xr(is,j)) {
	  printf("error: xr j=%d is=%d ir=%d\n",j,is,is); fflush(stdout);
	  assert(rest.xr(is,j)==sd.xr(is,j));}

      if(extracted.xr(0,j)!=sd.xr(is,j)) {
	  printf("error: xr j=%d is=%d ie=%d\n",j,is,0); fflush(stdout);
	  assert(extracted.xr(0,j)==sd.xr(is,j));}
      ir=is; is=ipt+1;
      for(; is<npts; is++, ir++)
	if(rest.xr(ir,j)!=sd.xr(is,j)) {
	  printf("error: xr j=%d is=%d ir=%d\n",j,is,ir); fflush(stdout);
	  assert(rest.xr(ir,j)==sd.xr(is,j));}
    }

    for(int j=0; j<nout; j++){
      int ir, is;
      for(is=0; is<ipt; is++)
	if(rest.y(is,j)!=sd.y(is,j)) {
	  printf("error: y j=%d is=%d ir=%d\n",j,is,is); fflush(stdout);
	  assert(rest.y(is,j)==sd.y(is,j));}

      if(extracted.y(0,j)!=sd.y(is,j)) {
	  printf("error: y j=%d is=%d ie=%d\n",j,is,0); fflush(stdout);
	  assert(extracted.y(0,j)==sd.y(is,j));}

      ir=is; is=ipt+1;
      for(; is<npts; is++, ir++)
	if(rest.y(ir,j)!=sd.y(is,j)) {
	  printf("error: y j=%d is=%d ir=%d\n",j,is,ir); fflush(stdout);
	  assert(rest.y(ir,j)==sd.y(is,j));}
    }
    */

    nkm::KrigingModel km(rest, km_params);
    km.create();

    temp_double=km.evaluate(extracted.xr)-extracted.y(0,extracted.getJOut());
    
    press_score_km+=temp_double*temp_double;
  }
  
   
  for(int ipt=0; ipt<npts; ipt++) {

    sd.extractPoints(rest,extracted,ipt);

    nkm::LinearRegressionModel lrm(rest);
    //printf("\nNpoly=%d rms=%g\n",lrm.getNPoly(),lrm.getRMS());

    temp_double=lrm.evaluate(extracted.xr)-extracted.y(0,extracted.getJOut());

    press_score_lrm+=temp_double*temp_double;
  }
  

  //I think Press is mean(squared(error)) (not RMS error)
  press_score_km =press_score_km/npts;
  press_score_lrm=press_score_lrm/npts;
  

  //kmrms=0.0/0.0;
  //printf("rms_km = %g\npress_score_km =%g\n",
  // kmrms,press_score_km *temp_double*temp_double);
  printf("rms_km = %g\npress_score_km =%g\npress_score_lrm=%g\n",
	 kmrms,press_score_km,press_score_lrm);
  
  return;
}
