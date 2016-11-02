#include "NKM_SurfMat.hpp"
#include "NKM_SurfData.hpp"
#include "NKM_KrigingModel.hpp"
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <cstdlib>


#define __DONT_COMPILE__
using std::cout;
using std::endl;
using std::string;
using namespace std;
using std::ostringstream;



int valgrind_test(int itest);
void timing_tests();
void nightly_test();



void hack();
void check_matrix();
void compare_sample_designs();
void compare_sample_designs_pav(int nvarsr);
void nested_Krig_vs_GEK_herbie_smooth_herbie_2D_4D_8D();
void GPAIS_build_and_eval_mean_adjvar(int Nvarsr, 
				      std::string& model_type,
				      std::string& optimization_method,
				      std::string& corr_func_family,
				      std::string& corr_func_param);
void generic_build_and_eval_mean_adjvar(int Nvarsr, 
					std::string& model_type,
					std::string& optimization_method,
					std::string& corr_func_family,
					std::string& corr_func_param,
					std::string& handle_ill_cond1,
					std::string& handle_ill_cond2);



void gen_sample_design_by_pivoted_cholesky() {
  int NumGuesses=100;
  int Npts=2048;
  int NptsGuess=4096;
  int Ndim=2;
  std::string filename="unit_hypercube_nested_design_2D_2048pts.txt";
  int imod=104395303; //a large prime number
  double dmod=static_cast<double>(imod);
  double L=0.155*std::pow(1.0/static_cast<double>(Npts),1.0/static_cast<double>(Ndim));
  double negtheta=-0.5/(L*L);
  double dtemp;

  nkm::MtxDbl x(NptsGuess,Ndim);
  nkm::MtxDbl xfinal(Npts,Ndim);
  nkm::MtxDbl R(NptsGuess,NptsGuess);
  for(int idim=0; idim<Ndim; ++idim) {
    x(0,idim)=0.5;
    for(int i=1; i<NptsGuess; ++i)
      x(i,idim)=static_cast<double>(std::rand()%imod)/dmod;
  }
  for(int j=0; j<NptsGuess-1; ++j)
    for(int i=j+1; i<NptsGuess; ++i) {
      dtemp=x(i,0)-x(j,0);
      R(i,j)=dtemp*dtemp;
    }
  for(int idim=1; idim<Ndim-1; ++idim)
    for(int j=0; j<NptsGuess-1; ++j)
      for(int i=j+1; i<NptsGuess; ++i) {
	dtemp=x(i,idim)-x(j,idim);
	R(i,j)+=dtemp*dtemp;
      }      
  for(int j=0; j<NptsGuess; ++j) {
    R(j,j)=1.0;
    for(int i=j+1; i<NptsGuess; ++i) {
      dtemp=x(i,Ndim)-x(j,Ndim);
      R(j,i)=R(i,j)=std::exp(negtheta*(R(i,j)+dtemp*dtemp));
    }
  }
  int info=0;
  char uplo='B';
  nkm::MtxInt ipiv(NptsGuess,1);  
  int ld_R=R.getNRowsAct();
  double min_allowed_rcond=std::pow(2.0,-40.0);
  int rank=-Npts;

  PIVOTCHOL_F77(&uplo, &NptsGuess, R.ptr(0,0), &ld_R,
    		ipiv.ptr(0,0), &rank, &min_allowed_rcond, &info); 
  printf("iLoop=0 rank=%d/%d\n",rank,NptsGuess);
  if(rank>Npts)
    rank=Npts;
  for(int idim=0; idim<Ndim; ++idim)
    for(int i=0; i<rank; ++i)
      xfinal(i,idim)=x(ipiv(i,0)-1,idim);

  for(int iLoop=1; iLoop<NumGuesses; ++iLoop) {

    for(int idim=0; idim<Ndim; ++idim) {
      //don't do i=0 because it doesn't change
      for(int i=1; i<rank; ++i)
	x(i,idim)=xfinal(i,idim);
      for(int i=rank; i<NptsGuess; ++i)
	x(i,idim)=static_cast<double>(std::rand()%imod)/dmod;
    }
    
    for(int j=0; j<NptsGuess-1; ++j)
      for(int i=j+1; i<NptsGuess; ++i) {
	dtemp=x(i,0)-x(j,0);
	R(i,j)=dtemp*dtemp;
      }
    for(int idim=1; idim<Ndim-1; ++idim)
      for(int j=0; j<NptsGuess-1; ++j)
	for(int i=j+1; i<NptsGuess; ++i) {
	  dtemp=x(i,idim)-x(j,idim);
	  R(i,j)+=dtemp*dtemp;
	}      
    for(int j=0; j<NptsGuess; ++j) {
      R(j,j)=1.0;
      for(int i=j+1; i<NptsGuess; ++i) {
	dtemp=x(i,Ndim)-x(j,Ndim);
	R(j,i)=R(i,j)=std::exp(negtheta*(R(i,j)+dtemp*dtemp));
      }
    }
    
    info=0; rank=-Npts;
    PIVOTCHOL_F77(&uplo, &NptsGuess, R.ptr(0,0), &ld_R,
		  ipiv.ptr(0,0), &rank, &min_allowed_rcond, &info); 
    int itemp=rank;
    if(rank>Npts)
      rank=Npts;
    for(int idim=0; idim<Ndim; ++idim)
      for(int i=0; i<rank; ++i)
	xfinal(i,idim)=x(ipiv(i,0)-1,idim);        
    printf("iLoop=%d done rank=%d/%d\n",iLoop,itemp,NptsGuess);
  }
  
  FILE* fpout=fopen(filename.c_str(),"w");
  for(int i=0; i<Npts; ++i) {
    fprintf(fpout,"%14.12f",xfinal(i,0));
    for(int idim=1; idim<Ndim; ++idim)
      fprintf(fpout," %14.12f",xfinal(i,idim));
    fprintf(fpout,"\n");
  }
  fclose(fpout);

  return;
}

//small fast tests to run valgrind on to catch memory errors and other debugging
int valgrind_test(int itest){
  std::cout << "********************************************************************************\n"
	    << "Running Valgrind/debugging test # " << itest << "\n"
	    << "********************************************************************************"
	    << std::endl;

  int ndim=2;
  int ndim_int=0;
  int nout=3;
  int iout=0;
  int sd_der_order=1;
  int ncol_skip=0;
  std::string sdbuildfilename;
  std::string sdevalfilename="grad_validate2d_10.spd";
  std::map< std::string, std::string> km_params;  
  km_params["order"]="4";
  //km_params["matern"]="2.5";
  if(itest<6) {
    sdbuildfilename="grad_validate2d_10.spd";
    std::cout << "test function = 2D Rosenbrock\n"
	      << "build file=\"" << sdbuildfilename << "\"\n" 
	      << "eval file=\"" << sdevalfilename << "\"" 
	      << std::endl;


    if(itest>=3) { //tests 3, 4, 5 are GEK
      km_params["derivative_order"]="1";
      std::cout << "Model Type is Gradient Enhanced Kriging" << std::endl;
    }
    else //tests 0, 1, 2 are Kriging
      std::cout << "Model Type is Kriging" << std::endl;

    if((itest==1)||(itest==4)) { //find a nugget if it is needed
      km_params["find_nugget"]="1";
      std::cout << "Handling ill conditioning by adding a nugget (it "
		<< "might not be activated)" << std::endl;
    }
    else if((itest==2)||(itest==5)) {//prescribed nugget
      km_params["nugget"]="0.001";
      std::cout << "Adding a prescribed nugget" << std::endl;
    }
    else{
      //test 0 and test 3 use pivoted cholesky to handle ill conditioning
      std::cout << "Handling ill conditioning with pivoted Cholesky (it " 
		<< "might not be activated)" << std::endl;
    }
  }
  else if(itest<10) {    
    sdbuildfilename="grad_validate2d_100.spd";
    std::cout << "test function = 2D Rosenbrock\n"
	      << "build file=\"" << sdbuildfilename << "\"\n" 
	      << "eval file=\"" << sdevalfilename << "\"" 
	      << std::endl;
    
    if(itest>=8) { //tests 8 and 9 are GEK
      km_params["derivative_order"]="1";
      std::cout << "Model Type is Gradient Enhanced Kriging" << std::endl;
    }
    else //tests 6 and 7 are Kriging
      std::cout << "Model Type is Kriging" << std::endl;     

    if((itest==7)||(itest==9)) { //find a nugget if it is needed
      km_params["find_nugget"]="1";
      std::cout << "Handling ill conditioning by adding a nugget (this should\n"
		<< "get activated for the Gaussian correlation function)"
		<< std::endl;
    }
    else{
      //tests 6 and 8 use pivoted cholesky
      std::cout << "Handling ill conditioning with pivoted Cholesky (this\n" 
		<< "should get activated for the Gaussian correlation function)"
		<< std::endl;
    }
  }
  else{
    std::cerr << "unknown valgrind test number, itest=" << itest << std::endl;
    return 1;
  }
  std::cout << "********************************************************************************"
	    << std::endl;

  std::cerr << "reading evaluation surfdata" << std::endl;
  nkm::SurfData sde(sdevalfilename,
		    ndim,ndim_int,nout,iout,sd_der_order,ncol_skip);
  std::cerr << "reading build surfdata" << std::endl;
  nkm::SurfData sdb(sdbuildfilename,
		    ndim,ndim_int,nout,iout,sd_der_order,ncol_skip);
  std::cerr << "calling model constructor" << std::endl;  
  nkm::KrigingModel km(sdb, km_params);   
  std::cerr << "calling create()" << std::endl;  
  km.create();
  std::cerr << "declaring ye, d1ye, d2ye, vye" << std::endl;
  nkm::MtxDbl ye(1,10), d1ye(2,10), d2ye(3,10), vye(1,10);
  std::cerr << "evaluating model adjusted mean, ye" << std::endl;
  km.evaluate(ye,sde.xr);
  std::cerr << "storing adjusted mean in evaluation surfdata" << std::endl;  
  sde.y.putRows(ye,iout); 
  std::cerr << "evaluating first derivative of y, d1ye" << std::endl;
  km.evaluate_d1y(d1ye,sde.xr);
  std::cerr << "storing d1ye in evaluation surfdata" << std::endl;
  sde.putDerY(d1ye,1); 
  std::cerr << "evaluating second derivative of y, d2ye" << std::endl;
  km.evaluate_d2y(d2ye,sde.xr);
  std::cerr << "storing d2ye in evaluation surfdata" << std::endl;
  sde.putDerY(d2ye,2); 
  std::cerr << "writing evaluation surfdata" << std::endl;
  sde.write("valgrind_test_output.spd");
  std::cerr << "evaluating adjusted variance" << std::endl;
  km.eval_variance(vye,sde.xr);
  return 0;
}




void timing_tests() {  
  std::map< std::string, std::string> km_params;
  struct timeval tv;  

  nkm::MtxDbl yeval(1,10000);  
  int ndim;
  string sdbuildfilename;
  string sdevalfilename;
  string problemname;
  string modelname;
  int nout;
  int iout;
  
  std::cout << "********************************************************************************\n"
	    << "********************************************************************************\n"
	    << "Running Timing/Performance Tests\n"
	    << "********************************************************************************"
	    << std::endl;

  
  km_params["optimization_method"]="global";
  for(int itest=0; itest<6; ++itest) {
    if(itest<4) {
      nout=3;
      ndim=2;
      sdbuildfilename="grad_validate2d_500.spd";
      sdevalfilename ="grad_validate2d_10K.spd";
      km_params["lower_bounds"]="-2.0 -2.0";
      km_params["upper_bounds"]="2.0 2.0";
      if(itest<2) {
	iout=0;
	problemname="2D Rosenbrock";
      }
      else{
	iout=2;
	problemname="2D Herbie";
      }
      if((itest==0)||(itest==2)) {
	km_params["derivative_order"]="0";
	modelname="Kriging";
      }
      else {
	km_params["derivative_order"]="1";
	modelname="Gradient Enhanced Kriging";
      }
    }
    else if(itest<6) {
      sdevalfilename="grad_paviani10d_10K.spd";
      nout=1;
      iout=0;
      ndim=10;
      problemname="10D Paviani";      
      km_params["lower_bounds"]=
	" 2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0";
      km_params["upper_bounds"]=
	"10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0";
      if(itest==4) {
	sdbuildfilename="grad_paviani10d_2500.spd";    
	km_params["derivative_order"]="0";
	modelname="Kriging";
      }
      else {
	sdbuildfilename="grad_paviani10d_500.spd";    
	km_params["derivative_order"]="1";
	modelname="Gradient Enhanced Kriging";
      }
    }

    std::cout << 
"********************************************************************************\n"
	      << "Building a " << modelname << " model for " << problemname
	      << "\nfrom data file \"" << sdbuildfilename << "\" using "
	      << km_params["optimization_method"] << " optimization" 
	      << std::endl;

    int sd_der_order=1;
    int ncol_skip=0;
    int ndim_int=0;

    //time the reading of eval surfdata
    gettimeofday(&tv, NULL);
    long int sde_read_start_sec=tv.tv_sec;
    long int sde_read_start_usec=tv.tv_usec;
    nkm::SurfData sde(sdevalfilename,
		      ndim,ndim_int,nout,iout,sd_der_order,ncol_skip);
    gettimeofday(&tv, NULL);
    long int sde_read_stop_sec=tv.tv_sec;
    long int sde_read_stop_usec=tv.tv_usec;
    double sde_read_time=
      static_cast<double>(sde_read_stop_sec -sde_read_start_sec)+
      static_cast<double>(sde_read_stop_usec-sde_read_start_usec)/1000000.0;

    //time the reading of build surfdata
    gettimeofday(&tv, NULL);
    long int sdb_read_start_sec=tv.tv_sec;
    long int sdb_read_start_usec=tv.tv_usec;
    nkm::SurfData sdb(sdbuildfilename,
		     ndim,ndim_int,nout,iout,sd_der_order,ncol_skip);
    gettimeofday(&tv, NULL);
    long int sdb_read_stop_sec=tv.tv_sec;
    long int sdb_read_stop_usec=tv.tv_usec;
    double sdb_read_time=
      static_cast<double>(sdb_read_stop_sec -sdb_read_start_sec)+
      static_cast<double>(sdb_read_stop_usec-sdb_read_start_usec)/1000000.0;

    //time the model constructor (processing of km params and initialize 
    //some things)
    gettimeofday(&tv, NULL);
    long int km_construct_start_sec =tv.tv_sec;
    long int km_construct_start_usec=tv.tv_usec;
    nkm::KrigingModel km(sdb, km_params); 
    gettimeofday(&tv, NULL);
    long int km_construct_stop_sec =tv.tv_sec;
    long int km_construct_stop_usec=tv.tv_usec;
    double km_construct_time=
      static_cast<double>(km_construct_stop_sec -km_construct_start_sec)+
      static_cast<double>(km_construct_stop_usec-km_construct_start_usec)/
      1000000.0;

    //time the model creation (the optimization to find correlation lengths)
    gettimeofday(&tv, NULL);
    long int km_create_start_sec =tv.tv_sec;
    long int km_create_start_usec=tv.tv_usec;
    km.create();
    gettimeofday(&tv, NULL);
    long int km_create_stop_sec =tv.tv_sec;
    long int km_create_stop_usec=tv.tv_usec;
    double km_create_time=
      static_cast<double>(km_create_stop_sec -km_create_start_sec)+
      static_cast<double>(km_create_stop_usec-km_create_start_usec)/1000000.0;

    //time the evaluation of the model adjuseted mean at 10K points
    gettimeofday(&tv, NULL);
    long int km_eval_start_sec =tv.tv_sec;
    long int km_eval_start_usec=tv.tv_usec;
    km.evaluate(yeval,sde.xr);
    gettimeofday(&tv, NULL);
    long int km_eval_stop_sec =tv.tv_sec;
    long int km_eval_stop_usec=tv.tv_usec;
    double km_eval_time=
      static_cast<double>(km_eval_stop_sec -km_eval_start_sec)+
      static_cast<double>(km_eval_stop_usec-km_eval_start_usec)/1000000.0;

    //might want to also time the putting of the evalutated means into the eval
    //surfdata
    sde.y.putRows(yeval,iout); 
    //sde.putDerY(yeval,0); //could also do this
    //sde.putDerY(yeval,0,iout); //could also do this

    //time the writing of eval surfdata
    gettimeofday(&tv, NULL);
    long int sde_write_start_sec =tv.tv_sec;
    long int sde_write_start_usec=tv.tv_usec;
    sde.write("timing_test_output.spd");
    gettimeofday(&tv, NULL);
    long int sde_write_stop_sec =tv.tv_sec;
    long int sde_write_stop_usec=tv.tv_usec;
    double sde_write_time=
      static_cast<double>(sde_write_stop_sec -sde_write_start_sec)+
      static_cast<double>(sde_write_stop_usec-sde_write_start_usec)/1000000.0;

    std::cout << "When building a " << modelname << " model for " << problemname
	      << "\nfrom data file \"" << sdbuildfilename << "\" using "
	      << km_params["optimization_method"] << " optimization\n" 
	      << "reading the evaluation surfdata (10K points) took " 
	      << sde_read_time << " seconds\n"
	      << "reading the build surfdata took " 
	      << sdb_read_time << " seconds\n"
	      << "model constructor took " 
	      << km_construct_time << " seconds\n"
	      << "model create() took " 
	      << km_create_time << " seconds\n"
	      << "model evaluation at 10K points took "
	      << km_eval_time << " seconds\n"
	      << "writing the evaluation surfdata (10K points) took " 
	      << sde_write_time << " seconds" 
	      << std::endl;
  }
  std::cout << "********************************************************************************"
	    << std::endl;
  return;
}

int main(int argc, char* argv[])
{
  
  if(argc < 2) {
    std::cerr << "you need to specify what test you want to run" << std::endl;
    return 1;
  }
  std::string testname;
  {  
    std::ostringstream oss;
    oss << argv[1];
    testname=oss.str();
  }
  if(testname=="timing")
    timing_tests();
  else if(testname=="nightly")
    nightly_test();
  else if(testname=="valgrind") {
    int itest=0;
    if(argc>=3) {
      std::stringstream ss;
      ss << argv[2];
      ss >> itest;
    }
    return valgrind_test(itest);
  }
  else if(testname=="GPAIS") {
    //command line input is like this
    //e.g. <executable_name> GPAIS 6D GEK global matern infinity
    //e.g. <executable_name> GPAIS 6D gek global powered_exponential 2
    //
    //e.g. <executable_name> GPAIS 2D kriging local matern 0.5
    //e.g. <executable_name> GPAIS 2D Krig local powered_exponential 1
    if(argc!=7) {
      std::cerr << "some examples of GPAIS commmand line input are:\n"
		<< argv[0] << " GPAIS 6D GEK global matern infinity\n"
		<< argv[0] << " GPAIS 6D gek global powered_exponential 2\n"
		<< argv[0] << " GPAIS 2D kriging local matern 0.5\n"
		<< argv[0] << " GPAIS 2D Krig local powered_exponential 1\n"
		<< argv[0] << " GPAIS 2D KM none matern 1.5\n"
		<< argv[0] << " GPAIS 2D krig none matern 1.5\n"
		<< "example 1 is equivalent to example 2\n"
		<< "example 3 is equivalent to example 4\n"
		<< "example 5 is equivalent to example 6\n" 
		<< "examples 3 through 6 require a file named \"corrlen.txt\""
		<< " to be present\n"
		<< "files named \"buildfile.spd;\" and \"evalfile_in.spd\" "
		<< "must always be present"
		<< std::endl;
      return 1;
    }
    int Nvarsr;
    { //{ if for scoping
      std::stringstream ss;
      ss << argv[2];
      ss >> Nvarsr; //this should take the integer before the "D"
    }
    std::string model_type;
    { //{ is for scoping
      std::ostringstream oss;
      oss << argv[3];
      model_type = oss.str();
    }
    std::string optimization_method;
    { //{ is for scoping
      std::ostringstream oss;
      oss << argv[4];
      optimization_method = oss.str();
    }
    std::string corr_func_family;
    { //{ is for scoping
      std::ostringstream oss;
      oss << argv[5];
      corr_func_family = oss.str();
    }
    std::string corr_func_param;
    { //{ is for scoping
      std::ostringstream oss;
      oss << argv[6];
      corr_func_param = oss.str();
    }
    GPAIS_build_and_eval_mean_adjvar(Nvarsr, model_type, optimization_method,
				     corr_func_family, corr_func_param);
  }
  else if(testname=="generic"){
    std::string handle_ill_cond1="pivot_chol";
    std::string handle_ill_cond2="nothing";

    //command line input is like this
    //e.g. <executable_name> generic 6D GEK global matern infinity 
    //e.g. <executable_name> generic 6D gek global powered_exponential 2 nugget 0.00000001
    //
    //e.g. <executable_name> generic 2D kriging local matern 0.5 find_nugget 0
    //e.g. <executable_name> generic 2D Krig local powered_exponential 1 find_nugget 1
    if((argc!=7)&&(argc!=9)) {
      std::cerr << "some examples of generic commmand line input are:\n"
		<< argv[0] << " generic 6D GEK global matern infinity\n"
		<< argv[0] << " generic 6D gek global powered_exponential 2\n"
		<< argv[0] << " generic 2D kriging local matern 0.5\n"
		<< argv[0] << " generic 2D Krig local powered_exponential 1\n"
		<< argv[0] << " generic 2D KM none matern 1.5\n"
		<< argv[0] << " generic 2D krig none matern 1.5\n"
		<< argv[0] << " generic 2D krig global powered_exponential 2 nugget 0.00000001\n"
		<< argv[0] << " generic 2D krig global powered_exponential 2 find_nugget 0\n"
		<< argv[0] << " generic 2D krig global powered_exponential 2 find_nugget 1\n"
		<< "example 1 is equivalent to example 2\n"
		<< "example 3 is equivalent to example 4\n"
		<< "example 5 is equivalent to example 6\n" 
		<< "examples 3 through 6 require a file named \"corrlen.txt\""
		<< " to be present\n"
		<< "files named \"buildfile.spd;\" and \"evalfile_in.spd\" "
		<< "must always be present"
		<< std::endl;
      return 1;
    }
    int Nvarsr;
    { //{ if for scoping
      std::stringstream ss;
      ss << argv[2];
      ss >> Nvarsr; //this should take the integer before the "D"
    }
    std::string model_type;
    { //{ is for scoping
      std::ostringstream oss;
      oss << argv[3];
      model_type = oss.str();
    }
    std::string optimization_method;
    { //{ is for scoping
      std::ostringstream oss;
      oss << argv[4];
      optimization_method = oss.str();
    }
    std::string corr_func_family;
    { //{ is for scoping
      std::ostringstream oss;
      oss << argv[5];
      corr_func_family = oss.str();
    }
    std::string corr_func_param;
    { //{ is for scoping
      std::ostringstream oss;
      oss << argv[6];
      corr_func_param = oss.str();
    }
    if(argc==9) {
      std::ostringstream oss1;
      oss1 << argv[7];
      handle_ill_cond1 = oss1.str();
      std::ostringstream oss2;
      oss2 << argv[8];
      handle_ill_cond2 = oss2.str();
    }

    generic_build_and_eval_mean_adjvar(Nvarsr, model_type, optimization_method,
				       corr_func_family, corr_func_param,
				       handle_ill_cond1, handle_ill_cond2);
  }
  else{
    std::cerr << "unknown test name" << std::endl;
    return 1;
  }

  return 0;
}



void check_matrix() 
{
  nkm::MtxDbl I10(10,10), I5;

  I10.identity();
  I5.copy(I10);
  I5.resize(5,5);

  nkm::MtxDbl a(10,10),  A(10,10), AChol(10,10), Ainv(10,10), IA(10,10), IA2(10,10), IA3(10,10);
  
  int nrows=10, ncols=10;

  for(int j=0; j<nrows; ++j)
    for(int i=0; i<ncols; ++i) 
      a(i,j)=static_cast<double>(std::rand());

  nkm::matrix_mult(A,a,a,0.0,1.0,'N','T');
  AChol.copy(A);
  int info;
  double rcond;
  nkm::Chol_fact(AChol,info,rcond);
  Ainv.copy(AChol);
  nkm::inverse_after_Chol_fact(Ainv);
  nkm::matrix_mult(IA,A,Ainv);
  nkm::matrix_mult(IA2,Ainv,A);
  nkm::solve_after_Chol_fact(IA3,AChol,A);

  nkm::MtxDbl B(15,15), BChol(15,15), Binv(15,15), Binv2(15,15), IB(15,15), IB2(15,15), IB3(15,15), IB4(14,14);

  B.copy(A);
  BChol.copy(B);
  nkm::Chol_fact(BChol,info,rcond);  
  Binv.copy(BChol);
  nkm::inverse_after_Chol_fact(Binv);
  nkm::matrix_mult(IB,B,Binv);
  nkm::matrix_mult(IB2,Binv,B);
  nkm::solve_after_Chol_fact(IB3,BChol,B);
  nkm::solve_after_Chol_fact(Binv2,BChol,I10);
  nkm::matrix_mult(IB4,Binv2,B,0.0,1.0,'N','T');

  FILE *fp=fopen("test_mat.txt","w");
  nrows=I10.getNRows();
  ncols=I10.getNCols();
  fprintf(fp,"I10 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,I10.getNRowsAct(),I10.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",I10(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",I10(i,j));
  }

  nrows=I5.getNRows();
  ncols=I5.getNCols();
  fprintf(fp,"\n\nI5 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,I5.getNRowsAct(),I5.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",I5(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",I5(i,j));
  }


  nrows=a.getNRows();
  ncols=a.getNCols();
  fprintf(fp,"\n\na (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,a.getNRowsAct(),a.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",a(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",a(i,j));
  }

  nrows=A.getNRows();
  ncols=A.getNCols();
  fprintf(fp,"\n\nA (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,A.getNRowsAct(),A.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",A(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",A(i,j));
  }

  
  int j;
  nrows=AChol.getNRows();
  ncols=AChol.getNCols();
  fprintf(fp,"\n\nAChol (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,AChol.getNRowsAct(),AChol.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",AChol(i,0));
    for(j=1; j<=i; ++j)
      fprintf(fp," %12.6g",AChol(i,j));
    for(; j<ncols; ++j)
      fprintf(fp," %12.6g",0.0);
  }
  

  
  nrows=Ainv.getNRows();
  ncols=Ainv.getNCols();
  fprintf(fp,"\n\nAinv (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Ainv.getNRowsAct(),Ainv.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Ainv(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Ainv(i,j));
  }

  nrows=IA.getNRows();
  ncols=IA.getNCols();
  fprintf(fp,"\n\nIA (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IA.getNRowsAct(),IA.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IA(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IA(i,j));
  }


  nrows=IA2.getNRows();
  ncols=IA2.getNCols();
  fprintf(fp,"\n\nIA2 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IA2.getNRowsAct(),IA2.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IA2(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IA2(i,j));
  }


  nrows=IA3.getNRows();
  ncols=IA3.getNCols();
  fprintf(fp,"\n\nIA3 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IA3.getNRowsAct(),IA3.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IA3(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IA3(i,j));
  }


  nrows=B.getNRows();
  ncols=B.getNCols();
  fprintf(fp,"\n\nB (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,B.getNRowsAct(),B.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",B(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",B(i,j));
  }

  nrows=BChol.getNRows();
  ncols=BChol.getNCols();
  fprintf(fp,"\n\nBChol (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,BChol.getNRowsAct(),BChol.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",BChol(i,0));
    for(j=1; j<=i; ++j)
      fprintf(fp," %12.6g",BChol(i,j));
    for(; j<ncols; ++j)
      fprintf(fp," %12.6g",0.0);
  }
  
  nrows=Binv.getNRows();
  ncols=Binv.getNCols();
  fprintf(fp,"\n\nBinv (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Binv.getNRowsAct(),Binv.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Binv(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Binv(i,j));
  }

  nrows=Binv2.getNRows();
  ncols=Binv2.getNCols();
  fprintf(fp,"\n\nBinv2 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Binv2.getNRowsAct(),Binv2.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Binv2(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Binv2(i,j));
  }

  nrows=IB.getNRows();
  ncols=IB.getNCols();
  fprintf(fp,"\n\nIB (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IB.getNRowsAct(),IB.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IB(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IB(i,j));
  }

  nrows=IB2.getNRows();
  ncols=IB2.getNCols();
  fprintf(fp,"\n\nIB2 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IB2.getNRowsAct(),IB2.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IB2(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IB2(i,j));
  }

  nrows=IB3.getNRows();
  ncols=IB3.getNCols();
  fprintf(fp,"\n\nIB3 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IB3.getNRowsAct(),IB3.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IB3(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IB3(i,j));
  }


  nrows=IB4.getNRows();
  ncols=IB4.getNCols();
  fprintf(fp,"\n\nIB4 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IB4.getNRowsAct(),IB4.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IB4(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IB4(i,j));
  }


  nkm::MtxDbl C(5,5),CChol(5,5), Cinv(5,5), IC(5,5);

  nrows=ncols=5;
  for(int j=0; j<ncols; ++j) 
    for(int i=0; i<nrows; ++i)
      C(i,j)=A(i,j);
  //for(int k=0; k<nrows*ncols; ++k)
  //C(k)=A(k);

  CChol.copy(C);
  nkm::Chol_fact(CChol,info,rcond);  
  Cinv.copy(CChol);
  nkm::inverse_after_Chol_fact(Cinv);
  nkm::matrix_mult(IC,C,Cinv);

  nrows=C.getNRows();
  ncols=C.getNCols();
  fprintf(fp,"\n\nC (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,C.getNRowsAct(),C.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",C(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",C(i,j));
  }

  
  nrows=CChol.getNRows();
  ncols=CChol.getNCols();
  fprintf(fp,"\n\nCChol (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,CChol.getNRowsAct(),CChol.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",CChol(i,0));
    for(j=1; j<=i; ++j)
      fprintf(fp," %12.6g",CChol(i,j));
    for(; j<ncols; ++j)
      fprintf(fp," %12.6g",0.0);
  }
    
  nrows=Cinv.getNRows();
  ncols=Cinv.getNCols();
  fprintf(fp,"\n\nCinv (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Cinv.getNRowsAct(),Cinv.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Cinv(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Cinv(i,j));
  }

  nrows=IC.getNRows();
  ncols=IC.getNCols();
  fprintf(fp,"\n\nIC (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IC.getNRowsAct(),IC.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IC(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IC(i,j));
  }



  nkm::MtxDbl D(15,15), DChol(15,15), Dinv(15,15), ID(15,15);

  D.copy(B); D.resize(5,5);

  DChol.copy(D);
  nkm::Chol_fact(DChol,info,rcond);  
  Dinv.copy(DChol);
  nkm::inverse_after_Chol_fact(Dinv);
  nkm::matrix_mult(ID,D,Dinv);


  nrows=D.getNRows();
  ncols=D.getNCols();
  fprintf(fp,"\n\nD (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,D.getNRowsAct(),D.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",D(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",D(i,j));
  }

  
  nrows=DChol.getNRows();
  ncols=DChol.getNCols();
  fprintf(fp,"\n\nDChol (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,DChol.getNRowsAct(),DChol.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",DChol(i,0));
    for(j=1; j<=i; ++j)
      fprintf(fp," %12.6g",DChol(i,j));
    for(; j<ncols; ++j)
      fprintf(fp," %12.6g",0.0);
  }
    
  nrows=Dinv.getNRows();
  ncols=Dinv.getNCols();
  fprintf(fp,"\n\nDinv (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Dinv.getNRowsAct(),Dinv.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Dinv(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Dinv(i,j));
  }

  nrows=ID.getNRows();
  ncols=ID.getNCols();
  fprintf(fp,"\n\nID (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,ID.getNRowsAct(),ID.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",ID(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",ID(i,j));
  }

  nkm::MtxDbl E1(13,13), E2(13,13), E3(15,25);
  E1.newSize(10,10);
  nrows=ncols=10;
  for(int j=0; j<ncols; ++j)
    for(int i=0; i<nrows; ++i)
      E1(i,j)=B(i,j);
  //for(int ij=0; ij<nrows*ncols; ++ij)
  //E1(ij)=B(ij);

  nrows=E1.getNRows();
  ncols=E1.getNCols();
  fprintf(fp,"\n\nE1 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,E1.getNRowsAct(),E1.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",E1(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",E1(i,j));
  }

  E1.reshape(20,5);
  nrows=E1.getNRows();
  ncols=E1.getNCols();
  fprintf(fp,"\n\nE1 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,E1.getNRowsAct(),E1.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",E1(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",E1(i,j));
  }
  
  E2.copy(E1);
  E2.reshape(10,10);

  nrows=E2.getNRows();
  ncols=E2.getNCols();
  fprintf(fp,"\n\nE2 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,E2.getNRowsAct(),E2.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",E2(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",E2(i,j));
  }

  E3.copy(B);
  E3.reshape(5,20);

  nrows=E3.getNRows();
  ncols=E3.getNCols();
  fprintf(fp,"\n\nE3 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,E3.getNRowsAct(),E3.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",E3(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",E3(i,j));
  }



  fclose(fp);

  return;
}


#ifndef __DONT_COMPILE__

void nested_Krig_vs_GEK_herbie_smooth_herbie_2D_4D_8D(){
  //string build_herbie_2D="gradHerbie_NestedLHS_2D_1024pts.spd";
  //string build_smooth_2D="gradSmoothHerbie_NestedLHS_2D_1024pts.spd";
  string build_herbie_2D="gradHerbie_NestedLHS_2D_2048pts.spd";
  string build_smooth_2D="gradSmoothHerbie_NestedLHS_2D_2048pts.spd";
  //string build_herbie_2D="gradHerbie_PivotChol_2D_1024pts.spd";
  //string build_smooth_2D="gradSmoothHerbie_PivotChol_2D_1024pts.spd";
  string valid_herbie_2D="gradHerbie_NestedLHS_2D_16384pts.spd";
  string valid_smooth_2D="gradSmoothHerbie_NestedLHS_2D_16384pts.spd";

  //string build_herbie_4D="gradHerbie_NestedLHS_4D_1024pts.spd";
  //string build_smooth_4D="gradSmoothHerbie_NestedLHS_4D_1024pts.spd";
  string build_herbie_4D="gradHerbie_NestedLHS_4D_4096pts.spd";
  string build_smooth_4D="gradSmoothHerbie_NestedLHS_4D_4096pts.spd";
  //string build_herbie_4D="gradHerbie_PivotChol_4D_2048pts.spd";
  //string build_smooth_4D="gradSmoothHerbie_PivotChol_4D_2048pts.spd";
  string valid_herbie_4D="gradHerbie_NestedLHS_4D_16384pts.spd";
  string valid_smooth_4D="gradSmoothHerbie_NestedLHS_4D_16384pts.spd";

  //string build_herbie_8D="gradHerbie_NestedLHS_8D_512pts.spd";
  //string build_smooth_8D="gradSmoothHerbie_NestedLHS_8D_512pts.spd";
  string build_herbie_8D="gradHerbie_NestedLHS_8D_2048pts.spd";
  string build_smooth_8D="gradSmoothHerbie_NestedLHS_8D_2048pts.spd";
  //string build_herbie_8D="gradHerbie_PivotChol_8D_1024pts.spd";
  //string build_smooth_8D="gradSmoothHerbie_PivotChol_8D_1024pts.spd";
  string valid_herbie_8D="gradHerbie_NestedLHS_8D_16384pts.spd";
  string valid_smooth_8D="gradSmoothHerbie_NestedLHS_8D_16384pts.spd";

  string build_herbie_filename, valid_herbie_filename;
  string build_smooth_filename, valid_smooth_filename;

  FILE *fpout1=fopen("GradKrigingPaperHerbieEffectOfDimensionStudyTableNestedLHSPivotCholKrigR.txt","w");
  FILE *fpout2=fopen("GradKrigingPaperSmoothHerbieEffectOfDimensionStudyTableNestedLHSPivotCholKrigR.txt","w");

  std::map< std::string, std::string> herbie_krig_params;    
  std::map< std::string, std::string> herbie_GEK_params;    
  std::map< std::string, std::string> smooth_krig_params;    
  std::map< std::string, std::string> smooth_GEK_params;    
  herbie_krig_params["order"] = "2";
  herbie_krig_params["reduced_polynomial"]=nkm::toString<bool>(true);
  herbie_krig_params["optimization_method"]="none";

  for(int ndimpow=1; ndimpow<=3; ++ndimpow) { //loop over the number of
    //dimensions for the effect of dimension study
    
    int Ndim=static_cast<int> (std::pow(2.0,static_cast<double>(ndimpow)));
    switch(Ndim){
    case 2:
      herbie_krig_params["lower_bounds"]="-2.0 -2.0";
      herbie_krig_params["upper_bounds"]="2.0 2.0";
      herbie_GEK_params=herbie_krig_params;
      herbie_GEK_params["derivative_order"]="1";
      smooth_krig_params=herbie_krig_params;
      smooth_GEK_params =herbie_GEK_params;

      build_herbie_filename=build_herbie_2D;
      build_smooth_filename=build_smooth_2D;
      valid_herbie_filename=valid_herbie_2D;
      valid_smooth_filename=valid_smooth_2D;
      
      break;
    case 4:
      herbie_krig_params["lower_bounds"]="-2.0 -2.0 -2.0 -2.0";
      herbie_krig_params["upper_bounds"]="2.0 2.0 2.0 2.0";
      herbie_GEK_params=herbie_krig_params;
      herbie_GEK_params["derivative_order"]="1";
      smooth_krig_params=herbie_krig_params;
      smooth_GEK_params =herbie_GEK_params;

      build_herbie_filename=build_herbie_4D;
      build_smooth_filename=build_smooth_4D;
      valid_herbie_filename=valid_herbie_4D;
      valid_smooth_filename=valid_smooth_4D;

      break;
    case 8:
      herbie_krig_params["lower_bounds"]="-2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0";
      herbie_krig_params["upper_bounds"]="2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0";
      herbie_GEK_params=herbie_krig_params;
      herbie_GEK_params["derivative_order"]="1";
      smooth_krig_params=herbie_krig_params;
      smooth_GEK_params =herbie_GEK_params;

      build_herbie_filename=build_herbie_8D;
      build_smooth_filename=build_smooth_8D;
      valid_herbie_filename=valid_herbie_8D;
      valid_smooth_filename=valid_smooth_8D;

      break;
    default:
      std::cerr << "Error: haven't coded for the " << Ndim << " dimensional case." << std::endl;
      assert(false);
    } 
    nkm::SurfData sd_build_herbie(build_herbie_filename, Ndim, 0, 1, 0, 1, 0);
    nkm::SurfData sd_build_smooth(build_smooth_filename, Ndim, 0, 1, 0, 1, 0);
    nkm::SurfData sd_valid_herbie(valid_herbie_filename, Ndim, 0, 1, 0, 1, 0);
    nkm::SurfData sd_valid_smooth(valid_smooth_filename, Ndim, 0, 1, 0, 1, 0);
    nkm::SurfData sd_build_temp;
    
    
    int NptsBuild=sd_build_herbie.getNPts();
    int NptsValid=sd_valid_herbie.getNPts();
    assert((NptsBuild==sd_build_smooth.getNPts())&&
	   (NptsValid==sd_valid_smooth.getNPts()));
    int Nref=static_cast<int>(std::log(static_cast<double>(NptsBuild)/static_cast<double>(2*Ndim))/std::log(2.0));
    //if(Nref>2) Nref=2; //fast test for debug
    nkm::MtxDbl yeval(NptsValid,1);
    nkm::MtxInt ipts(NptsBuild,1);
    nkm::MtxDbl error_metric(Nref+1,9);  error_metric.zero();
    for(int i=0; i<NptsBuild; ++i)
      ipts(i,0)=i;
    for(int iref=0; iref<=Nref; ++iref) { //nested sample design loop
      int NptsThis=2*Ndim*static_cast<int>(std::pow(2.0,static_cast<double>(iref)));
      if(true) {
	//make this run faster by feeding it the correlation lengths generated from the nested LHS design the first time I rank it (so can quickly calculate the mean absolute error, originally I only calculated the RMSE)
	switch(Ndim){
	case 2:
	  switch(NptsThis){
	  case 4:
	    herbie_GEK_params["correlation_lengths"] ="0.500044 0.500132";
	    smooth_GEK_params["correlation_lengths"] ="0.500044 0.500132";
	    break;
	  case 8:
	    herbie_krig_params["correlation_lengths"]="1.26045  0.353585";
	    herbie_GEK_params["correlation_lengths"] ="0.490945 0.399543";
	    smooth_krig_params["correlation_lengths"]="2.01378  0.353585";
	    smooth_GEK_params["correlation_lengths"] ="0.896517 0.918406";
	    break;
	  case 16:
	    herbie_krig_params["correlation_lengths"]="2.11395  0.481163";
	    herbie_GEK_params["correlation_lengths"] ="0.269972 0.327727";
	    smooth_krig_params["correlation_lengths"]="0.250022 0.535061";
	    smooth_GEK_params["correlation_lengths"] ="0.901693 0.916258";
	    break;
	  case 32:
	    herbie_krig_params["correlation_lengths"]="0.337192 0.405882";
	    herbie_GEK_params["correlation_lengths"] ="0.364354 0.349339";
	    smooth_krig_params["correlation_lengths"]="0.801452 0.818999";
	    smooth_GEK_params["correlation_lengths"] ="0.954075 0.964208";
	    break;
	  case 64:
	    herbie_krig_params["correlation_lengths"]="0.362288 0.386812";
	    herbie_GEK_params["correlation_lengths"] ="0.394655 0.396955";
	    smooth_krig_params["correlation_lengths"]="0.942504 0.946496";
	    smooth_GEK_params["correlation_lengths"] ="0.708602 0.843246";
	    break;
	  case 128:
	    herbie_krig_params["correlation_lengths"]="0.341277 0.345084";
	    herbie_GEK_params["correlation_lengths"] ="0.393593 0.378305";
	    smooth_krig_params["correlation_lengths"]="0.713824 0.793784";
	    smooth_GEK_params["correlation_lengths"] ="0.681042 0.641354";
	    break;
	  case 256:
	    herbie_krig_params["correlation_lengths"]="0.427906 0.428056";
	    herbie_GEK_params["correlation_lengths"] ="0.343611 0.345492";
	    smooth_krig_params["correlation_lengths"]="0.409262 0.631571";
	    smooth_GEK_params["correlation_lengths"] ="0.662672 0.702565";
	    break;
	  case 512:
	    herbie_krig_params["correlation_lengths"]="0.332171 0.335345";
	    herbie_GEK_params["correlation_lengths"] ="0.328333 0.323684";
	    smooth_krig_params["correlation_lengths"]="0.223043 0.532758";
	    smooth_GEK_params["correlation_lengths"] ="0.686387 0.673577";
	    break;
	  case 1024:
	    herbie_krig_params["correlation_lengths"]="0.182501 0.255595";
	    herbie_GEK_params["correlation_lengths"] ="0.295556 0.261237";
	    smooth_krig_params["correlation_lengths"]="0.107725 0.473198";
	    smooth_GEK_params["correlation_lengths"] ="0.636192 0.639449";
	    break;
	  case 2048:
	    herbie_krig_params["correlation_lengths"]="0.176487 0.118631";
	    smooth_krig_params["correlation_lengths"]="0.176487 0.118631";
	    break;
	  default:
	    std::cerr << "NptsThis is not a known size" << std::endl;
	    assert(false);
	  } //switch(NptsThis)
	  break;
	case 4:
	  switch(NptsThis){
	  case 8:
	    herbie_GEK_params["correlation_lengths"] ="0.986542 1.25723 0.752366  1.20457";
	    smooth_GEK_params["correlation_lengths"] ="1.20457  1.10578 0.763173  1.05946";
	    break;
	  case 16:
	    herbie_krig_params["correlation_lengths"]="0.632662 0.503578 15.6613  1.30938";
	    herbie_GEK_params["correlation_lengths"] ="0.817832 0.903696 0.606163 0.817832";
	    smooth_krig_params["correlation_lengths"]="0.510812 0.510812 1.57612  0.719313";
	    smooth_GEK_params["correlation_lengths"] ="0.903696 0.817832 0.750758 0.817832";
	    break;
	  case 32:
	    herbie_krig_params["correlation_lengths"]="0.524469 13.1696  13.1696  0.579533";
	    herbie_GEK_params["correlation_lengths"] ="0.687712 0.770831 0.579533 0.658908";
	    smooth_krig_params["correlation_lengths"]="0.50972  4.32949  13.3587  0.55526";
	    smooth_GEK_params["correlation_lengths"] ="0.888995 0.816084 0.816084 0.839698";
	    break;
	  case 64:
	    herbie_krig_params["correlation_lengths"]="0.530866 1.75907  0.473623 0.508631";
	    herbie_GEK_params["correlation_lengths"] ="0.473623 0.578295 0.603575 0.629961";
	    smooth_krig_params["correlation_lengths"]="0.603575 0.376989 11.2333  0.508631";
	    smooth_GEK_params["correlation_lengths"] ="0.887095 0.912763 0.887095 0.925875";
	    break;
	  case 128:
	    herbie_krig_params["correlation_lengths"]="0.745955 0.937167 0.317008 0.602285";
	    herbie_GEK_params["correlation_lengths"] ="0.446403 0.43385  0.446403 0.507544";
	    smooth_krig_params["correlation_lengths"]="0.602285 1.61135  0.360427 0.593756";
	    smooth_GEK_params["correlation_lengths"] ="0.964285 0.950629 0.964285 0.964285";
	    break;
	  case 256:
	    herbie_krig_params["correlation_lengths"]="0.776901 0.683312 0.683312 0.799381";
	    herbie_GEK_params["correlation_lengths"] ="0.408916 0.408916 0.445449 0.458339";
	    smooth_krig_params["correlation_lengths"]="0.788061 0.776901 0.776901 0.810863";
	    smooth_GEK_params["correlation_lengths"] ="0.962224 0.962224 0.962224 0.962224";
	    break;
	  case 512:
	    herbie_krig_params["correlation_lengths"]="0.62593  0.574595 0.653293 0.672196";
	    herbie_GEK_params["correlation_lengths"] ="0.408042 0.390952 0.408042 0.396567";
	    smooth_krig_params["correlation_lengths"]="0.919951 0.919951 0.919951 0.906924";
	    smooth_GEK_params["correlation_lengths"] ="1.00214  0.960168 0.98795  1.00214";
	    break;
	  case 1024:
	    herbie_krig_params["correlation_lengths"]="0.549352 0.462937 0.504297 0.565248";
	    herbie_GEK_params["correlation_lengths"] ="0.373776 0.373776 0.379145 0.373776";
	    smooth_krig_params["correlation_lengths"]="0.958116 0.944548 0.958116 0.958116";
	    smooth_GEK_params["correlation_lengths"] ="0.842697 0.842697 0.842697 0.89217";
	    break;
	  case 2048:
	    herbie_krig_params["correlation_lengths"]="0.482142 0.357356 0.412136 0.482142";
	    herbie_GEK_params["correlation_lengths"] ="0.4426   0.436332 0.4426   0.4426";
	    smooth_krig_params["correlation_lengths"]="0.956068 0.956068 0.956068 0.956068";
	    smooth_GEK_params["correlation_lengths"] ="0.840896 0.840896 0.840896 0.840896";
	    break;
	  case 4096:
	    herbie_krig_params["correlation_lengths"]="0.39969  0.38845  0.37218  0.405432";
	    smooth_krig_params["correlation_lengths"]="0.875781 0.875781 0.875781 0.851153";
	    break;
	  default:
	    std::cerr << "NptsThis is not a known size" << std::endl;
	    assert(false);
	  } //switch(NptsThis)
	  break;
	case 8:
	  switch(NptsThis){
	  case 16:
	    herbie_GEK_params["correlation_lengths"] ="1.25992  1.25992  1.62868  1.25992  1.25992  1.25992  1.85175  1.25992";
	    smooth_GEK_params["correlation_lengths"] ="0.753977 1.25992  1.25992  0.857244 1.25992  1.25992  1.25992  1.25992";
	    break;
	  case 32:
	    herbie_krig_params["correlation_lengths"]="0.786096 17.1154  0.6914   17.1154  17.1154  17.1154  1.15535  17.1154";
	    herbie_GEK_params["correlation_lengths"] ="0.786096 0.786096 0.786096 0.786096 0.6914   0.786096 0.786096 0.786096";
	    smooth_krig_params["correlation_lengths"]="1.69806  17.1154  0.6914   1.69806  17.1154  17.1154  0.786096 1.15535";
	    smooth_GEK_params["correlation_lengths"] ="0.786096 1.15535  1.15535  0.786096 1.15535  1.15535  0.893762 0.786096";
	    break;
	  case 64:
	    herbie_krig_params["correlation_lengths"]="0.634017 4.94358  4.94358  2.28857  15.6949  1.05946  0.720853 15.6949";
	    herbie_GEK_params["correlation_lengths"] ="0.720853 0.720853 0.720853 0.720853 0.720853 0.720853 0.720853 0.819584";
	    smooth_krig_params["correlation_lengths"]="0.634017 4.94358  4.94358  4.94358  15.6949  1.05946  0.720853 15.6949";
	    smooth_GEK_params["correlation_lengths"] ="0.720853 1.05946  1.05946  1.05946  1.05946  1.05946  0.931836 1.05946";
	    break;
	  case 128:
	    herbie_krig_params["correlation_lengths"]="0.661025 0.581396 1.42789  1.42789  0.661025 1.42789  0.661025 14.3923";
	    herbie_GEK_params["correlation_lengths"] ="0.661025 0.661025 0.661025 0.661025 0.661025 0.661025 0.661025 0.751561";
	    smooth_krig_params["correlation_lengths"]="0.661025 0.661025 0.661025 1.42789  0.661025 0.661025 0.661025 0.581396";
	    smooth_GEK_params["correlation_lengths"] ="0.971532 0.971532 0.971532 0.971532 0.971532 0.971532 0.854497 0.971532";
	    break;
	  case 256:
	    herbie_krig_params["correlation_lengths"]="0.606163 1.30938  1.30938  0.606163 0.606163 1.30938  0.606163 0.533142";
	    herbie_GEK_params["correlation_lengths"] ="0.890899 0.890899 0.890899 0.890899 0.890899 0.890899 0.783578 0.890899";
	    smooth_krig_params["correlation_lengths"]="1.30938  0.533142 1.30938  0.606163 0.606163 0.606163 13.1978  0.606163";
	    smooth_GEK_params["correlation_lengths"] ="0.890899 0.890899 0.890899 1.01292  0.890899 0.890899 0.890899 0.890899";
	    break;
	  case 512:
	    herbie_krig_params["correlation_lengths"]="0.488894 1.20071  1.20071  0.555854 1.20071  1.20071  0.555854 0.555854";
	    herbie_GEK_params["correlation_lengths"] ="0.816958 0.816958 0.816958 0.816958 0.816958 0.816958 0.718544 0.816958";
	    smooth_krig_params["correlation_lengths"]="0.555854 0.555854 1.20071  1.20071  1.20071  1.20071  0.555854 0.488894";
	    smooth_GEK_params["correlation_lengths"] ="0.816958 0.816958 0.816958 0.928851 0.816958 0.816958 0.816958 0.816958";
	    break;
	  case 1024:
	    herbie_krig_params["correlation_lengths"]="1.25186  0.50972  3.49564  0.50972  1.10106  0.50972  1.10106  0.50972";
	    herbie_GEK_params["correlation_lengths"] ="0.749154 0.749154 0.749154 0.749154 0.749154 0.749154 0.749154 0.749154";
	    smooth_krig_params["correlation_lengths"]="2.37841  0.50972  5.13766  0.50972  1.25186  0.749154 0.749154 0.50972";
	    smooth_GEK_params["correlation_lengths"] ="0.749154 0.749154 0.749154 0.85176  0.749154 0.749154 0.749154 0.749154";
	    break;
	  case 2048:
	    herbie_krig_params["correlation_lengths"]="0.686977 0.686977 0.531434 0.467416 4.71125  3.20551  0.686977 0.686977";
	    smooth_krig_params["correlation_lengths"]="2.18102  0.467416 0.531434 0.467416 3.20551  2.18102  0.686977 0.686977";
	    break;
	  default:
	    std::cerr << "NptsThis is not a known size" << std::endl;
	    assert(false);
	  } //switch(NptsThis) 
	  break;
	default:
	  std::cerr << "Ndim is not a known number of dimensions" << std::endl;
	  assert(false);
	} //switch{Ndim) 
      }


      error_metric(iref,0)=static_cast<double>(NptsThis);
      ipts.resize(NptsThis,1); //relies on actual and apparent sizes of the matrix 
      //class being different and that resize() doesn't copy or overwrite or shrink 
      //or enlarge unless it needs a bigger size than it actually has OR the user 
      //"forces" it to resize, neither of these 2 cases are true in the original 
      //implementation.of this function

      {//limit km and gkm for herbie to this scope

	sd_build_herbie.getPoints(sd_build_temp,ipts);
	printf("Herbie: Ndim=%d Npts=%4d/%-4d",Ndim,NptsThis,NptsBuild);
	fprintf(fpout1,"Herbie: Ndim=%d Npts=%4d/%-4d",Ndim,NptsThis,NptsBuild);

	if(iref>0) { //if have enough equations for a reduced quadratic trend
	  nkm::KrigingModel km(sd_build_temp,herbie_krig_params); km.create();
	  km.evaluate(yeval,sd_valid_herbie.xr);
	  for(int i=0; i<NptsValid; ++i) {
	    double tmpdbl=yeval(i,0)-sd_valid_herbie.y(i,0);
	    error_metric(iref,1)+=std::fabs(tmpdbl);
	    error_metric(iref,2)+=std::pow(tmpdbl,2.0);
	  }
	  error_metric(iref,1)/=static_cast<double>(NptsValid);
	  error_metric(iref,2)=std::sqrt(error_metric(iref,2)/static_cast<double>(NptsValid));
	  printf(" Krig_MAE=%12.6g Krig_RMSE=%12.6g",error_metric(iref,1),error_metric(iref,2));
	  fprintf(fpout1," Krig_MAE=%12.6g Krig_RMSE=%12.6g",error_metric(iref,1),error_metric(iref,2));
	}
	else{
	  printf(" Krig_MAE=NaN          Krig_RMSE=NaN         ");
	  fprintf(fpout1," Krig_MAE=NaN          Krig_RMSE=NaN         ");
	}
	fflush(fpout1);
	if(iref<Nref) { //Npts*(1+Ndim) equations makes for a BIG correlation matrix 
	  //(slow emulator construction) and I don't need the largest Npts for the
	  //Gradient Enhanced Kriging Paper
	  nkm::KrigingModel gkm(sd_build_temp,herbie_GEK_params); gkm.create();
	  gkm.evaluate(yeval,sd_valid_herbie.xr);
	  for(int i=0; i<NptsValid; ++i) {
	    double tmpdbl=yeval(i,0)-sd_valid_herbie.y(i,0);
	    error_metric(iref,3)+=std::fabs(tmpdbl);
	    error_metric(iref,4)+=std::pow(tmpdbl,2.0);
	  }
	  error_metric(iref,3)/=static_cast<double>(NptsValid);	
	  error_metric(iref,4)=std::sqrt(error_metric(iref,4)/static_cast<double>(NptsValid));	
	  printf(" GEK_MAE=%12.6g GEK_RMSE=%12.6g %5d/%-5d\n",
		 error_metric(iref,3),error_metric(iref,4),
		 gkm.getNumEqnKeep(),gkm.getNumEqnAvail());
	  fprintf(fpout1," GEK_MAE=%12.6g GEK_RMSE=%12.6g %5d/%-5d\n",
		  error_metric(iref,3),error_metric(iref,4),
		  gkm.getNumEqnKeep(),gkm.getNumEqnAvail());
	}
	else{
	  printf(" GEK_MAE=NaN          GEK_RMSE=NaN            NaN/NaN  \n");
	  fprintf(fpout1," GEK_MAE=NaN          GEK_RMSE=NaN            NaN/NaN  \n");
	}	
	fflush(fpout1);
      } //end herbie scope

      {//limit km and gkm for SMOOTH herbie to this scope
	sd_build_smooth.getPoints(sd_build_temp,ipts);
	printf("Smooth: Ndim=%d Npts=%4d/%-4d",Ndim,NptsThis,NptsBuild);
	fprintf(fpout2,"Smooth: Ndim=%d Npts=%4d/%-4d",Ndim,NptsThis,NptsBuild);
	
	if(iref>0) { //if have enough equations for a reduced quadratic trend
	  nkm::KrigingModel km(sd_build_temp,smooth_krig_params); km.create();
	  km.evaluate(yeval,sd_valid_smooth.xr);
	  for(int i=0; i<NptsValid; ++i) {
	    double tmpdbl=yeval(i,0)-sd_valid_smooth.y(i,0);
	    error_metric(iref,5)+=std::fabs(tmpdbl);
	    error_metric(iref,6)+=std::pow(tmpdbl,2.0);
	  }
	  error_metric(iref,5)/=static_cast<double>(NptsValid);
	  error_metric(iref,6)=std::sqrt(error_metric(iref,6)/static_cast<double>(NptsValid));
	  printf(" Krig_MAE=%12.6g Krig_RMSE=%12.6g",error_metric(iref,5),error_metric(iref,6));
	  fprintf(fpout2," Krig_MAE=%12.6g Krig_RMSE=%12.6g",error_metric(iref,5),error_metric(iref,6));
	}
	else{
	  printf(" Krig_MAE=NaN          Krig_RMSE=NaN         ");
	  fprintf(fpout2," Krig_MAE=NaN          Krig_RMSE=NaN         ");
	}
	fflush(fpout2);	
	if(iref<Nref) { //Npts*(1+Ndim) equations makes for a BIG correlation matrix 
	  //(slow emulator construction) and I don't need the largest Npts for the
	  //Gradient Enhanced Kriging Paper

	  nkm::KrigingModel gkm(sd_build_temp,smooth_GEK_params); gkm.create();
	  gkm.evaluate(yeval,sd_valid_smooth.xr);
	  for(int i=0; i<NptsValid; ++i) {
	    double tmpdbl=yeval(i,0)-sd_valid_smooth.y(i,0);
	    error_metric(iref,7)+=std::fabs(tmpdbl);
	    error_metric(iref,8)+=std::pow(tmpdbl,2.0);
	  }
	  error_metric(iref,7)/=static_cast<double>(NptsValid);	
	  error_metric(iref,8)=std::sqrt(error_metric(iref,8)/static_cast<double>(NptsValid));	
	  printf(" GEK_MAE=%12.6g GEK_RMSE=%12.6g %5d/%-5d\n",
		 error_metric(iref,7),error_metric(iref,8),
		 gkm.getNumEqnKeep(),gkm.getNumEqnAvail());
	  fprintf(fpout2," GEK_MAE=%12.6g GEK_RMSE=%12.6g %5d/%-5d\n",
		  error_metric(iref,7),error_metric(iref,8),
		  gkm.getNumEqnKeep(),gkm.getNumEqnAvail());
	}
	else{
	  printf(" GEK_MAE=NaN          GEK_RMSE=NaN            NaN/NaN  \n");
	  fprintf(fpout2," GEK_MAE=NaN          GEK_RMSE=NaN            NaN/NaN  \n");
	}	
	fflush(fpout2);
      } //end SMOOTH herbie scope
    } //for(int iref=0; iref<=Nref; ++iref)
    
  } //for(int ndimpow=1; ndimpow<3; ++ndimpow)

  fclose(fpout1);
  fclose(fpout2);
  return;
} //end of function

void compare_sample_designs() {
  
  string build_filename ="build_file.spd";
  nkm::SurfData sd2dbuild(build_filename , 2, 0, 3, 0, 1, 0);
  string valid_filename ="valid_file.spd";
  nkm::SurfData sd2dvalid(valid_filename , 2, 0, 3, 0, 1, 0);
  FILE* fpout=fopen("compare_out.txt","w");
  
  nkm::MtxDbl yeval(16384,1);
  int iout; //the 0th output column is Rosenbrock    
  double rmse;


  std::map< std::string, std::string> km_params;
  km_params["lower_bounds"]="-2.0 -2.0";
  km_params["upper_bounds"]="2.0 2.0";

  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);

  iout=0;
  sd2dbuild.setIOut(iout);
  sd2dvalid.setIOut(iout);
  nkm::KrigingModel kmros( sd2dbuild, km_params); kmros.create();

  //evaluate error the rosenbrock kriging model at 2^14=16384 validation points
  kmros.evaluate(yeval,sd2dvalid.xr);
  rmse=0.0;
  for(int i=0; i<16384; ++i)
    rmse+=std::pow(yeval(i,0)-sd2dvalid.y(iout,i),2);
  rmse=std::sqrt(rmse/16384.0);
  fprintf(fpout,"%22.16g\n",rmse);


  iout=1;
  sd2dbuild.setIOut(iout);
  sd2dvalid.setIOut(iout);
  nkm::KrigingModel kmshu( sd2dbuild, km_params); kmshu.create();

  //evaluate error the shubert kriging model at 2^14=16384 validation points
  kmshu.evaluate(yeval,sd2dvalid.xr);
  rmse=0.0;
  for(int i=0; i<16384; ++i)
    rmse+=std::pow(yeval(i,0)-sd2dvalid.y(iout,i),2);
  rmse=std::sqrt(rmse/16384.0);
  fprintf(fpout,"%22.16g\n",rmse);


  iout=2;
  sd2dbuild.setIOut(iout);
  sd2dvalid.setIOut(iout);
  nkm::KrigingModel kmherb( sd2dbuild, km_params); kmherb.create();

  //evaluate error the herbie kriging model at 2^14=16384 validation points
  kmherb.evaluate(yeval,sd2dvalid.xr);
  rmse=0.0;
  for(int i=0; i<16384; ++i)
    rmse+=std::pow(yeval(i,0)-sd2dvalid.y(iout,i),2);
  rmse=std::sqrt(rmse/16384.0);
  fprintf(fpout,"%22.16g\n",rmse);

  fclose(fpout);


  return;
}

void compare_sample_designs_pav(int nvarsr) {
  
  string build_filename ="build_file.spd";
  nkm::SurfData sd_pav_build(build_filename , nvarsr, 0, 1, 0, 1, 0);
  string valid_filename ="valid_file.spd";
  nkm::SurfData sd_pav_valid(valid_filename , nvarsr, 0, 1, 0, 1, 0);
  FILE* fpout=fopen("compare_out.txt","w");

  nkm::MtxDbl yeval(16384,1);
  double rmse;

  std::map< std::string, std::string> km_params;
  {  
    ostringstream os;  
    os << "2.0";
    for(int i=1; i<nvarsr; ++i)
      os << " 2.0";
    km_params["lower_bounds"]=os.str();
  }
  {  
    ostringstream os;  
    os << "10.0";
    for(int i=1; i<nvarsr; ++i)
      os << " 10.0";
    km_params["upper_bounds"]=os.str();
  }

  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  
  nkm::KrigingModel km( sd_pav_build, km_params); km.create();

  //evaluate error the rosenbrock kriging model at 2^14=16384 validation points
  km.evaluate(yeval,sd_pav_valid.xr);
  rmse=0.0;
  for(int i=0; i<16384; ++i)
    rmse+=std::pow(yeval(i,0)-sd_pav_valid.y(i,0),2);
  rmse=std::sqrt(rmse/16384.0);
  fprintf(fpout,"%22.16g\n",rmse);


  fclose(fpout);


  return;
}
#endif

std::string mtxdbl_2_string(nkm::MtxDbl& md) { 
  std::ostringstream oss;
  oss.precision(16);
  oss << md(0,0);
  for(int i=1; i<md.getNRows(); ++i)
    oss << " " << md(i,0);
  for(int j=1; j<md.getNCols(); ++j)
    for(int i=0; i<md.getNElems(); ++i)
      oss << " " << md(i,j);
  return oss.str();
}


void generic_build_and_eval_mean_adjvar(int Nvarsr, 
					std::string& model_type,
					std::string& optimization_method,
					std::string& corr_func_family,
					std::string& corr_func_param,
					std::string& handle_ill_cond1,
					std::string& handle_ill_cond2) {
  std::map< std::string, std::string> km_params;
  km_params[corr_func_family]=corr_func_param;
  km_params["optimization_method"]=optimization_method;
  int der_order;
  //make the model type uppercase
  std::transform(model_type.begin(), model_type.end(), model_type.begin(), ::toupper);
  if(model_type=="GEK") {    
    km_params["derivative_order"]="1";
    der_order=1;
  }
  else{
    km_params["derivative_order"]="0";
    der_order=0;
  }

  //km_params["order"] = "0"; //polynomial trend order = 0 is most robust which
  //is good for pathelogical problems
  
  //pivot_chol flag isn't used so it will be ignored.
  km_params[handle_ill_cond1]=handle_ill_cond2;
    
  if((optimization_method=="local")||
     (optimization_method=="none")) {
    //*local is local starting from a set of specified correlation lengths
    //*none is use the specified set of correlation lengths without 
    // optimization, i.e. to reproduce a model
    std::ifstream infile("corrlen.txt",ios::in);
    std::string corrlen_str;
    getline(infile,corrlen_str);
    km_params["correlation_lengths"]=corrlen_str;
  }

  string buildfile="buildfile.spd";
  string evalfile_in ="evalfile_in.spd";
  nkm::SurfData sdbuild(buildfile, Nvarsr, 0, 1, 0, der_order, 0);
  nkm::SurfData sdeval(evalfile_in, Nvarsr, 0, 1, 0, 0, 0);  
  int Neval=sdeval.getNPts();
  nkm::KrigingModel km(sdbuild,km_params); km.create();
  nkm::MtxDbl y(1,Neval);
  nkm::MtxDbl vary(1,Neval);
  km.evaluate(y,sdeval.xr);
  km.eval_variance(vary,sdeval.xr);
  FILE* fp=fopen("evalfile.out","w");
  for(int ipt=0; ipt<Neval; ++ipt) {
    for(int ivar=0; ivar<Nvarsr; ++ivar)
      fprintf(fp,"%22.16g ",sdeval.xr(ivar,ipt));
    fprintf(fp,"%22.16g %22.16g\n",y(0,ipt),vary(0,ipt));
  }
  fclose(fp);
  return;
}




//used for quick develop/test of matlab implementation of GPAIS
void GPAIS_build_and_eval_mean_adjvar(int Nvarsr, 
				      std::string& model_type,
				      std::string& optimization_method,
				      std::string& corr_func_family,
				      std::string& corr_func_param) {
  std::map< std::string, std::string> km_params;
  km_params[corr_func_family]=corr_func_param;
  km_params["optimization_method"]=optimization_method;
  int der_order;
  //make the model type uppercase
  std::transform(model_type.begin(), model_type.end(), model_type.begin(), ::toupper);
  if(model_type=="GEK") {    
    km_params["derivative_order"]="1";
    der_order=1;
  }
  else{
    km_params["derivative_order"]="0";
    der_order=0;
  }

  km_params["order"] = "0"; //polynomial trend order = 0 is most robust which
  //is good for pathelogical problems
  km_params["find_nugget"]="1"; //use a nugget rather than pivoted Cholesky
  //to handle ill conditioning, this is critical for GPAIS to work well

  if((optimization_method=="local")||
     (optimization_method=="none")) {
    //*local is local starting from a set of specified correlation lengths
    //*none is use the specified set of correlation lengths without 
    // optimization, i.e. to reproduce a model
    std::ifstream infile("corrlen.txt",ios::in);
    std::string corrlen_str;
    getline(infile,corrlen_str);
    km_params["correlation_lengths"]=corrlen_str;
  }

  string buildfile="buildfile.spd";
  string evalfile_in ="evalfile_in.spd";
  nkm::SurfData sdbuild(buildfile, Nvarsr, 0, 1, 0, der_order, 0);
  nkm::SurfData sdeval(evalfile_in, Nvarsr, 0, 1, 0, 0, 0);  
  int Neval=sdeval.getNPts();
  nkm::KrigingModel km(sdbuild,km_params); km.create();
  nkm::MtxDbl y(1,Neval);
  nkm::MtxDbl vary(1,Neval);
  km.evaluate(y,sdeval.xr);
  km.eval_variance(vary,sdeval.xr);
  FILE* fp=fopen("evalfile.out","w");
  for(int ipt=0; ipt<Neval; ++ipt) {
    for(int ivar=0; ivar<Nvarsr; ++ivar)
      fprintf(fp,"%22.16g ",sdeval.xr(ivar,ipt));
    fprintf(fp,"%22.16g %22.16g\n",y(0,ipt),vary(0,ipt));
  }
  fclose(fp);
  return;
}

void nightly_test()
{
  //filenames for 2D surfdata
  string validate2d_10 ="grad_validate2d_10.spd";
  string validate2d_100="grad_validate2d_100.spd";
  string validate2d_500="grad_validate2d_500.spd"; //GEK skips to save time
  string validate2d_10K="grad_validate2d_10K.spd";
  nkm::SurfData sd2d10( validate2d_10 , 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d500(validate2d_500, 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d100(validate2d_100, 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d10K(validate2d_10K, 2, 0, 3, 0, 1, 0);

  //filenames for 10D surfdata
  string paviani10d_50  ="grad_paviani10d_50.spd";
  string paviani10d_500 ="grad_paviani10d_500.spd"; //GEK skips to save time
  string paviani10d_10K ="grad_paviani10d_10K.spd";
  nkm::SurfData sdpav50( paviani10d_50 , 10, 0, 1, 0, 1, 0);
  nkm::SurfData sdpav500(paviani10d_500, 10, 0, 1, 0, 1, 0);
  nkm::SurfData sdpav10K(paviani10d_10K, 10, 0, 1, 0, 1, 0);

  nkm::MtxDbl yeval10( 1,   10); //only for 2D
  nkm::MtxDbl yeval50( 1,   50); //only for paviani
  nkm::MtxDbl yeval100(1,  100); //only for 2D
  nkm::MtxDbl yeval500(1,  500);
  nkm::MtxDbl yeval10K(1,10000);

  //for GEK will evaluate the derivative and compare them to what's
  //in the surfdata that GEK was built from, derivative reproduction
  //checking is only implemented for the 2D GEK test functions, it's 
  //not implemented for the 10D paviani test function
  nkm::MtxDbl d1y2d10(    2, 10);
  nkm::MtxDbl d1y2d100(   2,100);
  nkm::MtxDbl d1ysd2d10(  2, 10); //sd
  nkm::MtxDbl d1ysd2d100( 2,100); //sd

  std::map< std::string, std::string> km_params;

  std::cout << "********************************************************************************\n"
	    << "********************************************************************************\n"
	    << "Running \"nightly\" Kriging Model tests\n" 
	    << "********************************************************************************\n"
	    << "********************************************************************************\n"
	    << std::endl;
  { //the { is for scope we don't want km500, iout, error_stats, or the boost 
    //serialization variablesto exist afterwards
    //km_params["order"]="3";

    std::cout << "Running Kriging \"all defaults\" (no specified settings) "
	      << "test on 2D Rosenbrock with 500 points" << std::endl;
    int iout=0; //Rosenbrock
    sd2d500.setIOut(iout);
    sd2d10K.setIOut(iout);
    nkm::MtxDbl error_stats(3,4); error_stats.zero();
    nkm::KrigingModel km500(sd2d500, km_params); km500.create();
    //evaluate error the 500 pt kriging model at build points
    km500.evaluate(yeval500,sd2d500.xr);
    for(int j=0; j<500; ++j)
      error_stats(2,0)+=std::pow(yeval500(0,j)-sd2d500.y(iout,j),2);
    error_stats(2,1)=std::sqrt(error_stats(2,0)/500.0);
    
    //evaluate error the 500 pt kriging model at 10K new points
    km500.evaluate(yeval10K,sd2d10K.xr);
    for(int j=0; j<10000; ++j)
      error_stats(2,2)+=std::pow(yeval10K(0,j)-sd2d10K.y(iout,j),2);
    error_stats(2,3)=std::sqrt(error_stats(2,2)/10000.0);

    std::cout << "# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points"
	      << "\n" << setw(12) << 500
	      << ", " << setw(19) << setprecision(7) << error_stats(2,0)
	      << ", " << setw(20) << setprecision(7) << error_stats(2,1)
	      << ", " << setw(17) << setprecision(7) << error_stats(2,2)
	      << ", " << setw(18) << setprecision(7) 
	      << error_stats(2,3) << std::endl;

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
    std::cout << "Testing Kriging Save and Reload: ";
    nkm::MtxDbl yRL10K(1,10000); //compare reload eval to original eval
    nkm::SurfPackModel *kmOrigPtr=&km500;
    {
      std::ofstream nkm_km_ofstream("km.sav");
      
      boost::archive::text_oarchive output_archive(nkm_km_ofstream);
      output_archive << kmOrigPtr;
    }
    nkm::SurfPackModel *kmReload;
    {
      std::ifstream nkm_km_ifstream("km.sav");
      boost::archive::text_iarchive input_archive(nkm_km_ifstream);
      input_archive >> kmReload;
    }
    kmReload->evaluate(yRL10K,sd2d10K.xr);  
    delete kmReload; 

    bool ifdiff=false;
    for(int j=0; j<10000; ++j) {
      ifdiff=!(yeval10K(0,j)==yRL10K(0,j));
      if(ifdiff) break;
    }
    if(ifdiff==false)
      std::cout << "Passed" << std::endl;
    else
      std::cout << "Failed" << std::endl;
#endif
      std::cout << "********************************************************************************" << std::endl;
  } //end of default settings Kriging test
  { //the { is for scope we don't want der_order, oss, km100, iout, error_stats,
    //d1error_stats or the boost serialization variables to exist afterwards
    //km_params["order"]="3";


    std::cout << "Running Gradient Enhanced Kriging (GEK) \"all defaults\" "
	      << "(no specified settings\nother than derivative_order) test "
	      << "on 2D Rosenbrock with 100 points" << std::endl;
    
    int der_order=1;
    std::ostringstream oss;
    oss << der_order;
    km_params["derivative_order"]=oss.str();
    int iout=0; //Rosenbrock
    sd2d500.setIOut(iout);
    sd2d10K.setIOut(iout);
    nkm::MtxDbl error_stats(3,4); error_stats.zero();
    nkm::MtxDbl d1error_stats(2,4); d1error_stats.zero();
    nkm::KrigingModel km100(sd2d100, km_params); km100.create();
    //evaluate error the 100 pt kriging model at build points
    km100.evaluate(yeval100,sd2d100.xr);
    for(int j=0; j<100; ++j)
      error_stats(1,0)+=std::pow(yeval100(0,j)-sd2d100.y(iout,j),2);
    error_stats(1,1)=std::sqrt(error_stats(1,0)/100);
    
    //evaluate error the 100 pt kriging model at 10K points
    km100.evaluate(yeval10K,sd2d10K.xr);
    for(int j=0; j<10000; ++j)
      error_stats(1,2)+=std::pow(yeval10K(0,j)-sd2d10K.y(iout,j),2);
    error_stats(1,3)=std::sqrt(error_stats(1,2)/10000.0); 

    km100.evaluate_d1y(d1y2d100,sd2d100.xr);
    for(int j=0; j<100; ++j) {
      d1error_stats(1,0)+=std::pow(d1y2d100(0,j)-d1ysd2d100(0,j),2);
      d1error_stats(1,2)+=std::pow(d1y2d100(1,j)-d1ysd2d100(1,j),2);
    }
    d1error_stats(1,1)=std::sqrt(d1error_stats(1,0)/100.0);
    d1error_stats(1,3)=std::sqrt(d1error_stats(1,2)/100.0);
    std::cout << "# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points"
	      << "\n" << setw(12) << 100 
	      << ", " << setw(19) << setprecision(7) << error_stats(1,0)
	      << ", " << setw(20) << setprecision(7) << error_stats(1,1)
	      << ", " << setw(17) << setprecision(7) << error_stats(1,2)
	      << ", " << setw(18) << setprecision(7) << error_stats(1,3) 
	      << "\n" 
	      << "# of samples, d1y0 SSE at build points, d1y0 RMSE at build points, d1y1 SSE at build points, d1y1 RMSE at build points"
	      << "\n" << setw(12) << 100
	      << ", " << setw(24) << setprecision(7) 
	      << d1error_stats(1,0)
	      << ", " << setw(25) << setprecision(7)
	      << d1error_stats(1,1)
		      << ", " << setw(24) << setprecision(7)
	      << d1error_stats(1,2) 
	      << ", " << setw(25) << setprecision(7)
	      << d1error_stats(1,3) << std::endl;
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
    std::cout << "Testing GEK Save and Reload: ";
    nkm::MtxDbl yRL10K(1,10000); //compare reload eval to original eval
    nkm::MtxDbl d1yRL2d100( 2,100); //compare reload eval to orginal eval
    nkm::SurfPackModel *kmOrigPtr=&km100;
    {
      std::ofstream nkm_km_ofstream("km.sav");
      
      boost::archive::text_oarchive output_archive(nkm_km_ofstream);
      output_archive << kmOrigPtr;
    }
    nkm::SurfPackModel *kmReload;
    {
      std::ifstream nkm_km_ifstream("km.sav");
      boost::archive::text_iarchive input_archive(nkm_km_ifstream);
      input_archive >> kmReload;
    }
    kmReload->evaluate(yRL10K,sd2d10K.xr);  
    kmReload->evaluate_d1y(d1yRL2d100,sd2d100.xr);
    delete kmReload; 

    bool ifdiff=false;
    for(int j=0; j<10000; ++j) {
      ifdiff=!(yeval10K(0,j)==yRL10K(0,j));
      if(ifdiff) break;
    }
    if(ifdiff==false)
      for(int j=0; j<100; ++j) {
	ifdiff=!((d1y2d100(0,j)==d1yRL2d100(0,j))&&
		 (d1y2d100(1,j)==d1yRL2d100(1,j)));
	if(ifdiff) break;
      }
    if(ifdiff==false)
      std::cout << "Passed" << std::endl;
    else
      std::cout << "Failed" << std::endl;
#endif
      std::cout << "********************************************************************************" << std::endl;
  } //end of GEK default (not set) parameters test

  //now we want to test a whole bunch of different options

  //we're going to test 5 correlation functions
  //0 Gaussian (Matern infinity)
  //1 Matern 5/2
  //2 Matern 3/2
  //3 Exponential (Matern 1/2)
  //4 Powered Exponential with power=1.5;
  //the first three of these are available for GEK all 5 of these are 
  //available for Kriging
  nkm::MtxInt num_corr_func_tests(2,1);
  num_corr_func_tests(0,0)=5; //Kriging
  num_corr_func_tests(1,0)=3; //GEK

  std::vector<std::string> test_corr_func_family(5), test_corr_func_param(5), 
    model_name(2), output_name_2D(3), handle_ill_cond(3);
  test_corr_func_family[0]="matern"; test_corr_func_param[0]="infinity";
  test_corr_func_family[1]="matern"; test_corr_func_param[1]="2.5";
  test_corr_func_family[2]="matern"; test_corr_func_param[2]="1.5";
  test_corr_func_family[3]="matern"; test_corr_func_param[3]="0.5";
  test_corr_func_family[4]="powered_exponential"; 
  test_corr_func_param[4]="1.5";

  model_name[0]="Kriging";
  model_name[1]="Gradient Enhanced Kriging";

  output_name_2D[0]="Rosenbrock";
  output_name_2D[1]="Shubert";
  output_name_2D[2]="Herbie";

  handle_ill_cond[0]="via pivoted Cholesky";
  handle_ill_cond[1]="by adding a nugget found from rcond(R)";
  handle_ill_cond[2]="by adding a nugget where rcond(R)==0.0 is assumed";


  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  km_params["optimization_method"]="global_local";

  for(int der_order=0; der_order<=1; ++der_order) {
    { //for scope
      std::ostringstream oss;
      oss << der_order;
      km_params["derivative_order"]=oss.str();
    } //for scope
    //do 2D test functions
    km_params["lower_bounds"]="-2.0 -2.0";
    km_params["upper_bounds"]="2.0 2.0";
    //loop over test functions
    for(int iout=0; iout<3; ++iout) {
      //km_params["order"]="3";

      sd2d10.setIOut(iout);
      sd2d100.setIOut(iout);
      sd2d500.setIOut(iout);
      sd2d10K.setIOut(iout);
      if(der_order==1) {
	sd2d10.getDerY( d1ysd2d10 ,1,iout);
	sd2d100.getDerY(d1ysd2d100,1,iout);
      }

      //loop over correlation functions
      for(int icorrfunc=0; icorrfunc<num_corr_func_tests(der_order,0); 
	  ++icorrfunc) {
	
	km_params[test_corr_func_family[icorrfunc]]=
	  test_corr_func_param[icorrfunc];

	for(int icond=0; icond<3; ++icond) {
	  if(icond==1)
	    km_params["find_nugget"]="1";
	  else if(icond==2)
	    km_params["find_nugget"]="0";


	  
	  nkm::MtxDbl error_stats(3,4); error_stats.zero();
	  nkm::MtxDbl d1error_stats(2,4); d1error_stats.zero();
	  nkm::KrigingModel km10( sd2d10 , km_params);       
	  std::cout << "****************************************************************************\n"
		    << model_name[der_order] << " on the 2D " 
		    << output_name_2D[iout] << " test function\nUsing the " 
		    << km10.get_corr_func() << " correlation function\n" 
		    << "With correlation lengths found by the "
		    << km_params["optimization_method"] << " optimization method\n"
		    << "Handling ill-conditioning " << handle_ill_cond[icond] 
		    << "\n"
		    << "****************************************************************************"
		    << std::endl;
	  km10.create();
	  
	  //evaluate error the 10 pt kriging model at build points  
	  km10.evaluate(yeval10,sd2d10.xr);
	  for(int j=0; j<10; ++j)
	    error_stats(0,0)+=std::pow(yeval10(0,j)-sd2d10.y(iout,j),2);
	  error_stats(0,1)=std::sqrt(error_stats(0,0)/10.0);
	  
	  //evaluate error the 10 pt kriging model at 10K new points
	  km10.evaluate(yeval10K,sd2d10K.xr);
	  for(int j=0; j<10000; ++j)
	    error_stats(0,2)+=std::pow(yeval10K(0,j)-sd2d10K.y(iout,j),2);
	  error_stats(0,3)=std::sqrt(error_stats(0,2)/10000.0);

	  nkm::KrigingModel km100(sd2d100, km_params); km100.create();
	  //evaluate error the 100 pt kriging model at build points
	  km100.evaluate(yeval100,sd2d100.xr);
	  for(int j=0; j<100; ++j)
	    error_stats(1,0)+=std::pow(yeval100(0,j)-sd2d100.y(iout,j),2);
	  error_stats(1,1)=std::sqrt(error_stats(1,0)/100);
	  
	  //evaluate error the 100 pt kriging model at 10K points
	  km100.evaluate(yeval10K,sd2d10K.xr);
	  for(int j=0; j<10000; ++j)
	    error_stats(1,2)+=std::pow(yeval10K(0,j)-sd2d10K.y(iout,j),2);
	  error_stats(1,3)=std::sqrt(error_stats(1,2)/10000.0);
	  
	  if(der_order==0) {
	    nkm::KrigingModel km500(sd2d500, km_params); km500.create();
	    //evaluate error the 500 pt kriging model at build points
	    km500.evaluate(yeval500,sd2d500.xr);
	    for(int j=0; j<500; ++j)
	      error_stats(2,0)+=std::pow(yeval500(0,j)-sd2d500.y(iout,j),2);
	    error_stats(2,1)=std::sqrt(error_stats(2,0)/500.0);
	    
	    //evaluate error the 500 pt kriging model at 10K new points
	    km500.evaluate(yeval10K,sd2d10K.xr);
	    for(int j=0; j<10000; ++j)
	      error_stats(2,2)+=std::pow(yeval10K(0,j)-sd2d10K.y(iout,j),2);
	    error_stats(2,3)=std::sqrt(error_stats(2,2)/10000.0);
	  } else{
	    km10.evaluate_d1y(d1y2d10,sd2d10.xr);
	    for(int j=0; j<10; ++j) {
	      d1error_stats(0,0)+=std::pow(d1y2d10(0,j)-d1ysd2d10(0,j),2);
	      d1error_stats(0,2)+=std::pow(d1y2d10(1,j)-d1ysd2d10(1,j),2);
	    }
	    d1error_stats(0,1)=std::sqrt(d1error_stats(0,0)/10.0);
	    d1error_stats(0,3)=std::sqrt(d1error_stats(0,2)/10.0);
	    km100.evaluate_d1y(d1y2d100,sd2d100.xr);
	    for(int j=0; j<100; ++j) {
	      d1error_stats(1,0)+=std::pow(d1y2d100(0,j)-d1ysd2d100(0,j),2);
	      d1error_stats(1,2)+=std::pow(d1y2d100(1,j)-d1ysd2d100(1,j),2);
	    }
	    d1error_stats(1,1)=std::sqrt(d1error_stats(1,0)/100.0);
	    d1error_stats(1,3)=std::sqrt(d1error_stats(1,2)/100.0);
	  }
	  	    
	  std::cout << "# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points";
	  std::cout << "\n" << setw(12) << 10 
		    << ", " << setw(19) << setprecision(7) << error_stats(0,0)
		    << ", " << setw(20) << setprecision(7) << error_stats(0,1)
		    << ", " << setw(17) << setprecision(7) << error_stats(0,2)
		    << ", " << setw(18) << setprecision(7) << error_stats(0,3);
	  std::cout << "\n" << setw(12) << 100 
		    << ", " << setw(19) << setprecision(7) << error_stats(1,0)
		    << ", " << setw(20) << setprecision(7) << error_stats(1,1)
		    << ", " << setw(17) << setprecision(7) << error_stats(1,2)
		    << ", " << setw(18) << setprecision(7) << error_stats(1,3);
	  if(der_order==0)
	    std::cout << "\n" << setw(12) << 500
		      << ", " << setw(19) << setprecision(7) << error_stats(2,0)
		      << ", " << setw(20) << setprecision(7) << error_stats(2,1)
		      << ", " << setw(17) << setprecision(7) << error_stats(2,2)
		      << ", " << setw(18) << setprecision(7) 
		      << error_stats(2,3);
	  std::cout << std::endl;
	  if(der_order>=1) {
	    std::cout << "# of samples, d1y0 SSE at build points, d1y0 RMSE at build points, d1y1 SSE at build points, d1y1 RMSE at build points";
	    std::cout << "\n" << setw(12) << 10 
		      << ", " << setw(24) << setprecision(7) 
		      << d1error_stats(0,0)
		      << ", " << setw(25) << setprecision(7)
		      << d1error_stats(0,1)
		      << ", " << setw(24) << setprecision(7)
		      << d1error_stats(0,2) 
		      << ", " << setw(25) << setprecision(7)
		      << d1error_stats(0,3);
	    std::cout << "\n" << setw(12) << 100
		      << ", " << setw(24) << setprecision(7) 
		      << d1error_stats(1,0)
		      << ", " << setw(25) << setprecision(7)
		      << d1error_stats(1,1)
		      << ", " << setw(24) << setprecision(7)
		      << d1error_stats(1,2) 
		      << ", " << setw(25) << setprecision(7)
		      << d1error_stats(1,3) << std::endl;
	  }


	  //fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,
	  if((icond==1)||(icond==2))
	    km_params.erase("find_nugget");	  
	}//icond for 2D: Pivoted Cholesky, add a nugget
	
	//can't specify both Matern and Powered Exponential or will get an error
	km_params.erase(test_corr_func_family[icorrfunc]);
      }//icorrfunc
      
    } //iout: rosenbrock, shubert, herbie
    
    //now do the paviani 10D test function
    km_params["lower_bounds"]=
      " 2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0";
    km_params["upper_bounds"]=
      "10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0";
    
    for(int icorrfunc=0; icorrfunc<num_corr_func_tests(der_order,0); 
	++icorrfunc) {
      
      km_params[test_corr_func_family[icorrfunc]]=
	test_corr_func_param[icorrfunc];

      for(int icond=0; icond<3; ++icond) {
	if(icond==1)
	  km_params["find_nugget"]="1";
	else if(icond==2)
	  km_params["find_nugget"]="0";

	nkm::MtxDbl error_stats(2,4); error_stats.zero();
	nkm::KrigingModel km50( sdpav50 , km_params); 
	std::cout << "****************************************************************************\n"
		  << model_name[der_order] << " on the 10D Paviani" 
		  << " test function\nUsing the " 
		  << km50.get_corr_func() << " correlation function\n" 
		  << "With correlation lengths found by the "
		  << km_params["optimization_method"] << "optimization method\n"
		  << "Handling ill-conditioning " << handle_ill_cond[icond] 
		  << "\n"
		  << "****************************************************************************"
		  << std::endl;
	km50.create();
	
	//evaluate error the 50 pt kriging model at build points  
	km50.evaluate(yeval50,sdpav50.xr);
	for(int j=0; j<50; ++j)
	  error_stats(0,0)+=std::pow(yeval50(0,j)-sdpav50.y(0,j),2);
	error_stats(0,1)=std::sqrt(error_stats(0,0)/50.0);
	
	//evaluate error the 50 pt kriging model at 10K new points
	km50.evaluate(yeval10K,sdpav10K.xr);
	for(int j=0; j<10000; ++j)
	  error_stats(0,2)+=std::pow(yeval10K(0,j)-sdpav10K.y(0,j),2);
	error_stats(0,3)=std::sqrt(error_stats(0,2)/10000.0);
	
	if(der_order==0) {
	  nkm::KrigingModel km500(sdpav500, km_params); km500.create();
	  
	  //evaluate error the 500 pt kriging model at build points
	  km500.evaluate(yeval500,sdpav500.xr);
	  for(int j=0; j<500; ++j)
	    error_stats(1,0)+=std::pow(yeval500(0,j)-sdpav500.y(0,j),2);
	  error_stats(1,1)=std::sqrt(error_stats(1,0)/500.0);
	  
	  //evaluate error the 500 pt kriging model at 10K new points
	  km500.evaluate(yeval10K,sdpav10K.xr);
	  for(int j=0; j<10000; ++j)
	    error_stats(1,2)+=std::pow(yeval10K(0,j)-sdpav10K.y(0,j),2);
	  error_stats(1,3)=std::sqrt(error_stats(1,2)/10000.0);
	}
	
	std::cout << "# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points";
	std::cout << "\n" << setw(12) << 10 
		  << ", " << setw(19) << setprecision(7) << error_stats(0,0)
		  << ", " << setw(20) << setprecision(7) << error_stats(0,1)
		  << ", " << setw(17) << setprecision(7) << error_stats(0,2)
		  << ", " << setw(18) << setprecision(7) << error_stats(0,3);
	if(der_order==0)
	  std::cout << "\n" << setw(12) << 500
		    << ", " << setw(19) << setprecision(7) << error_stats(1,0)
		    << ", " << setw(20) << setprecision(7) << error_stats(1,1)
		    << ", " << setw(17) << setprecision(7) << error_stats(1,2)
		    << ", " << setw(18) << setprecision(7) << error_stats(1,3);
	std::cout << std::endl;
	
	if((icond==1)||(icond==2)) {
	  //if use nugget then it doesn't use pivoted Cholesky
	  km_params.erase("find_nugget");	  
	}
      } //icond for Paviani: pivoted Cholesky, add a nugget
      //can't specify both Matern and Powered Exponential or will get an error
      km_params.erase(test_corr_func_family[icorrfunc]);
    } //icorrfunc for Paviani
  } //der_order: Kriging, GEK

  sd2d10.clear();
  sd2d500.clear();
  sd2d100.clear();
  sd2d10K.clear();

  sdpav50.clear();
  sdpav500.clear();
  sdpav10K.clear();

  yeval10.clear();
  yeval50.clear();
  yeval100.clear();
  yeval500.clear();

  return;
}


