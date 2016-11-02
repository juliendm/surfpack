#ifdef HAVE_CONFIG_H
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#endif
#include "surfpack_system_headers.h"

#include "surfpack.h"
#include "SurfData.h"
#include "SurfpackInterface.h"
#include "AxesBounds.h"
#include "KrigingModel.h"
#include "LinearRegressionModel.h"

#if !defined(HAVE_GETTIMEOFDAY) && (defined(_MSC_VER) || defined(__MINGW32__))
#include <windows.h>
#endif

using std::accumulate;
using std::cout;
using std::endl;
using std::setw;
using std::string;
using std::sort;
using std::vector;
using SurfpackInterface::CreateAxes;
using SurfpackInterface::CreateSample;
using SurfpackInterface::CreateSurface;
using SurfpackInterface::Save;

// Modified from http://mywebpage.netscape.com/yongweiwu/timeval.h.txt
#if !defined(HAVE_GETTIMEOFDAY) && (defined(_MSC_VER) || defined(__MINGW32__))
int gettimeofday (struct timeval *tv, void* tz)
{
  union {
    __int64 ns100; /*time since 1 Jan 1601 in 100ns units */
    FILETIME ft;
  } now;

  GetSystemTimeAsFileTime (&now.ft);
  tv->tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
  tv->tv_sec = (long) ((now.ns100 - 116444736000000000LL) / 10000000LL);
  return (0);
}
#endif

double time_difference(struct timeval& starttime, struct timeval& endtime)
{
  return (((double)endtime.tv_sec+(1.0e-06)*endtime.tv_usec) -
	((double)starttime.tv_sec+(1.0e-06)*starttime.tv_usec)) ;

}

int main(int argc, char** argv)
{
  AxesBounds* ab = CreateAxes(string("-2 2 | -2 2"));
  SurfpackModelFactory* smf = new KrigingModelFactory;
  vector<string> functions;
  vector<double> thetas(2,1.0);
  vector<unsigned> setsizes;
  vector<double> onepoint(1);
  const int num_trials = 1;
  vector<double> times_one_size(num_trials);
  vector<double> responses(2);
  double timeneeded;
  SurfData krigtimes;
  for (unsigned i = 30; i <= 50; i += 5) setsizes.push_back(i);
  setsizes.push_back(100);
  //for (unsigned i = 500; i < 1000; i += 25) setsizes.push_back(i);
  //for (unsigned i = 1000; i <= 2500; i += 100) setsizes.push_back(i);
  functions.push_back(string("rosenbrock"));
  struct timeval t1;
  struct timeval t2;
  smf->add("correlations","1.0 1.0");
  for (unsigned setsize = 0; setsize < setsizes.size(); setsize++) {
    cout << setw(8) << setsizes[setsize] ;
    for (unsigned trial = 0; trial < num_trials; trial++) {
      SurfData* sd = CreateSample(ab,setsizes[setsize]);
      SurfpackInterface::Evaluate(sd,functions);
      gettimeofday(&t1,NULL);
      SurfpackModel* km = smf->Build(*sd); 
      gettimeofday(&t2,NULL);
      timeneeded = time_difference(t1,t2); 
      times_one_size[trial] = timeneeded;
      cout <<  setw(15) << timeneeded;
      delete sd; sd = 0;
      delete km; km = 0;
    }
    sort(times_one_size.begin(),times_one_size.end());
    double avgtime = times_one_size[times_one_size.size()/2];
    double meantime = accumulate(times_one_size.begin(),times_one_size.end(),0.0)/times_one_size.size();
    onepoint[0] = setsizes[setsize];
    responses[0] = avgtime;
    responses[1] = meantime;
    krigtimes.addPoint(SurfPoint(onepoint,responses));
    cout << " avg: " << setw(15) 
         << avgtime << endl;
  }
  Save(&krigtimes,string("one_trial_3050to5000.spd"));
  SurfpackModelFactory* pmf = new LinearRegressionModelFactory;
  pmf->add("order","3");
  SurfpackModel* pm = pmf->Build(krigtimes); 

  Save(pm,string("poly3_krigtimes.sps"));
  // cleanup
  delete ab; ab = 0;
  delete smf; smf = 0;
  return 0;
}
