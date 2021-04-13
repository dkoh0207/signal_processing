#ifndef __SIGPROC_MULTITHREADING_FREQUENCYFILTERS_CXX__
#define __SIGPROC_MULTITHREADING_FREQUENCYFILTERS_CXX__

#include "FrequencyFilters.h"


// sigproc_multithreading::FrequencyFilters::FrequencyFilters(
//   const int numTimeSteps)
// {
//   fTimeVec.resize(numTimeSteps, 0.);
//   fFrequencyVec.resize(numTimeSteps, std::complex<float>(0., 0.));
  
//   fForwardPlan = fftwf_plan_dft_r2c_1d(
//     numTimeSteps, 
//     fTimeVec.data(), 
//     reinterpret_cast<fftwf_complex*>(fFrequencyVec.data()), 
//     FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

//   fInversePlan = fftwf_plan_dft_c2r_1d(
//     numTimeSteps, 
//     reinterpret_cast<fftwf_complex*>(fFrequencyVec.data()), 
//     fTimeVec.data(), 
//     FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

//   return;
// }


// sigproc_multithreading::FrequencyFilters::~FrequencyFilters()
// {
//   fftwf_destroy_plan(fForwardPlan);
//   fftwf_destroy_plan(fInversePlan);

//   fTimeVec.clear();
//   fFrequencyVec.clear();

//   return;
// }


#endif
