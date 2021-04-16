/**
 * \file FrequencyFilters.h
 *
 * \ingroup sigproc_multithreading
 * 
 * \brief Class def header for a class FrequencyFilters
 *
 * @author koh0207
 */

/** \addtogroup sigproc_multithreading

    @{*/
#ifndef __SIGPROC_MULTITHREADING_FREQUENCYFILTERS_H__
#define __SIGPROC_MULTITHREADING_FREQUENCYFILTERS_H__

#include <complex>
#include "SigprocParallelDefs.h"
#include "fftw3.h"

namespace sigproc_multithreading {

  /**
     \class FrequencyFilters
     User defined class FrequencyFilters ... these comments are used to generate
     doxygen documentation!
  */

  class FrequencyFilters{
    
    public:

      using TimeVec = ConcurrentVector<float>;
      using FrequencyVec = ConcurrentVector<std::complex<float>>;
    
      /// Default constructor
      FrequencyFilters(const unsigned int numTimeSteps = 4096){}

      void highPassButterWorth1D(const ConcurrentVector<float>&,
                                 const size_t,
                                 const unsigned int,
                                 ConcurrentVector<float>&);
      
      /// Default destructor
      ~FrequencyFilters(){}
    
    private:

      TimeVec fTimeVec;
      FrequencyVec fFrequencyVec;

      float samplingRate = 2500000.0;

      fftwf_plan fForwardPlan;
      fftwf_plan fInversePlan;
    
  };
}

#endif
/** @} */ // end of doxygen group 

