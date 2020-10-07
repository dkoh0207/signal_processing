/**
 * \file FrequencyFilters.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class FrequencyFilters
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_FREQUENCYFILTERS_H__
#define __SIGPROC_TOOLS_FREQUENCYFILTERS_H__

#include <vector>
#include <complex>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>
#include <iterator>

#include "fftw3.h"

namespace sigproc_tools {

  /**
     \class FrequencyFilters
     User defined class FrequencyFilters ... these comments are used to generate
     doxygen documentation!
  */

  template <class T> using Waveform = std::vector<T>;
  template <class T> using Waveform2D = std::vector<std::vector<T>>;

  class FrequencyFilters1D{
    
    public:

      using TimeVec = std::vector<float>;
      using FrequencyVec = std::vector<std::complex<float>>;

      /// Default constructor
      FrequencyFilters1D(const int numTimeSteps = 4096);

      // void lowPassHard1D(const Waveform<float>&,
      //                   const float,
      //                   Waveform<float>&);

      void lowPassButterWorth1D(const Waveform<float>&,
                                const size_t,
                                const unsigned int,
                                Waveform<float>&);

      void lowPassGaussian1D(const Waveform<float>&,
                             const size_t,
                             Waveform<float>&);

      void highPassHard1D(const Waveform<float>&,
                          const size_t,
                          Waveform<float>&);

      void highPassButterWorth1D(const Waveform<float>&,
                                 const size_t,
                                 const unsigned int,
                                 Waveform<float>&);

      void highPassGaussian1D(const Waveform<float>&,
                              const size_t,
                              Waveform<float>&);
                              
      void forwardFFT(TimeVec&, FrequencyVec&) const;
      void inverseFFT(FrequencyVec&, TimeVec&) const;

      void filterImage(const Waveform2D<float>&,
                       const size_t,
                       Waveform2D<float>&,
                       const unsigned int = 4,
                       const int mode = 0);
              
      /// Default destructor
      ~FrequencyFilters1D();
    
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

