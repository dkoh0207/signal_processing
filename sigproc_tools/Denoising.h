/**
 * \file Denoising.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class Denoising
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_DENOISING_H__
#define __SIGPROC_TOOLS_DENOISING_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "Morph1D.h"
#include "Morph2D.h"

namespace sigproc_tools {

  /**
     \class Denoising
     User defined class Denoising ... these comments are used to generate
     doxygen documentation!
  */
  class Denoising{
    
    public:
      
      /// Default constructor
      Denoising(){}


      void getSelectVals(
        const std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        std::vector<std::vector<bool>>&,
        const unsigned int,
        const float);

      void getSelectVals(
        const std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        std::vector<std::vector<bool>>&,
        const unsigned int,
        const float);

      void getSelectVals(
        const std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        std::vector<std::vector<bool>>&,
        const unsigned int,
        const float);


      void removeCoherentNoise1D(
        std::vector<std::vector<short> >&,
        const std::vector<std::vector<short> >&,
        std::vector<std::vector<short> >&,
        std::vector<std::vector<short> >&,
        std::vector<std::vector<bool> >&,
        std::vector<std::vector<short> >&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float 
      );

      void removeCoherentNoise1D(
        std::vector<std::vector<float> >&,
        const std::vector<std::vector<float> >&,
        std::vector<std::vector<float> >&,
        std::vector<std::vector<float> >&,
        std::vector<std::vector<bool> >&,
        std::vector<std::vector<float> >&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float 
      );

      void removeCoherentNoise1D(
        std::vector<std::vector<double> >&, 
        const std::vector<std::vector<double> >&,
        std::vector<std::vector<double> >&,
        std::vector<std::vector<double> >&,
        std::vector<std::vector<bool> >&,
        std::vector<std::vector<double> >&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float 
      );


      void removeCoherentNoise2D(
        std::vector<std::vector<short> >&,
        const std::vector<std::vector<short> >&, 
        std::vector<std::vector<short> >&,
        std::vector<std::vector<short> >&,
        std::vector<std::vector<bool> >&,
        std::vector<std::vector<short> >&,
        const char, 
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      void removeCoherentNoise2D(
        std::vector<std::vector<float> >&,
        const std::vector<std::vector<float> >&, 
        std::vector<std::vector<float> >&,
        std::vector<std::vector<float> >&,
        std::vector<std::vector<bool> >&,
        std::vector<std::vector<float> >&,
        const char, 
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      void removeCoherentNoise2D(
        std::vector<std::vector<double> >&, 
        const std::vector<std::vector<double> >&, 
        std::vector<std::vector<double> >&,
        std::vector<std::vector<double> >&,
        std::vector<std::vector<bool> >&,
        std::vector<std::vector<double> >&,
        const char, 
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

    
    /// Default destructor
    ~Denoising(){}

    private:

      template <typename T>
      void getSelectVals(
        const std::vector<std::vector<T>>& waveforms,
        const std::vector<std::vector<T>>& morphedWaveforms,
        std::vector<std::vector<bool>>& selectVals,
        const unsigned int window,
        const float thresholdFactor
      );


      template <typename T>
      void removeCoherentNoise1D(
        std::vector<std::vector<T> >& waveLessCoherent, 
        const std::vector<std::vector<T> >& filteredWaveforms, 
        std::vector<std::vector<T> >& morphedWaveforms, 
        std::vector<std::vector<T> >& intrinsicRMS,
        std::vector<std::vector<bool> >& selectVals,
        std::vector<std::vector<T> >& correctedMedians,
        const char filterName='d', 
        const unsigned int grouping=64,
        const unsigned int structuringElement=5,
        const unsigned int window=0,
        const float thresholdFactor=2.5);


      template <typename T>
      void removeCoherentNoise2D(
        std::vector<std::vector<T> >& waveLessCoherent, 
        const std::vector<std::vector<T> >& filteredWaveforms,
        std::vector<std::vector<T> >& morphedWaveforms, 
        std::vector<std::vector<T> >& intrinsicRMS,
        std::vector<std::vector<bool> >& selectVals,
        std::vector<std::vector<T> >& correctedMedians,
        const char filterName='g',
        const unsigned int grouping=64, 
        const unsigned int structuringElementx=5,
        const unsigned int structuringElementy=20,
        const unsigned int window=0,
        const float thresholdFactor=2.5);
    
  };
}

#endif
/** @} */ // end of doxygen group 

