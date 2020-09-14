/**
 * \file Morph1DFast.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class Morph1DFast
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_MORPH1DFAST_H__
#define __SIGPROC_TOOLS_MORPH1DFAST_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include <limits>

namespace sigproc_tools {

  /**
     \class Morph1DFast
     User defined class Morph1DFast ... these comments are used to generate
     doxygen documentation!
  */

  template <class T> using Waveform = std::vector<T>;

  class Morph1DFast{
    
    public:
      
      /// Default constructor
      Morph1DFast(){}

      void getDilation(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getDilation(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getDilation(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;


      void getErosion(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getErosion(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getErosion(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;


      void getGradient(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getGradient(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getGradient(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;

      void getDilationTest(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getDilationTest(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getDilationTest(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;

      // void getMedian(const Waveform<short>&,
      //                 const unsigned int,
      //                 Waveform<short>&) const;

      // void getMedian(const Waveform<float>&,
      //                 const unsigned int,
      //                 Waveform<float>&) const;

      // void getMedian(const Waveform<double>&,
      //                 const unsigned int,
      //                 Waveform<double>&) const;
    
    /// Default destructor
      ~Morph1DFast(){}

    private:

      template <typename T> 
      void getDilation(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& dilationVec) const;

      template <typename T> 
      void getErosion(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& erosionVec) const;

      template <typename T> 
      void getGradient(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& gradientVec) const;

      template <typename T> 
      void getDilationTest(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& dilationVec) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 

