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
#include <assert.h>

namespace sigproc_tools {

  /**
     \class Morph1DFast
     User defined class Morph1DFast ... these comments are used to generate
     doxygen documentation!
  */

  template <class T> using Waveform = std::vector<T>;
  template <class T> using Array2D = std::vector<std::vector<T>>;

  class Morph1DFast{
    
    public:
      
      /// Default constructor
      Morph1DFast(){}

      void getDilation(const Waveform<bool>&,
                       const unsigned int,
                       Waveform<bool>&) const;

      void getDilation(const Waveform<short>&,
                       const unsigned int,
                       Waveform<short>&) const;

      void getDilation(const Waveform<float>&,
                       const unsigned int,
                       Waveform<float>&) const;

      void getDilation(const Waveform<double>&,
                       const unsigned int,
                       Waveform<double>&) const;

      // Column Major

      void getDilation(const Array2D<float>&,
                       const unsigned int,
                       Array2D<float>&,
                       const unsigned int) const;


      void getErosion(const Waveform<bool>&,
                      const unsigned int,
                      Waveform<bool>&) const;

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


      void getOpening(const Waveform<bool>&,
                      const unsigned int,
                      Waveform<bool>&) const;

      void getOpening(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getOpening(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getOpening(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;


      void getClosing(const Waveform<bool>&,
                      const unsigned int,
                      Waveform<bool>&) const;            

      void getClosing(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getClosing(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getClosing(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;

    
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
      void getOpening(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& openingVec) const;

      template <typename T>
      void getClosing(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& closingVec) const;

      template <typename T> 
      void getGradient(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& gradientVec) const;
    
      template <typename T> 
      void getDilation(
        const Array2D<T>& inputArray2D,
        const unsigned int structuringElementy,
        Array2D<T>& dilation2D,
        const unsigned int columnNum) const;

  };
}

#endif
/** @} */ // end of doxygen group 

