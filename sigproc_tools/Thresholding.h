/**
 * \file LocalThresholding.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class LocalThresholding
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_LOCALTHRESHOLDING_H__
#define __SIGPROC_TOOLS_LOCALTHRESHOLDING_H__

#include <vector>
#include <cmath>
#include <algorithm> 
#include <functional>
#include <iostream>

namespace sigproc_tools {

  /**
     \class Thresholding
     User defined class LocalThresholding ... these comments are used to generate
     doxygen documentation!
  */

  template <class T> using Array2D = std::vector<std::vector<T>>;

  class Thresholding{
    
    public:
      
      /// Default constructor
      Thresholding(){}

      void Niblack(const Array2D<float>& waveform2D,
                  Array2D<bool>& outputROI,
                  const float k,
                  const int sx,
                  const int sy) const;

      void Sauvola(const Array2D<float>& waveform2D,
                  Array2D<bool>& outputROI,
                  const float k,
                  const float R,
                  const int sx,
                  const int sy) const;

      void Otsu(const Array2D<float>& waveform2D,
                Array2D<bool>& outputBinary) const;

      void globalMean(const Array2D<float>& waveform2D,
                      Array2D<bool>& outputBinary,
                      float k=2.0) const;

      void computeIntegralImage(const Array2D<float>& waveform2D,
                                Array2D<float>& integral2D) const;
      
      /// Default destructor
      ~Thresholding(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

