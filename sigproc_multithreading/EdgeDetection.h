/**
 * \file EdgeDetection.h
 *
 * \ingroup sigproc_multithreading
 * 
 * \brief Class def header for a class EdgeDetection
 *
 * @author koh0207
 */

/** \addtogroup sigproc_multithreading

    @{*/
#ifndef __SIGPROC_MULTITHREADING_EDGEDETECTION_H__
#define __SIGPROC_MULTITHREADING_EDGEDETECTION_H__

#include "SigprocParallelDefs.h"
#include "assert.h"

namespace sigproc_multithreading {

  /**
     \class EdgeDetection
     User defined class EdgeDetection ... these comments are used to generate
     doxygen documentation!
  */
  class EdgeDetection{
    
  public:
    
    /// Default constructor
    EdgeDetection(){}


    // Separable Implementation of Sobel Filtering

    void SobelX(const ConcurrentArray2D<float>& input2D,
                ConcurrentArray2D<float>& gradientX) const;

    void SobelXRow(const ConcurrentVector<float>& inputRow,
                   ConcurrentVector<float>& outputRow) const;

    void SobelXCol(const ConcurrentArray2D<float>& input2D,
                   ConcurrentArray2D<float>& outputRow) const;

    void SobelY(const ConcurrentArray2D<float>& input2D,
                ConcurrentArray2D<float>& gradientY) const;

    void SobelYRow(const ConcurrentVector<float>& inputRow,
                   ConcurrentVector<float>& outputRow) const;

    void SobelYCol(const ConcurrentArray2D<float>& input2D,
                   ConcurrentArray2D<float>& gradientY) const;

    void EdgeNMSInterpolation(const ConcurrentArray2D<float>& gradient2D,
                              const ConcurrentArray2D<float>& gradX,
                              const ConcurrentArray2D<float>& gradY,
                              const ConcurrentArray2D<float>& degrees2D,
                              ConcurrentArray2D<float>& output2D) const;

    void DoubleThresholding(
      const ConcurrentArray2D<float>& doneNMS2D,
      ConcurrentArray2D<bool>& binary2D,
      ConcurrentVector<int>& strongEdgeRows,
      ConcurrentVector<int>& strongEdgeCols,
      ConcurrentVector<int>& weakEdgeRows,
      ConcurrentVector<int>& weakEdgeCols,
      float lowThreshold,
      float highThreshold
    ) const;

    void HysteresisThresholding(
      const ConcurrentArray2D<bool>& binary2D,
      const ConcurrentVector<int>& strongEdgeRows,
      const ConcurrentVector<int>& strongEdgeCols,
      const ConcurrentVector<int>& weakEdgeRows,
      const ConcurrentVector<int>& weakEdgeCols,
      ConcurrentArray2D<bool>& output2D
    ) const;

    // Combine Double + Hysteresis for optimal performance
    void HysteresisThresholdingCombined(
      const ConcurrentArray2D<float>& doneNMS2D,
      ConcurrentArray2D<bool>& output2D,
      float lowThreshold,
      float highThreshold) const;

    void Canny(
      const ConcurrentArray2D<float>& waveLessCoherent,
      ConcurrentArray2D<bool>& output2D,
      const unsigned int sx,
      const unsigned int sy,
      const float sigma_x,
      const float sigma_y,
      const float sigma_r,
      const float lowThreshold,
      const float highThreshold,
      const char mode
    ) const;
    
    /// Default destructor
    ~EdgeDetection(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

