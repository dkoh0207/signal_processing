#ifndef __SIGPROC_TOOLS_EDGEDETECTION_CXX__
#define __SIGPROC_TOOLS_EDGEDETECTION_CXX__

#include "EdgeDetection.h"


void sigproc_tools::EdgeDetection::Convolve2D(
  const Array2D<float>& input2D,
  Array2D<float>& output2D,
  const Array2D<float>& kernel) const
{
  // Input kernel must be normalized. 

  int numChannels = input2D.size();
  int numTicks = input2D.at(0).size();

  int kernelX = kernel.size() / 2;
  int kernelY = kernel.at(0).size() / 2;

  int kernelWidth = kernel.size();
  int kernelHeight = kernel.at(0).size();

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {

      for (int m=0; m<kernelWidth; ++m) {
        for (int n=0; n<kernelHeight; ++n) {
          
          int offsetX = m - kernelX;
          int offsetY = n - kernelY;

          int ix = i + offsetX;
          int iy = j + offsetY;

          if (0 <= ix && ix < numChannels && 0 <= iy && iy < numTicks) {
            output2D[i][j] += kernel[m][n] * input2D[ix][iy];
          }
        }
      }
    }
  }
  return;
}


void sigproc_tools::EdgeDetection::SobelX(
  const Array2D<float>& input2D,
  Array2D<float>& gradient) const
{
  Array2D<float> kernel(3);
  for (auto& v : kernel) {
    v.resize(3);
  }
  // Define SobelX
  kernel[0][0] = 1.0;
  kernel[0][1] = 0.0;
  kernel[0][2] = -1.0;
  kernel[1][0] = 2.0;
  kernel[1][1] = 0.0;
  kernel[1][2] = -2.0;
  kernel[2][0] = 1.0;
  kernel[2][1] = 0.0;
  kernel[2][2] = -1.0;

  Convolve2D(input2D, gradient, kernel);
  return;
}


void sigproc_tools::EdgeDetection::SobelY(
  const Array2D<float>& input2D,
  Array2D<float>& gradient) const
{
  Array2D<float> kernel(3);
  for (auto& v : kernel) {
    v.resize(3);
  }

  // Define SobelY
  kernel[0][0] = 1.0;
  kernel[0][1] = 2.0;
  kernel[0][2] = 1.0;
  kernel[1][0] = 0.0;
  kernel[1][1] = 0.0;
  kernel[1][2] = 0.0;
  kernel[2][0] = -1.0;
  kernel[2][1] = -2.0;
  kernel[2][2] = -1.0;

  Convolve2D(input2D, gradient, kernel);
  return;
}


void sigproc_tools::EdgeDetection::Sobel(
  const Array2D<float>& input2D,
  Array2D<float>& gradient,
  Array2D<float>& direction) const
{
  int numChannels = input2D.size();
  int numTicks = input2D.at(0).size();

  Array2D<float> sobelX;
  Array2D<float> sobelY;

  sobelX.resize(numChannels);
  for (auto& v : sobelX) {
    v.resize(numTicks);
  }

  sobelY.resize(numChannels);
  for (auto& v : sobelY) {
    v.resize(numTicks);
  }

  SobelX(input2D, sobelX);
  SobelY(input2D, sobelY);

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      float g = sqrt(sobelX[i][j]*sobelX[i][j] + sobelY[i][j]*sobelY[i][j]);
      gradient[i][j] = g;
      float gradDir = atan2(sobelY[i][j], sobelX[i][j]);
      direction[i][j] = gradDir * 180.0 / M_PI;
    }
  }
  return;
}


void sigproc_tools::EdgeDetection::EdgeNMS(
  const Array2D<float>& gradient2D,
  const Array2D<float>& degrees2D,
  Array2D<float>& output2D
) const
{
  int numChannels = gradient2D.size();
  int numTicks = gradient2D.at(0).size();

  // Compensate for boundary
  for (int i=1; i<numChannels-1; ++i) {
    for (int j=1; j<numTicks-1; ++j) {
      if ((112.5 <= degrees2D[i][j] && degrees2D[i][j] < 157.5)
        || (-67.5 <= degrees2D[i][j] && degrees2D[i][j] < -22.5)) {
        // Upper Diag Left - Lower Diag Right
        float v1 = gradient2D[i-1][j+1];
        float v2 = gradient2D[i+1][j-1];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2) {
          output2D[i][j] = gradient2D[i][j];
        } 
        else {
          output2D[i][j] = 0.0;
        }
      }
      else if ((22.5 <= degrees2D[i][j] && degrees2D[i][j] < 67.5) 
        || (-157.5 <= degrees2D[i][j] && degrees2D[i][j] < -112.5)) {
        // Upper Diag Right - Lower Diag Left
        float v1 = gradient2D[i+1][j+1];
        float v2 = gradient2D[i-1][j-1];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2) {
          output2D[i][j] = gradient2D[i][j];
        } 
        else {
          output2D[i][j] = 0.0;
        }
      }
      else if ((67.5 <= degrees2D[i][j] && degrees2D[i][j] < 112.5)
        || (-112.5 <= degrees2D[i][j] && degrees2D[i][j] < -67.5)) {
        // Up-Down
        float v1 = gradient2D[i][j-1];
        float v2 = gradient2D[i][j+1];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2) {
          output2D[i][j] = gradient2D[i][j];
        } 
        else {
          output2D[i][j] = 0.0;
        }
      }
      else{
        // Left-Right
        float v1 = gradient2D[i+1][j];
        float v2 = gradient2D[i-1][j];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2) {
          output2D[i][j] = gradient2D[i][j];
        } 
        else {
          output2D[i][j] = 0.0;
        }
      }
    }
  }

  return;
}


void sigproc_tools::EdgeDetection::EdgeNMSInterpolation(
  const Array2D<float>& gradient2D,
  const Array2D<float>& gradX,
  const Array2D<float>& gradY,
  const Array2D<float>& degrees2D,
  Array2D<float>& output2D) const
{
  int numChannels = gradient2D.size();
  int numTicks = gradient2D.at(0).size();

  // Implementation was partly adapted from:
  // https://github.com/JustinLiang/ComputerVisionProjects/blob/master/CannyEdgeDetector/CannyEdgeDetector.m

  // Compensate for boundary
  for (int i=1; i<numChannels-1; ++i) {
    for (int j=1; j<numTicks-1; ++j) {
      if ((0 <= degrees2D[i][j] && degrees2D[i][j] < 45) 
        || (-180 <= degrees2D[i][j] && -135 > degrees2D[i][j])) {
        float v1 = gradient2D[i+1][j+1];
        float v2 = gradient2D[i+1][j];
        float v3 = gradient2D[i-1][j];
        float v4 = gradient2D[i-1][j-1];
        float G = abs(gradY[i][j] / gradX[i][j]);
        if ((gradient2D[i][j] >= v1 + (v2 - v1)*G) && 
            (gradient2D[i][j] >= v4 + G*(v3 - v4))) {
          output2D[i][j] = G;
        } else {
          output2D[i][j] = 0;
        }
      }
      else if ((-45 <= degrees2D[i][j] && degrees2D[i][j] < 0) 
        || (135 <= degrees2D[i][j] && 180 >= degrees2D[i][j])) {
        float v1 = gradient2D[i+1][j-1];
        float v2 = gradient2D[i+1][j];
        float v3 = gradient2D[i-1][j];
        float v4 = gradient2D[i-1][j+1];
        float G = abs(gradY[i][j] / gradX[i][j]);
        if ((gradient2D[i][j] >= v1 + (v2 - v1)*G) && 
            (gradient2D[i][j] >= v4 + G*(v3 - v4))) {
          output2D[i][j] = G;
        } else {
          output2D[i][j] = 0;
        }
      }
      else if ((45 <= degrees2D[i][j] && degrees2D[i][j] < 90) 
        || (-135 <= degrees2D[i][j] && -90 > degrees2D[i][j])) {
        float v1 = gradient2D[i+1][j+1];
        float v2 = gradient2D[i][j+1];
        float v3 = gradient2D[i][j-1];
        float v4 = gradient2D[i-1][j-1];
        float G = abs(gradX[i][j] / gradY[i][j]);
        if ((gradient2D[i][j] >= v1 + (v2 - v1)*G) && 
            (gradient2D[i][j] >= v4 + G*(v3 - v4))) {
          output2D[i][j] = G;
        } else {
          output2D[i][j] = 0;
        }
      }
      else {
        float v1 = gradient2D[i-1][j+1];
        float v2 = gradient2D[i][j+1];
        float v3 = gradient2D[i][j-1];
        float v4 = gradient2D[i+1][j-1];
        float G = abs(gradX[i][j] / gradY[i][j]);
        if ((gradient2D[i][j] >= v1 + (v2 - v1)*G) && 
            (gradient2D[i][j] >= v4 + G*(v3 - v4))) {
          output2D[i][j] = G;
        } else {
          output2D[i][j] = 0;
        }
      }
    }
  }

  return;
}


void sigproc_tools::EdgeDetection::DoubleThresholding(
  const Array2D<float>& doneNMS2D,
  Array2D<bool>& binary2D,
  std::vector<int>& strongEdgeRows,
  std::vector<int>& strongEdgeCols,
  std::vector<int>& weakEdgeRows,
  std::vector<int>& weakEdgeCols,
  float lowThreshold,
  float highThreshold) const
{
  int numChannels = doneNMS2D.size();
  int numTicks = doneNMS2D.at(0).size();

  // Implementation was partly adapted from:
  // https://github.com/JustinLiang/ComputerVisionProjects/blob/master/CannyEdgeDetector/CannyEdgeDetector.m

  strongEdgeRows.reserve(numChannels * numTicks);
  strongEdgeCols.reserve(numChannels * numTicks);
  weakEdgeRows.reserve(numChannels * numTicks);
  weakEdgeCols.reserve(numChannels * numTicks);
  // Enumerator for strong / weak edges

  // Compensate for boundary
  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      if (doneNMS2D[i][j] >= highThreshold) {
        binary2D[i][j] = true;
        strongEdgeRows.push_back(i);
        strongEdgeCols.push_back(j);
      }
      else if (doneNMS2D[i][j] < highThreshold && doneNMS2D[i][j] >= lowThreshold) {
        weakEdgeRows.push_back(i);
        weakEdgeCols.push_back(j);
      }
      else {
        binary2D[i][j] = false;
      }
    }
  }
  return;
}


void sigproc_tools::EdgeDetection::HysteresisThresholding(
  const Array2D<bool>& binary2D,
  const std::vector<int>& strongEdgeRows,
  const std::vector<int>& strongEdgeCols,
  const std::vector<int>& weakEdgeRows,
  const std::vector<int>& weakEdgeCols,
  Array2D<bool>& output2D) const
{
  int numChannels = binary2D.size();
  int numTicks = binary2D.at(0).size();

  bool converged = false;

  sigproc_tools::Morph2DFast morph2d;

  // Initialize fixed point iteration buffer
  Array2D<bool> tempBuffer(numChannels);
  for (auto& v : tempBuffer) {
    v.resize(numTicks);
  }

  for (size_t i=0; i<strongEdgeRows.size(); ++i) {
    tempBuffer[strongEdgeRows[i]][strongEdgeCols[i]] = true;
  }

  assert(strongEdgeRows.size() == strongEdgeCols.size());
  assert(weakEdgeRows.size() == weakEdgeCols.size());

  // Construct strong edge / weak edge 2D binary arrays from indices
  // REMARK: binary2D already contains strong edges. 
  // Array2D<bool> strongEdges2D(numChannels);
  // for (auto& v : strongEdge2D) {
  //   v.resize(numTicks);
  // }
  Array2D<bool> weakEdges2D(numChannels);
  for (auto& v : weakEdges2D) {
    v.resize(numTicks);
  }

  for (size_t i=0; i<weakEdgeRows.size(); ++i) {
    weakEdges2D[weakEdgeRows[i]][weakEdgeCols[i]] = true;
  }

  // This fixed point iteration always converges. 
  while (!converged) {
    // Dilation + compute overlap
    morph2d.getDilation(tempBuffer, 3, 3, output2D);
    for (int i=0; i<numChannels; ++i) {
      for (int j=0; j<numTicks; ++j) {
        // Pixels have to be true in both dilated binary image and weak edges
        output2D[i][j] = (output2D[i][j] && weakEdges2D[i][j]);
      }
    }
    converged = true;
    // Check array equality
    for (int i=0; i<numChannels; ++i) {
      for (int j=0; j<numTicks; ++j) {
        if (tempBuffer[i][j] != output2D[i][j]) {
          converged = false;
          goto nextIter; // goto is used to avoid outer loop break
        }
      }
    }
    nextIter:
    // Update buffer
    for (int i=0; i<numChannels; ++i) {
      for (int j=0; j<numTicks; ++j) {
        tempBuffer[i][j] = output2D[i][j];
      }
    }
  }
  return;
}

#endif
