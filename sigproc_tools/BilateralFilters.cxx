#ifndef __SIGPROC_TOOLS_BILATERALFILTERS_CXX__
#define __SIGPROC_TOOLS_BILATERALFILTERS_CXX__

#include "BilateralFilters.h"


void sigproc_tools::BilateralFilters::bilateral(
  const Array2D<float>& input2D,
  Array2D<float>& output2D,
  const unsigned int sx,
  const unsigned int sy,
  float sigma_x,
  float sigma_y,
  float sigma_r) const
{

  int numChannels = input2D.size();
  int numTicks = input2D.at(0).size();
  int xHalfWindowSize(sx / 2);
  int yHalfWindowSize(sy / 2);

  output2D.resize(numChannels);
  for (auto& v : output2D) {
    v.resize(numTicks);
  }

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {

      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      int lowerBoundx = std::max(lbx, 0);
      int upperBoundx = std::min(ubx, (int) numChannels);
      int lowerBoundy = std::max(lby, 0);
      int upperBoundy = std::min(uby, (int) numTicks);

      float newPixel = 0.0;
      float normFactor = 0.0;

      for (int ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (int iy=lowerBoundy; iy<upperBoundy; ++iy) {
          float w = 1.0;
          w = exp(- (ix-i)*(ix-i) / (2.0*sigma_x*sigma_x) - (iy-j)*(iy-j) / (2.0*sigma_y*sigma_y));
          w *= exp(-(input2D[i][j] - input2D[ix][iy])*(input2D[i][j] - input2D[ix][iy]) / (2.0 * sigma_r * sigma_r));
          newPixel += w * input2D[ix][iy];
          normFactor += w;
        }
      }

      newPixel = newPixel / normFactor;
      output2D[i][j] = newPixel;
    }
  }
  return;
}


#endif
