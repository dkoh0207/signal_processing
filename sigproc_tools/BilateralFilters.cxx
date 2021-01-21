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
  /*
  Brute force implementation of Bilateral Filtering
  WARNING: O(N^2) complexity, algorithm may be slow (~17s for 1050 x 4098)
  */
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

void sigproc_tools::BilateralFilters::Gaussian(
  const Array2D<float>& input2D,
  const float sigma,
  Array2D<float>& output2D) const
{
  Gaussian<float>(input2D, sigma, output2D);
  return;
}

void sigproc_tools::BilateralFilters::Gaussian(
  const Array2D<std::complex<float>>& input2D,
  const float sigma,
  Array2D<std::complex<float>>& output2D) const
{
  Gaussian<std::complex<float>>(input2D, sigma, output2D);
  return;
}

template <typename T>
void sigproc_tools::BilateralFilters::Gaussian(
  const Array2D<T>& input2D,
  const float sigma,
  Array2D<T>& output2D) const
{
  /*
  The constant time O(1) recursive gaussian filter
  See Young and Vliet, Recursive Implementation of the Gaussian Filter
  */
  assert(sigma >= 0.5);

  size_t numChannels = input2D.size();
  size_t numTicks = input2D.at(0).size();

  assert(output2D.size() == numChannels);
  assert(output2D.at(0).size() == numTicks);
  // Initialize Recursive Gaussian Parameters
  // Compute q-value
  float q = 0.0;
  if (sigma >= 2.5) {
    q = 0.98711 * sigma  - 0.96330;
  } 
  else {
    q = 3.97156-4.14554 * sqrt(1 - 0.26891 * sigma);
  }

  // Compute linear filter coefficients
  float b0 = 1.57825 + 2.44413*q + 1.4281*std::pow(q, 2.0) + 0.422205*std::pow(q, 3.0);
  float b1 = (2.44413*q) + (2.85619*std::pow(q, 2.0)) + (1.26661*std::pow(q, 3.0));
  float b2 = -(1.4281*std::pow(q, 2.0)) - (1.26661*std::pow(q, 3.0));
  float b3 = 0.422205*std::pow(q, 3.0);
  float B = 1 - (b1 + b2 + b3) / b0;

  std::vector<float> filterCoeffs = {b0, b1, b2, b3, B};
  
  // Gaussian filter is separable, so 1D gaussian applied row/column in seq.
  for (size_t i=0; i<numChannels; ++i) {
    Gaussian1D(input2D[i], filterCoeffs, output2D[i]);
  }
  for (size_t j=0; j<numTicks; ++j) {
    std::vector<T> columnIn(numChannels);
    std::vector<T> columnOut(numChannels);
    for (size_t i=0; i<numChannels; ++i) {
      columnIn[i] = input2D[i][j];
    }
    Gaussian1D(columnIn, filterCoeffs, columnOut);
    for (size_t i=0; i<numChannels; ++i) {
      output2D[i][j] = columnOut[i];
    }
  }

  return;
}


void sigproc_tools::BilateralFilters::Gaussian1D(
  const std::vector<float>& input1D,
  const std::vector<float>& filterCoeffs,
  std::vector<float>& output1D) const
{
  Gaussian1D<float>(input1D, filterCoeffs, output1D);
  return;
}

void sigproc_tools::BilateralFilters::Gaussian1D(
  const std::vector<std::complex<float>>& input1D,
  const std::vector<float>& filterCoeffs,
  std::vector<std::complex<float>>& output1D) const
{
  Gaussian1D<std::complex<float>>(input1D, filterCoeffs, output1D);
  return;
}

template <typename T>
void sigproc_tools::BilateralFilters::Gaussian1D(
  const std::vector<T>& input1D,
  const std::vector<float>& filterCoeffs,
  std::vector<T>& output1D) const
{
  assert(filterCoeffs.size() == 5);
  assert(output1D.size() == input1D.size());

  size_t N = input1D.size();

  float b0 = filterCoeffs[0];
  float b1 = filterCoeffs[1];
  float b2 = filterCoeffs[2];
  float b3 = filterCoeffs[3];
  float B = filterCoeffs[4];

  std::vector<T> w(N+6); // Padding left 3 + right 3
  w[0] = input1D[0];
  w[1] = input1D[0];
  w[2] = input1D[0];
  w[N] = input1D[0];
  w[N+1] = input1D[0];
  w[N+2] = input1D[0];

  std::vector<T> c(N+6); // Padding left 3 + right 3
  c[0] = input1D[0];
  c[1] = input1D[0];
  c[2] = input1D[0];
  c[N] = input1D[0];
  c[N+1] = input1D[0];
  c[N+2] = input1D[0];

  // Do forward pass
  for (size_t i=3; i<N+3; ++i) {
    w[i] = B*input1D[i] + (b1*w[i-1] + b2*w[i-2] + b3*w[i-3]) / b0;
  }
  // Do backward pass
  for (size_t i=N+2; i>2; --i) {
    c[i] = B*w[i] + (b1*c[i+1] + b2*c[i+2] + b3*c[i+3]) / b0;
  }
  
  for (size_t i=0; i<N; ++i) {
    output1D[i] = c[i+3];
  }
  return;
}


void sigproc_tools::BilateralFilters::Discretize(
  const Array2D<float>& input2D,
  const unsigned int maxPixelValue,
  Array2D<short>& output2D) const
{
  Discretize<short>(input2D, maxPixelValue, output2D);
  return;
}

void sigproc_tools::BilateralFilters::Discretize(
  const Array2D<float>& input2D,
  const unsigned int maxPixelValue,
  Array2D<float>& output2D) const
{
  Discretize<float>(input2D, maxPixelValue, output2D);
  return;
}

template <typename T>
void sigproc_tools::BilateralFilters::Discretize(
  const Array2D<float>& input2D,
  const unsigned int maxPixelValue,
  Array2D<T>& output2D) const
{

  size_t numChannels = input2D.size();
  size_t numTicks = input2D.at(0).size();
  assert(output2D.size() == numChannels);
  assert(output2D.at(0).size() == numTicks);
  assert(maxPixelValue > 0);

  sigproc_tools::MiscUtils utils;
  float lb = utils.computeMinimum(input2D);
  float ub = utils.computeMaximum(input2D);

  float m = ((float) maxPixelValue) / (ub - lb);

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<numTicks; ++j) {
      T newVal = (T) std::round(m * (input2D[i][j] - lb));
      output2D[i][j] = newVal;
    }
  }
  return;
}


void sigproc_tools::BilateralFilters::CSUBilateral(
  const Array2D<float>& input2D,
  Array2D<float>& output2D,
  const unsigned int T,
  const float sigma_s,
  const float sigma_r) const
{

  size_t numChannels = input2D.size();
  size_t numTicks = input2D.at(0).size();
  std::complex<float> I(0, 1);

  assert(output2D.size() == numChannels);
  assert(output2D.at(0).size() == numTicks);

  float gamma = M_PI / (2.0 * ((float) T));
  float rho = gamma * sigma_r;
  size_t N = 3;
  if (sigma_r < (1.0 / gamma)) {
    N = (size_t) std::ceil( std::pow(gamma * sigma_r, -2.0));
  }
  
  Array3D<std::complex<float>> H(N+1);
  Array3D<std::complex<float>> G(N+1);
  Array3D<std::complex<float>> D(N+1);

  Array3D<std::complex<float>> Hbar(N+1);
  Array3D<std::complex<float>> Gbar(N+1);

  for (size_t n=0; n<=N; ++n) {
    H[n].resize(numChannels);
    G[n].resize(numChannels);
    D[n].resize(numChannels);
    Hbar[n].resize(numChannels);
    Gbar[n].resize(numChannels);
    for (size_t i=0; i<numChannels; ++i) {
      H[n][i].resize(numTicks);
      G[n][i].resize(numTicks);
      D[n][i].resize(numTicks);
      Hbar[n][i].resize(numTicks);
      Gbar[n][i].resize(numTicks);
    }
  }

  sigproc_tools::MiscUtils utils;

  float w = gamma / (rho * sqrt((float) N));
  float c0 = std::pow(2.0, -(float)N);

  // Set up auxiliary images H, G, D;
  for (size_t n=0; n<=N; ++n) {
    float binom = (float) utils.nChoosek(N, n);
    for (size_t i=0; i<numChannels; ++i) {
      for (size_t j=0; j<numTicks; ++j) {
        float q = (float) (2 * (int) n - (int) N);
        std::complex<float> z = q * I * w * input2D[i][j];
        H[n][i][j] = std::exp(z);
        G[n][i][j] = H[n][i][j] * input2D[i][j];
        D[n][i][j] = c0 * binom * std::exp(-z);
      }
    }
  }

  // Run O(1) Gaussian Filtering
  for (size_t n=0; n<=N; ++n) {
    Gaussian(H[n], sigma_s, Hbar[n]);
    Gaussian(G[n], sigma_s, Gbar[n]);
  }

  // Compute output
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<numTicks; ++j) {
      std::complex<float> summedGbarX;
      std::complex<float> summedHbarX;
      for (size_t n=0; n<=N; ++n) {
        summedGbarX += D[n][i][j] * Gbar[n][i][j];
        summedHbarX += D[n][i][j] * Hbar[n][i][j];
      }
      output2D[i][j] = summedGbarX.real() / summedHbarX.real();
    }
  }
  return;

}


#endif
