#ifndef __SIGPROC_TOOLS_WAVELETTRANSFORM_CXX__
#define __SIGPROC_TOOLS_WAVELETTRANSFORM_CXX__

#include "WaveletTransform.h"

void sigproc_tools::WaveletTransform::wt1(
  std::vector<float>& wf,
  Daubechies4& wavelet,
  const int levels,
  const int isign) const
{
  int nn, n=wf.size();
  int minRes = n / std::pow(2, levels);
  if (n < 4) return;

  if (isign >= 0) {
    wavelet.condition(wf, n, 1);
    for (nn=n; nn>=minRes; nn>>=1) wavelet.transform(wf, nn, isign);
  } else {
    for (nn=minRes; nn<=n; nn<<=1) wavelet.transform(wf, nn, isign);
    wavelet.condition(wf, n, -1);
  }
}


void sigproc_tools::WaveletTransform::wt2(
  const std::vector<std::vector<float>>& inputWaveform,
  std::vector<std::vector<float>>& transform,
  const size_t numChannels,
  const size_t nTicks,
  const int isign,
  Daubechies4& wavelet,
  const int levelChannels,
  const int levelTicks) const
{
  size_t ntot = numChannels * nTicks;
  if (ntot & (ntot-1)) throw("All lengths must be powers of 2 in wtn");

  transform.resize(inputWaveform.size());
  for (size_t i=0; i<numChannels; ++i) {
    transform[i].resize(inputWaveform.at(0).size());
  }

  std::vector<float> tickWksp(nTicks);
  std::vector<float> channelWksp(numChannels);

  for (size_t i=0; i<numChannels; ++i) {
    transform[i] = inputWaveform[i];
    sigproc_tools::WaveletTransform::wt1(transform[i], wavelet, levelTicks, isign);
  }
  for (size_t i=0; i<nTicks; ++i) {
    std::vector<float> tempWf(numChannels);
    for (size_t j=0; j<numChannels; ++j) {
      tempWf[j] = transform[j][i];
    }
    sigproc_tools::WaveletTransform::wt1(tempWf, wavelet, levelChannels, isign);
    for (size_t j=0; j<numChannels; ++j) {
      transform[j][i] = tempWf[j];
    }
  }
}

void sigproc_tools::WaveletTransform::estimateLevelNoise(
  const std::vector<std::vector<float>>& transform,
  const int levels,
  std::vector<float>& noiseEstimates) const
{
  int numChannels = (int) transform.size();
  int nTicks = (int) transform.at(0).size();

  noiseEstimates.reserve(levels);
  // std::cout << "1 noiseEstimates.size(): " << noiseEstimates.size() << std::endl;

  sigproc_tools::MiscUtils utils;

  // Compute noise at level j by MAD of detail coefficients
  for (int l=1; l< levels+1; ++l) {
    // std::cout << "Level: " << l << std::endl;
    size_t detailChannels = numChannels >> l;
    size_t detailTicks = nTicks >> l;
    // std::cout << "detailChannels: " << detailChannels << std::endl;
    // std::cout << "detailTicks: " << detailTicks << std::endl;
    std::vector<float> noiseWksp;
    noiseWksp.reserve(detailChannels * detailTicks * 3);
    // std::cout << noiseWksp.size() << std::endl;
    // std::cout << detailChannels * detailTicks * 3 << std::endl;
    for (size_t i=detailChannels; i < detailChannels*2; ++i) {
      for (size_t j=detailTicks; j < detailTicks*2; ++j) {
        noiseWksp.push_back(transform[i][j]);
      }
    }
    for (size_t i=0; i < detailChannels; ++i) {
      for (size_t j=detailTicks; j < detailTicks*2; ++j) {
        noiseWksp.push_back(transform[i][j]);
      }
    }
    for (size_t i=detailChannels; i < detailChannels*2; ++i) {
      for (size_t j=0; j < detailTicks; ++j) {
        noiseWksp.push_back(transform[i][j]);
      }
    }
    // std::cout << noiseWksp.size() << std::endl;
    float noiseJ = 0.0;
    noiseJ = utils.estimateMAD(noiseWksp);
    // std::cout << noiseJ << std::endl;
    noiseEstimates.push_back(noiseJ);
    // std::cout << "2 noiseEstimates.size(): " << noiseEstimates.size() << std::endl;
  }
  // std::cout << "3 noiseEstimates.size(): " << noiseEstimates.size() << std::endl;
}

void sigproc_tools::WaveletTransform::WaveletWienerShrink(
  std::vector<std::vector<float>>& transform,
  const int levels,
  const std::vector<float>& noiseEstimates) const
{
  /*
  Wavelet Domain Wiener Filtering

  For a more rigorous implementation, one must use a different wavelet basis
  in estimating the noise variances. 
  */
  int numChannels = (int) transform.size();
  int nTicks = (int) transform.at(0).size();

  if ((size_t) levels != noiseEstimates.size()) {
    std::cout << "Number of levels " << levels << 
      " must agree with length of noise estimates " << 
        noiseEstimates.size() << std::endl;
  }

  sigproc_tools::MiscUtils utils;

  // Compute noise at level j by MAD of detail coefficients
  for (int l=1; l< levels+1; ++l) {
    float noiseJ = noiseEstimates[l-1];
    size_t detailChannels = numChannels >> l;
    size_t detailTicks = nTicks >> l;

    for (size_t i=detailChannels; i < detailChannels*2; ++i) {
      for (size_t j=detailTicks; j < detailTicks*2; ++j) {
        transform[i][j] = transform[i][j] * (std::pow(transform[i][j], 2.0) / 
          (std::pow(transform[i][j], 2.0) + std::pow(noiseJ, 2.0)));
      }
    }
    for (size_t i=0; i < detailChannels; ++i) {
      for (size_t j=detailTicks; j < detailTicks*2; ++j) {
        transform[i][j] = transform[i][j] * (std::pow(transform[i][j], 2.0) / 
          (std::pow(transform[i][j], 2.0) + std::pow(noiseJ, 2.0)));
      }
    }
    for (size_t i=detailChannels; i < detailChannels*2; ++i) {
      for (size_t j=0; j < detailTicks; ++j) {
        transform[i][j] = transform[i][j] * (std::pow(transform[i][j], 2.0) / 
          (std::pow(transform[i][j], 2.0) + std::pow(noiseJ, 2.0)));
      }
    }
  }
}



#endif
