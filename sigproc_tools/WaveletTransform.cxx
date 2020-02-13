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
    for (nn=n; nn>=minRes; nn>>=1) wavelet.transformInterval(wf, nn, isign);
  } else {
    for (nn=minRes; nn<=n; nn<<=1) wavelet.transformInterval(wf, nn, isign);
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
  const int levels) const
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
    sigproc_tools::WaveletTransform::wt1(transform[i], wavelet, levels, isign);
  }
  for (size_t i=0; i<nTicks; ++i) {
    std::vector<float> tempWf(numChannels);
    for (size_t j=0; j<numChannels; ++j) {
      tempWf[j] = transform[j][i];
    }
    sigproc_tools::WaveletTransform::wt1(tempWf, wavelet, levels, isign);
    for (size_t j=0; j<numChannels; ++j) {
      transform[j][i] = tempWf[j];
    }
  }
}

#endif
