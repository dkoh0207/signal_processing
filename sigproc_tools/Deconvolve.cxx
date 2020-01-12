#ifndef __SIGPROC_TOOLS_DECONVOLVE_CXX__
#define __SIGPROC_TOOLS_DECONVOLVE_CXX__

#include "Deconvolve.h"

#endif

void sigproc_tools::Deconvolve::filterLee(
  std::vector<std::vector<short>>& deconvolvedWaveform,
  const std::vector<std::vector<short>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy)
{
  filterLee<short>(
    deconvolvedWaveform, waveLessCoherent, sx, sy, noiseVar);
}

void sigproc_tools::Deconvolve::filterLee(
  std::vector<std::vector<float>>& deconvolvedWaveform,
  const std::vector<std::vector<float>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy)
{
  filterLee<float>(
    deconvolvedWaveform, waveLessCoherent, sx, sy, noiseVar);
}

void sigproc_tools::Deconvolve::filterLee(
  std::vector<std::vector<double>>& deconvolvedWaveform,
  const std::vector<std::vector<double>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy)
{
  filterLee<double>(
    deconvolvedWaveform, waveLessCoherent, sx, sy, noiseVar);
}

template <typename T>
void sigproc_tools::Deconvolve::filterLee(
  std::vector<std::vector<T>>& deconvolvedWaveform,
  const std::vector<std::vector<T>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy)
{
  size_t numChannels = waveLessCoherent.size();
  size_t nTicks = waveLessCoherent.at(0).size();
  int xHalfWindowSize(sx / 2);
  int yHalfWindowSize(sy / 2);

  deconvolvedWaveform.resize(numChannels);
  for (size_t i=0; i<numChannels; ++i) {
    deconvolvedWaveform[i].resize(nTicks);
  }

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, apply a adaptive local wiener filter.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      size_t lowerBoundx = std::max(lbx, 0);
      size_t upperBoundx = std::min(ubx, (int) numChannels);
      size_t lowerBoundy = std::max(lby, 0);
      size_t upperBoundy = std::min(uby, (int) nTicks);
      std::vector<T> x;
      x.reserve(sx * sy);
      std::vector<T> xsq;
      xsq.reserve(sx * sy);
      for (size_t ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (size_t iy=lowerBoundy; iy<upperBoundy; ++iy) {
          x.push_back(waveLessCoherent[ix][iy]);
          xsq.push_back(waveLessCoherent[ix][iy] * waveLessCoherent[ix][iy]);
        }
      }
      T localMean = std::accumulate(x.begin(), x.end(), 0.0) / x.size() ;
      T localSquare = std::accumulate(xsq.begin(), xsq.end(), 0.0) / x.size() ;
      T localVar = localSquare - localMean * localMean;
      if (noiseVar > localVar) {
        deconvolvedWaveform[i][j] = localMean;
      } else {
        deconvolvedWaveform[i][j] = localMean + (1 - noiseVar / localVar) * 
          (waveLessCoherent[i][j] - localMean);
      }
    }
  }
  return;
}


void sigproc_tools::Deconvolve::filterLeeMedian(
  std::vector<std::vector<short>>& deconvolvedWaveform,
  const std::vector<std::vector<short>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy)
{
  filterLeeMedian<short>(
    deconvolvedWaveform, waveLessCoherent, noiseVar, sx, sy);
}

void sigproc_tools::Deconvolve::filterLeeMedian(
  std::vector<std::vector<float>>& deconvolvedWaveform,
  const std::vector<std::vector<float>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy)
{
  filterLeeMedian<float>(
    deconvolvedWaveform, waveLessCoherent, noiseVar, sx, sy);
}

void sigproc_tools::Deconvolve::filterLeeMedian(
  std::vector<std::vector<double>>& deconvolvedWaveform,
  const std::vector<std::vector<double>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy)
{
  filterLeeMedian<double>(
    deconvolvedWaveform, waveLessCoherent, noiseVar, sx, sy);
}

template <typename T>
void sigproc_tools::Deconvolve::filterLeeMedian(
  std::vector<std::vector<T>>& deconvolvedWaveform,
  const std::vector<std::vector<T>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy)
{
  size_t numChannels = waveLessCoherent.size();
  size_t nTicks = waveLessCoherent.at(0).size();
  int xHalfWindowSize(sx / 2);
  int yHalfWindowSize(sy / 2);

  deconvolvedWaveform.resize(numChannels);
  for (size_t i=0; i<numChannels; ++i) {
    deconvolvedWaveform[i].resize(nTicks);
  }

  sigproc_tools::MiscUtils utils;

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, apply a adaptive local wiener filter.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      size_t lowerBoundx = std::max(lbx, 0);
      size_t upperBoundx = std::min(ubx, (int) numChannels);
      size_t lowerBoundy = std::max(lby, 0);
      size_t upperBoundy = std::min(uby, (int) nTicks);
      std::vector<T> x;
      x.reserve(sx * sy);
      std::vector<T> xsq;
      xsq.reserve(sx * sy);
      for (size_t ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (size_t iy=lowerBoundy; iy<upperBoundy; ++iy) {
          x.push_back(waveLessCoherent[ix][iy]);
          xsq.push_back(waveLessCoherent[ix][iy] * waveLessCoherent[ix][iy]);
        }
      }
      T localMean = std::accumulate(x.begin(), x.end(), 0.0) / x.size() ;
      T localSquare = std::accumulate(xsq.begin(), xsq.end(), 0.0) / x.size() ;
      T localVar = localSquare - localMean * localMean;
      T localMedian = utils.computeMedian(x);
      if (noiseVar > localVar) {
        deconvolvedWaveform[i][j] = localMedian;
      } else {
        deconvolvedWaveform[i][j] = localMedian + (1 - noiseVar / localVar) * 
          (waveLessCoherent[i][j] - localMedian);
      }
    }
  }
  return;
}


void sigproc_tools::Deconvolve::filterLeeEnhanced(
  std::vector<std::vector<short>>& deconvolvedWaveform,
  const std::vector<std::vector<short>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy,
  const float a,
  const float epsilon)
{
  filterLeeEnhanced<short>(
    deconvolvedWaveform, waveLessCoherent, sx, sy, noiseVar, a, epsilon);
}

void sigproc_tools::Deconvolve::filterLeeEnhanced(
  std::vector<std::vector<float>>& deconvolvedWaveform,
  const std::vector<std::vector<float>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy,
  const float a,
  const float epsilon)
{
  filterLeeEnhanced<float>(
    deconvolvedWaveform, waveLessCoherent, sx, sy, noiseVar, a, epsilon);
}

void sigproc_tools::Deconvolve::filterLeeEnhanced(
  std::vector<std::vector<double>>& deconvolvedWaveform,
  const std::vector<std::vector<double>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy,
  const float a,
  const float epsilon)
{
  filterLeeEnhanced<double>(
    deconvolvedWaveform, waveLessCoherent, sx, sy, noiseVar, a, epsilon);
}

template <typename T>
void sigproc_tools::Deconvolve::filterLeeEnhanced(
  std::vector<std::vector<T>>& deconvolvedWaveform,
  const std::vector<std::vector<T>>& waveLessCoherent,
  const float noiseVar,
  const unsigned int sx,
  const unsigned int sy,
  const float a,
  const float epsilon)
{
  size_t numChannels = waveLessCoherent.size();
  size_t nTicks = waveLessCoherent.at(0).size();
  int xHalfWindowSize(sx / 2);
  int yHalfWindowSize(sy / 2);

  deconvolvedWaveform.resize(numChannels);
  for (size_t i=0; i<numChannels; ++i) {
    deconvolvedWaveform[i].resize(nTicks);
  }

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, apply a adaptive local wiener filter.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      size_t lowerBoundx = std::max(lbx, 0);
      size_t upperBoundx = std::min(ubx, (int) numChannels);
      size_t lowerBoundy = std::max(lby, 0);
      size_t upperBoundy = std::min(uby, (int) nTicks);
      std::vector<T> x;
      x.reserve(sx * sy);
      std::vector<T> xsq;
      xsq.reserve(sx * sy);
      std::vector<float> weight;
      weight.reserve(sx * sy);
      for (size_t ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (size_t iy=lowerBoundy; iy<upperBoundy; ++iy) {
          float eps = std::pow(epsilon * std::sqrt(noiseVar), 2);
          float fsq = std::pow(waveLessCoherent[i][j] - waveLessCoherent[ix][iy], 2.0);
          float w = 1.0 / (1.0 + a * std::max(eps, fsq));
          x.push_back(waveLessCoherent[ix][iy]);
          xsq.push_back(waveLessCoherent[ix][iy] * waveLessCoherent[ix][iy]);
          weight.push_back(w);
        }
      }
      float normWeight = std::accumulate(weight.begin(), weight.end(), 0.0);
      for (auto& w : weight) {
        w = w / normWeight;
      }
      T localMean = std::inner_product(
        x.begin(), x.end(), weight.begin(), 0.0) / x.size();
      T localSquare = std::inner_product(
        xsq.begin(), xsq.end(), weight.begin(), 0.0) / x.size();
      T localVar = localSquare - localMean * localMean;
      if (noiseVar > localVar) {
        deconvolvedWaveform[i][j] = localMean;
      } else {
        deconvolvedWaveform[i][j] = localMean + (1 - noiseVar / localVar) * 
          (waveLessCoherent[i][j] - localMean);
      }
    }
  }
  return;
}


template <typename T>
void inverseFilter(
  std::vector<std::vector<T>>& deconvolvedWaveform,
  const std::vector<std::vector<T>>& waveLessCoherent,
  const std::vector<std::vector<T>>& responseFunction)
{
  return;
}