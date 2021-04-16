#ifndef __SIGPROC_TOOLS_DENOISING_CXX__
#define __SIGPROC_TOOLS_DENOISING_CXX__

#include "Denoising.h"


void sigproc_tools::Denoising::getSelectVals(
  const std::vector<std::vector<short>>& waveforms,
  const std::vector<std::vector<short>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const unsigned int window,
  const float thresholdFactor)
{
  getSelectVals<short>(waveforms, morphedWaveforms, selectVals,
    roi, window, thresholdFactor);
}

void sigproc_tools::Denoising::getSelectVals(
  const std::vector<std::vector<float>>& waveforms,
  const std::vector<std::vector<float>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const unsigned int window,
  const float thresholdFactor)
{
  getSelectVals<float>(waveforms, morphedWaveforms, selectVals,
    roi, window, thresholdFactor);
}

void sigproc_tools::Denoising::getSelectVals(
  const std::vector<std::vector<double>>& waveforms,
  const std::vector<std::vector<double>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const unsigned int window,
  const float thresholdFactor)
{
  getSelectVals<double>(waveforms, morphedWaveforms, selectVals,
    roi, window, thresholdFactor);
}

template <typename T>
void sigproc_tools::Denoising::getSelectVals(
  const std::vector<std::vector<T>>& waveforms,
  const std::vector<std::vector<T>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const unsigned int window,
  const float thresholdFactor)
{
  auto numChannels = waveforms.size();
  auto nTicks = waveforms.at(0).size();

  for (size_t i=0; i<numChannels; ++i) {
    T median = 0.0;
    std::vector<T> localVec = morphedWaveforms[i];
    if (localVec.size() % 2 == 0) {
      const auto m1 = localVec.begin() + localVec.size() / 2 - 1;
      const auto m2 = localVec.begin() + localVec.size() / 2;
      std::nth_element(localVec.begin(), m1, localVec.end());
      const auto e1 = *m1;
      std::nth_element(localVec.begin(), m2, localVec.end());
      const auto e2 = *m2;
      median = (e1 + e2) / 2.0;
    } else {
      const auto m = localVec.begin() + localVec.size() / 2;
      std::nth_element(localVec.begin(), m, localVec.end());
      median = *m;
    }
    std::vector<T> baseVec;
    baseVec.resize(localVec.size());
    for (size_t j=0; j<baseVec.size(); ++j) {
      baseVec[j] = morphedWaveforms[i][j] - median;
    }

    float rms;
    rms = std::sqrt(std::inner_product(baseVec.begin(),
      baseVec.end(), baseVec.begin(), 0.) / float(baseVec.size()));
    float threshold;
    threshold = thresholdFactor * rms;

    for (size_t j=0; j<nTicks; ++j) {
      if (abs(morphedWaveforms[i][j]) > threshold) {
        // Check Bounds
        selectVals[i][j] = true;
        int lb = j - (int) window;
        int ub = j + (int) window + 1;
        size_t lowerBound = std::max(lb, 0);
        size_t upperBound = std::min(ub, (int) nTicks);
        for (size_t k=lowerBound; k<upperBound; ++k) {
          roi[i][k] = true;
        }
      } else {
        selectVals[i][j] = false;
      }
    }
  }
  return;
}


void sigproc_tools::Denoising::removeCoherentNoise1D(
  std::vector<std::vector<short>>& waveLessCoherent,
  const std::vector<std::vector<short>>& filteredWaveforms,
  std::vector<std::vector<short>>& morphedWaveforms,
  std::vector<std::vector<short>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  std::vector<std::vector<short>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise1D<short>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms,
    intrinsicRMS, selectVals, roi, correctedMedians,
    filterName, grouping, structuringElement, window, thresholdFactor);
  return;
}

void sigproc_tools::Denoising::removeCoherentNoise1D(
  std::vector<std::vector<float>>& waveLessCoherent,
  const std::vector<std::vector<float>>& filteredWaveforms,
  std::vector<std::vector<float>>& morphedWaveforms,
  std::vector<std::vector<float>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  std::vector<std::vector<float>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise1D<float>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms,
    intrinsicRMS, selectVals, roi, correctedMedians,
    filterName, grouping, structuringElement, window, thresholdFactor);
  return;
}

void sigproc_tools::Denoising::removeCoherentNoise1D(
  std::vector<std::vector<double>>& waveLessCoherent,
  const std::vector<std::vector<double>>& filteredWaveforms,
  std::vector<std::vector<double>>& morphedWaveforms,
  std::vector<std::vector<double>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  std::vector<std::vector<double>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise1D<double>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms,
    intrinsicRMS, selectVals, roi, correctedMedians,
    filterName, grouping, structuringElement, window, thresholdFactor);
  return;
}

template <typename T>
void sigproc_tools::Denoising::removeCoherentNoise1D(
  std::vector<std::vector<T>>& waveLessCoherent,
  const std::vector<std::vector<T>>& filteredWaveforms,
  std::vector<std::vector<T>>& morphedWaveforms,
  std::vector<std::vector<T>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  std::vector<std::vector<T>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  auto numChannels = filteredWaveforms.size();
  auto nTicks = filteredWaveforms.at(0).size();
  auto nGroups = numChannels / grouping;

  // Coherent noise subtracted denoised waveforms
  waveLessCoherent.resize(filteredWaveforms.size());
  for (auto& v : waveLessCoherent) {
    v.resize(filteredWaveforms.at(0).size());
  }

  // Waveform with morphological filter applied
  morphedWaveforms.resize(filteredWaveforms.size());
  for (auto& v : morphedWaveforms) {
    v.resize(filteredWaveforms.at(0).size());
  }

  // Regions to protect waveform from coherent noise subtraction.
  selectVals.resize(filteredWaveforms.size());
  for (auto& v : selectVals) {
    v.resize(filteredWaveforms.at(0).size());
  }

  roi.resize(filteredWaveforms.size());
  for (auto& v : roi) {
    v.resize(filteredWaveforms.at(0).size());
  }

  correctedMedians.resize(nGroups);
  for (auto& v : correctedMedians) {
    v.resize(nTicks);
  }

  intrinsicRMS.resize(nGroups);
  for (auto& v : intrinsicRMS) {
    v.resize(nTicks);
  }

  sigproc_tools::Morph1D denoiser;

  switch (filterName) {
    case 'd':
      for (size_t i=0; i<numChannels; ++i) {
        denoiser.getDilation(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    case 'e':
      for (size_t i=0; i<numChannels; ++i) {
        denoiser.getErosion(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    case 'a':
      for (size_t i=0; i<numChannels; ++i) {
        denoiser.getAverage(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    case 'g':
      for (size_t i=0; i<numChannels; ++i) {
        denoiser.getGradient(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    default:
      for (size_t i=0; i<numChannels; ++i) {
        denoiser.getDilation(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
  }

  getSelectVals(filteredWaveforms, morphedWaveforms,
    selectVals, roi, window, thresholdFactor);

  for (size_t i=0; i<nTicks; ++i) {
    for (size_t j=0; j<nGroups; ++j) {
      size_t group_start = j * grouping;
      size_t group_end = (j+1) * grouping;
      // Compute median.
      std::vector<T> v;
      for (size_t c=group_start; c<group_end; ++c) {
        if (!selectVals[c][i]) {
          v.push_back(filteredWaveforms[c][i]);
        }
      }
      T median = (T) 0;
      if (v.size() > 0) {
        if (v.size() % 2 == 0) {
          const auto m1 = v.begin() + v.size() / 2 - 1;
          const auto m2 = v.begin() + v.size() / 2;
          std::nth_element(v.begin(), m1, v.end());
          const auto e1 = *m1;
          std::nth_element(v.begin(), m2, v.end());
          const auto e2 = *m2;
          median = (e1 + e2) / 2.0;
        } else {
          const auto m = v.begin() + v.size() / 2;
          std::nth_element(v.begin(), m, v.end());
          median = *m;
        }
      }
      correctedMedians[j][i] = median;
      for (auto k=group_start; k<group_end; ++k) {
        if (!selectVals[k][i]) {
          waveLessCoherent[k][i] = filteredWaveforms[k][i] - median;
        } else {
          waveLessCoherent[k][i] = filteredWaveforms[k][i];
        }
      }
    }
  }

  T rms = (T) 0;
  for (size_t i=0; i<nGroups; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      std::vector<T> v;
      for (size_t k=i*grouping; k<(i+1)*grouping; ++k) {
        v.push_back(waveLessCoherent[k][j]);
      }
      rms = std::sqrt(
        std::inner_product(
          v.begin(), v.end(), v.begin(), 0.) / T(v.size()));
      intrinsicRMS[i][j] = (T) rms;
    }
  }
  return;
}


void sigproc_tools::Denoising::removeCoherentNoise2D(
  std::vector<std::vector<short>>& waveLessCoherent,
  const std::vector<std::vector<short>>& filteredWaveforms,
  std::vector<std::vector<short>>& morphedWaveforms,
  std::vector<std::vector<short>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  std::vector<std::vector<short>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise2D<short>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms,
    intrinsicRMS, selectVals, roi, correctedMedians,
    filterName, grouping, structuringElementx, structuringElementy,
    window, thresholdFactor);
  return;
}

void sigproc_tools::Denoising::removeCoherentNoise2D(
  std::vector<std::vector<float>>& waveLessCoherent,
  const std::vector<std::vector<float>>& filteredWaveforms,
  std::vector<std::vector<float>>& morphedWaveforms,
  std::vector<std::vector<float>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  std::vector<std::vector<float>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise2D<float>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms,
    intrinsicRMS, selectVals, roi, correctedMedians,
    filterName, grouping, structuringElementx, structuringElementy,
    window, thresholdFactor);
  return;
}

void sigproc_tools::Denoising::removeCoherentNoise2D(
  std::vector<std::vector<double>>& waveLessCoherent,
  const std::vector<std::vector<double>>& filteredWaveforms,
  std::vector<std::vector<double>>& morphedWaveforms,
  std::vector<std::vector<double>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  std::vector<std::vector<double>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise2D<double>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms,
    intrinsicRMS, selectVals, roi, correctedMedians,
    filterName, grouping, structuringElementx, structuringElementy,
    window, thresholdFactor);
  return;
}

template <typename T>
void sigproc_tools::Denoising::removeCoherentNoise2D(
  std::vector<std::vector<T>>& waveLessCoherent,
  const std::vector<std::vector<T>>& filteredWaveforms,
  std::vector<std::vector<T>>& morphedWaveforms,
  std::vector<std::vector<T>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  std::vector<std::vector<T>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor)
{
  auto numChannels = filteredWaveforms.size();
  auto nTicks = filteredWaveforms.at(0).size();
  auto nGroups = numChannels / grouping;

  // Coherent noise subtracted denoised waveforms
  waveLessCoherent.resize(filteredWaveforms.size());
  for (auto& v : waveLessCoherent) {
    v.resize(filteredWaveforms.at(0).size());
  }

  // Regions to protect waveform from coherent noise subtraction.
  selectVals.resize(filteredWaveforms.size());
  for (auto& v : selectVals) {
    v.resize(filteredWaveforms.at(0).size());
  }

  roi.resize(filteredWaveforms.size());
  for (auto& v : roi) {
    v.resize(filteredWaveforms.at(0).size());
  }

  correctedMedians.resize(nGroups);
  for (auto& v : correctedMedians) {
    v.resize(nTicks);
  }

  intrinsicRMS.resize(nGroups);
  for (auto& v : intrinsicRMS) {
    v.resize(nTicks);
  }

  sigproc_tools::Morph2D denoiser;

  std::vector<std::vector<T>> dilation;
  std::vector<std::vector<T>> erosion;
  std::vector<std::vector<T>> average;
  std::vector<std::vector<T>> gradient;

  denoiser.getFilter2D(filteredWaveforms, structuringElementx,
    structuringElementy, dilation, erosion, average, gradient);

  switch (filterName) {
    case 'd':
      getSelectVals(filteredWaveforms, dilation,
        selectVals, roi, window, thresholdFactor);
      morphedWaveforms = dilation;
      break;
    case 'e':
      getSelectVals(filteredWaveforms, erosion,
        selectVals, roi, window, thresholdFactor);
      morphedWaveforms = erosion;
      break;
    case 'a':
      getSelectVals(filteredWaveforms, average,
        selectVals, roi, window, thresholdFactor);
      morphedWaveforms = average;
      break;
    case 'g':
      getSelectVals(filteredWaveforms, gradient,
        selectVals, roi, window, thresholdFactor);
      morphedWaveforms = gradient;
      break;
    default:
      getSelectVals(filteredWaveforms, gradient,
        selectVals, roi, window, thresholdFactor);
      morphedWaveforms = gradient;
      break;
  }

  for (size_t i=0; i<nTicks; ++i) {
    for (size_t j=0; j<nGroups; ++j) {
      size_t group_start = j * grouping;
      size_t group_end = (j+1) * grouping;
      // Compute median.
      std::vector<T> v;
      for (size_t c=group_start; c<group_end; ++c) {
        if (!selectVals[c][i]) {
          v.push_back(filteredWaveforms[c][i]);
        }
      }
      T median = (T) 0.0;
      if (v.size() > 0) {
        if (v.size() % 2 == 0) {
          const auto m1 = v.begin() + v.size() / 2 - 1;
          const auto m2 = v.begin() + v.size() / 2;
          std::nth_element(v.begin(), m1, v.end());
          const auto e1 = *m1;
          std::nth_element(v.begin(), m2, v.end());
          const auto e2 = *m2;
          median = (e1 + e2) / 2.0;
        } else {
          const auto m = v.begin() + v.size() / 2;
          std::nth_element(v.begin(), m, v.end());
          median = *m;
        }
      }
      correctedMedians[j][i] = median;
      for (size_t k=group_start; k<group_end; ++k) {
        if (!selectVals[k][i]) {
          waveLessCoherent[k][i] = filteredWaveforms[k][i] - median;
        } else {
          waveLessCoherent[k][i] = filteredWaveforms[k][i];
        }
      }
    }
  }

  float rms = 0.0;
  for (size_t i=0; i<nGroups; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      std::vector<T> v;
      for (size_t k=i*grouping; k<(i+1)*grouping; ++k) {
        v.push_back(waveLessCoherent[k][j]);
      }
      rms = std::sqrt(
        std::inner_product(
          v.begin(), v.end(), v.begin(), 0.) / T(v.size()));
      intrinsicRMS[i][j] = (T) rms;
    }
  }
  return;
}


void sigproc_tools::Denoising::subtractPedestals(
  std::vector<std::vector<float>>& wf2d,
  const std::vector<std::vector<bool>>& selectVals)
{
  size_t numChannels = wf2d.size();
  size_t nTicks = wf2d.at(0).size();

  std::vector<float> medians;
  medians.reserve(numChannels);

  sigproc_tools::MiscUtils utils;

  for (size_t i=0; i<numChannels; ++i) {
    std::vector<float> v;
    v.reserve(nTicks);
    for (size_t j=0; j<nTicks; ++j) {
      if (!selectVals[i][j]) {
        v.push_back(wf2d[i][j]);
      }
    }
    float median = utils.computeMedian(v);
    medians.push_back(median);
  }

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      if (!selectVals[i][j]) {
        wf2d[i][j] = wf2d[i][j] - medians[i];
      }
    }
  }
}


#endif
