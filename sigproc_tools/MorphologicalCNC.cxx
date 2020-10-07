#ifndef __SIGPROC_TOOLS_MORPHOLOGICALCNC_CXX__
#define __SIGPROC_TOOLS_MORPHOLOGICALCNC_CXX__

#include "MorphologicalCNC.h"

void sigproc_tools::MorphologicalCNC::getSelectVals(
  const std::vector<std::vector<short>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const unsigned int window,
  const float thresholdFactor) const
{
  getSelectVals<short>(morphedWaveforms, selectVals,
    roi, window, thresholdFactor);
}

void sigproc_tools::MorphologicalCNC::getSelectVals(
  const std::vector<std::vector<float>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const unsigned int window,
  const float thresholdFactor) const
{
  getSelectVals<float>(morphedWaveforms, selectVals,
    roi, window, thresholdFactor);
}

void sigproc_tools::MorphologicalCNC::getSelectVals(
  const std::vector<std::vector<double>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const unsigned int window,
  const float thresholdFactor) const
{
  getSelectVals<double>(morphedWaveforms, selectVals,
    roi, window, thresholdFactor);
}

template <typename T>
void sigproc_tools::MorphologicalCNC::getSelectVals(
  const std::vector<std::vector<T>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const unsigned int window,
  const float thresholdFactor) const
{
  auto numChannels = morphedWaveforms.size();
  auto nTicks = morphedWaveforms.at(0).size();
  sigproc_tools::MiscUtils utils;

  for (size_t i=0; i<numChannels; ++i) {
    float median = 0.0;
    std::vector<T> localVec = morphedWaveforms[i];
    median = utils.computeMedian(localVec);
    std::vector<T> absoluteDeviation(localVec.size());
    for (size_t j=0; j<absoluteDeviation.size(); ++j) {
      absoluteDeviation[j] = std::abs(morphedWaveforms[i][j] - median);
    }

    float mad = utils.computeMedian(absoluteDeviation);
    float threshold;
    threshold = thresholdFactor * mad;

    for (size_t j=0; j<nTicks; ++j) {
      if (std::abs(morphedWaveforms[i][j] - median) > threshold) {
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


void sigproc_tools::MorphologicalCNC::denoiseMorph1D(
  std::vector<std::vector<short> >& waveLessCoherent,
  const std::vector<std::vector<short> >& fullEvent,
  std::vector<std::vector<bool> >& selectVals,
  std::vector<std::vector<bool> >& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor) const
{
  denoiseMorph1D<short>(waveLessCoherent, fullEvent, selectVals, roi,
    filterName, grouping, structuringElement, window, thresholdFactor);
}

void sigproc_tools::MorphologicalCNC::denoiseMorph1D(
  std::vector<std::vector<float> >& waveLessCoherent,
  const std::vector<std::vector<float> >& fullEvent,
  std::vector<std::vector<bool> >& selectVals,
  std::vector<std::vector<bool> >& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor) const
{
  denoiseMorph1D<float>(waveLessCoherent, fullEvent, selectVals, roi,
    filterName, grouping, structuringElement, window, thresholdFactor);
}

void sigproc_tools::MorphologicalCNC::denoiseMorph1D(
  std::vector<std::vector<double> >& waveLessCoherent,
  const std::vector<std::vector<double> >& fullEvent,
  std::vector<std::vector<bool> >& selectVals,
  std::vector<std::vector<bool> >& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor) const
{
  denoiseMorph1D<double>(waveLessCoherent, fullEvent, selectVals, roi,
    filterName, grouping, structuringElement, window, thresholdFactor);
}


template <typename T>
void sigproc_tools::MorphologicalCNC::denoiseMorph1D(
  std::vector<std::vector<T> >& waveLessCoherent,
  const std::vector<std::vector<T> >& fullEvent,
  std::vector<std::vector<bool> >& selectVals,
  std::vector<std::vector<bool> >& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor) const
{
  size_t numChannels = fullEvent.size();
  size_t nTicks = fullEvent.at(0).size();
  size_t nGroups = numChannels / grouping;

  sigproc_tools::Morph1DFast filter;
  sigproc_tools::MiscUtils utils;
  std::vector<std::vector<T>> morphedWaveforms;
  morphedWaveforms.resize(numChannels);
  for (auto& v : morphedWaveforms) {
    v.resize(nTicks);
  }

  switch (filterName) {
    case 'd':
      for (size_t i=0; i<numChannels; ++i) {
        filter.getDilation(fullEvent[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    case 'e':
      for (size_t i=0; i<numChannels; ++i) {
        filter.getErosion(fullEvent[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    case 'g':
      for (size_t i=0; i<numChannels; ++i) {
        filter.getGradient(fullEvent[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    default:
      for (size_t i=0; i<numChannels; ++i) {
        filter.getDilation(fullEvent[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
  }

  getSelectVals(morphedWaveforms, selectVals, roi, window, thresholdFactor);

  for (size_t i=0; i<nTicks; ++i) {
    for (size_t j=0; j<nGroups; ++j) {
      size_t group_start = j * grouping;
      size_t group_end = (j+1) * grouping;
      // Compute median.
      std::vector<T> v;
      v.reserve(grouping);
      for (size_t c=group_start; c<group_end; ++c) {
        if (!selectVals[c][i]) {
          v.push_back(fullEvent[c][i]);
        }
      }
      float median = 0.0;
      if (v.size() > 0) {
        median = utils.computeMedian(v);
      }
      for (auto k=group_start; k<group_end; ++k) {
        if (!selectVals[k][i]) {
          waveLessCoherent[k][i] = fullEvent[k][i] - median;
        } else {
          waveLessCoherent[k][i] = fullEvent[k][i];
        }
      }
    }
  }
  return;
}


void sigproc_tools::MorphologicalCNC::denoiseMorph2D(
  std::vector<std::vector<short>>& waveLessCoherent,
  const std::vector<std::vector<short>>& fullEvent,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor) const
{
  denoiseMorph2D<short>(
    waveLessCoherent, fullEvent, selectVals, roi,
    filterName, grouping, structuringElementx, structuringElementy,
    window, thresholdFactor);
  return;
}

void sigproc_tools::MorphologicalCNC::denoiseMorph2D(
  std::vector<std::vector<float>>& waveLessCoherent,
  const std::vector<std::vector<float>>& fullEvent,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor) const
{
  denoiseMorph2D<float>(
    waveLessCoherent, fullEvent, selectVals, roi,
    filterName, grouping, structuringElementx, structuringElementy,
    window, thresholdFactor);
  return;
}

void sigproc_tools::MorphologicalCNC::denoiseMorph2D(
  std::vector<std::vector<double>>& waveLessCoherent,
  const std::vector<std::vector<double>>& fullEvent,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor) const
{
  denoiseMorph2D<double>(
    waveLessCoherent, fullEvent, selectVals, roi,
    filterName, grouping, structuringElementx, structuringElementy,
    window, thresholdFactor);
  return;
}

template <typename T>
void sigproc_tools::MorphologicalCNC::denoiseMorph2D(
  std::vector<std::vector<T>>& waveLessCoherent,
  const std::vector<std::vector<T>>& fullEvent,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<bool>>& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor) const
{
  auto numChannels = fullEvent.size();
  auto nTicks = fullEvent.at(0).size();
  auto nGroups = numChannels / grouping;

  sigproc_tools::Morph2DFast filter;
  sigproc_tools::MiscUtils utils;
  std::vector<std::vector<T>> morphedWaveforms;
  morphedWaveforms.resize(numChannels);
  for (auto& v : morphedWaveforms) {
    v.resize(nTicks);
  }

  switch (filterName) {
    case 'd':
      filter.getDilation(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, roi, window, thresholdFactor);
      break;
    case 'e':
      filter.getErosion(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, roi, window, thresholdFactor);
      break;
    case 'g':
      filter.getGradient(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, roi, window, thresholdFactor);
      break;
    default:
      filter.getDilation(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, roi, window, thresholdFactor);
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
          v.push_back(fullEvent[c][i]);
        }
      }
      float median = 0.0;
      if (v.size() > 0) {
        median = utils.computeMedian(v);
      }
      for (size_t k=group_start; k<group_end; ++k) {
        if (!selectVals[k][i]) {
          waveLessCoherent[k][i] = fullEvent[k][i] - median;
        } else {
          waveLessCoherent[k][i] = fullEvent[k][i];
        }
      }
    }
  }
  return;
}

void sigproc_tools::MorphologicalCNC::denoiseHough2D(
  Array2D<short>& waveLessCoherent,
  const Array2D<short>& fullEvent,
  Array2D<bool>& selectVals,
  Array2D<bool>& refinedSelectVals,
  Array2D<bool>& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor,
  const size_t thetaSteps,
  const unsigned int houghThreshold,
  const unsigned int nmsWindowSize,
  const unsigned int angleWindow, 
  const unsigned int dilationX,
  const unsigned int dilationY,
  const unsigned int maxLines,
  const float eps) const
{
  denoiseHough2D<short>(
    waveLessCoherent, 
    fullEvent,
    selectVals,
    refinedSelectVals,
    roi,
    filterName,
    grouping,
    structuringElementx,
    structuringElementy,
    window,
    thresholdFactor,
    thetaSteps,
    houghThreshold,
    nmsWindowSize,
    angleWindow,
    dilationX,
    dilationY,
    maxLines,
    eps
  );
}

void sigproc_tools::MorphologicalCNC::denoiseHough2D(
  Array2D<float>& waveLessCoherent,
  const Array2D<float>& fullEvent,
  Array2D<bool>& selectVals,
  Array2D<bool>& refinedSelectVals,
  Array2D<bool>& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor,
  const size_t thetaSteps,
  const unsigned int houghThreshold,
  const unsigned int nmsWindowSize,
  const unsigned int angleWindow, 
  const unsigned int dilationX,
  const unsigned int dilationY,
  const unsigned int maxLines,
  const float eps) const
{
  denoiseHough2D<float>(
    waveLessCoherent, 
    fullEvent,
    selectVals,
    refinedSelectVals,
    roi,
    filterName,
    grouping,
    structuringElementx,
    structuringElementy,
    window,
    thresholdFactor,
    thetaSteps,
    houghThreshold,
    nmsWindowSize,
    angleWindow,
    dilationX,
    dilationY,
    maxLines,
    eps
  );
}

void sigproc_tools::MorphologicalCNC::denoiseHough2D(
  Array2D<double>& waveLessCoherent,
  const Array2D<double>& fullEvent,
  Array2D<bool>& selectVals,
  Array2D<bool>& refinedSelectVals,
  Array2D<bool>& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor,
  const size_t thetaSteps,
  const unsigned int houghThreshold,
  const unsigned int nmsWindowSize,
  const unsigned int angleWindow, 
  const unsigned int dilationX,
  const unsigned int dilationY,
  const unsigned int maxLines,
  const float eps) const
{
  denoiseHough2D<double>(
    waveLessCoherent, 
    fullEvent,
    selectVals,
    refinedSelectVals,
    roi,
    filterName,
    grouping,
    structuringElementx,
    structuringElementy,
    window,
    thresholdFactor,
    thetaSteps,
    houghThreshold,
    nmsWindowSize,
    angleWindow,
    dilationX,
    dilationY,
    maxLines,
    eps
  );
}

template <typename T>
void sigproc_tools::MorphologicalCNC::denoiseHough2D(
  Array2D<T>& waveLessCoherent,
  const Array2D<T>& fullEvent,
  Array2D<bool>& selectVals,
  Array2D<bool>& refinedSelectVals,
  Array2D<bool>& roi,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor,
  const size_t thetaSteps,
  const unsigned int houghThreshold,
  const unsigned int nmsWindowSize,
  const unsigned int angleWindow, 
  const unsigned int dilationX,
  const unsigned int dilationY,
  const unsigned int maxLines,
  const float eps) const
{
  auto numChannels = fullEvent.size();
  auto nTicks = fullEvent.at(0).size();
  auto nGroups = numChannels / grouping;

  sigproc_tools::Morph2DFast filter;
  sigproc_tools::MiscUtils utils;
  sigproc_tools::LineDetection lineModule;
  std::vector<std::vector<T>> morphedWaveforms;
  morphedWaveforms.resize(numChannels);
  for (auto& v : morphedWaveforms) {
    v.resize(nTicks);
  }

  switch (filterName) {
    case 'd':
      filter.getDilation(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, roi, window, thresholdFactor);
      break;
    case 'e':
      filter.getErosion(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, roi, window, thresholdFactor);
      break;
    case 'g':
      filter.getGradient(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, roi, window, thresholdFactor);
      break;
    default:
      filter.getDilation(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, roi, window, thresholdFactor);
      break;
  }

  refinedSelectVals.resize(selectVals.size());
  for (auto& v : refinedSelectVals) {
    v.resize(selectVals.at(0).size());
  }

  lineModule.refineSelectVals(
    selectVals,
    refinedSelectVals,
    thetaSteps,
    houghThreshold,
    angleWindow,
    maxLines,
    nmsWindowSize,
    dilationX,
    dilationY,
    eps);

  for (size_t i=0; i<nTicks; ++i) {
    for (size_t j=0; j<nGroups; ++j) {
      size_t group_start = j * grouping;
      size_t group_end = (j+1) * grouping;
      // Compute median.
      std::vector<T> v;
      for (size_t c=group_start; c<group_end; ++c) {
        if (!refinedSelectVals[c][i]) {
          v.push_back(fullEvent[c][i]);
        }
      }
      float median = 0.0;
      if (v.size() > 0) {
        median = utils.computeMedian(v);
      }
      for (size_t k=group_start; k<group_end; ++k) {
        if (!refinedSelectVals[k][i]) {
          waveLessCoherent[k][i] = fullEvent[k][i] - median;
        } else {
          waveLessCoherent[k][i] = fullEvent[k][i];
        }
      }
    }
  }
  return;
}
#endif
