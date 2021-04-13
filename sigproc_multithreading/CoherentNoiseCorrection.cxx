#ifndef __SIGPROC_MULTITHREADING_COHERENTNOISECORRECTION_CXX__
#define __SIGPROC_MULTITHREADING_COHERENTNOISECORRECTION_CXX__

#include "CoherentNoiseCorrection.h"


float sigproc_multithreading::CoherentNoiseCorrection::computeMedian(
  const ConcurrentVector<float>& inputVector) const
{
  float median = 0.0;
  ConcurrentVector<float> localVec = inputVector;
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
  return median;
}


void sigproc_multithreading::CoherentNoiseCorrection::getSelectVals(
  const ConcurrentArray2D<float>& morphedWaveforms,
  ConcurrentArray2D<bool>& selectVals,
  const unsigned int window,
  const float thresholdFactor) const
{
  size_t numChannels = morphedWaveforms.size();
  size_t nTicks = morphedWaveforms.at(0).size();

  for (size_t i=0; i<numChannels; ++i) {
    float median = 0.0;
    median = computeMedian(morphedWaveforms[i]);
    ConcurrentVector<float> absoluteDeviation(nTicks);
    for (size_t j=0; j<nTicks; ++j) {
      absoluteDeviation[j] = std::abs(morphedWaveforms[i][j] - median);
    }

    float mad = computeMedian(absoluteDeviation);
    float threshold = thresholdFactor * mad;

    for (size_t j=0; j<nTicks; ++j) {
      if (std::abs(morphedWaveforms[i][j] - median) > threshold) {
        selectVals[i][j] = true;
      } else {
        selectVals[i][j] = false;
      }
    }
  }
  return;
}


void sigproc_multithreading::CoherentNoiseCorrection::denoiseMorph2D(
  ConcurrentArray2D<float>& waveLessCoherent,
  ConcurrentArray2D<float>& morphedWaveforms,
  const ConcurrentArray2D<float>& fullEvent,
  ConcurrentArray2D<bool>& selectVals,
  const char filterName,
  const unsigned int grouping,
  const unsigned int groupingOffset,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor) const
{
  size_t numChannels = fullEvent.size();
  size_t nTicks = fullEvent.at(0).size();
  assert(groupingOffset < numChannels);
  size_t nGroups = ((int) numChannels - (int) groupingOffset) / grouping;

  sigproc_multithreading::Morph2DFast filter;

  switch (filterName) {
    case 'd':
      filter.getDilation(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, window, thresholdFactor);
      break;
    case 'e':
      filter.getErosion(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, window, thresholdFactor);
      break;
    case 'g':
      filter.getGradient(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, window, thresholdFactor);
      break;
    default:
      filter.getDilation(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms,
        selectVals, window, thresholdFactor);
      break;
  }

  for (size_t i=0; i<nTicks; ++i) {
    for (size_t j=0; j<nGroups; ++j) {
      size_t group_start = j * grouping + (size_t) groupingOffset;
      size_t group_end = (j+1) * grouping + (size_t) groupingOffset;
      // Compute median.
      ConcurrentVector<float> v;
      for (size_t c=group_start; c<group_end; ++c) {
        if (!selectVals[c][i]) {
          v.push_back(fullEvent[c][i]);
        }
      }
      float median = 0.0;
      if (v.size() > 0) median = computeMedian(v);
      for (size_t k=group_start; k<group_end; ++k) {
        if (!selectVals[k][i]) {
          waveLessCoherent[k][i] = fullEvent[k][i] - median;
        } else {
          waveLessCoherent[k][i] = fullEvent[k][i];
        }
      }
    }
  }

  // Compensate for offset in channel groupings
  if (groupingOffset > 0) {
    for (size_t i=0; i<nTicks; ++i) {
      ConcurrentVector<float> v;
      for (size_t c=0; c<groupingOffset; ++c) {
        if (!selectVals[c][i]) {
          v.push_back(fullEvent[c][i]);
        }
      }
      float median = 0.0;
      if (v.size() > 0) {
        median = computeMedian(v);
      }
      for (size_t k=0; k<groupingOffset; ++k) {
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


void sigproc_multithreading::CoherentNoiseCorrection::denoiseHough2D(
  ConcurrentArray2D<float>& waveLessCoherent,
  ConcurrentArray2D<float>& morphedWaveforms,
  const ConcurrentArray2D<float>& fullEvent,
  ConcurrentArray2D<bool>& selectVals,
  ConcurrentArray2D<bool>& refinedSelectVals,
  const char filterName,
  const unsigned int grouping,
  const unsigned int groupingOffset,
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
  size_t numChannels = fullEvent.size();
  size_t nTicks = fullEvent.at(0).size();
  assert(groupingOffset < numChannels);
  size_t nGroups = ((int) numChannels - (int) groupingOffset) / grouping;

  sigproc_multithreading::Morph2DFast filter;
  // sigproc_multithreading::LineDetection lineModule;

  switch (filterName) {
    case 'd':
      filter.getDilation(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms, selectVals, window, thresholdFactor);
      break;
    case 'e':
      filter.getErosion(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms, selectVals, window, thresholdFactor);
      break;
    case 'g':
      filter.getGradient(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms, selectVals, window, thresholdFactor);
      break;
    default:
      filter.getDilation(fullEvent,
        structuringElementx, structuringElementy, morphedWaveforms);
      getSelectVals(morphedWaveforms, selectVals, window, thresholdFactor);
      break;
  }

  // lineModule.refineSelectVals(
  //   selectVals,
  //   refinedSelectVals,
  //   thetaSteps,
  //   houghThreshold,
  //   angleWindow,
  //   maxLines,
  //   nmsWindowSize,
  //   dilationX,
  //   dilationY,
  //   eps);

  for (size_t i=0; i<refinedSelectVals.size(); ++i) {
    for (size_t j=0; j<refinedSelectVals.at(0).size(); ++j) {
      refinedSelectVals[i][j] = refinedSelectVals[i][j] && selectVals[i][j];
    }
  }

  for (size_t i=0; i<nTicks; ++i) {
    for (size_t j=0; j<nGroups; ++j) {
      size_t group_start = j * grouping + (size_t) groupingOffset;
      size_t group_end = (j+1) * grouping + (size_t) groupingOffset;
      // Compute median.
      ConcurrentVector<float> v;
      for (size_t c=group_start; c<group_end; ++c) {
        if (!refinedSelectVals[c][i]) {
          v.push_back(fullEvent[c][i]);
        }
      }
      float median = 0.0;
      if (v.size() > 0) {
        median = computeMedian(v);
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

  // Compensate for offset in channel groupings
  if (groupingOffset > 0) {
    for (size_t i=0; i<nTicks; ++i) {
      ConcurrentVector<float> v;
      for (size_t c=0; c<groupingOffset; ++c) {
        if (!refinedSelectVals[c][i]) {
          v.push_back(fullEvent[c][i]);
        }
      }
      float median = 0.0;
      if (v.size() > 0) {
        median = computeMedian(v);
      }
      for (size_t k=0; k<groupingOffset; ++k) {
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
