#ifndef __SIGPROC_TOOLS_FINDROI2D_CXX__
#define __SIGPROC_TOOLS_FINDROI2D_CXX__

#include "FindROI2D.h"

void sigproc_tools::FindROI2D::applyChainFilter(
  const Array2D<float>& waveform2D,
  Array2D<float>& fullEvent,
  Array2D<bool>& outputROI,
  Array2D<float>& waveLessCoherent,
  Array2D<float>& morphedWaveform2D,
  Array2D<float>& finalErosion2D,
  // Default Parameters
  size_t FREQUENCY_THRESHOLD = 30,
  size_t FREQUENCY_FILTER_SMOOTHNESS_ORDER = 2,
  size_t FREQUENCY_FILTER_MODE = 0,

  char MORPHOLOGICAL_FILTER_NAME = 'e',
  const unsigned int CHANNEL_GROUPING = 32,
  const unsigned int STRUCTURING_ELEMENT_X = 7,
  const unsigned int STRUCTURING_ELEMENT_Y = 20,
  const unsigned int ROI_EXPAND_WINDOW_SIZE = 10,
  const float MORPHOLOGICAL_THRESHOLD_FACTOR = 2.5,

  const size_t THETASTEPS = 360,
  const unsigned int HOUGH_THRESHOLD = 400,
  const unsigned int NMS_WINDOW_SIZE = 20,
  const unsigned int ANGLE_WINDOW = 20,

  // float NOISE_VARIANCE = 20.0;
  const unsigned int ADFILTER_SX = 5,
  const unsigned int ADFILTER_SY = 15,

  const float GLOBAL_THRESHOLDING_FACTOR = 2.0) const
{
  // All input arrays must have the same dimensions as waveform2D

  int numChannels = waveform2D.size();
  int numTicks = waveform2D.at(0).size();

  fullEvent.resize(numChannels);
  for (auto& v : fullEvent) {
    v.resize(numTicks);
  }

  sigproc_tools::MiscUtils utils;

  // 1. Remove Pedestals
  for (int i=0; i<numChannels; ++i) {
    float median = utils.computeMedian(waveform2D[i]);
    for (int j=0; j<numTicks; ++j) {
      fullEvent[i][j] = waveform2D[i][j] - median;
    }
  }

  // 2. Buffer for intermediate computations
  Array2D<float> buffer(numChannels);
  Array2D<bool> selectVals(numChannels);
  Array2D<bool> rois(numChannels);
  Array2D<bool> refinedSelectVals(numChannels);

  for (auto& v : buffer) {
    v.resize(numTicks);
  }

  for (auto& v : selectVals) {
    v.resize(numTicks);
  }

  for (auto& v : rois) {
    v.resize(numTicks);
  }

  for (auto& v : refinedSelectVals) {
    v.resize(numTicks);
  }

  // 3. Apply frequency high pass filters
  sigproc_tools::FrequencyFilters1D freqFilt;
  freqFilt.filterImage(
    fullEvent,
    FREQUENCY_THRESHOLD,
    buffer,
    FREQUENCY_FILTER_SMOOTHNESS_ORDER,
    FREQUENCY_FILTER_MODE);
  
  // 4. Run Coherent Noise Correction
  sigproc_tools::MorphologicalCNC denoiser;

  buffer.resize(numChannels);
  for (auto& v : buffer) {
    v.resize(numTicks);
  }

  denoiser.denoiseHough2D(
    waveLessCoherent,
    morphedWaveform2D,
    buffer,
    selectVals,
    refinedSelectVals,
    rois,
    MORPHOLOGICAL_FILTER_NAME,
    CHANNEL_GROUPING,
    STRUCTURING_ELEMENT_X,
    STRUCTURING_ELEMENT_Y,
    ROI_EXPAND_WINDOW_SIZE,
    MORPHOLOGICAL_THRESHOLD_FACTOR,
    THETASTEPS,
    HOUGH_THRESHOLD,
    NMS_WINDOW_SIZE,
    ANGLE_WINDOW);

  // 5. Run Adaptive (Incoherent) noise filtering

  sigproc_tools::AdaptiveWiener adFilter;
  
  buffer.resize(numChannels);
  for (auto& v : buffer) {
    v.resize(numTicks);
  }

  adFilter.MMWFStar(
    buffer,
    waveLessCoherent,
    ADFILTER_SX,
    ADFILTER_SY
  );

  finalErosion2D.resize(numChannels);
  for (auto& v : finalErosion2D) {
    v.resize(numTicks);
  }

  sigproc_tools::Morph2DFast morph2D;

  buffer.resize(numChannels);
  for (auto& v : buffer) {
    v.resize(numTicks);
  }

  if (MORPHOLOGICAL_FILTER_NAME == 'e') {
    morph2D.getErosion(
      buffer,
      STRUCTURING_ELEMENT_X,
      STRUCTURING_ELEMENT_Y,
      finalErosion2D
    );
  }
  else if (MORPHOLOGICAL_FILTER_NAME == 'd') {
    morph2D.getDilation(
      buffer,
      STRUCTURING_ELEMENT_X,
      STRUCTURING_ELEMENT_Y,
      finalErosion2D
    );
  }
  else if (MORPHOLOGICAL_FILTER_NAME == 'g') {
    morph2D.getGradient(
      buffer,
      STRUCTURING_ELEMENT_X,
      STRUCTURING_ELEMENT_Y,
      finalErosion2D
    );
  }
  else {
    morph2D.getGradient(
      buffer,
      STRUCTURING_ELEMENT_X,
      STRUCTURING_ELEMENT_Y,
      finalErosion2D
    );
  }

  // Take absolute value
  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      if (finalErosion2D[i][j] < 0) {
        finalErosion2D[i][j] = -finalErosion2D[i][j];
      }
    }
  }

  // 6. Run Thresholding
  sigproc_tools::Thresholding thresholder;
  thresholder.globalMean(
    finalErosion2D,
    outputROI,
    GLOBAL_THRESHOLDING_FACTOR
  );

  // 7. Refine ROI




  return;
}

#endif
