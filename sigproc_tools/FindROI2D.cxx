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
  const unsigned int STRUCTURING_ELEMENT_Y = 28,
  const unsigned int ROI_EXPAND_WINDOW_SIZE = 10,
  const float MORPHOLOGICAL_THRESHOLD_FACTOR = 2.5,

  const size_t THETASTEPS = 360,
  const unsigned int HOUGH_THRESHOLD = 999999,
  const unsigned int NMS_WINDOW_SIZE = 20,
  const unsigned int ANGLE_WINDOW = 20,

  // float NOISE_VARIANCE = 20.0;
  const unsigned int ADFILTER_SX = 5,
  const unsigned int ADFILTER_SY = 20,

  const unsigned int BINARY_CLOSING_SX = 31,
  const unsigned int BINARY_CLOSING_SY = 31,

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
  std::cout << "1" << std::endl;
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
  std::cout << "2" << std::endl;
  sigproc_tools::AdaptiveWiener adFilter;
  std::cout << "3" << std::endl;
  buffer.resize(numChannels);
  for (auto& v : buffer) {
    v.resize(numTicks);
  }
  std::cout << "4" << std::endl;
  adFilter.MMWFStar(
    buffer,
    waveLessCoherent,
    ADFILTER_SX,
    ADFILTER_SY
  );
  std::cout << "5" << std::endl;
  finalErosion2D.resize(numChannels);
  for (auto& v : finalErosion2D) {
    v.resize(numTicks);
  }
  std::cout << "6" << std::endl;
  sigproc_tools::Morph2DFast morph2D;

  buffer.resize(numChannels);
  for (auto& v : buffer) {
    v.resize(numTicks);
  }
  std::cout << "7" << std::endl;
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
  std::cout << "8" << std::endl;
  // Take absolute value
  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      if (finalErosion2D[i][j] < 0) {
        finalErosion2D[i][j] = -finalErosion2D[i][j];
      }
    }
  }
  std::cout << "9" << std::endl;
  // 6. Run Thresholding
  sigproc_tools::Thresholding thresholder;
  thresholder.globalMean(
    finalErosion2D,
    rois,
    GLOBAL_THRESHOLDING_FACTOR
  );
  std::cout << "10" << std::endl;
  morph2D.getClosing(rois, BINARY_CLOSING_SX, BINARY_CLOSING_SY, outputROI);
  // 7. Refine ROI
  return;
}


// void sigproc_tools::FindROI2D::applyCannyFilter(
//   const Array2D<float>& waveform2D,
//   Array2D<float>& fullEvent,
//   Array2D<bool>& outputROI,
//   Array2D<float>& waveLessCoherent,
//   Array2D<float>& morphedWaveform2D,
//   // Default Parameters
//   size_t FREQUENCY_THRESHOLD = 30,
//   size_t FREQUENCY_FILTER_SMOOTHNESS_ORDER = 2,
//   size_t FREQUENCY_FILTER_MODE = 0,

//   // Coherent Noise Correction Parameters
//   char MORPHOLOGICAL_FILTER_NAME = 'e',
//   const unsigned int CHANNEL_GROUPING = 32,
//   const unsigned int CHANNEL_GROUPING_OFFSET = 0,
//   const unsigned int STRUCTURING_ELEMENT_X = 7,
//   const unsigned int STRUCTURING_ELEMENT_Y = 7,
//   const unsigned int ROI_EXPAND_WINDOW_SIZE = 10,
//   const float MORPHOLOGICAL_THRESHOLD_FACTOR = 2.5,

//   // Hough Transform Parameters
//   const size_t THETASTEPS = 360,
//   const unsigned int HOUGH_THRESHOLD = 400,
//   const unsigned int NMS_WINDOW_SIZE = 20,
//   const unsigned int ANGLE_WINDOW = 20,

//   // float NOISE_VARIANCE = 20.0;
//   const unsigned int ADFILTER_SX = 7,
//   const unsigned int ADFILTER_SY = 7,
//   const float sigma_x = 5.0, 
//   const float sigma_y = 5.0, 
//   const float sigma_r = 30.0, 
//   const float lowThreshold = 3.0,
//   const float highThreshold = 15.0,

//   const unsigned int BINARY_CLOSING_SX = 13,
//   const unsigned int BINARY_CLOSING_SY = 13) const
// {
//   // All input arrays must have the same dimensions as waveform2D

//   size_t numChannels = waveform2D.size();
//   size_t numTicks = waveform2D.at(0).size();

//   Denoising basicDenoise;
//   FrequencyFilters1D freqFilter;
//   MorphologicalCNC CNC;
//   Morph2DFast morph;
//   EdgeDetection edges;
//   BilateralFilters bilateral;

//   Array2D<float> buffer(numChannels, std::vector<float>(numTicks));

//   Array2D<bool> binaryBuffer(numChannels, std::vector<bool>(numTicks));

//   Array2D<float> sobelX(numChannels, std::vector<float>(numTicks));

//   Array2D<float> sobelY(numChannels, std::vector<float>(numTicks));

//   Array2D<float> gradient(numChannels, std::vector<float>(numTicks));

//   Array2D<float> direction(numChannels, std::vector<float>(numTicks));

//   for (size_t i=0; i<numChannels; ++i) {
//     for (size_t j=0; j<numTicks; ++j) {
//       buffer[i][j] = waveform2D[i][j];
//     }
//   }

//   basicDenoise.subtractPedestals(buffer);

//   freqFilter.filterImage(buffer, FREQUENCY_THRESHOLD, waveLessCoherent,
//                          FREQUENCY_FILTER_SMOOTHNESS_ORDER, 0);

//   CNC.simpleCNC(waveLessCoherent, buffer, CHANNEL_GROUPING, CHANNEL_GROUPING_OFFSET);

//   edges.Sobel(buffer, sobelX, sobelY, gradient, direction);

//   bilateral.directional(buffer, direction, waveLessCoherent, 
//     STRUCTURING_ELEMENT_X, STRUCTURING_ELEMENT_Y, 
//     sigma_x, sigma_y, sigma_r, THETASTEPS);

//   morph.getDilation(waveLessCoherent, ADFILTER_SX, ADFILTER_SY, buffer);

//   edges.Sobel(buffer, sobelX, sobelY, gradient, direction);

//   edges.EdgeNMSInterpolation(gradient, sobelX, sobelY, direction, buffer);

//   edges.HTFastLowMem(buffer, lowThreshold, highThreshold, binaryBuffer);

//   morph.getDilation(binaryBuffer, BINARY_CLOSING_SX, BINARY_CLOSING_SY, outputROI);

//   return;
// }

#endif
