#ifndef __SIGPROC_MULTITHREADING_ROIFINDER_CXX__
#define __SIGPROC_MULTITHREADING_ROIFINDER_CXX__

#include "ROIFinder.h"


void sigproc_multithreading::ROIFinder1::applyChain(
    const Array2D<float>& input2D,
    Array2D<bool>& output2D,
    size_t FREQUENCY_THRESHOLD,
    size_t FREQUENCY_FILTER_SMOOTHNESS_ORDER,
    char MORPHOLOGICAL_FILTER_NAME,
    const unsigned int CHANNEL_GROUPING,
    const unsigned int CHANNEL_GROUPING_OFFSET,
    const unsigned int STRUCTURING_ELEMENT_X,
    const unsigned int STRUCTURING_ELEMENT_Y,
    const float MORPHOLOGICAL_THRESHOLD_FACTOR,

    const size_t THETASTEPS,
    const unsigned int HOUGH_THRESHOLD,
    const unsigned int NMS_WINDOW_SIZE,
    const unsigned int ANGLE_WINDOW,

    const unsigned int ADFILTER_SX,
    const unsigned int ADFILTER_SY,
    const float sigma_x,
    const float sigma_y,
    const float sigma_r,
    const float lowThreshold,
    const float highThreshold,

    const unsigned int BINARY_CLOSING_SX,
    const unsigned int BINARY_CLOSING_SY) const
{
    int numChannels = input2D.size();
    int numTicks = input2D.at(0).size();

    Array2D<float> buffer(numChannels, VectorFloat(numTicks));
    Array2D<float> smoothedWf(numChannels, VectorFloat(numTicks));

    Array2D<bool> selectVals(numChannels, VectorBool(numTicks));
    Array2D<bool> refinedSelectVals(numChannels, VectorBool(numTicks));
    Array2D<float> waveLessCoherent(numChannels, VectorFloat(numTicks));

    Array2D<float> gradX(numChannels, VectorFloat(numTicks));
    Array2D<float> gradY(numChannels, VectorFloat(numTicks));
    Array2D<float> direction(numChannels, VectorFloat(numTicks));

    sigproc_tools::MiscUtils utils;

    // 1. Remove Pedestals
    
    for (int i=0; i<numChannels; ++i) {
        float median = utils.computeMedian(input2D[i]);
        for (int j=0; j<numTicks; ++j) {
            buffer[i][j] = input2D[i][j] - median;
        }
    }

    // 2. Coherent Noise Correction

    sigproc_tools::MorphologicalCNC denoiser;

    denoiser.denoiseRestrictedHough2D(
        waveLessCoherent, buffer, selectVals, refinedSelectVals, 
        MORPHOLOGICAL_FILTER_NAME, CHANNEL_GROUPING, CHANNEL_GROUPING_OFFSET,
        STRUCTURING_ELEMENT_X, STRUCTURING_ELEMENT_Y, MORPHOLOGICAL_THRESHOLD_FACTOR,
        THETASTEPS, HOUGH_THRESHOLD, NMS_WINDOW_SIZE, NMS_WINDOW_SIZE, 15, 15, ANGLE_WINDOW);

    // 3. High Pass Filtering

    sigproc_tools::FrequencyFilters1D freqFilt;
    freqFilt.filterImage(
        waveLessCoherent,
        FREQUENCY_THRESHOLD,
        buffer,
        FREQUENCY_FILTER_SMOOTHNESS_ORDER, 0);

    // 4. Sobel Filtering

    sigproc_tools::EdgeDetection edgeUtils;
    edgeUtils.SobelRads(buffer, gradX, gradY, waveLessCoherent, direction);

    // 5. Bilateral Filtering
    sigproc_tools::BilateralFilters smoother;
    smoother.directional(buffer, direction, smoothedWf,
        ADFILTER_SX, ADFILTER_SY, sigma_x, sigma_y, sigma_r, 360);

    // 6. Enhance signal by morphing 
    sigproc_tools::Morph2DFast morph2D;
    if (MORPHOLOGICAL_FILTER_NAME == 'e') {
        morph2D.getErosion(smoothedWf, ADFILTER_SX, ADFILTER_SY, buffer);
    }
    else if (MORPHOLOGICAL_FILTER_NAME == 'd') {
        morph2D.getDilation(smoothedWf, ADFILTER_SX, ADFILTER_SY, buffer);
    } 
    else if (MORPHOLOGICAL_FILTER_NAME == 'g') {
        morph2D.getGradient(smoothedWf, ADFILTER_SX, ADFILTER_SY, buffer);
    }
    else {
        morph2D.getGradient(smoothedWf, ADFILTER_SX, ADFILTER_SY, buffer);
    }

    // 7. Second Sobel Filtering
    edgeUtils.Sobel(buffer, gradX, gradY, waveLessCoherent, direction);
    edgeUtils.EdgeNMSInterpolation(waveLessCoherent, gradX, gradY, direction, smoothedWf);

    // 8. Thresholding
    edgeUtils.HTFastLowMem(smoothedWf, lowThreshold, highThreshold, selectVals);

    // 9. Morph Binary
    morph2D.getDilation(selectVals, BINARY_CLOSING_SX, BINARY_CLOSING_SY, output2D);

    return;
}

#endif
