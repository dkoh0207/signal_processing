#ifndef __SIGPROC_MULTITHREADING_COHERENTNOISECORRECTION_CXX__
#define __SIGPROC_MULTITHREADING_COHERENTNOISECORRECTION_CXX__

#include "CoherentNoiseCorrection.h"

using namespace sigproc_multithreading;

short sigproc_multithreading::CoherentNoiseCorrection::computeMedian(
    const ConcurrentVector<short> &cVector) const
{
    short median = computeMedian<short>(cVector);
    return median;
}

float sigproc_multithreading::CoherentNoiseCorrection::computeMedian(
    const ConcurrentVector<float> &cVector) const
{
    float median = computeMedian<float>(cVector);
    return median;
}

double sigproc_multithreading::CoherentNoiseCorrection::computeMedian(
    const ConcurrentVector<double> &cVector) const
{
    double median = computeMedian<double>(cVector);
    return median;
}

template <typename T>
T sigproc_multithreading::CoherentNoiseCorrection::computeMedian(
    const ConcurrentVector<T> &cVector) const
{
    T median = (T) 0;

    if (cVector.size() < 1) return median;

    ConcurrentVector<T> localVec = cVector;
    const auto m = localVec.begin() + localVec.size() / 2;
    std::nth_element(localVec.begin(), m, localVec.end());
    median = *m;
    return median;
}

void sigproc_multithreading::CoherentNoiseCorrection::getSelectVals(
    const Array2D<short> &morphedWaveforms,
    Array2D<bool> &selectVals,
    const short threshold) const
{
    getSelectVals<short>(morphedWaveforms, selectVals, threshold);
    return;
}

void sigproc_multithreading::CoherentNoiseCorrection::getSelectVals(
    const Array2D<float> &morphedWaveforms,
    Array2D<bool> &selectVals,
    const float threshold) const
{
    getSelectVals<float>(morphedWaveforms, selectVals, threshold);
    return;
}

void sigproc_multithreading::CoherentNoiseCorrection::getSelectVals(
    const Array2D<double> &morphedWaveforms,
    Array2D<bool> &selectVals,
    const double threshold) const
{
    getSelectVals<double>(morphedWaveforms, selectVals, threshold);
    return;
}

template <typename T>
void sigproc_multithreading::CoherentNoiseCorrection::getSelectVals(
    const Array2D<T> &morphedWaveforms,
    Array2D<bool> &selectVals,
    const T threshold) const
{
  size_t numChannels = morphedWaveforms.size();
  size_t nTicks = morphedWaveforms.at(0).size();

  for (size_t i=0; i<numChannels; ++i) {

    for (size_t j=0; j<nTicks; ++j) {

        if (std::abs(morphedWaveforms[i][j]) > threshold) {
            selectVals[i][j] = true;
        }
        else selectVals[i][j] = false;
    }
  }
  return;
}

void sigproc_multithreading::CoherentNoiseCorrection::ParallelDenoiseMorph2D(
    const Array2D<short> &inputArray2D,
    Array2D<bool> &selectVals,
    const char filterName,
    Array2D<short> &waveLessCoherent,
    const size_t fStructuringElementx,
    const size_t fStructuringElementy,
    const size_t grouping,
    const size_t groupingOffset,
    const short threshold) const
{
    ParallelDenoiseMorph2D<short>(
        inputArray2D,
        selectVals,
        filterName,
        waveLessCoherent,
        fStructuringElementx,
        fStructuringElementy, 
        grouping,
        groupingOffset,
        threshold
    );
    return;
}

void sigproc_multithreading::CoherentNoiseCorrection::ParallelDenoiseMorph2D(
    const Array2D<float> &inputArray2D,
    Array2D<bool> &selectVals,
    const char filterName,
    Array2D<float> &waveLessCoherent,
    const size_t fStructuringElementx,
    const size_t fStructuringElementy,
    const size_t grouping,
    const size_t groupingOffset,
    const float threshold) const
{
    ParallelDenoiseMorph2D<float>(
        inputArray2D,
        selectVals,
        filterName,
        waveLessCoherent,
        fStructuringElementx,
        fStructuringElementy, 
        grouping,
        groupingOffset,
        threshold
    );
    return;
}

void sigproc_multithreading::CoherentNoiseCorrection::ParallelDenoiseMorph2D(
    const Array2D<double> &inputArray2D,
    Array2D<bool> &selectVals,
    const char filterName,
    Array2D<double> &waveLessCoherent,
    const size_t fStructuringElementx,
    const size_t fStructuringElementy,
    const size_t grouping,
    const size_t groupingOffset,
    const double threshold) const
{
    ParallelDenoiseMorph2D<double>(
        inputArray2D,
        selectVals,
        filterName,
        waveLessCoherent,
        fStructuringElementx,
        fStructuringElementy, 
        grouping,
        groupingOffset,
        threshold
    );
    return;
}

template <typename T>
void sigproc_multithreading::CoherentNoiseCorrection::ParallelDenoiseMorph2D(
    const Array2D<T> &inputArray2D,
    Array2D<bool> &selectVals,
    const char filterName,
    Array2D<T> &waveLessCoherent,
    const size_t fStructuringElementx,
    const size_t fStructuringElementy,
    const size_t grouping,
    const size_t groupingOffset,
    const T threshold) const
{
    // Compute Morphological Filter
    const size_t numChannels = inputArray2D.size();
    const size_t numTicks = inputArray2D.at(0).size();
    const size_t nGroups = (size_t) ((int) numChannels - 
                                     (int) groupingOffset) / grouping;

    Array2D<T> buffer(numChannels, Vector<T>(numTicks));

    sigproc_multithreading::Morph2DFast morph2D;

    switch (filterName) {
        case 'd':
            morph2D.getDilation(inputArray2D,
                fStructuringElementx, fStructuringElementy, buffer);
            break;
        case 'e':
            morph2D.getErosion(inputArray2D,
                fStructuringElementx, fStructuringElementy, buffer);
            break;
        case 'g':
            morph2D.getGradient(inputArray2D,
                fStructuringElementx, fStructuringElementy, buffer);
            break;
        default:
            morph2D.getDilation(inputArray2D,
                fStructuringElementx, fStructuringElementy, buffer);
        break;
    }

    getSelectVals(buffer, selectVals, threshold);

    tbb::parallel_for( (size_t) 0, numTicks, (size_t) 1,
        [this, &grouping, &groupingOffset, &nGroups, 
         &selectVals, &inputArray2D, &waveLessCoherent](size_t i) 
        {
            for (size_t j=0; j<nGroups; ++j) {

                size_t group_start = j * grouping + (size_t) groupingOffset;
                size_t group_end = (j+1) * grouping + (size_t) groupingOffset;

                // Compute median. Need a concurrent vector due to push_back
                ConcurrentVector<T> v;

                for (size_t c=group_start; c<group_end; ++c) {

                    if (!selectVals[c][i]) v.push_back(inputArray2D[c][i]);
                }

                float median = computeMedian(v);

                for (size_t k=group_start; k<group_end; ++k) {

                    if (!selectVals[k][i]) {
                        waveLessCoherent[k][i] = inputArray2D[k][i] - median;
                    } else {
                        waveLessCoherent[k][i] = inputArray2D[k][i];
                    }
                }
            }
        }
    );

    if (groupingOffset > 0) {
        for (size_t i=0; i<numTicks; ++i) {

            ConcurrentVector<T> v;

            for (size_t c=0; c<groupingOffset; ++c) {
                if (!selectVals[c][i]) v.push_back(inputArray2D[c][i]);
            }

            float median = computeMedian(v);

            for (size_t k=0; k<groupingOffset; ++k) {

                if (!selectVals[k][i]) {
                    waveLessCoherent[k][i] = inputArray2D[k][i] - median;
                } else {
                    waveLessCoherent[k][i] = inputArray2D[k][i];
                }
            }
        }
    }

    return;
}


#endif
