#ifndef __ALGORITHMS_WAVEFORMUTILS_CXX__
#define __ALGORITHMS_WAVEFORMUTILS_CXX__

#include "WaveformUtils.h"

void algorithms::WaveformUtils::getWaveformParams(const std::vector<short>& waveform,
                                                 float& mean,
                                                 float& median,
                                                 float& mode,
                                                 float& skewness,
                                                 float& rms)
{
  getWaveformParams<short>(waveform, mean, median, mode, skewness, rms);
  return;
}

void algorithms::WaveformUtils::getWaveformParams(const std::vector<float>& waveform,
                                                 float& mean,
                                                 float& median,
                                                 float& mode,
                                                 float& skewness,
                                                 float& rms)
{
  getWaveformParams<float>(waveform, mean, median, mode, skewness, rms);
  return;
}

void algorithms::WaveformUtils::getWaveformParams(const std::vector<double>& waveform,
                                                 float& mean,
                                                 float& median,
                                                 float& mode,
                                                 float& skewness,
                                                 float& rms)
{
  getWaveformParams<double>(waveform, mean, median, mode, skewness, rms);
  return;
}


template <typename T> void algorithms::WaveformUtils::getWaveformParams(const std::vector<T>& waveform,
                                                 float& mean,
                                                 float& median,
                                                 float& mode,
                                                 float& skewness,
                                                 float& rms)
{
  /*
  Calculate waveform parameters for a given 1D waveform.

  INPUTS:
    - waveform: 1D RawDigit Waveform.
  */
  typename std::vector<T> localWaveform = waveform;
  float realMean(float(std::accumulate(localWaveform.begin(),localWaveform.end(),0))/float(localWaveform.size()));
  if (localWaveform.size() % 2 == 0) {
    const auto m1 = localWaveform.begin() + localWaveform.size() / 2 - 1;
    const auto m2 = localWaveform.begin() + localWaveform.size() / 2;
    std::nth_element(localWaveform.begin(), m1, localWaveform.end());
    std::nth_element(localWaveform.begin(), m2, localWaveform.end());
    median = (*m1 + *m2) / 2.0;
  } else {
    median = localWaveform[localWaveform.size()/2];
  }
  mean   = realMean;
  std::vector<float> adcLessPedVec;
  adcLessPedVec.resize(localWaveform.size());
  std::transform(localWaveform.begin(),localWaveform.end(),adcLessPedVec.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));
  rms      = std::sqrt(std::inner_product(adcLessPedVec.begin(), adcLessPedVec.end(), adcLessPedVec.begin(), 0.) / float(adcLessPedVec.size()));
  skewness = 3. * float(realMean - median) / rms;
  return;
}



void algorithms::WaveformUtils::getDilation(const Waveform<short>& waveform,
                                            const unsigned int     structuringElement,
                                            Waveform<short>& dilationVec) const
{
  getDilation<short>(waveform, structuringElement, dilationVec);
  return;
}

void algorithms::WaveformUtils::getDilation(const Waveform<float>& waveform,
                                            const unsigned int     structuringElement,
                                            Waveform<float>& dilationVec) const
{
  getDilation<float>(waveform, structuringElement, dilationVec);
  return;
}


void algorithms::WaveformUtils::getDilation(const Waveform<double>& waveform,
                                            const unsigned int     structuringElement,
                                            Waveform<double>& dilationVec) const
{
  getDilation<double>(waveform, structuringElement, dilationVec);
  return;
}

template <typename T> void algorithms::WaveformUtils::getDilation(const Waveform<T>& inputWaveform,
                                                                  const unsigned int     structuringElement,
                                                                  Waveform<T>&       dilationVec) const
{
    // Set the window size
    int halfWindowSize(structuringElement/2);
    // Initialize min and max elements

    std::pair<typename Waveform<T>::const_iterator,typename Waveform<T>::const_iterator> minMaxItr =
            std::minmax_element(inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);

    typename Waveform<T>::const_iterator minElementItr = minMaxItr.first;
    typename Waveform<T>::const_iterator maxElementItr = minMaxItr.second;

    // Initialize the erosion and dilation vectors
    dilationVec.resize(inputWaveform.size());

    // Now loop through remaining elements and complete the vectors
    typename Waveform<T>::iterator maxItr = dilationVec.begin();

    for (typename Waveform<T>::const_iterator inputItr = inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
        {
            if (std::distance(minElementItr,inputItr) >= halfWindowSize)
                minElementItr = std::min_element(inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) < *minElementItr)
                minElementItr = inputItr + halfWindowSize;
            if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
                maxElementItr = std::max_element(inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) > *maxElementItr)
                maxElementItr = inputItr + halfWindowSize;
        }
        // Update the vectors
        *maxItr++ = *maxElementItr;
    }
    return;
}


void algorithms::WaveformUtils::getErosionDilationAverageDifference(const Waveform<short>& waveform,
                                                        const unsigned int     structuringElement,
                                                        Waveform<short>&       erosionVec,
                                                        Waveform<short>&       dilationVec,
                                                        Waveform<short>&       averageVec,
                                                        Waveform<short>&       differenceVec) const
{
    getErosionDilationAverageDifference<short>(waveform, structuringElement, erosionVec, dilationVec, averageVec, differenceVec);
    return;
}

void algorithms::WaveformUtils::getErosionDilationAverageDifference(const Waveform<float>& waveform,
                                                        const unsigned int     structuringElement,
                                                        Waveform<float>&       erosionVec,
                                                        Waveform<float>&       dilationVec,
                                                        Waveform<float>&       averageVec,
                                                        Waveform<float>&       differenceVec) const
{
    getErosionDilationAverageDifference<float>(waveform, structuringElement, erosionVec, dilationVec, averageVec, differenceVec);
    return;
}

void algorithms::WaveformUtils::getErosionDilationAverageDifference(const Waveform<double>& waveform,
                                                        const unsigned int     structuringElement,
                                                        Waveform<double>&       erosionVec,
                                                        Waveform<double>&       dilationVec,
                                                        Waveform<double>&       averageVec,
                                                        Waveform<double>&       differenceVec) const
{
    getErosionDilationAverageDifference<double>(waveform, structuringElement, erosionVec, dilationVec, averageVec, differenceVec);
    return;
}

template <typename T> void algorithms::WaveformUtils::getErosionDilationAverageDifference(const Waveform<T>& inputWaveform,
                                                                              const unsigned int     structuringElement,
                                                                              Waveform<T>&       erosionVec,
                                                                              Waveform<T>&       dilationVec,
                                                                              Waveform<T>&       averageVec,
                                                                              Waveform<T>&       differenceVec) const
{
    // Set the window size
    int halfWindowSize(structuringElement/2);
    // Initialize min and max elements

    std::pair<typename Waveform<T>::const_iterator,typename Waveform<T>::const_iterator> minMaxItr =
            std::minmax_element(inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);

    typename Waveform<T>::const_iterator minElementItr = minMaxItr.first;
    typename Waveform<T>::const_iterator maxElementItr = minMaxItr.second;

    // Initialize the erosion and dilation vectors
    erosionVec.resize(inputWaveform.size());
    dilationVec.resize(inputWaveform.size());
    averageVec.resize(inputWaveform.size());
    differenceVec.resize(inputWaveform.size());

    // Now loop through remaining elements and complete the vectors
    typename Waveform<T>::iterator minItr = erosionVec.begin();
    typename Waveform<T>::iterator maxItr = dilationVec.begin();
    typename Waveform<T>::iterator aveItr = averageVec.begin();
    typename Waveform<T>::iterator difItr = differenceVec.begin();

    for (typename Waveform<T>::const_iterator inputItr = inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
        {
            if (std::distance(minElementItr,inputItr) >= halfWindowSize)
                minElementItr = std::min_element(inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) < *minElementItr)
                minElementItr = inputItr + halfWindowSize;
            if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
                maxElementItr = std::max_element(inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) > *maxElementItr)
                maxElementItr = inputItr + halfWindowSize;
        }
        // Update the vectors
        *minItr++ = *minElementItr;
        *maxItr++ = *maxElementItr;
        *aveItr++ = 0.5 * (*maxElementItr + *minElementItr);
        *difItr++ = *maxElementItr - *minElementItr;
    }
    return;
}


void algorithms::WaveformUtils::getOpeningAndClosing(const Waveform<short>& erosionVec,
                                         const Waveform<short>& dilationVec,
                                         const unsigned int     structuringElement,
                                         Waveform<short>&       openingVec,
                                         Waveform<short>&       closingVec) const
{
    getOpeningAndClosing<short>(erosionVec, dilationVec, structuringElement, openingVec, closingVec);
    return;
}

void algorithms::WaveformUtils::getOpeningAndClosing(const Waveform<float>& erosionVec,
                                         const Waveform<float>& dilationVec,
                                         const unsigned int     structuringElement,
                                         Waveform<float>&       openingVec,
                                         Waveform<float>&       closingVec) const
{
    getOpeningAndClosing<float>(erosionVec, dilationVec, structuringElement, openingVec, closingVec);
    return;
}

void algorithms::WaveformUtils::getOpeningAndClosing(const Waveform<double>& erosionVec,
                                         const Waveform<double>& dilationVec,
                                         const unsigned int     structuringElement,
                                         Waveform<double>&       openingVec,
                                         Waveform<double>&       closingVec) const
{
    getOpeningAndClosing<double>(erosionVec, dilationVec, structuringElement, openingVec, closingVec);
    return;
}

template <typename T> void algorithms::WaveformUtils::getOpeningAndClosing(const Waveform<T>& erosionVec,
                                                               const Waveform<T>& dilationVec,
                                                               const unsigned int     structuringElement,
                                                               Waveform<T>&       openingVec,
                                                               Waveform<T>&       closingVec)  const
{
    // Set the window size
    int halfWindowSize(structuringElement/2);
    // Start with the opening, here we get the max element in the input erosion vector
    typename Waveform<T>::const_iterator maxElementItr = std::max_element(erosionVec.begin(),erosionVec.begin()+halfWindowSize);
    // Initialize the opening vector
    openingVec.resize(erosionVec.size());
    // Now loop through remaining elements and complete the vectors
    typename Waveform<T>::iterator maxItr = openingVec.begin();
    for (typename Waveform<T>::const_iterator inputItr = erosionVec.begin(); inputItr != erosionVec.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,erosionVec.end()) > halfWindowSize)
        {
            if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
                maxElementItr = std::max_element(inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) > *maxElementItr)
                maxElementItr = inputItr + halfWindowSize;
        }
        // Update the vectors
        *maxItr++ = *maxElementItr;
    }
    // Now go with the closling, here we get the min element in the input dilation vector
    typename Waveform<T>::const_iterator minElementItr = std::min_element(dilationVec.begin(),dilationVec.begin()+halfWindowSize);
    // Initialize the opening and closing vectors
    closingVec.resize(dilationVec.size());
    // Now loop through remaining elements and complete the vectors
    typename Waveform<T>::iterator minItr = closingVec.begin();
    for (typename Waveform<T>::const_iterator inputItr = dilationVec.begin(); inputItr != dilationVec.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,dilationVec.end()) > halfWindowSize)
        {
            if (std::distance(minElementItr,inputItr) >= halfWindowSize)
                minElementItr = std::min_element(inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) < *minElementItr)
                minElementItr = inputItr + halfWindowSize;
        }
        // Update the vectors
        *minItr++ = *minElementItr;
    }
    return;
}



void algorithms::WaveformUtils::getMorph2D(const std::vector<std::vector<float> >& waveforms,
                                            const unsigned int grouping,
                                            const unsigned int nTicks,
                                            const unsigned int structuringElementx,
                                            const unsigned int structuringElementy,
                                            std::vector<std::vector<float> >& dilation2D,
                                            std::vector<std::vector<float> >& erosion2D,
                                            std::vector<std::vector<float> >& gradient2D,
                                            const unsigned int window)
{
  WaveformUtils wUtils;
  auto numChannels = waveforms.size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  dilation2D.resize(waveforms.size());
  erosion2D.resize(waveforms.size());
  gradient2D.resize(waveforms.size());

  for (auto i=0; i<waveforms.size(); ++i) {
    dilation2D[i].resize(waveforms.at(0).size());
    erosion2D[i].resize(waveforms.at(0).size());
    gradient2D[i].resize(waveforms.at(0).size());
  }

  for (auto i=0; i<numChannels; ++i) {
    for (auto j=0; j<nTicks; ++j) {
      // For each center pixel, do 2D morphological filtering.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      int lowerBoundx = std::max(lbx, 0);
      int upperBoundx = std::min(ubx, (int) numChannels);
      int lowerBoundy = std::max(lby, 0);
      int upperBoundy = std::min(uby, (int) nTicks);
      std::vector<float> v;
      for (auto ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (auto iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v.push_back(waveforms[ix][iy]);
        }
      }
      float erosion = *std::min_element(v.begin(), v.end());
      float dilation = *std::max_element(v.begin(), v.end());
      float gradient = dilation - erosion;
      dilation2D[i][j] = dilation;
      erosion2D[i][j] = erosion;
      gradient2D[i][j] = gradient;
    }
  }
  return;
}




#endif   
