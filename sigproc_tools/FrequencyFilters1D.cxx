#ifndef __SIGPROC_TOOLS_FREQUENCYFILTERS_CXX__
#define __SIGPROC_TOOLS_FREQUENCYFILTERS_CXX__

#include "FrequencyFilters1D.h"

sigproc_tools::FrequencyFilters1D::FrequencyFilters1D(
  const int numTimeSteps)
{
  fTimeVec.resize(numTimeSteps, 0.);
  fFrequencyVec.resize(numTimeSteps, std::complex<float>(0., 0.));
  
  fForwardPlan = fftwf_plan_dft_r2c_1d(
    numTimeSteps, 
    fTimeVec.data(), 
    reinterpret_cast<fftwf_complex*>(fFrequencyVec.data()), 
    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

  fInversePlan = fftwf_plan_dft_c2r_1d(
    numTimeSteps, 
    reinterpret_cast<fftwf_complex*>(fFrequencyVec.data()), 
    fTimeVec.data(), 
    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

  return;
}


sigproc_tools::FrequencyFilters1D::~FrequencyFilters1D()
{
  fftwf_destroy_plan(fForwardPlan);
  fftwf_destroy_plan(fInversePlan);

  fTimeVec.clear();
  fFrequencyVec.clear();

  return;
}


void sigproc_tools::FrequencyFilters1D::lowPassButterWorth1D(
  const Waveform<float>& inputWaveform1D,
  const size_t threshold,
  const unsigned int order,
  Waveform<float>& outputWaveform1D)
{
  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    fTimeVec[i] = inputWaveform1D[i];
  }

  forwardFFT(fTimeVec, fFrequencyVec);

  size_t freqVecSize = fFrequencyVec.size();
  size_t timeVecSize = fTimeVec.size();

  if (threshold < 1 or threshold > freqVecSize / 2) 
    throw std::runtime_error(
      "FrequencyFilters: Thresholding frequency bin must be at least one and less than half of total bin count.");

  std::vector<float> filterFunc;
  filterFunc.resize(freqVecSize);

  float d0 = (float) threshold;

  for (size_t i=0; i<filterFunc.size()/2; ++i) {
    float filterValue;
    float d = (float) i;
    filterValue = 1.0 / (1.0 + std::pow( d / d0, 2 * order));
    filterFunc[i] = filterValue;
    filterFunc[freqVecSize-i-1] = filterValue;
  }

  // Ad Hoc compensation for odd length frequency arrays.
  filterFunc[freqVecSize/2] = filterFunc[freqVecSize/2-1];

  outputWaveform1D.resize(timeVecSize);

  for (size_t i=0; i<threshold; ++i) {
    fFrequencyVec[i] = fFrequencyVec[i] * filterFunc[i];
  }

  inverseFFT(fFrequencyVec, fTimeVec);

  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    outputWaveform1D[i] = fTimeVec[i];
  }

  return;
}


void sigproc_tools::FrequencyFilters1D::lowPassGaussian1D(
  const Waveform<float>& inputWaveform1D,
  const size_t threshold,
  Waveform<float>& outputWaveform1D)
{
  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    fTimeVec[i] = inputWaveform1D[i];
  }

  forwardFFT(fTimeVec, fFrequencyVec);

  size_t freqVecSize = fFrequencyVec.size();
  size_t timeVecSize = fTimeVec.size();

  if (threshold < 1 or threshold > freqVecSize / 2) 
    throw std::runtime_error(
      "FrequencyFilters: Thresholding frequency bin must be at least one and less than half of total bin count.");

  std::vector<float> filterFunc;
  filterFunc.resize(freqVecSize);

  float d0 = (float) threshold;

  for (size_t i=0; i<filterFunc.size()/2; ++i) {
    float filterValue;
    float d = (float) i;
    filterValue = std::exp(-1.0 * std::pow(d, 2) / (2.0 * std::pow(d0, 2)));
    filterFunc[i] = filterValue;
    filterFunc[freqVecSize-i-1] = filterValue;
  }

  // Ad Hoc compensation for odd length frequency arrays.
  filterFunc[freqVecSize/2] = filterFunc[freqVecSize/2-1];

  outputWaveform1D.resize(timeVecSize);

  for (size_t i=0; i<threshold; ++i) {
    fFrequencyVec[i] = fFrequencyVec[i] * filterFunc[i];
  }

  inverseFFT(fFrequencyVec, fTimeVec);

  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    outputWaveform1D[i] = fTimeVec[i];
  }

  return;
}


void sigproc_tools::FrequencyFilters1D::highPassHard1D(
  const Waveform<float>& inputWaveform1D,
  const size_t threshold,
  Waveform<float>& outputWaveform1D)
{
  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    fTimeVec[i] = inputWaveform1D[i];
  }

  forwardFFT(fTimeVec, fFrequencyVec);

  size_t freqVecSize = fFrequencyVec.size();
  size_t timeVecSize = fTimeVec.size();

  outputWaveform1D.resize(timeVecSize);

  if (threshold < 1 or threshold > freqVecSize / 2) 
    throw std::runtime_error(
      "FrequencyFilters: Thresholding frequency bin must be at least one and less than half of total bin count.");

  for (size_t i=0; i<threshold; ++i) {
    fFrequencyVec[i] = std::complex<float>(0., 0.);
    fFrequencyVec[freqVecSize-i] = std::complex<float>(0., 0.);
  }

  inverseFFT(fFrequencyVec, fTimeVec);

  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    outputWaveform1D[i] = fTimeVec[i];
  }

  return;
}


void sigproc_tools::FrequencyFilters1D::highPassButterWorth1D(
  const Waveform<float>& inputWaveform1D,
  const size_t threshold,
  const unsigned int order,
  Waveform<float>& outputWaveform1D)
{
  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    fTimeVec[i] = inputWaveform1D[i];
  }

  forwardFFT(fTimeVec, fFrequencyVec);

  size_t freqVecSize = fFrequencyVec.size();
  size_t timeVecSize = fTimeVec.size();

  if (threshold < 1 or threshold > freqVecSize / 2) 
    throw std::runtime_error(
      "FrequencyFilters: Thresholding frequency bin must be at least one and less than half of total bin count.");

  std::vector<float> filterFunc;
  filterFunc.resize(freqVecSize);

  float d0 = (float) threshold;

  for (size_t i=0; i<filterFunc.size()/2; ++i) {
    float filterValue;
    float d = (float) i;
    filterValue = 1.0 - 1.0 / (1.0 + std::pow( d / d0, 2 * order));
    filterFunc[i] = filterValue;
    filterFunc[freqVecSize-i-1] = filterValue;
  }

  // Ad Hoc compensation for odd length frequency arrays.
  filterFunc[freqVecSize/2] = filterFunc[freqVecSize/2-1];

  outputWaveform1D.resize(timeVecSize);

  for (size_t i=0; i<threshold; ++i) {
    fFrequencyVec[i] = fFrequencyVec[i] * filterFunc[i];
  }

  inverseFFT(fFrequencyVec, fTimeVec);

  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    outputWaveform1D[i] = fTimeVec[i];
  }

  return;
}


void sigproc_tools::FrequencyFilters1D::highPassGaussian1D(
  const Waveform<float>& inputWaveform1D,
  const size_t threshold,
  Waveform<float>& outputWaveform1D)
{
  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    fTimeVec[i] = inputWaveform1D[i];
  }

  forwardFFT(fTimeVec, fFrequencyVec);

  size_t freqVecSize = fFrequencyVec.size();
  size_t timeVecSize = fTimeVec.size();

  if (threshold < 1 or threshold > freqVecSize / 2) 
    throw std::runtime_error(
      "FrequencyFilters: Thresholding frequency bin must be at least one and less than half of total bin count.");

  std::vector<float> filterFunc;
  filterFunc.resize(freqVecSize);

  float d0 = (float) threshold;

  for (size_t i=0; i<filterFunc.size()/2; ++i) {
    float filterValue;
    float d = (float) i;
    filterValue = 1.0 - std::exp(-1.0 * std::pow(d, 2) / (2.0 * std::pow(d0, 2)));
    filterFunc[i] = filterValue;
    filterFunc[freqVecSize-i-1] = filterValue;
  }

  // Ad Hoc compensation for odd length frequency arrays.
  filterFunc[freqVecSize/2] = filterFunc[freqVecSize/2-1];

  outputWaveform1D.resize(timeVecSize);

  for (size_t i=0; i<threshold; ++i) {
    fFrequencyVec[i] = fFrequencyVec[i] * filterFunc[i];
  }

  inverseFFT(fFrequencyVec, fTimeVec);

  for (size_t i=0; i<inputWaveform1D.size(); ++i) {
    outputWaveform1D[i] = fTimeVec[i];
  }

  return;
}


void sigproc_tools::FrequencyFilters1D::forwardFFT(
  TimeVec& timeVec, 
  FrequencyVec& frequencyVec) const
{

    if (timeVec.size() != fTimeVec.size()) 
      throw std::runtime_error(
        "ICARUSFFT: Input time vector size does not match expected");

    frequencyVec.resize(timeVec.size());
    fftwf_execute_dft_r2c(
      fForwardPlan, 
      timeVec.data(), 
      reinterpret_cast<fftwf_complex*>(frequencyVec.data()));
    // Reflect the output frequency vector

    size_t vecSize    = timeVec.size();
    size_t nyquistBin = vecSize/2 + 1;

    for(size_t idx = nyquistBin; idx < timeVec.size(); idx++)
      frequencyVec[idx] = std::conj(frequencyVec[vecSize - idx]);

    return;

}


void sigproc_tools::FrequencyFilters1D::inverseFFT(
  FrequencyVec& frequencyVec, 
  TimeVec& timeVec) const
{
    if (frequencyVec.size() < fFrequencyVec.size()) 
        throw std::runtime_error(
          "ICARUSFFT: Input frequency vector size does not match expected");

    timeVec.resize(frequencyVec.size());

    fftwf_execute_dft_c2r(
      fInversePlan, 
      reinterpret_cast<fftwf_complex*>(frequencyVec.data()), 
      timeVec.data());

    // Now normalize

    float normFactor = 1. / float(timeVec.size());
    std::transform(
      timeVec.begin(),
      timeVec.end(),
      timeVec.begin(),
      std::bind1st(std::multiplies<float>(),
      normFactor));

    return;
}


void sigproc_tools::FrequencyFilters1D::filterImage(
  const Waveform2D<float>& inputWaveform2D,
  const size_t threshold,
  Waveform2D<float>& outputWaveform2D,
  const unsigned int order,
  const int mode)
{
  size_t numChannels = inputWaveform2D.size();
  size_t nTicks = inputWaveform2D.at(0).size();

  outputWaveform2D.resize(numChannels);
  for (auto& v : outputWaveform2D) {
    v.resize(nTicks);
  }

  switch (mode) {
    case 0:
      for (size_t i=0; i<numChannels; ++i) {
        highPassButterWorth1D(inputWaveform2D[i], threshold, order, outputWaveform2D[i]);
      }
      break;
    case 1:
      for (size_t i=0; i<numChannels; ++i) {
        highPassHard1D(inputWaveform2D[i], threshold, outputWaveform2D[i]);
      }
      break;
    case 2:
      for (size_t i=0; i<numChannels; ++i) {
        highPassGaussian1D(inputWaveform2D[i], threshold, outputWaveform2D[i]);
      }
      break;
    case 3:
      for (size_t i=0; i<numChannels; ++i) {
        lowPassButterWorth1D(inputWaveform2D[i], threshold, order, outputWaveform2D[i]);
      }
      break;
    case 4:
      for (size_t i=0; i<numChannels; ++i) {
        lowPassGaussian1D(inputWaveform2D[i], threshold, outputWaveform2D[i]);
      }
      break;
    default:
      for (size_t i=0; i<numChannels; ++i) {
        highPassButterWorth1D(inputWaveform2D[i], threshold, order, outputWaveform2D[i]);
      }
      break;
  }
  return;
}




#endif
