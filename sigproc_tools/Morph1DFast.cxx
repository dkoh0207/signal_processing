#ifndef __SIGPROC_TOOLS_MORPH1DFAST_CXX__
#define __SIGPROC_TOOLS_MORPH1DFAST_CXX__

#include "Morph1DFast.h"



void sigproc_tools::Morph1DFast::getDilation(
  const Waveform<bool>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<bool>& dilationVec) const
{
  size_t N = inputWaveform.size();
  size_t k = (size_t) structuringElement;

  assert(dilationVec.size() == N);

  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  size_t bufferSize = N + 2 * (k/2) + (k - (N % k));
  size_t windowSize = k/2;
  size_t paddingSize = k - (N % k);
  std::vector<bool> suffixArr(bufferSize);
  std::vector<bool> prefixArr(bufferSize);

  // Padding Operations on Buffers
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = false;
    prefixArr[i] = false;
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = false;
    prefixArr[i] = false;
  }

  // Compute Prefix and Suffix Buffers
  for (size_t i=0; i<N+paddingSize; ++i) {
    if (i % k == 0) {
      prefixArr[i+windowSize] = inputWaveform[i];
    } else {
      prefixArr[i+windowSize] = (prefixArr[i+windowSize-1] || inputWaveform[i]);
    }
  }

  for (size_t i=N+paddingSize; i!=0; --i) {
    if (i > N) {
      // Compensate for divisibility padding (must be -inf)
      continue;
    }
    else if (i % k == 0) {
      suffixArr[i+windowSize-1] = inputWaveform[i-1];
    } 
    else {
      suffixArr[i+windowSize-1] = (suffixArr[i+windowSize] || inputWaveform[i-1]);
    }
  }

  int prefixIndex = 0;
  int suffixIndex = 0;

  for (size_t i=0; i<N; ++i) {
    prefixIndex = i + 2 * windowSize;
    suffixIndex = i;
    dilationVec[i] = (prefixArr[prefixIndex] || suffixArr[suffixIndex]);
  }
  return;
}

void sigproc_tools::Morph1DFast::getDilation(
  const Waveform<short>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<short>& dilationVec) const
{
  getDilation<short>(inputWaveform, structuringElement, dilationVec);
  return;
}

void sigproc_tools::Morph1DFast::getDilation(
  const Waveform<float>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<float>& dilationVec) const
{
  getDilation<float>(inputWaveform, structuringElement, dilationVec);
  return;
}

void sigproc_tools::Morph1DFast::getDilation(
  const Waveform<double>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<double>& dilationVec) const
{
  getDilation<double>(inputWaveform, structuringElement, dilationVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getDilation(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& dilationVec) const
{
  size_t N = inputWaveform.size();
  assert(dilationVec.size() == N);
  size_t k = (size_t) structuringElement;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  size_t bufferSize = N + 2 * (k/2) + (k - (N % k));
  size_t windowSize = k/2;
  size_t paddingSize = k - (N % k);
  std::vector<T> suffixArr(bufferSize);
  std::vector<T> prefixArr(bufferSize);

  // Padding Operations on Buffers
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::min();
    prefixArr[i] = std::numeric_limits<T>::min();
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::min();
    prefixArr[i] = std::numeric_limits<T>::min();
  }

  // Compute Prefix and Suffix Buffers
  for (size_t i=0; i<N+paddingSize; ++i) {
    if (i % k == 0) {
      prefixArr[i+windowSize] = inputWaveform[i];
    } else {
      prefixArr[i+windowSize] = std::max(prefixArr[i+windowSize-1], inputWaveform[i]);
    }
  }

  for (size_t i=N+paddingSize; i!=0; --i) {
    if (i > N) {
      // Compensate for divisibility padding (must be -inf)
      continue;
    }
    else if (i % k == 0) {
      suffixArr[i+windowSize-1] = inputWaveform[i-1];
    } 
    else {
      suffixArr[i+windowSize-1] = std::max(suffixArr[i+windowSize], inputWaveform[i-1]);
    }
  }

  int prefixIndex = 0;
  int suffixIndex = 0;

  for (size_t i=0; i<N; ++i) {
    prefixIndex = i + 2 * windowSize;
    suffixIndex = i;
    dilationVec[i] = std::max(prefixArr[prefixIndex],
      suffixArr[suffixIndex]);
  }
  return;
}


void sigproc_tools::Morph1DFast::getErosion(
  const Waveform<bool>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<bool>& erosionVec) const
{
  size_t N = inputWaveform.size();
  assert(erosionVec.size() == N);
  size_t k = (size_t) structuringElement;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  size_t bufferSize = N + 2 * (k/2) + (k - (N % k));
  size_t windowSize = k/2;
  size_t paddingSize = k - (N % k);
  std::vector<bool> suffixArr(bufferSize);
  std::vector<bool> prefixArr(bufferSize);

  // Padding Operations on Buffers
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = true;
    prefixArr[i] = true;
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = true;
    prefixArr[i] = true;
  }

  // Compute Prefix and Suffix Buffers
  for (size_t i=0; i<N+paddingSize; ++i) {
    if (i % k == 0) {
      prefixArr[i+windowSize] = inputWaveform[i];
    } else {
      prefixArr[i+windowSize] = (prefixArr[i+windowSize-1] && inputWaveform[i]);
    }
  }

  for (size_t i=N+paddingSize; i!=0; --i) {
    if (i > N) {
      // Compensate for divisibility padding (must be -inf)
      continue;
    }
    else if (i % k == 0) {
      suffixArr[i+windowSize-1] = inputWaveform[i-1];
    } 
    else {
      suffixArr[i+windowSize-1] = (suffixArr[i+windowSize] && inputWaveform[i-1]);
    }
  }

  int prefixIndex = 0;
  int suffixIndex = 0;

  for (size_t i=0; i<N; ++i) {
    prefixIndex = i + 2 * windowSize;
    suffixIndex = i;
    erosionVec[i] = (prefixArr[prefixIndex] && suffixArr[suffixIndex]);
  }
  return;
}

void sigproc_tools::Morph1DFast::getErosion(
  const Waveform<short>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<short>& erosionVec) const
{
  getErosion<short>(inputWaveform, structuringElement, erosionVec);
  return;
}

void sigproc_tools::Morph1DFast::getErosion(
  const Waveform<float>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<float>& erosionVec) const
{
  getErosion<float>(inputWaveform, structuringElement, erosionVec);
  return;
}

void sigproc_tools::Morph1DFast::getErosion(
  const Waveform<double>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<double>& erosionVec) const
{
  getErosion<double>(inputWaveform, structuringElement, erosionVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getErosion(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& erosionVec) const
{
  size_t N = inputWaveform.size();
  assert(erosionVec.size() == N);
  size_t k = (size_t) structuringElement;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  size_t bufferSize = N + 2 * (k/2) + (k - (N % k));
  size_t windowSize = k/2;
  size_t paddingSize = k - (N % k);
  std::vector<T> suffixArr(bufferSize);
  std::vector<T> prefixArr(bufferSize);

  // Padding Operations on Buffers
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::max();
    prefixArr[i] = std::numeric_limits<T>::max();
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::max();
    prefixArr[i] = std::numeric_limits<T>::max();
  }

  // Compute Prefix and Suffix Buffers
  for (size_t i=0; i<N+paddingSize; ++i) {
    if (i % k == 0) {
      prefixArr[i+windowSize] = inputWaveform[i];
    } else {
      prefixArr[i+windowSize] = std::min(prefixArr[i+windowSize-1], inputWaveform[i]);
    }
  }

  for (size_t i=N+paddingSize; i!=0; --i) {
    if (i > N) {
      // Compensate for divisibility padding (must be -inf)
      continue;
    }
    else if (i % k == 0) {
      suffixArr[i+windowSize-1] = inputWaveform[i-1];
    } 
    else {
      suffixArr[i+windowSize-1] = std::min(suffixArr[i+windowSize], inputWaveform[i-1]);
    }
  }

  int prefixIndex = 0;
  int suffixIndex = 0;

  for (size_t i=0; i<N; ++i) {
    prefixIndex = i + 2 * windowSize;
    suffixIndex = i;
    erosionVec[i] = std::min(prefixArr[prefixIndex],
      suffixArr[suffixIndex]);
  }
  return;
}


void sigproc_tools::Morph1DFast::getGradient(
  const Waveform<short>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<short>& gradientVec) const
{
  getGradient<short>(inputWaveform, structuringElement, gradientVec);
  return;
}

void sigproc_tools::Morph1DFast::getGradient(
  const Waveform<float>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<float>& gradientVec) const
{
  getGradient<float>(inputWaveform, structuringElement, gradientVec);
  return;
}

void sigproc_tools::Morph1DFast::getGradient(
  const Waveform<double>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<double>& gradientVec) const
{
  getGradient<double>(inputWaveform, structuringElement, gradientVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getGradient(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& gradientVec) const
{
  size_t N = inputWaveform.size();
  assert(gradientVec.size() == N);
  std::vector<T> erosionVec(N);
  std::vector<T> dilationVec(N);
  getDilation<T>(inputWaveform, structuringElement, dilationVec);
  getErosion<T>(inputWaveform, structuringElement, erosionVec);

  for (size_t i=0; i<N; ++i) {
    gradientVec[i] = dilationVec[i] - erosionVec[i];
  }
  return;
}


void sigproc_tools::Morph1DFast::getOpening(
  const Waveform<bool>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<bool>& openingVec) const
{
  size_t N = inputWaveform.size();
  assert(openingVec.size() == N);
  std::vector<bool> tempVec(N);
  getErosion(inputWaveform, structuringElement, tempVec);
  getDilation(tempVec, structuringElement, openingVec);
  return;
}

void sigproc_tools::Morph1DFast::getOpening(
  const Waveform<short>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<short>& openingVec) const
{
  getOpening<short>(inputWaveform, structuringElement, openingVec);
  return;
}

void sigproc_tools::Morph1DFast::getOpening(
  const Waveform<float>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<float>& openingVec) const
{
  getOpening<float>(inputWaveform, structuringElement, openingVec);
  return;
}

void sigproc_tools::Morph1DFast::getOpening(
  const Waveform<double>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<double>& openingVec) const
{
  getOpening<double>(inputWaveform, structuringElement, openingVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getOpening(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& openingVec) const
{
  size_t N = inputWaveform.size();
  assert(openingVec.size() == N);
  std::vector<T> tempVec(N);
  getErosion<T>(inputWaveform, structuringElement, tempVec);
  getDilation<T>(tempVec, structuringElement, openingVec);
  return;
}


void sigproc_tools::Morph1DFast::getClosing(
  const Waveform<bool>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<bool>& closingVec) const
{
  size_t N = inputWaveform.size();
  assert(closingVec.size() == N);
  std::vector<bool> tempVec(N);
  getDilation(inputWaveform, structuringElement, tempVec);
  getErosion(tempVec, structuringElement, closingVec);
  return;
}

void sigproc_tools::Morph1DFast::getClosing(
  const Waveform<short>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<short>& closingVec) const
{
  getClosing<short>(inputWaveform, structuringElement, closingVec);
  return;
}

void sigproc_tools::Morph1DFast::getClosing(
  const Waveform<float>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<float>& closingVec) const
{
  getClosing<float>(inputWaveform, structuringElement, closingVec);
  return;
}

void sigproc_tools::Morph1DFast::getClosing(
  const Waveform<double>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<double>& closingVec) const
{
  getClosing<double>(inputWaveform, structuringElement, closingVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getClosing(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& closingVec) const
{
  size_t N = inputWaveform.size();
  assert(closingVec.size() == N);
  std::vector<T> tempVec(N);
  getDilation<T>(inputWaveform, structuringElement, tempVec);
  getErosion<T>(tempVec, structuringElement, closingVec);
  return;
}

#endif
