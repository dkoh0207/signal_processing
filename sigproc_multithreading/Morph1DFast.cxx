#ifndef __SIGPROC_MULTITHREADING_MORPH1DFAST_CXX__
#define __SIGPROC_MULTITHREADING_MORPH1DFAST_CXX__

#include "Morph1DFast.h"


void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentVector<bool>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<bool>& dilationVec) const
{
  const size_t N = inputVector.size();
  const size_t k = (size_t) structuringElement;

  assert(dilationVec.size() == N);
  
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  const size_t windowSize = k/2;
  const size_t paddingSize = (k - (N % k)) % k;
  const size_t bufferSize = N + 2 * windowSize + paddingSize;
  ConcurrentVector<bool> suffixArr(bufferSize);
  ConcurrentVector<bool> prefixArr(bufferSize);

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
  tbb::parallel_invoke(
    [&prefixArr, &inputVector, &windowSize, &paddingSize, &N, &k]() {
      for (size_t i=0; i<N+paddingSize; ++i) {
        if (i % k == 0) {
          prefixArr[i+windowSize] = inputVector[i];
        } else if ((i % k == 0) && (i < N)) {
          prefixArr[i+windowSize] = (prefixArr[i+windowSize-1] || inputVector[i]);
        } else {
          continue;
        }
      }
    },
    [&suffixArr, &inputVector, &windowSize, &paddingSize, &N, &k]() {
      for (size_t i=N+paddingSize; i!=0; --i) {
        if (i > N) {
          // Compensate for divisibility padding (must be -inf)
          continue;
        }
        else if (i % k == 0) {
          suffixArr[i+windowSize-1] = inputVector[i-1];
        } 
        else {
          suffixArr[i+windowSize-1] = (suffixArr[i+windowSize] || inputVector[i-1]);
        }
      }
    }
  );

  size_t prefixIndex = 0;
  size_t suffixIndex = 0;

  tbb::parallel_for( (size_t) windowSize, N+windowSize, (size_t) 1, 
    [&prefixIndex, &suffixIndex, &dilationVec, 
     &prefixArr, &suffixArr, &windowSize](size_t i) 
    {
      prefixIndex = i + windowSize;
      suffixIndex = i - windowSize;
      dilationVec[i-windowSize] = (prefixArr[prefixIndex] || suffixArr[suffixIndex]);
    }
  );
  return;
}

void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentArray2D<bool>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<bool>& dilation2D,
  const unsigned int columnNum) const
{
  const size_t N = inputArray2D.size();
  const size_t k = (size_t) structuringElementx;
  assert(columnNum < inputArray2D.at(0).size());
  assert(dilation2D.size() == N);

  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  const size_t windowSize = k/2;
  const size_t paddingSize = (k - (N % k)) % k;
  const size_t bufferSize = N + 2 * windowSize + paddingSize;
  ConcurrentVector<bool> suffixArr(bufferSize);
  ConcurrentVector<bool> prefixArr(bufferSize);

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
  tbb::parallel_invoke(
    [&prefixArr, &inputArray2D, &windowSize, &paddingSize, &columnNum, &N, &k]() {
      for (size_t i=0; i<N+paddingSize; ++i) {
        if (i % k == 0) {
          prefixArr[i+windowSize] = inputArray2D[i][columnNum];
        } else {
          prefixArr[i+windowSize] = (prefixArr[i+windowSize-1] || inputArray2D[i][columnNum]);
        }
      }
    },
    [&suffixArr, &inputArray2D, &windowSize, &paddingSize, &columnNum, &N, &k]() {
      for (size_t i=N+paddingSize; i!=0; --i) {
        if (i > N) {
          // Compensate for divisibility padding (must be -inf)
          continue;
        }
        else if (i % k == 0) {
          suffixArr[i+windowSize-1] = inputArray2D[i-1][columnNum];
        } 
        else {
          suffixArr[i+windowSize-1] = (suffixArr[i+windowSize] || inputArray2D[i-1][columnNum]);
        }
      }
    }
  );

  size_t prefixIndex = 0;
  size_t suffixIndex = 0;

  tbb::parallel_for( (size_t) windowSize, N+windowSize, (size_t) 1, 
    [&prefixIndex, &suffixIndex, &dilation2D, 
     &prefixArr, &suffixArr, &windowSize, &columnNum](size_t i) 
    {
      prefixIndex = i + windowSize;
      suffixIndex = i - windowSize;
      dilation2D[i-windowSize][columnNum] = (prefixArr[prefixIndex] || suffixArr[suffixIndex]);
    }
  );
  return;
}

void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentVector<short>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<short>& dilationVec) const
{
  getDilation<short>(inputVector, structuringElement, dilationVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentVector<float>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<float>& dilationVec) const
{
  getDilation<float>(inputVector, structuringElement, dilationVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentVector<double>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<double>& dilationVec) const
{
  getDilation<double>(inputVector, structuringElement, dilationVec);
  return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentVector<T>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<T>& dilationVec) const
{
  const size_t N = inputVector.size();
  assert(dilationVec.size() == N);
  const size_t k = (size_t) structuringElement;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  const size_t windowSize = k/2;
  const size_t paddingSize = (k - (N % k)) % k;
  const size_t bufferSize = N + 2 * windowSize + paddingSize;
  ConcurrentVector<T> suffixArr(bufferSize);
  ConcurrentVector<T> prefixArr(bufferSize);

  // Padding Operations on Buffers
  // This could be done with parallel_for and parallel_invoke, yet as they run
  // through not many entries anyway it may not be worth it. 
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::min();
    prefixArr[i] = std::numeric_limits<T>::min();
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::min();
    prefixArr[i] = std::numeric_limits<T>::min();
  }

  // Prefix computation is independent from suffix computation,
  // so we use tbb::parallel_invoke. 
  // However, each prefix/suffix computation loop must respect the order of
  // element, so we cannot use tbb::parallel_for here. 

  tbb::parallel_invoke(
    [&suffixArr, &inputVector, &windowSize, &paddingSize, &N, &k]() {
      for (size_t i=N+paddingSize; i!=0; --i) {
        if (i > N) {
          // Compensate for divisibility padding (must be -inf)
          continue;
        }
        else if (i % k == 0) {
          suffixArr[i+windowSize-1] = inputVector[i-1];
        } 
        else {
          suffixArr[i+windowSize-1] = std::max(suffixArr[i+windowSize], 
                                               inputVector[i-1]);
        }
      }
    },
    [&prefixArr, &inputVector, &windowSize, &paddingSize, &N, &k]() {
      for (size_t i=0; i<N+paddingSize; ++i) {
        if (i % k == 0) {
          prefixArr[i+windowSize] = inputVector[i];
        }
        else if ((i % k == 0) && (i < N)) {
          prefixArr[i+windowSize] = std::max(prefixArr[i+windowSize-1], 
                                             inputVector[i]);
        }
        else {
          continue;
        }
      }
    }
  );

  // Compute Prefix and Suffix Buffers

  size_t prefixIndex = 0;
  size_t suffixIndex = 0;

  tbb::parallel_for((size_t) windowSize, N+windowSize, (size_t) 1, 
    [&prefixIndex, &suffixIndex, &dilationVec, 
     &prefixArr, &suffixArr, &windowSize](size_t i) {
      prefixIndex = i + windowSize;
      suffixIndex = i - windowSize;
      dilationVec[i-windowSize] = std::max(prefixArr[prefixIndex],
        suffixArr[suffixIndex]);
    }
  );
  return;
}






void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentVector<bool>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<bool>& erosionVec) const
{
  const size_t N = inputVector.size();
  assert(erosionVec.size() == N);
  const size_t k = (size_t) structuringElement;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  const size_t windowSize = k/2;
  const size_t paddingSize = (k - (N % k)) % k;
  const size_t bufferSize = N + 2 * windowSize + paddingSize;
  ConcurrentVector<bool> suffixArr(bufferSize);
  ConcurrentVector<bool> prefixArr(bufferSize);

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
  tbb::parallel_invoke(
    [&prefixArr, &inputVector, &windowSize, &paddingSize, &N, &k]() {
      for (size_t i=0; i<N+paddingSize; ++i) {
        if (i % k == 0) {
          prefixArr[i+windowSize] = inputVector[i];
        } 
        else if ((i % k == 0) && (i < N)) {
          prefixArr[i+windowSize] = (prefixArr[i+windowSize-1] && inputVector[i]);
        }
        else {
          continue;
        }
      }
    },
    [&suffixArr, &inputVector, &windowSize, &paddingSize, &N, &k]() {
      for (size_t i=N+paddingSize; i!=0; --i) {
        if (i > N) {
          // Compensate for divisibility padding (must be -inf)
          continue;
        }
        else if (i % k == 0) {
          suffixArr[i+windowSize-1] = inputVector[i-1];
        } 
        else {
          suffixArr[i+windowSize-1] = (suffixArr[i+windowSize] && inputVector[i-1]);
        }
      }
    }
  );

  int prefixIndex = 0;
  int suffixIndex = 0;

  tbb::parallel_for( (size_t) windowSize, N+windowSize, (size_t) 1, 
    [&prefixIndex, &suffixIndex, &erosionVec, 
     &prefixArr, &suffixArr, &windowSize](size_t i) 
    {
      prefixIndex = i + windowSize;
      suffixIndex = i - windowSize;
      erosionVec[i-windowSize] = (prefixArr[prefixIndex] && suffixArr[suffixIndex]);
    }
  );
  return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<bool>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<bool>& erosion2D,
  const unsigned int columnNum) const
{
  const size_t N = inputArray2D.size();
  assert(erosion2D.size() == N);
  assert(columnNum < erosion2D.at(0).size());
  const size_t k = (size_t) structuringElementx;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  const size_t windowSize = k/2;
  const size_t paddingSize = (k - (N % k)) % k;
  const size_t bufferSize = N + 2 * windowSize + paddingSize;
  ConcurrentVector<bool> suffixArr(bufferSize);
  ConcurrentVector<bool> prefixArr(bufferSize);

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
  tbb::parallel_invoke(
    [&prefixArr, &inputArray2D, &windowSize, &paddingSize, &columnNum, &N, &k]() {
      for (size_t i=0; i<N+paddingSize; ++i) {
        if (i % k == 0) {
          prefixArr[i+windowSize] = inputArray2D[i][columnNum];
        }
        else if ((i % k == 0) && (i < N)) {
          prefixArr[i+windowSize] = (prefixArr[i+windowSize-1] && inputArray2D[i][columnNum]);
        }
        else{
          continue;
        }
      }
    },
    [&suffixArr, &inputArray2D, &windowSize, &paddingSize, &columnNum, &N, &k]() {
      for (size_t i=N+paddingSize; i!=0; --i) {
        if (i > N) {
          // Compensate for divisibility padding (must be -inf)
          continue;
        }
        else if (i % k == 0) {
          suffixArr[i+windowSize-1] = inputArray2D[i-1][columnNum];
        } 
        else {
          suffixArr[i+windowSize-1] = (suffixArr[i+windowSize] && inputArray2D[i-1][columnNum]);
        }
      }
    }
  );

  int prefixIndex = 0;
  int suffixIndex = 0;

  tbb::parallel_for( (size_t) windowSize, N+windowSize, (size_t) 1, 
    [&prefixIndex, &suffixIndex, &erosion2D, 
     &prefixArr, &suffixArr, &windowSize, &columnNum](size_t i) 
    {
      prefixIndex = i + windowSize;
      suffixIndex = i - windowSize;
      erosion2D[i-windowSize][columnNum] = (prefixArr[prefixIndex] && suffixArr[suffixIndex]);
    }
  );
  return;
}


void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentVector<short>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<short>& erosionVec) const
{
  getErosion<short>(inputVector, structuringElement, erosionVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentVector<float>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<float>& erosionVec) const
{
  getErosion<float>(inputVector, structuringElement, erosionVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentVector<double>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<double>& erosionVec) const
{
  getErosion<double>(inputVector, structuringElement, erosionVec);
  return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentVector<T>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<T>& erosionVec) const
{
  const size_t N = inputVector.size();
  assert(erosionVec.size() == N);
  const size_t k = (size_t) structuringElement;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  const size_t windowSize = k/2;
  const size_t paddingSize = (k - (N % k)) % k;
  const size_t bufferSize = N + 2 * windowSize + paddingSize;
  ConcurrentVector<T> suffixArr(bufferSize);
  ConcurrentVector<T> prefixArr(bufferSize);

  // Padding Operations on Buffers
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::max();
    prefixArr[i] = std::numeric_limits<T>::max();
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::max();
    prefixArr[i] = std::numeric_limits<T>::max();
  }

  tbb::parallel_invoke(
    [&suffixArr, &inputVector, &windowSize, &paddingSize, &N, &k]() {
      for (size_t i=N+paddingSize; i!=0; --i) {
        if (i > N) {
          // Compensate for divisibility padding (must be -inf)
          continue;
        }
        else if (i % k == 0) {
          suffixArr[i+windowSize-1] = inputVector[i-1];
        } 
        else {
          suffixArr[i+windowSize-1] = std::min(suffixArr[i+windowSize], 
                                               inputVector[i-1]);
        }
      }
    },
    [&prefixArr, &inputVector, &windowSize, &paddingSize, &N, &k]() {
      for (size_t i=0; i<N+paddingSize; ++i) {
        if (i % k == 0) {
          prefixArr[i+windowSize] = inputVector[i];
        }
        else if ((i % k == 0) && (i < N)) {
          prefixArr[i+windowSize] = std::min(prefixArr[i+windowSize-1], 
                                             inputVector[i]);
        }
        else {
          continue;
        }
      }
    }
  );

  int prefixIndex = 0;
  int suffixIndex = 0;

  tbb::parallel_for((size_t) windowSize, N+windowSize, (size_t) 1, 
    [&prefixIndex, &suffixIndex, &erosionVec, 
     &prefixArr, &suffixArr, &windowSize](size_t i) {
      prefixIndex = i + windowSize;
      suffixIndex = i - windowSize;
      erosionVec[i-windowSize] = std::min(prefixArr[prefixIndex],
        suffixArr[suffixIndex]);
    }
  );
  return;
}


void sigproc_multithreading::Morph1DFast::getGradient(
  const ConcurrentVector<short>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<short>& gradientVec) const
{
  getGradient<short>(inputVector, structuringElement, gradientVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getGradient(
  const ConcurrentVector<float>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<float>& gradientVec) const
{
  getGradient<float>(inputVector, structuringElement, gradientVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getGradient(
  const ConcurrentVector<double>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<double>& gradientVec) const
{
  getGradient<double>(inputVector, structuringElement, gradientVec);
  return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getGradient(
  const ConcurrentVector<T>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<T>& gradientVec) const
{
  size_t N = inputVector.size();
  assert(gradientVec.size() == N);

  ConcurrentVector<T> erosionVec(N);
  ConcurrentVector<T> dilationVec(N);

  tbb::parallel_invoke(
    [&]() {
      getErosion<T>(inputVector, structuringElement, erosionVec);
    },
    [&]() {
      getDilation<T>(inputVector, structuringElement, dilationVec);
    }
  );


  // Can be made faster by SIMD Vectorization with PSTL
  for (size_t i=0; i<N; ++i) {
    gradientVec[i] = dilationVec[i] - erosionVec[i];
  }
  return;
}


void sigproc_multithreading::Morph1DFast::getOpening(
  const ConcurrentVector<bool>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<bool>& openingVec) const
{
  size_t N = inputVector.size();
  assert(openingVec.size() == N);
  ConcurrentVector<bool> tempVec(N);
  getErosion(inputVector, structuringElement, tempVec);
  getDilation(tempVec, structuringElement, openingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getOpening(
  const ConcurrentVector<short>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<short>& openingVec) const
{
  getOpening<short>(inputVector, structuringElement, openingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getOpening(
  const ConcurrentVector<float>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<float>& openingVec) const
{
  getOpening<float>(inputVector, structuringElement, openingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getOpening(
  const ConcurrentVector<double>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<double>& openingVec) const
{
  getOpening<double>(inputVector, structuringElement, openingVec);
  return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getOpening(
  const ConcurrentVector<T>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<T>& openingVec) const
{
  size_t N = inputVector.size();
  assert(openingVec.size() == N);
  ConcurrentVector<T> tempVec(N);
  getErosion<T>(inputVector, structuringElement, tempVec);
  getDilation<T>(tempVec, structuringElement, openingVec);
  return;
}


void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<bool>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<bool>& closingVec) const
{
  size_t N = inputVector.size();
  assert(closingVec.size() == N);
  ConcurrentVector<bool> tempVec(N);
  getDilation(inputVector, structuringElement, tempVec);
  getErosion(tempVec, structuringElement, closingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<short>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<short>& closingVec) const
{
  getClosing<short>(inputVector, structuringElement, closingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<float>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<float>& closingVec) const
{
  getClosing<float>(inputVector, structuringElement, closingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<double>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<double>& closingVec) const
{
  getClosing<double>(inputVector, structuringElement, closingVec);
  return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<T>& inputVector,
  const unsigned int structuringElement,
  ConcurrentVector<T>& closingVec) const
{
  size_t N = inputVector.size();
  assert(closingVec.size() == N);
  ConcurrentVector<T> tempVec(N);
  getDilation<T>(inputVector, structuringElement, tempVec);
  getErosion<T>(tempVec, structuringElement, closingVec);
  return;
}


// Column Major Operations
// We implement column major versions of dilation and erosion, as 
// N-dimensional morphological filters are all derivable from ND 
// dilation and ND erosion filters.
void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentArray2D<short>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<short>& dilation2D,
  const unsigned int columnNum) const
{
  getDilation<short>(inputArray2D, structuringElementx, dilation2D, columnNum);
  return;
}

void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentArray2D<float>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<float>& dilation2D,
  const unsigned int columnNum) const
{
  getDilation<float>(inputArray2D, structuringElementx, dilation2D, columnNum);
  return;
}

void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentArray2D<double>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<double>& dilation2D,
  const unsigned int columnNum) const
{
  getDilation<double>(inputArray2D, structuringElementx, dilation2D, columnNum);
  return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentArray2D<T>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<T>& dilation2D,
  const unsigned int columnNum) const
{
  const size_t N = inputArray2D.size();
  assert(dilation2D.size() == N);
  assert(columnNum < inputArray2D.at(0).size());
  const size_t k = (size_t) structuringElementx;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  const size_t windowSize = k/2;
  const size_t paddingSize = (k - (N % k)) % k;
  const size_t bufferSize = N + 2 * windowSize + paddingSize;
  ConcurrentVector<T> suffixArr(bufferSize);
  ConcurrentVector<T> prefixArr(bufferSize);

  // Padding Operations on Buffers
  // This could be done with parallel_for and parallel_invoke, yet as they run
  // through not many entries anyway it may not be worth it. 
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::min();
    prefixArr[i] = std::numeric_limits<T>::min();
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::min();
    prefixArr[i] = std::numeric_limits<T>::min();
  }

  // Prefix computation is independent from suffix computation,
  // so we use tbb::parallel_invoke. 
  // However, each prefix/suffix computation loop must respect the order of
  // element, so we cannot use tbb::parallel_for here. 

  for (size_t i=N+paddingSize; i!=0; --i) {
    if (i > N) {
      // Compensate for divisibility padding (must be -inf)
      continue;
    }
    else if (i % k == 0) {
      suffixArr[i+windowSize-1] = inputArray2D[i-1][columnNum];
    } 
    else {
      suffixArr[i+windowSize-1] = std::max(suffixArr[i+windowSize], 
                                            inputArray2D[i-1][columnNum]);
    }
  }

  for (size_t i=0; i<N+paddingSize; ++i) {
    if (i % k == 0) {
      prefixArr[i+windowSize] = inputArray2D[i][columnNum];
    }
    else if ((i % k == 0) && (i < N)) {
      prefixArr[i+windowSize] = std::max(prefixArr[i+windowSize-1], 
                                          inputArray2D[i][columnNum]);
    }
    else {
      continue;
    }
  }

  // Compute Prefix and Suffix Buffers

  size_t prefixIndex = 0;
  size_t suffixIndex = 0;

  for (size_t i=windowSize; i<N+windowSize; ++i) {
    prefixIndex = i + windowSize;
    suffixIndex = i - windowSize;
    dilation2D[i-windowSize][columnNum] = std::max(prefixArr[prefixIndex],
      suffixArr[suffixIndex]);
  }
  return;
}


void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<short>& inputArray2D,
  const unsigned int structuringElementy,
  ConcurrentArray2D<short>& dilation2D,
  const unsigned int columnNum) const
{
  getErosion<short>(inputArray2D, structuringElementy, dilation2D, columnNum);
  return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<float>& inputArray2D,
  const unsigned int structuringElementy,
  ConcurrentArray2D<float>& dilation2D,
  const unsigned int columnNum) const
{
  getErosion<float>(inputArray2D, structuringElementy, dilation2D, columnNum);
  return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<double>& inputArray2D,
  const unsigned int structuringElementy,
  ConcurrentArray2D<double>& dilation2D,
  const unsigned int columnNum) const
{
  getErosion<double>(inputArray2D, structuringElementy, dilation2D, columnNum);
  return;
}


template <typename T>
void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<T>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<T>& erosion2D,
  const unsigned int columnNum) const
{
  const size_t N = inputArray2D.size();
  assert(erosion2D.size() == N);
  assert(columnNum < inputArray2D.at(0).size());
  const size_t k = (size_t) structuringElementx;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  const size_t windowSize = k/2;
  const size_t paddingSize = (k - (N % k)) % k;
  const size_t bufferSize = N + 2 * windowSize + paddingSize;
  ConcurrentVector<T> suffixArr(bufferSize);
  ConcurrentVector<T> prefixArr(bufferSize);

  // Padding Operations on Buffers
  // This could be done with parallel_for and parallel_invoke, yet as they run
  // through not many entries anyway it may not be worth it. 
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::max();
    prefixArr[i] = std::numeric_limits<T>::max();
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::max();
    prefixArr[i] = std::numeric_limits<T>::max();
  }

  // Prefix computation is independent from suffix computation,
  // so we use tbb::parallel_invoke. 
  // However, each prefix/suffix computation loop must respect the order of
  // element, so we cannot use tbb::parallel_for here. 

  tbb::parallel_invoke(
    [&suffixArr, &inputArray2D, &windowSize, &paddingSize, &columnNum, &N, &k]() {
      for (size_t i=N+paddingSize; i!=0; --i) {
        if (i > N) {
          // Compensate for divisibility padding (must be -inf)
          continue;
        }
        else if (i % k == 0) {
          suffixArr[i+windowSize-1] = inputArray2D[i-1][columnNum];
        } 
        else {
          suffixArr[i+windowSize-1] = std::min(suffixArr[i+windowSize], 
                                               inputArray2D[i-1][columnNum]);
        }
      }
    },
    [&prefixArr, &inputArray2D, &windowSize, &paddingSize, &columnNum, &N, &k]() {
      for (size_t i=0; i<N+paddingSize; ++i) {
        if (i % k == 0) {
          prefixArr[i+windowSize] = inputArray2D[i][columnNum];
        }
        else if ((i % k == 0) && (i < N)) {
          prefixArr[i+windowSize] = std::min(prefixArr[i+windowSize-1], 
                                             inputArray2D[i][columnNum]);
        }
        else {
          continue;
        }
      }
    }
  );

  // Compute Prefix and Suffix Buffers

  size_t prefixIndex = 0;
  size_t suffixIndex = 0;

  tbb::parallel_for((size_t) windowSize, N+windowSize, (size_t) 1, 
    [&prefixIndex, &suffixIndex, &erosion2D, 
     &prefixArr, &suffixArr, &windowSize, &columnNum](size_t i) {
      prefixIndex = i + windowSize;
      suffixIndex = i - windowSize;
      erosion2D[i-windowSize][columnNum] = std::min(prefixArr[prefixIndex],
        suffixArr[suffixIndex]);
    }
  );
  return;
}

#endif
