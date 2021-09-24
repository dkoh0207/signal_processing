#ifndef __SIGPROC_MULTITHREADING_MORPH1DFAST_CXX__
#define __SIGPROC_MULTITHREADING_MORPH1DFAST_CXX__

#include "Morph1DFast.h"


void sigproc_multithreading::Morph1DFast::getDilation(
    const ConcurrentVector<bool>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<bool>& dilationVec) const
{
     /*
      Module for 1D Dilation Filter - special handling for bool arrays

      INPUTS:
      - waveform: 1D Pedestal Corrected Waveform.
      - fStructuringElement: Size of moving window
  
      MODIFIES:
      - dilationVec: Returned Dilation ConcurrentVector.
    */
    if (dilationVec.size() != inputWaveform.size())
    {
        std::cout << "Dilation1D: output dilation vector not same size as "
            "input waveform array" << std::endl;
        return;
    }

    const size_t N = inputWaveform.size();

    if (N <= fStructuringElement) 
    {
        std::cout << "Dilation1D: Input array size " << N << 
            " must be greater than structuring element size " << 
            fStructuringElement << std::endl;
        return;
    }
    const size_t windowSize  = fStructuringElement/2;
    const size_t paddingSize = (fStructuringElement - 
                               (N % fStructuringElement)) % fStructuringElement;
    const size_t bufferSize  = N + 2 * windowSize + paddingSize;

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

    // Compute Prefix (h_x) and Suffix (g_x)

    int kint = (int) fStructuringElement;

    tbb::parallel_invoke(
        [&prefixArr, &inputWaveform, &kint, &windowSize, &N]()
        {
            for (int i=N-1; i>-1; --i) {
                if ((i+1) % kint == 0) {
                    prefixArr[i+windowSize] = inputWaveform[i];
                } 
                else if ((i+1) % kint != 0) {
                    prefixArr[i+windowSize] = (prefixArr[i+1+windowSize] || 
                                               inputWaveform[i]);
                }
                else continue;
            }
        },
        [&suffixArr, &inputWaveform, &kint, &windowSize, &N]()
        {
            for (size_t i=0; i<N; ++i) {

                if (i > N) {
                    // Compensate for divisibility padding (must be -inf)
                    continue;
                }
                else if (i % kint == 0) {
                    suffixArr[i+windowSize] = inputWaveform[i];
                } 
                else {
                    suffixArr[i+windowSize] = (suffixArr[i-1+windowSize] || 
                                               inputWaveform[i]);
                }
            }
        }
    );

    int prefixIndex = 0;
    int suffixIndex = 0;

    tbb::parallel_for( (size_t) 0, N, (size_t) 1, 
        [&prefixIndex, &suffixIndex, &dilationVec, 
        &prefixArr, &suffixArr, &windowSize](size_t i)
        {
            suffixIndex = i + 2 * windowSize;
            prefixIndex = i;
            dilationVec[i] = (prefixArr[prefixIndex] || suffixArr[suffixIndex]);
        });

    return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getDilation(
    const ConcurrentVector<T>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<T>& dilationVec) const
{
     /*
      Module for 1D Dilation Filter - special handling for bool arrays

      INPUTS:
      - waveform: 1D Pedestal Corrected Waveform.
      - fStructuringElement: Size of moving window
  
      MODIFIES:
      - dilationVec: Returned Dilation ConcurrentVector.
    */
    if (dilationVec.size() != inputWaveform.size())
    {
        std::cout << "Dilation1D: output dilation vector not same size as "
            "input waveform array" << std::endl;
        return;
    }

    const size_t N = inputWaveform.size();

    if (N <= fStructuringElement) 
    {
        std::cout << "Dilation1D: Input array size " << N << 
            " must be greater than structuring element size " << 
            fStructuringElement << std::endl;
        return;
    }
    const size_t windowSize  = fStructuringElement/2;
    const size_t paddingSize = (fStructuringElement - 
                               (N % fStructuringElement)) % fStructuringElement;
    const size_t bufferSize  = N + 2 * windowSize + paddingSize;

    ConcurrentVector<T> suffixArr(bufferSize);
    ConcurrentVector<T> prefixArr(bufferSize);

    // Padding Operations on Buffers
    for (size_t i=0; i<windowSize; ++i) {
        suffixArr[i] = std::numeric_limits<T>::min();
        prefixArr[i] = std::numeric_limits<T>::min();
    }

    for (size_t i=N+windowSize; i<bufferSize; ++i) {
        suffixArr[i] = std::numeric_limits<T>::min();
        prefixArr[i] = std::numeric_limits<T>::min();
    }


    // Compute Prefix (h_x) and Suffix (g_x)

    int kint = (int) fStructuringElement;

    tbb::parallel_invoke(
        [&prefixArr, &inputWaveform, &kint, &windowSize, &N]()
        {
            for (int i=N-1; i>-1; --i) {
                if ((i+1) % kint == 0) {
                    prefixArr[i+windowSize] = inputWaveform[i];
                } 
                else if ((i+1) % kint != 0) {
                    prefixArr[i+windowSize] = (prefixArr[i+1+windowSize] || 
                                               inputWaveform[i]);
                }
                else continue;
            }
        },
        [&suffixArr, &inputWaveform, &kint, &windowSize, &N]()
        {
            for (size_t i=0; i<N; ++i) {

                if (i > N) {
                    // Compensate for divisibility padding (must be -inf)
                    continue;
                }
                else if (i % kint == 0) {
                    suffixArr[i+windowSize] = inputWaveform[i];
                } 
                else {
                    suffixArr[i+windowSize] = std::max(
                        suffixArr[i-1+windowSize], inputWaveform[i]);
                }
            }
        }
    );

    int prefixIndex = 0;
    int suffixIndex = 0;

    tbb::parallel_for( (size_t) 0, N, (size_t) 1, 
        [&prefixIndex, &suffixIndex, &dilationVec, 
        &prefixArr, &suffixArr, &windowSize](size_t i)
        {
            suffixIndex = i + 2 * windowSize;
            prefixIndex = i;
            dilationVec[i] = std::max(prefixArr[prefixIndex], 
                                      suffixArr[suffixIndex]);
        });

    return;
}

void sigproc_multithreading::Morph1DFast::getDilation(
    const ConcurrentVector<short>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<short>& dilationVec) const
{
    getDilation<short>(inputWaveform, fStructuringElement, dilationVec);
    return;
}

void sigproc_multithreading::Morph1DFast::getDilation(
    const ConcurrentVector<float>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<float>& dilationVec) const
{
    getDilation<float>(inputWaveform, fStructuringElement, dilationVec);
    return;
}

void sigproc_multithreading::Morph1DFast::getDilation(
    const ConcurrentVector<double>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<double>& dilationVec) const
{
    getDilation<double>(inputWaveform, fStructuringElement, dilationVec);
    return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
    const ConcurrentVector<bool>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<bool>& erosionVec) const
{
     /*
      Module for 1D Erosion Filter - special handling for bool arrays

      INPUTS:
      - waveform: 1D Pedestal Corrected Waveform.
      - fStructuringElement: Size of moving window
  
      MODIFIES:
      - erosionVec: Returned Erosion ConcurrentVector.
    */
    if (erosionVec.size() != inputWaveform.size())
    {
        std::cout << "Dilation1D: output dilation vector not same size as "
            "input waveform array" << std::endl;
        return;
    }

    const size_t N = inputWaveform.size();

    if (N <= fStructuringElement) 
    {
        std::cout << "Dilation1D: Input array size " << N << 
            " must be greater than structuring element size " << 
            fStructuringElement << std::endl;
        return;
    }
    const size_t windowSize  = fStructuringElement/2;
    const size_t paddingSize = (fStructuringElement - 
                               (N % fStructuringElement)) % fStructuringElement;
    const size_t bufferSize  = N + 2 * windowSize + paddingSize;

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

    // Compute Prefix (h_x) and Suffix (g_x)

    int kint = (int) fStructuringElement;

    tbb::parallel_invoke(
        [&prefixArr, &inputWaveform, &kint, &windowSize, &N]()
        {
            for (int i=N-1; i>-1; --i) {
                if ((i+1) % kint == 0) {
                    prefixArr[i+windowSize] = inputWaveform[i];
                } 
                else if ((i+1) % kint != 0) {
                    prefixArr[i+windowSize] = (prefixArr[i+1+windowSize] && 
                                               inputWaveform[i]);
                }
                else continue;
            }
        },
        [&suffixArr, &inputWaveform, &kint, &windowSize, &N]()
        {
            for (size_t i=0; i<N; ++i) {

                if (i > N) {
                    // Compensate for divisibility padding (must be -inf)
                    continue;
                }
                else if (i % kint == 0) {
                    suffixArr[i+windowSize] = inputWaveform[i];
                } 
                else {
                    suffixArr[i+windowSize] = (suffixArr[i-1+windowSize] && 
                                               inputWaveform[i]);
                }
            }
        }
    );

    int prefixIndex = 0;
    int suffixIndex = 0;

    tbb::parallel_for( (size_t) 0, N, (size_t) 1, 
        [&prefixIndex, &suffixIndex, &erosionVec, 
        &prefixArr, &suffixArr, &windowSize](size_t i)
        {
            suffixIndex = i + 2 * windowSize;
            prefixIndex = i;
            erosionVec[i] = (prefixArr[prefixIndex] && suffixArr[suffixIndex]);
        });

    return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getErosion(
    const ConcurrentVector<T>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<T>& erosionVec) const
{
    if (erosionVec.size() != inputWaveform.size())
    {
        std::cout << "Erosion1D: output erosion vector not same size as "
            "input waveform array" << std::endl;
        return;
    }

    const size_t N = inputWaveform.size();

    if (N <= fStructuringElement) 
    {
        std::cout << "Erosion1D: Input array size " << N << 
            " must be greater than structuring element size " << 
            fStructuringElement << std::endl;
        return;
    }

    const size_t windowSize  = fStructuringElement/2;
    const size_t paddingSize = (fStructuringElement - 
                               (N % fStructuringElement)) % fStructuringElement;
    const size_t bufferSize  = N + 2 * windowSize + paddingSize;

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

    int kint = (int) fStructuringElement;
    // Compute Prefix and Suffix Buffers
    tbb::parallel_invoke(
        [&prefixArr, &inputWaveform, &kint, &windowSize, &N]()
        {
            for (int i=N-1; i>-1; --i) {
                if ((i+1) % kint == 0) {
                    prefixArr[i+windowSize] = inputWaveform[i];
                } 
                else if ((i+1) % kint != 0) {
                    prefixArr[i+windowSize] = std::min(
                        prefixArr[i+1+windowSize], inputWaveform[i]);
                }
                else continue;
            }
        },
        [&suffixArr, &inputWaveform, &kint, &windowSize, &N]()
        {
            for (size_t i=0; i<N; ++i) {

                if (i > N) {
                    // Compensate for divisibility padding (must be -inf)
                    continue;
                }
                else if (i % kint == 0) {
                    suffixArr[i+windowSize] = inputWaveform[i];
                } 
                else {
                    suffixArr[i+windowSize] = std::min(
                        suffixArr[i-1+windowSize], inputWaveform[i]);
                }
            }
        }
    );

    int prefixIndex = 0;
    int suffixIndex = 0;

    tbb::parallel_for( (size_t) 0, N, (size_t) 1, 
        [&prefixIndex, &suffixIndex, &erosionVec, 
        &prefixArr, &suffixArr, &windowSize](size_t i)
        {
            suffixIndex = i + 2 * windowSize;
            prefixIndex = i;
            erosionVec[i] = std::min(prefixArr[prefixIndex], 
                                     suffixArr[suffixIndex]);
        });
    
    return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
    const ConcurrentVector<short>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<short>& dilationVec) const
{
    getErosion<short>(inputWaveform, fStructuringElement, dilationVec);
    return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
    const ConcurrentVector<float>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<float>& dilationVec) const
{
    getErosion<float>(inputWaveform, fStructuringElement, dilationVec);
    return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
    const ConcurrentVector<double>& inputWaveform,
    const unsigned int fStructuringElement,
    ConcurrentVector<double>& dilationVec) const
{
    getErosion<double>(inputWaveform, fStructuringElement, dilationVec);
    return;
}

void sigproc_multithreading::Morph1DFast::getGradient(
  const ConcurrentVector<short>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<short>& gradientVec) const
{
    getGradient<short>(inputVector, fStructuringElement, gradientVec);
    return;
}

void sigproc_multithreading::Morph1DFast::getGradient(
  const ConcurrentVector<float>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<float>& gradientVec) const
{
  getGradient<float>(inputVector, fStructuringElement, gradientVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getGradient(
  const ConcurrentVector<double>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<double>& gradientVec) const
{
  getGradient<double>(inputVector, fStructuringElement, gradientVec);
  return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getGradient(
  const ConcurrentVector<T>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<T>& gradientVec) const
{
    size_t N = inputVector.size();
    assert(gradientVec.size() == N);

    ConcurrentVector<T> erosionVec(N);
    ConcurrentVector<T> dilationVec(N);

    tbb::parallel_invoke(
        [this, &inputVector, &fStructuringElement, &erosionVec]() 
        {
            getErosion<T>(inputVector, fStructuringElement, erosionVec);
        },
        [this, &inputVector, &fStructuringElement, &dilationVec]() 
        {
            getDilation<T>(inputVector, fStructuringElement, dilationVec);
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
  const unsigned int fStructuringElement,
  ConcurrentVector<bool>& openingVec) const
{
    size_t N = inputVector.size();
    assert(openingVec.size() == N);
    ConcurrentVector<bool> tempVec(N);
    getErosion(inputVector, fStructuringElement, tempVec);
    getDilation(tempVec, fStructuringElement, openingVec);
    return;
}

void sigproc_multithreading::Morph1DFast::getOpening(
  const ConcurrentVector<short>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<short>& openingVec) const
{
  getOpening<short>(inputVector, fStructuringElement, openingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getOpening(
  const ConcurrentVector<float>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<float>& openingVec) const
{
  getOpening<float>(inputVector, fStructuringElement, openingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getOpening(
  const ConcurrentVector<double>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<double>& openingVec) const
{
  getOpening<double>(inputVector, fStructuringElement, openingVec);
  return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getOpening(
  const ConcurrentVector<T>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<T>& openingVec) const
{
  size_t N = inputVector.size();
  assert(openingVec.size() == N);
  ConcurrentVector<T> tempVec(N);
  getErosion<T>(inputVector, fStructuringElement, tempVec);
  getDilation<T>(tempVec, fStructuringElement, openingVec);
  return;
}


void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<bool>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<bool>& closingVec) const
{
  size_t N = inputVector.size();
  assert(closingVec.size() == N);
  ConcurrentVector<bool> tempVec(N);
  getDilation(inputVector, fStructuringElement, tempVec);
  getErosion(tempVec, fStructuringElement, closingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<short>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<short>& closingVec) const
{
  getClosing<short>(inputVector, fStructuringElement, closingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<float>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<float>& closingVec) const
{
  getClosing<float>(inputVector, fStructuringElement, closingVec);
  return;
}

void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<double>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<double>& closingVec) const
{
  getClosing<double>(inputVector, fStructuringElement, closingVec);
  return;
}

template <typename T>
void sigproc_multithreading::Morph1DFast::getClosing(
  const ConcurrentVector<T>& inputVector,
  const unsigned int fStructuringElement,
  ConcurrentVector<T>& closingVec) const
{
  size_t N = inputVector.size();
  assert(closingVec.size() == N);
  ConcurrentVector<T> tempVec(N);
  getDilation<T>(inputVector, fStructuringElement, tempVec);
  getErosion<T>(tempVec, fStructuringElement, closingVec);
  return;
}


void sigproc_multithreading::Morph1DFast::getDilation(
  const ConcurrentArray2D<bool>& inputArray2D,
  const unsigned int fStructuringElementx,
  ConcurrentArray2D<bool>& dilation2D,
  const int columnNum) const
{
    if ((dilation2D.size() != inputArray2D.size()) || 
        (dilation2D.at(0).size() != inputArray2D.at(0).size()))
    {
        std::cout << "Dilation1D: output dilation vector not same size as "
            "input waveform array" << std::endl;
        return;
    }

    const size_t N = inputArray2D.size();

    if (N <= fStructuringElementx) 
    {
        std::cout << "Dilation1D: Input array size " << N << 
            " must be greater than structuring element size " << 
            fStructuringElementx << std::endl;
        return;
    }

    const size_t windowSize  = fStructuringElementx/2;
    const size_t paddingSize = (fStructuringElementx - 
                               (N % fStructuringElementx)) % fStructuringElementx;
    const size_t bufferSize  = N + 2 * windowSize + paddingSize;
    ConcurrentVector<bool> suffixArr(bufferSize);
    ConcurrentVector<bool> prefixArr(bufferSize);

    // Padding Operations on Buffers
    // This could be done with parallel_for and parallel_invoke, yet as they run
    // through not many entries anyway it may not be worth it. 
    for (size_t i=0; i<windowSize; ++i) {
        suffixArr[i] = false;
        prefixArr[i] = false;
    }

    for (size_t i=N+windowSize; i<bufferSize; ++i) {
        suffixArr[i] = false;
        prefixArr[i] = false;
    }

    // Prefix computation is independent from suffix computation,
    // so we use tbb::parallel_invoke. 
    // However, each prefix/suffix computation loop must respect the order of
    // element, so we cannot use tbb::parallel_for here. 

    int kint = (int) fStructuringElementx;

    tbb::parallel_invoke(
        [&columnNum, &windowSize, &inputArray2D, &N, &suffixArr, &kint]() {

            for (size_t i=0; i<N; ++i) {

                if (i > N) {
                    // Compensate for divisibility padding (must be -inf)
                    continue;
                }
                else if (i % kint == 0) {
                    suffixArr[i+windowSize] = inputArray2D[i][columnNum];
                } 
                else {
                    suffixArr[i+windowSize] = suffixArr[i-1+windowSize] || 
                                              inputArray2D[i][columnNum];
                }
            }
        },
        [&columnNum, &windowSize, &inputArray2D, &N, &prefixArr, &kint]() {

            for (int i=N-1; i>-1; --i) {

                if ((i+1) % kint == 0) {
                    prefixArr[i+windowSize] = inputArray2D[i][columnNum];
                } 
                else if ((i+1) % kint != 0) {
                    prefixArr[i+windowSize] = prefixArr[i+1+windowSize] ||
                                              inputArray2D[i][columnNum];
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

    tbb::parallel_for( (size_t) 0, N, (size_t) 1, 
        [&prefixIndex, &suffixIndex, &dilation2D, 
        &prefixArr, &suffixArr, &windowSize, &columnNum](size_t i)
        {
            suffixIndex = i + 2 * windowSize;
            prefixIndex = i;
            dilation2D[i][columnNum] = prefixArr[prefixIndex] || 
                                       suffixArr[suffixIndex];
        });

    return;
}

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
  const unsigned int fStructuringElementx,
  ConcurrentArray2D<T>& dilation2D,
  const int columnNum) const
{
    if ((dilation2D.size() != inputArray2D.size()) || 
        (dilation2D.at(0).size() != inputArray2D.at(0).size()))
    {
        std::cout << "Dilation1D: output dilation vector not same size as "
            "input waveform array" << std::endl;
        return;
    }

    const size_t N = inputArray2D.size();

    if (N <= fStructuringElementx) 
    {
        std::cout << "Dilation1D: Input array size " << N << 
            " must be greater than structuring element size " << 
            fStructuringElementx << std::endl;
        return;
    }

    const size_t windowSize  = fStructuringElementx/2;
    const size_t paddingSize = (fStructuringElementx - 
                               (N % fStructuringElementx)) % fStructuringElementx;
    const size_t bufferSize  = N + 2 * windowSize + paddingSize;
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

    int kint = (int) fStructuringElementx;

    tbb::parallel_invoke(
        [&columnNum, &windowSize, &inputArray2D, &N, &suffixArr, &kint]() {

            for (size_t i=0; i<N; ++i) {

                if (i > N) {
                    // Compensate for divisibility padding (must be -inf)
                    continue;
                }
                else if (i % kint == 0) {
                    suffixArr[i+windowSize] = inputArray2D[i][columnNum];
                } 
                else {
                    suffixArr[i+windowSize] = std::max(
                        suffixArr[i-1+windowSize], inputArray2D[i][columnNum]);
                }
            }
        },
        [&columnNum, &windowSize, &inputArray2D, &N, &prefixArr, &kint]() {

            for (int i=N-1; i>-1; --i) {

                if ((i+1) % kint == 0) {
                    prefixArr[i+windowSize] = inputArray2D[i][columnNum];
                } 
                else if ((i+1) % kint != 0) {
                    prefixArr[i+windowSize] = std::max(
                        prefixArr[i+1+windowSize], inputArray2D[i][columnNum]);
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

    tbb::parallel_for( (size_t) 0, N, (size_t) 1, 
        [&prefixIndex, &suffixIndex, &dilation2D, 
        &prefixArr, &suffixArr, &windowSize, &columnNum](size_t i)
        {
            suffixIndex = i + 2 * windowSize;
            prefixIndex = i;
            dilation2D[i][columnNum] = std::max(prefixArr[prefixIndex], 
                                                suffixArr[suffixIndex]);
        });

    return;
}


void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<bool>& inputArray2D,
  const unsigned int fStructuringElementx,
  ConcurrentArray2D<bool>& erosion2D,
  const int columnNum) const
{
    if ((erosion2D.size() != inputArray2D.size()) || 
        (erosion2D.at(0).size() != inputArray2D.at(0).size()))
    {
        std::cout << "Dilation1D: output dilation vector not same size as "
            "input waveform array" << std::endl;
        return;
    }

    const size_t N = inputArray2D.size();

    if (N <= fStructuringElementx) 
    {
        std::cout << "Dilation1D: Input array size " << N << 
            " must be greater than structuring element size " << 
            fStructuringElementx << std::endl;
        return;
    }

    const size_t windowSize  = fStructuringElementx/2;
    const size_t paddingSize = (fStructuringElementx - 
                               (N % fStructuringElementx)) % fStructuringElementx;
    const size_t bufferSize  = N + 2 * windowSize + paddingSize;
    ConcurrentVector<bool> suffixArr(bufferSize);
    ConcurrentVector<bool> prefixArr(bufferSize);

    // Padding Operations on Buffers
    // This could be done with parallel_for and parallel_invoke, yet as they run
    // through not many entries anyway it may not be worth it. 
    for (size_t i=0; i<windowSize; ++i) {
        suffixArr[i] = true;
        prefixArr[i] = true;
    }

    for (size_t i=N+windowSize; i<bufferSize; ++i) {
        suffixArr[i] = true;
        prefixArr[i] = true;
    }

    // Prefix computation is independent from suffix computation,
    // so we use tbb::parallel_invoke. 
    // However, each prefix/suffix computation loop must respect the order of
    // element, so we cannot use tbb::parallel_for here. 

    int kint = (int) fStructuringElementx;

    tbb::parallel_invoke(
        [&columnNum, &windowSize, &inputArray2D, &N, &suffixArr, &kint]() {

            for (size_t i=0; i<N; ++i) {

                if (i > N) {
                    // Compensate for divisibility padding (must be -inf)
                    continue;
                }
                else if (i % kint == 0) {
                    suffixArr[i+windowSize] = inputArray2D[i][columnNum];
                } 
                else {
                    suffixArr[i+windowSize] = suffixArr[i-1+windowSize] && 
                                              inputArray2D[i][columnNum];
                }
            }
        },
        [&columnNum, &windowSize, &inputArray2D, &N, &prefixArr, &kint]() {

            for (int i=N-1; i>-1; --i) {

                if ((i+1) % kint == 0) {
                    prefixArr[i+windowSize] = inputArray2D[i][columnNum];
                } 
                else if ((i+1) % kint != 0) {
                    prefixArr[i+windowSize] = prefixArr[i+1+windowSize] &&
                                              inputArray2D[i][columnNum];
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

    tbb::parallel_for( (size_t) 0, N, (size_t) 1, 
        [&prefixIndex, &suffixIndex, &erosion2D, 
        &prefixArr, &suffixArr, &windowSize, &columnNum](size_t i)
        {
            suffixIndex = i + 2 * windowSize;
            prefixIndex = i;
            erosion2D[i][columnNum] = prefixArr[prefixIndex] && 
                                      suffixArr[suffixIndex];
        });

    return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<short>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<short>& dilation2D,
  const unsigned int columnNum) const
{
  getErosion<short>(inputArray2D, structuringElementx, dilation2D, columnNum);
  return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<float>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<float>& dilation2D,
  const unsigned int columnNum) const
{
  getErosion<float>(inputArray2D, structuringElementx, dilation2D, columnNum);
  return;
}

void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<double>& inputArray2D,
  const unsigned int structuringElementx,
  ConcurrentArray2D<double>& dilation2D,
  const unsigned int columnNum) const
{
  getErosion<double>(inputArray2D, structuringElementx, dilation2D, columnNum);
  return;
}


template <typename T>
void sigproc_multithreading::Morph1DFast::getErosion(
  const ConcurrentArray2D<T>& inputArray2D,
  const unsigned int fStructuringElementx,
  ConcurrentArray2D<T>& erosion2D,
  const int columnNum) const
{
    if ( (erosion2D.size() != inputArray2D.size()) ||
         (erosion2D.at(0).size() != inputArray2D.at(0).size()) )
    {
        std::cout << "Erosion1D: output dilation vector not same size as "
            "input waveform array" << std::endl;
        return;
    }

    const size_t N = inputArray2D.size();

    if (N <= fStructuringElementx) 
    {
        std::cout << "Erosion1D: Input array size " << N << " must be "
            "greater than structuring element size " << 
            fStructuringElementx << std::endl;
        return;
    }

    const size_t windowSize  = fStructuringElementx/2;
    const size_t paddingSize = (fStructuringElementx - 
                               (N % fStructuringElementx)) % fStructuringElementx;
    const size_t bufferSize  = N + 2 * windowSize + paddingSize;

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

    int kint = (int) fStructuringElementx;

    tbb::parallel_invoke(
        [&suffixArr, &inputArray2D, &windowSize, 
         &paddingSize, &columnNum, &N, &kint]() {

            for (size_t i=0; i<N; ++i) {

                if (i > N) {
                    // Compensate for divisibility padding (must be -inf)
                    continue;
                }
                else if (i % kint == 0) {
                    suffixArr[i+windowSize] = inputArray2D[i][columnNum];
                } 
                else {
                    suffixArr[i+windowSize] = std::min(
                        suffixArr[i-1+windowSize], inputArray2D[i][columnNum]);
                }
            }
        },
        [&prefixArr, &inputArray2D, &windowSize, 
         &paddingSize, &columnNum, &N, &kint]() {

            for (int i=N-1; i>-1; --i) {

                if ((i+1) % kint == 0) {
                    prefixArr[i+windowSize] = inputArray2D[i][columnNum];
                } 
                else if ((i+1) % kint != 0) {
                    prefixArr[i+windowSize] = std::min(
                        prefixArr[i+1+windowSize], inputArray2D[i][columnNum]);
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
            suffixIndex = i + 2 * windowSize;
            prefixIndex = i;
            erosion2D[i][columnNum] = std::min(
                prefixArr[prefixIndex], suffixArr[suffixIndex]);
        }
    );
    return;
}

#endif
