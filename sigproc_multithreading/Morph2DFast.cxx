#ifndef __SIGPROC_MULTITHREADING_MORPH2DFAST_CXX__
#define __SIGPROC_MULTITHREADING_MORPH2DFAST_CXX__

#include "Morph2DFast.h"

using namespace sigproc_multithreading;

void sigproc_multithreading::Morph2DFast::getDilation(
  const ConcurrentArray2D<bool>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<bool>& dilation2D) const
{
  getDilation<bool>(inputArray2D, structuringElementx,
    structuringElementy, dilation2D);
  return;
}

void sigproc_multithreading::Morph2DFast::getDilation(
  const ConcurrentArray2D<float>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<float>& dilation2D) const
{
  getDilation<float>(inputArray2D, structuringElementx,
    structuringElementy, dilation2D);
  return;
}

template <typename T>
void sigproc_multithreading::Morph2DFast::getDilation(
  const ConcurrentArray2D<T>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<T>& dilation2D) const
{
  size_t numChannels = inputArray2D.size();
  size_t numTicks = inputArray2D.at(0).size();
  
  std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();
  std::chrono::high_resolution_clock::time_point allocateStart = funcStartTime;

  ConcurrentArray2D<T> buffer(numChannels, ConcurrentVector<T>(numTicks));

  std::chrono::high_resolution_clock::time_point allocateStop = std::chrono::high_resolution_clock::now();

  assert(dilation2D.size() == numChannels);
  assert(dilation2D.at(0).size() == numTicks);
  assert(buffer.size() == numChannels);
  assert(buffer.at(0).size() == numTicks);

  sigproc_multithreading::Morph1DFast fast1D;

  std::chrono::high_resolution_clock::time_point rowStart = allocateStop;

  tbb::parallel_for( (size_t) 0, numChannels, (size_t) 1,
    [&fast1D, &inputArray2D, &structuringElementy, &buffer](size_t i) {
      fast1D.getDilation(inputArray2D[i], structuringElementy, buffer[i]);
    });

  std::chrono::high_resolution_clock::time_point rowStop = std::chrono::high_resolution_clock::now();
  std::chrono::high_resolution_clock::time_point colStart = rowStop;
  // Column operations 
  tbb::parallel_for ( (size_t) 0, numTicks, (size_t) 1,
    [&fast1D, &buffer, &structuringElementx, &dilation2D](size_t j) {
      fast1D.getDilation(buffer, structuringElementx, dilation2D, j);
    });

  std::chrono::high_resolution_clock::time_point colStop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> allocateTime = std::chrono::duration_cast<std::chrono::duration<double>>(allocateStop - allocateStart);
  std::chrono::duration<double> rowTime = std::chrono::duration_cast<std::chrono::duration<double>>(rowStop - rowStart);
  std::chrono::duration<double> colTime = std::chrono::duration_cast<std::chrono::duration<double>>(colStop - colStart);

  std::cout << "*** Dilation ***  " << std::endl;
  std::cout << "Allocate: " << allocateTime.count() << ", Row: " << rowTime.count() << ", Col: " << colTime.count() << std::endl;

  return;
}


void sigproc_multithreading::Morph2DFast::getErosion(
  const ConcurrentArray2D<bool>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<bool>& erosion2D) const
{
  getErosion<bool>(inputArray2D, structuringElementx,
    structuringElementy, erosion2D);
  return;
}

void sigproc_multithreading::Morph2DFast::getErosion(
  const ConcurrentArray2D<float>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<float>& erosion2D) const
{
  getErosion<float>(inputArray2D, structuringElementx,
    structuringElementy, erosion2D);
  return;
}

template <typename T>
void sigproc_multithreading::Morph2DFast::getErosion(
  const ConcurrentArray2D<T>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<T>& erosion2D) const
{
  size_t numChannels = inputArray2D.size();
  size_t numTicks = inputArray2D.at(0).size();

  ConcurrentArray2D<T> buffer(numChannels, ConcurrentVector<T>(numTicks));

  assert(erosion2D.size() == numChannels);
  assert(erosion2D.at(0).size() == numTicks);
  assert(buffer.size() == numChannels);
  assert(buffer.at(0).size() == numTicks);

  sigproc_multithreading::Morph1DFast fast1D;

  tbb::parallel_for( (size_t) 0, numChannels, (size_t) 1,
    [&fast1D, &inputArray2D, &structuringElementy, &buffer](size_t i) {
      fast1D.getErosion(inputArray2D[i], structuringElementy, buffer[i]);
    });

  // Column operations 
  tbb::parallel_for ( (size_t) 0, numTicks, (size_t) 1,
    [&fast1D, &buffer, &structuringElementx, &erosion2D](size_t j) {
      fast1D.getErosion(buffer, structuringElementx, erosion2D, j);
    });

  return;
}


void sigproc_multithreading::Morph2DFast::getGradient(
  const ConcurrentArray2D<float>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<float>& gradient2D) const
{
  size_t numChannels = inputArray2D.size();
  size_t numTicks = inputArray2D.at(0).size();

  assert(gradient2D.size() == numChannels);
  assert(gradient2D.at(0).size() == numTicks);

  ConcurrentArray2D<float> dilation2D(numChannels, ConcurrentVector<float>(numTicks));
  ConcurrentArray2D<float> erosion2D(numChannels, ConcurrentVector<float>(numTicks));

  tbb::parallel_invoke(
    [&]() {
      getDilation<float>(inputArray2D, structuringElementx, structuringElementy, dilation2D);
    },
    [&]() {
      getErosion<float>(inputArray2D, structuringElementx, structuringElementy, erosion2D);
    }
  );

  // This is better done by adding a SIMD layer. 
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<numTicks; ++j) {
      gradient2D[i][j] = (dilation2D[i][j] - erosion2D[i][j]);
    }
  }
  return;
}


void sigproc_multithreading::Morph2DFast::getClosing(
  const ConcurrentArray2D<bool>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<bool>& closing2D) const
{
  getClosing<bool>(inputArray2D, structuringElementx,
    structuringElementy, closing2D);
  return;
}

void sigproc_multithreading::Morph2DFast::getClosing(
  const ConcurrentArray2D<float>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<float>& closing2D) const
{
  getClosing<float>(inputArray2D, structuringElementx,
    structuringElementy, closing2D);
  return;
}

template <typename T>
void sigproc_multithreading::Morph2DFast::getClosing(
  const ConcurrentArray2D<T>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<T>& closing2D) const
{
  size_t numChannels = inputArray2D.size();
  size_t numTicks = inputArray2D.at(0).size();

  assert(closing2D.size() == numChannels);
  assert(closing2D.at(0).size() == numTicks);

  ConcurrentArray2D<T> buffer(numChannels, ConcurrentVector<T>(numTicks));

  getDilation<T>(inputArray2D, structuringElementx, structuringElementy, buffer);
  getErosion<T>(buffer, structuringElementx, structuringElementy, closing2D);
  return;
}


void sigproc_multithreading::Morph2DFast::getOpening(
  const ConcurrentArray2D<bool>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<bool>& opening2D) const
{
  getOpening<bool>(inputArray2D, structuringElementx,
    structuringElementy, opening2D);
  return;
}


void sigproc_multithreading::Morph2DFast::getOpening(
  const ConcurrentArray2D<float>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<float>& opening2D) const
{
  getOpening<float>(inputArray2D, structuringElementx,
    structuringElementy, opening2D);
  return;
}

template <typename T>
void sigproc_multithreading::Morph2DFast::getOpening(
  const ConcurrentArray2D<T>& inputArray2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  ConcurrentArray2D<T>& opening2D) const
{
  size_t numChannels = inputArray2D.size();
  size_t numTicks = inputArray2D.at(0).size();

  assert(opening2D.size() == numChannels);
  assert(opening2D.at(0).size() == numTicks);

  ConcurrentArray2D<T> buffer(numChannels, ConcurrentVector<T>(numTicks));

  getErosion<T>(inputArray2D, structuringElementx, structuringElementy, buffer);
  getDilation<T>(buffer, structuringElementx, structuringElementy, opening2D);

  return;
}


template <typename T>
ConcurrentArray2D<T> sigproc_multithreading::Morph2DFast::converttoConcurrent(
  const Array2D<T>& inputSTL) const
{
  size_t numChannels = inputSTL.size();
  size_t numTicks = inputSTL.at(0).size();

  ConcurrentArray2D<T> outputTBB(numChannels);
  for (auto& v : outputTBB) {
    v.grow_to_at_least(numTicks);
  }

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<numTicks; ++j) {
      outputTBB[i][j] = inputSTL[i][j];
    }
  }

  return outputTBB;
}


ConcurrentArray2D<float> sigproc_multithreading::Morph2DFast::converttoConcurrent(
  const Array2D<float>& inputSTL) const
{
  ConcurrentArray2D<float> outputTBB = converttoConcurrent<float>(inputSTL);
  return outputTBB;
}

ConcurrentArray2D<bool> sigproc_multithreading::Morph2DFast::converttoConcurrent(
  const Array2D<bool>& inputSTL) const
{
  ConcurrentArray2D<bool> outputTBB = converttoConcurrent<bool>(inputSTL);
  return outputTBB;
}


template <typename T>
Array2D<T> sigproc_multithreading::Morph2DFast::converttoSTL(
  const ConcurrentArray2D<T>& inputTBB) const
{
  size_t numChannels = inputTBB.size();
  size_t numTicks = inputTBB.at(0).size();

  Array2D<T> outputSTL(numChannels);
  for (auto& v : outputSTL) {
    v.resize(numTicks);
  }

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<numTicks; ++j) {
      outputSTL[i][j] = inputTBB[i][j];
    }
  }

  return outputSTL;
}


Array2D<float> sigproc_multithreading::Morph2DFast::converttoSTL(
  const ConcurrentArray2D<float>& inputTBB) const
{
  Array2D<float> outputSTL = converttoSTL<float>(inputTBB);
  return outputSTL;
}

Array2D<bool> sigproc_multithreading::Morph2DFast::converttoSTL(
  const ConcurrentArray2D<bool>& inputTBB) const
{
  Array2D<bool> outputSTL = converttoSTL<bool>(inputTBB);
  return outputSTL;
}


#endif
