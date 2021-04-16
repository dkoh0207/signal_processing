#ifndef __SIGPROC_MULTITHREADING_MISCUTILS_CXX__
#define __SIGPROC_MULTITHREADING_MISCUTILS_CXX__

#include "MiscUtils.h"

using namespace sigproc_multithreading;

template <typename T>
ConcurrentArray2D<T> sigproc_multithreading::MiscUtils::converttoConcurrent(
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


ConcurrentArray2D<float> sigproc_multithreading::MiscUtils::converttoConcurrent(
  const Array2D<float>& inputSTL) const
{
  ConcurrentArray2D<float> outputTBB = converttoConcurrent<float>(inputSTL);
  return outputTBB;
}

ConcurrentArray2D<bool> sigproc_multithreading::MiscUtils::converttoConcurrent(
  const Array2D<bool>& inputSTL) const
{
  ConcurrentArray2D<bool> outputTBB = converttoConcurrent<bool>(inputSTL);
  return outputTBB;
}


template <typename T>
Array2D<T> sigproc_multithreading::MiscUtils::converttoSTL(
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


Array2D<float> sigproc_multithreading::MiscUtils::converttoSTL(
  const ConcurrentArray2D<float>& inputTBB) const
{
  Array2D<float> outputSTL = converttoSTL<float>(inputTBB);
  return outputSTL;
}

Array2D<bool> sigproc_multithreading::MiscUtils::converttoSTL(
  const ConcurrentArray2D<bool>& inputTBB) const
{
  Array2D<bool> outputSTL = converttoSTL<bool>(inputTBB);
  return outputSTL;
}


#endif
