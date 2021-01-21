#ifndef __SIGPROC_TOOLS_MORPH2DFAST_CXX__
#define __SIGPROC_TOOLS_MORPH2DFAST_CXX__

#include "Morph2DFast.h"


void sigproc_tools::Morph2DFast::getDilation(
  const std::vector<std::vector<bool> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<bool> >& dilation2D) const
{
  getDilation<bool>(waveform2D, structuringElementx,
    structuringElementy, dilation2D);
  return;
}

void sigproc_tools::Morph2DFast::getDilation(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& dilation2D) const
{
  getDilation<short>(waveform2D, structuringElementx,
    structuringElementy, dilation2D);
  return;
}

void sigproc_tools::Morph2DFast::getDilation(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& dilation2D) const
{
  getDilation<float>(waveform2D, structuringElementx,
    structuringElementy, dilation2D);
  return;
}

void sigproc_tools::Morph2DFast::getDilation(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& dilation2D) const
{
  getDilation<double>(waveform2D, structuringElementx,
    structuringElementy, dilation2D);
  return;
}

template <typename T>
void sigproc_tools::Morph2DFast::getDilation(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& dilation2D) const
{
  size_t numChannels = waveform2D.size();
  size_t nTicks = waveform2D.at(0).size();

  assert(dilation2D.size() == numChannels);
  assert(dilation2D.at(0).size() == nTicks);

  sigproc_tools::Morph1DFast fast1D;

  for (size_t i=0; i<numChannels; ++i) {
    fast1D.getDilation(waveform2D[i], structuringElementy, dilation2D[i]);
  }
  for (size_t j=0; j<nTicks; ++j) {
    std::vector<T> column(numChannels);
    std::vector<T> columnOut(numChannels);
    for (size_t i=0; i<numChannels; ++i) {
      column[i] = dilation2D[i][j];
    }
    fast1D.getDilation(column, structuringElementx, columnOut);
    for (size_t i=0; i<numChannels; ++i) {
      dilation2D[i][j] = columnOut[i];
    }
  }
  return;
}


void sigproc_tools::Morph2DFast::getErosion(
  const std::vector<std::vector<bool> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<bool> >& erosion2D) const
{
  getErosion<bool>(waveform2D, structuringElementx,
    structuringElementy, erosion2D);
  return;
}

void sigproc_tools::Morph2DFast::getErosion(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& erosion2D) const
{
  getErosion<short>(waveform2D, structuringElementx,
    structuringElementy, erosion2D);
  return;
}

void sigproc_tools::Morph2DFast::getErosion(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& erosion2D) const
{
  getErosion<float>(waveform2D, structuringElementx,
    structuringElementy, erosion2D);
  return;
}

void sigproc_tools::Morph2DFast::getErosion(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& erosion2D) const
{
  getErosion<double>(waveform2D, structuringElementx,
    structuringElementy, erosion2D);
  return;
}

template <typename T>
void sigproc_tools::Morph2DFast::getErosion(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& erosion2D) const
{
  size_t numChannels = waveform2D.size();
  size_t nTicks = waveform2D.at(0).size();

  assert(erosion2D.size() == numChannels);
  assert(erosion2D.at(0).size() == nTicks);

  sigproc_tools::Morph1DFast fast1D;

  for (size_t i=0; i<numChannels; ++i) {
    fast1D.getErosion(waveform2D[i], structuringElementy, erosion2D[i]);
  }
  for (size_t j=0; j<nTicks; ++j) {
    std::vector<T> column(numChannels);
    std::vector<T> columnOut(numChannels);
    for (size_t i=0; i<numChannels; ++i) {
      column[i] = erosion2D[i][j];
    }
    fast1D.getErosion(column, structuringElementx, columnOut);
    for (size_t i=0; i<numChannels; ++i) {
      erosion2D[i][j] = columnOut[i];
    }
  }
  return;
}


void sigproc_tools::Morph2DFast::getGradient(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& gradient2D) const
{
  getGradient<short>(waveform2D, structuringElementx,
    structuringElementy, gradient2D);
  return;
}

void sigproc_tools::Morph2DFast::getGradient(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& gradient2D) const
{
  getGradient<float>(waveform2D, structuringElementx,
    structuringElementy, gradient2D);
  return;
}

void sigproc_tools::Morph2DFast::getGradient(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& gradient2D) const
{
  getGradient<double>(waveform2D, structuringElementx,
    structuringElementy, gradient2D);
  return;
}

template <typename T>
void sigproc_tools::Morph2DFast::getGradient(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& gradient2D) const
{
  size_t numChannels = waveform2D.size();
  size_t nTicks = waveform2D.at(0).size();

  assert(gradient2D.size() == numChannels);
  assert(gradient2D.at(0).size() == nTicks);

  std::vector<std::vector<T> > dilation2D;
  dilation2D.resize(numChannels);
  for (auto& v : dilation2D) {
    v.resize(nTicks);
  }
  std::vector<std::vector<T> > erosion2D;
  erosion2D.resize(numChannels);
  for (auto& v : erosion2D) {
    v.resize(nTicks);
  }

  getDilation<T>(waveform2D, structuringElementx, structuringElementy, dilation2D);
  getErosion<T>(waveform2D, structuringElementx, structuringElementy, erosion2D);

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      gradient2D[i][j] = (dilation2D[i][j] - erosion2D[i][j]);
    }
  }
  return;
}


void sigproc_tools::Morph2DFast::getClosing(
  const std::vector<std::vector<bool> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<bool> >& closing2D) const
{
  getClosing<bool>(waveform2D, structuringElementx,
    structuringElementy, closing2D);
  return;
}

void sigproc_tools::Morph2DFast::getClosing(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& closing2D) const
{
  getClosing<short>(waveform2D, structuringElementx,
    structuringElementy, closing2D);
  return;
}

void sigproc_tools::Morph2DFast::getClosing(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& closing2D) const
{
  getClosing<float>(waveform2D, structuringElementx,
    structuringElementy, closing2D);
  return;
}

void sigproc_tools::Morph2DFast::getClosing(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& closing2D) const
{
  getClosing<double>(waveform2D, structuringElementx,
    structuringElementy, closing2D);
  return;
}

template <typename T>
void sigproc_tools::Morph2DFast::getClosing(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& closing2D) const
{
  size_t numChannels = waveform2D.size();
  size_t numTicks = waveform2D.at(0).size();

  assert(closing2D.size() == numChannels);
  assert(closing2D.at(0).size() == numTicks);

  std::vector<std::vector<T>> temp(numChannels);
  for (auto& v : temp) {
    v.resize(numTicks);
  }
  getDilation<T>(waveform2D, structuringElementx, structuringElementy, temp);
  getErosion<T>(temp, structuringElementx, structuringElementy, closing2D);
  return;
}


void sigproc_tools::Morph2DFast::getOpening(
  const std::vector<std::vector<bool> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<bool> >& opening2D) const
{
  getOpening<bool>(waveform2D, structuringElementx,
    structuringElementy, opening2D);
  return;
}

void sigproc_tools::Morph2DFast::getOpening(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& opening2D) const
{
  getOpening<short>(waveform2D, structuringElementx,
    structuringElementy, opening2D);
  return;
}

void sigproc_tools::Morph2DFast::getOpening(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& opening2D) const
{
  getOpening<float>(waveform2D, structuringElementx,
    structuringElementy, opening2D);
  return;
}

void sigproc_tools::Morph2DFast::getOpening(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& opening2D) const
{
  getOpening<double>(waveform2D, structuringElementx,
    structuringElementy, opening2D);
  return;
}

template <typename T>
void sigproc_tools::Morph2DFast::getOpening(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& opening2D) const
{
  size_t numChannels = waveform2D.size();
  size_t numTicks = waveform2D.at(0).size();

  assert(opening2D.size() == numChannels);
  assert(opening2D.at(0).size() == numTicks);

  std::vector<std::vector<T>> temp(numChannels);
  for (auto& v : temp) {
    v.resize(numTicks);
  }
  getErosion<T>(waveform2D, structuringElementx, structuringElementy, temp);
  getDilation<T>(temp, structuringElementx, structuringElementy, opening2D);

  return;
}

#endif
