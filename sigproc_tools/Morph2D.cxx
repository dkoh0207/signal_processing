#ifndef __SIGPROC_TOOLS_MORPH2D_CXX__
#define __SIGPROC_TOOLS_MORPH2D_CXX__

#include "Morph2D.h"


void sigproc_tools::Morph2D::getFilter2D(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& dilation2D,
  std::vector<std::vector<short> >& erosion2D,
  std::vector<std::vector<short> >& average2D,
  std::vector<std::vector<short> >& gradient2D) const
{
  getFilter2D<short>(waveform2D, structuringElementx, structuringElementy,
    dilation2D, erosion2D, average2D, gradient2D);
  return;
}

void sigproc_tools::Morph2D::getFilter2D(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& dilation2D,
  std::vector<std::vector<float> >& erosion2D,
  std::vector<std::vector<float> >& average2D,
  std::vector<std::vector<float> >& gradient2D) const
{
  getFilter2D<float>(waveform2D, structuringElementx, structuringElementy,
    dilation2D, erosion2D, average2D, gradient2D);
  return;
}

void sigproc_tools::Morph2D::getFilter2D(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& dilation2D,
  std::vector<std::vector<double> >& erosion2D,
  std::vector<std::vector<double> >& average2D,
  std::vector<std::vector<double> >& gradient2D) const
{
  getFilter2D<double>(waveform2D, structuringElementx, structuringElementy,
    dilation2D, erosion2D, average2D, gradient2D);
  return;
}

template <typename T>
void sigproc_tools::Morph2D::getFilter2D(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& dilation2D,
  std::vector<std::vector<T> >& erosion2D,
  std::vector<std::vector<T> >& average2D,
  std::vector<std::vector<T> >& gradient2D) const
{
  auto numChannels = waveform2D.size();
  auto nTicks = waveform2D.at(0).size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  dilation2D.resize(waveform2D.size());
  erosion2D.resize(waveform2D.size());
  average2D.resize(waveform2D.size());
  gradient2D.resize(waveform2D.size());

  for (auto i=0; i<waveform2D.size(); ++i) {
    dilation2D[i].resize(waveform2D.at(0).size());
    erosion2D[i].resize(waveform2D.at(0).size());
    average2D[i].resize(waveform2D.at(0).size());
    gradient2D[i].resize(waveform2D.at(0).size());
  }
  float dilation;
  float erosion;
  float gradient;
  float average;
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
      std::vector<T> v;
      for (auto ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (auto iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v.push_back(waveform2D[ix][iy]);
        }
      }
      dilation = *std::max_element(v.begin(), v.end());
      erosion = *std::min_element(v.begin(), v.end());
      average = 0.5 * (dilation + erosion);
      gradient = dilation - erosion;
      dilation2D[i][j] = dilation;
      erosion2D[i][j] = erosion;
      average2D[i][j] = average;
      gradient2D[i][j] = gradient;
    }
  }
  return;
}


void sigproc_tools::Morph2D::getGradient(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& gradient2D) const
{
  getGradient<short>(waveform2D, structuringElementx, 
    structuringElementy, gradient2D);
  return;
}

void sigproc_tools::Morph2D::getGradient(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& gradient2D) const
{
  getGradient<float>(waveform2D, structuringElementx, 
    structuringElementy, gradient2D);
  return;
}

void sigproc_tools::Morph2D::getGradient(
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
void sigproc_tools::Morph2D::getGradient(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& gradient2D) const
{
  auto numChannels = waveform2D.size();
  auto nTicks = waveform2D.at(0).size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  gradient2D.resize(waveform2D.size());

  for (auto i=0; i<waveform2D.size(); ++i) {
    gradient2D[i].resize(waveform2D.at(0).size());
  }
  float dilation;
  float erosion;
  float gradient;
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
      std::vector<T> v;
      for (auto ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (auto iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v.push_back(waveform2D[ix][iy]);
        }
      }
      dilation = *std::max_element(v.begin(), v.end());
      erosion = *std::min_element(v.begin(), v.end());
      gradient = dilation - erosion;
      gradient2D[i][j] = gradient;
    }
  }
  return;
}


void sigproc_tools::Morph2D::getOpeningAndClosing(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& opening2D,
  std::vector<std::vector<short> >& closing2D) const
{
  getOpeningAndClosing<short>(
    waveform2D, 
    structuringElementx, 
    structuringElementy, 
    opening2D,
    closing2D);
  return;
}

void sigproc_tools::Morph2D::getOpeningAndClosing(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& opening2D,
  std::vector<std::vector<float> >& closing2D) const
{
  getOpeningAndClosing<float>(
    waveform2D, 
    structuringElementx, 
    structuringElementy, 
    opening2D,
    closing2D);
  return;
}

void sigproc_tools::Morph2D::getOpeningAndClosing(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& opening2D,
  std::vector<std::vector<double> >& closing2D) const
{
  getOpeningAndClosing<double>(
    waveform2D, 
    structuringElementx, 
    structuringElementy, 
    opening2D,
    closing2D);
  return;
}


template <typename T>
void sigproc_tools::Morph2D::getOpeningAndClosing(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& opening2D,
  std::vector<std::vector<T> >& closing2D) const
{
  auto numChannels = waveform2D.size();
  auto nTicks = waveform2D.at(0).size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  opening2D.resize(waveform2D.size());
  closing2D.resize(waveform2D.size());

  std::vector<std::vector<T> > dilation2D;
  std::vector<std::vector<T> > erosion2D;
  std::vector<std::vector<T> > average2D;
  std::vector<std::vector<T> > gradient2D;

  getFilter2D(waveform2D, structuringElementx, structuringElementy,
              dilation2D, erosion2D, average2D, gradient2D);

  for (auto i=0; i<waveform2D.size(); ++i) {
    opening2D[i].resize(waveform2D.at(0).size());
    closing2D[i].resize(waveform2D.at(0).size());
  }
  float opening;
  float closing;
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
      std::vector<T> v1;
      std::vector<T> v2;
      for (auto ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (auto iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v1.push_back(dilation2D[ix][iy]);
          v2.push_back(erosion2D[ix][iy]);
        }
      }
      opening = *std::max_element(v2.begin(), v2.end());
      closing = *std::min_element(v1.begin(), v1.end());
      opening2D[i][j] = opening;
      closing2D[i][j] = closing;
    }
  }
  return;
}

#endif
