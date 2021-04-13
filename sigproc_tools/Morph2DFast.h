/**
 * \file Morph2DFast.h
 *
 * \ingroup sigproc_tools
 *
 * \brief Class def header for a class Morph2DFast
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_MORPH2DFAST_H__
#define __SIGPROC_TOOLS_MORPH2DFAST_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <assert.h>
#include <functional>
#include "Morph1DFast.h"

namespace sigproc_tools {

  /**
     \class Morph2DFast
     User defined class Morph2DFast ... these comments are used to generate
     doxygen documentation!
  */
  class Morph2DFast{

    public:

      /// Default constructor
      Morph2DFast(){}

        void getDilation(const std::vector<std::vector<bool> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<bool> >&) const;

        void getDilation(const std::vector<std::vector<short> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<short> >&) const;

        void getDilation(const std::vector<std::vector<float> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<float> >&) const;

        void getDilation(const std::vector<std::vector<double> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<double> >&) const;

        void getDilationFast(const std::vector<std::vector<float> >&,
                            const unsigned int,
                            const unsigned int,
                            std::vector<std::vector<float> >&) const;


        void getErosion(const std::vector<std::vector<bool> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<bool> >&) const;

        void getErosion(const std::vector<std::vector<short> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<short> >&) const;

        void getErosion(const std::vector<std::vector<float> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<float> >&) const;

        void getErosion(const std::vector<std::vector<double> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<double> >&) const;


        void getGradient(const std::vector<std::vector<short> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<short> >&) const;

        void getGradient(const std::vector<std::vector<float> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<float> >&) const;

        void getGradient(const std::vector<std::vector<double> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<double> >&) const;


        void updateMedian(short& medianVal,
                          unsigned int& numLowPixels,
                          std::vector<short>& hist,
                          unsigned int halfPixels,
                          unsigned int valRange) const;

        void updateMedian(int& medianVal,
                          unsigned int& numLowPixels,
                          std::vector<int>& hist,
                          unsigned int halfPixels,
                          unsigned int valRange) const;


        void getMedian(const std::vector<std::vector<short> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<short> >&) const;

        void getMedian(const std::vector<std::vector<int> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<int> >&) const;


        void getOpening(const std::vector<std::vector<bool> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<bool> >&) const;

        void getOpening(const std::vector<std::vector<short> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<short> >&) const;

        void getOpening(const std::vector<std::vector<float> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<float> >&) const;

        void getOpening(const std::vector<std::vector<double> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<double> >&) const;


        void getClosing(const std::vector<std::vector<bool> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<bool> >&) const;

        void getClosing(const std::vector<std::vector<short> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<short> >&) const;

        void getClosing(const std::vector<std::vector<float> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<float> >&) const;

        void getClosing(const std::vector<std::vector<double> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<double> >&) const;

        void contrastStretching(const std::vector<std::vector<float> >&,
                                const unsigned int,
                                const unsigned int,
                                std::vector<std::vector<float> >&,
                                const float) const;
      /// Default destructor
      ~Morph2DFast(){}


    private:

      template <typename T>
      void getDilation(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& dilation2D) const;

      template <typename T>
      void getErosion(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& erosion2D) const;

      template <typename T>
      void getGradient(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& gradient2D) const;

      template <typename T>
      void updateMedian(T& medianVal,
                        unsigned int& numLowPixels,
                        std::vector<T>& hist,
                        unsigned int halfPixels,
                        unsigned int valRange) const;

      template <typename T>
      void getMedian(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& median2D,
        const unsigned int pixelRange) const;

      template <typename T>
      void getClosing(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& closing2D) const;

      template <typename T>
      void getOpening(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& opening2D) const;

  };
}

#endif
/** @} */ // end of doxygen group
