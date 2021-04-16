/**
 * \file MiscUtils.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class MiscUtils
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_MISCUTILS_H__
#define __SIGPROC_TOOLS_MISCUTILS_H__

#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include <assert.h>

namespace sigproc_tools {

  /**
     \class MiscUtils
     Miscellaneous utility functions. 
  */

  template <class T> using Array2D = std::vector<std::vector<T>>;

  class MiscUtils{
    
    public:
      
      /// Default constructor
      MiscUtils(){}

      short computeMedian(const std::vector<short>& vec);
      float computeMedian(const std::vector<float>& vec);
      double computeMedian(const std::vector<double>& vec);

      short computeMaximum(const Array2D<short>& input2D);
      float computeMaximum(const Array2D<float>& input2D);
      double computeMaximum(const Array2D<double>& input2D);

      short computeMinimum(const Array2D<short>& input2D);
      float computeMinimum(const Array2D<float>& input2D);
      double computeMinimum(const Array2D<double>& input2D);

      float estimateNoiseVariance(
        const std::vector<std::vector<float>>& waveLessCoherent,
        const std::vector<std::vector<bool>>& selectVals);

      float estimateMAD(
        const std::vector<float>& wf);

      unsigned nChoosek(unsigned n, unsigned k) const;

      void drawIndextoImage(const std::vector<int>& rows,
        const std::vector<int>& cols,
        Array2D<bool>& output2D) const;

      void drawIndextoImage(const std::vector<size_t>& rows,
        const std::vector<size_t>& cols,
        Array2D<bool>& output2D) const;

      void drawIndextoImage(const std::vector<unsigned int>& rows,
        const std::vector<unsigned int>& cols,
        Array2D<bool>& output2D) const;

      void drawIndextoImage(const std::vector<short>& rows,
        const std::vector<short>& cols,
        Array2D<bool>& output2D) const;

      
      /// Default destructor
      ~MiscUtils(){}
    
    private:

      template <typename T>
      T computeMedian(const std::vector<T>& vec);

      template <typename T>
      T computeMaximum(const Array2D<T>& input2D);

      template <typename T>
      T computeMinimum(const Array2D<T>& input2D);

      template <typename T>
      void drawIndextoImage(const std::vector<T>& rows,
        const std::vector<T>& cols,
        Array2D<bool>& output2D) const;
      // template <typename T> T computeMedian(
      //   const std::vector<T>& waveform
      // );
  };
}

#endif
/** @} */ // end of doxygen group 

