/**
 * \file BilateralFilters.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class BilateralFilters
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_BILATERALFILTERS_H__
#define __SIGPROC_TOOLS_BILATERALFILTERS_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "MiscUtils.h"

namespace sigproc_tools {

  /**
     \class BilateralFilters
     User defined class BilateralFilters ... these comments are used to generate
     doxygen documentation!
  */

  template <class T> using Array2D = std::vector<std::vector<T>>;

  class BilateralFilters{
    
    public:
      
      /// Default constructor
      BilateralFilters(){}

      void bilateral(
        const Array2D<float>& input2D,
        Array2D<float>& output2D,
        const unsigned int sx,
        const unsigned int sy,
        float sigma_x,
        float sigma_y,
        float sigma_r
      ) const;

      void guidedBilateral(
        const Array2D<float>& input2D,
        Array2D<float>& output2D,
        const unsigned int sx,
        const unsigned int sy
      ) const;
      
      /// Default destructor
      ~BilateralFilters(){}
      
    };
}

#endif
/** @} */ // end of doxygen group 

