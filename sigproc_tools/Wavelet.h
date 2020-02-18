/**
 * \file Wavelet.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class Wavelet
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_WAVELET_H__
#define __SIGPROC_TOOLS_WAVELET_H__

#include <vector>
#include <cmath>

namespace sigproc_tools {

  /**
     \class Wavelet
     User defined class Wavelet ... these comments are used to generate
     doxygen documentation!
  */
  class Wavelet{
    
    public:

      // virtual void filt(
      //   std::vector<float>& inputWaveform,
      //   const size_t n, 
      //   const int isign) const = 0 {};

      virtual void condition(
        std::vector<float>& inputWaveform,
        const size_t n,
        const int isign) const = 0;

      virtual ~Wavelet(){};
    
  };

  class Daubechies4 {

    public: 

      std::vector<float> filters;
      unsigned int filterLength;

      Daubechies4(const unsigned int n=4);

      void transformPeriodic(
        std::vector<float>& inputWaveform,
        const size_t n,
        const int isign) const;

      void transform(
        std::vector<float>& inputWaveform,
        const size_t n,
        const int isign) const;

      void condition(
        std::vector<float>& inputWaveform,
        const size_t n,
        const int isign) const;

      ~Daubechies4(){};
  };

  class Haar {

    public:

      std::vector<float> filters;
      unsigned int filterLength;

      Haar(){};

      void transform(
        std::vector<float>& inputWaveform,
        const size_t n,
        const int isign) const;

      ~Haar(){};
  };
}

#endif
/** @} */ // end of doxygen group 

