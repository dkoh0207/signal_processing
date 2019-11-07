/**
 * \file WaveformUtils.h
 *
 * \ingroup algorithms
 * 
 * \brief Class def header for a class WaveformUtils
 *
 * @author koh0207
 */

/** \addtogroup algorithms

    @{*/
#ifndef __ALGORITHMS_WAVEFORMUTILS_H__
#define __ALGORITHMS_WAVEFORMUTILS_H__

#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>

namespace algorithms {

  /**
     \class WaveformUtils
     User defined class WaveformUtils ... these comments are used to generate
     doxygen documentation!
  */

  template <class T> using Waveform = std::vector<T>;

  class WaveformUtils{
    
  public:
    
    /// Default constructor
    WaveformUtils(){}

    void getWaveformParams(const std::vector<short>& waveform,
                                                 float& mean,
                                                 float& median,
                                                 float& mode,
                                                 float& skewness,
                                                 float& rms);
    
    void getErosionDilationAverageDifference(const Waveform<short>&,
                                             int,
                                             Waveform<short>&,
                                             Waveform<short>&,
                                             Waveform<short>&,
                                             Waveform<short>&)         const;

    void getErosionDilationAverageDifference(const Waveform<float>&,
                                             int,
                                             Waveform<float>&,
                                             Waveform<float>&,
                                             Waveform<float>&,
                                             Waveform<float>&)         const;

    void getErosionDilationAverageDifference(const Waveform<double>&,
                                             int,
                                             Waveform<double>&,
                                             Waveform<double>&,
                                             Waveform<double>&,
                                             Waveform<double>&)         const;

    void getOpeningAndClosing(const Waveform<short>&,  const Waveform<short>&,  int, Waveform<short>&,  Waveform<short>&)  const;
    void getOpeningAndClosing(const Waveform<float>&,  const Waveform<float>&,  int, Waveform<float>&,  Waveform<float>&)  const;
    void getOpeningAndClosing(const Waveform<double>&, const Waveform<double>&, int, Waveform<double>&, Waveform<double>&) const;

    /// Default destructor
    ~WaveformUtils(){}
    
  private:
    template <typename T> void triangleSmooth(const std::vector<T>&, std::vector<T>&, size_t = 0)                                        const;
    template <typename T> void medianSmooth(  const std::vector<T>&, std::vector<T>&, size_t = 3)                                        const;
    template <typename T> void getTruncatedMean(const std::vector<T>&, T&, int&)                                                         const;
    template <typename T> void getTruncatedMeanRMS(const std::vector<T>&, T, T&, T&, T&, int&)                                           const;
    template <typename T> void firstDerivative(const std::vector<T>&,  std::vector<T>&)                                                  const;
    template <typename T> void getErosionDilationAverageDifference(const Waveform<T>&,
                                                                   int,
                                                                   Waveform<T>&,
                                                                   Waveform<T>&,
                                                                   Waveform<T>&,
                                                                   Waveform<T>&) const;

    template <typename T> void getOpeningAndClosing(const Waveform<T>&,  const Waveform<T>&,  int, Waveform<T>&,  Waveform<T>&)  const;
  };
}

#endif
/** @} */ // end of doxygen group 

