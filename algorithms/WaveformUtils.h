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

    void getWaveformParams(const std::vector<short>&,
                                                 float&,
                                                 float&,
                                                 float&,
                                                 float&,
                                                 float&);

    void getWaveformParams(const std::vector<float>&,
                                                 float&,
                                                 float&,
                                                 float&,
                                                 float&,
                                                 float&);

    void getWaveformParams(const std::vector<double>&,
                                                 float&,
                                                 float&,
                                                 float&,
                                                 float&,
                                                 float&);
    
    void getDilation(const Waveform<short>&,
                     const unsigned int,
                     Waveform<short>&) const;

    void getDilation(const Waveform<float>&,
                     const unsigned int,
                     Waveform<float>&) const;

    void getDilation(const Waveform<double>&,
                     const unsigned int,
                     Waveform<double>&) const;
    
    void getErosionDilationAverageDifference(const Waveform<short>&,
                                             const unsigned int,
                                             Waveform<short>&,
                                             Waveform<short>&,
                                             Waveform<short>&,
                                             Waveform<short>&)         const;

    void getErosionDilationAverageDifference(const Waveform<float>&,
                                             const unsigned int,
                                             Waveform<float>&,
                                             Waveform<float>&,
                                             Waveform<float>&,
                                             Waveform<float>&)         const;

    void getErosionDilationAverageDifference(const Waveform<double>&,
                                             const unsigned int,
                                             Waveform<double>&,
                                             Waveform<double>&,
                                             Waveform<double>&,
                                             Waveform<double>&)         const;

    void getOpeningAndClosing(const Waveform<short>&,  const Waveform<short>&,  const unsigned int, Waveform<short>&,  Waveform<short>&)  const;
    void getOpeningAndClosing(const Waveform<float>&,  const Waveform<float>&,  const unsigned int, Waveform<float>&,  Waveform<float>&)  const;
    void getOpeningAndClosing(const Waveform<double>&, const Waveform<double>&, const unsigned int, Waveform<double>&, Waveform<double>&) const;

    /// Default destructor
    ~WaveformUtils(){}
    
  private:

    template <typename T> void getWaveformParams(const std::vector<T>& waveform,
                                                 float& mean,
                                                 float& median,
                                                 float& mode,
                                                 float& skewness,
                                                 float& rms);

    template <typename T> void getDilation(const Waveform<T>& inputWaveform,
                                           const unsigned int structuringElement,
                                           Waveform<T>& dilationVec) const;
    template <typename T> void getErosionDilationAverageDifference(const Waveform<T>&,
                                                                   const unsigned int,
                                                                   Waveform<T>&,
                                                                   Waveform<T>&,
                                                                   Waveform<T>&,
                                                                   Waveform<T>&) const;

    template <typename T> void getOpeningAndClosing(const Waveform<T>&,  const Waveform<T>&,  const unsigned int, Waveform<T>&,  Waveform<T>&)  const;
  };
}

#endif
/** @} */ // end of doxygen group 

