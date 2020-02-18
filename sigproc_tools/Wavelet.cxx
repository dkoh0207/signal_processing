#ifndef __SIGPROC_TOOLS_WAVELET_CXX__
#define __SIGPROC_TOOLS_WAVELET_CXX__

#include "Wavelet.h"

// Haar Wavelets

void sigproc_tools::Haar::transform(
  std::vector<float>& wf,
  const size_t n,
  const int isign) const
{
  /*
  From Numerical Recipes:
  Press, William H. and Teukolsky, Saul A. and Vetterling, William T. and Flannery, Brian P.
  Numerical Recipes 3rd Edition: The Art of Scientific Computing
  Chapter 13.11 Wavelet Transforms

  1D Wavelet Transform with Daubechies 4 filter with 
  periodic boundary conditions.

  INPUTS:
    - wf: input waveform to be transformed
    - isign: parameter for forward/inverse transform mode (i = 1 or -1)
    - wksp: transformed waveform. 
  */

  std::vector<float> wksp(n);
  const std::vector<float> C = this->filters;
  float invSqrt = 0.70710678118;

  size_t nh, i, j;
  nh = n >> 1;

  if (isign >= 0) {
    for (i=0, j=0; j<n-1; j+=2, ++i) {
        // Compute and store scaling coefficients
        wksp[i] = invSqrt * (wf[j] + wf[j+1]);
        // Compute and store detail coefficients
        wksp[i + nh] = invSqrt * (wf[j] - wf[j+1]);
      }
  } else {
    for (i=0, j=0; i<nh-1; ++i) {
      wksp[j++] = invSqrt * (wf[i] + wf[i+nh]);
      wksp[j++] = invSqrt * (wf[i] - wf[i+nh]);
    }
  }
  for (i=0; i<n; ++i) wf[i] = wksp[i];
  return;
}

// Daubechies Family of Wavelets

sigproc_tools::Daubechies4::Daubechies4(const unsigned int n)
{
  /*
  Initialize the Daubechies family of wavelet filters for different n.
  */
  if (n == 4) {
    // Daubechies 4 coefficient filter
    this->filters.resize(4);
    this->filters[0] = 0.4829629131445341;
    this->filters[1] = 0.8365163037378077;
    this->filters[2] = 0.2241438680420134;
    this->filters[3] = -0.1294095225512603;
  }
  // Include Db8, Db16, etc. 
  return;
}

void sigproc_tools::Daubechies4::transformPeriodic(
  std::vector<float>& wf,
  const size_t n,
  const int isign) const
{
  /*
  From Numerical Recipes:
  Press, William H. and Teukolsky, Saul A. and Vetterling, William T. and Flannery, Brian P.
  Numerical Recipes 3rd Edition: The Art of Scientific Computing
  Chapter 13.11 Wavelet Transforms

  1D Wavelet Transform with Daubechies 4 filter with 
  periodic boundary conditions.

  INPUTS:
    - wf: input waveform to be transformed
    - isign: parameter for forward/inverse transform mode (i = 1 or -1)
    - wksp: transformed waveform. 
  */

  if (n < 4) return;
  std::vector<float> wksp(n);
  const std::vector<float> C = this->filters;

  size_t nh, i, j;
  nh = n >> 1;

  if (isign >= 0) {
    for (i=0, j=0; j<n-3; j+=2, ++i) {
        // Compute and store scaling coefficients
        wksp[i] = C[0]*wf[j] + C[1]*wf[j+1] + C[2]*wf[j+2] + C[3]*wf[j+3];
        // Compute and store detail coefficients
        wksp[i + nh] = C[3]*wf[j]-C[2]*wf[j+1]+C[1]*wf[j+2]-C[0]*wf[j+3];
      }
      // Handle last two transforms
    wksp[i] = C[0]*wf[n-2]+C[1]*wf[n-1]+C[2]*wf[0]+C[3]*wf[1];
    wksp[i+nh] = C[3]*wf[n-2]-C[2]*wf[n-1]+C[1]*wf[0]-C[0]*wf[1];
  } else {
    wksp[0] = C[2]*wf[nh-1]+C[1]*wf[n-1]+C[0]*wf[0]+C[3]*wf[1];
    wksp[1] = C[3]*wf[nh-1]-C[0]*wf[n-1]+C[1]*wf[0]-C[0]*wf[1];
    for (i=0, j=2; i<nh-1; ++i) {
      wksp[j++] = C[2]*wf[i]+C[1]*wf[i+nh]+C[0]*wf[i+1]+C[3]*wf[i+nh+1];
      wksp[j++] = C[3]*wf[i]-C[0]*wf[i+nh]+C[1]*wf[i+1]-C[2]*wf[i+nh+1];
    }
  }
  for (i=0; i<n; ++i) wf[i] = wksp[i];
  return;
}

void sigproc_tools::Daubechies4::transform(
  std::vector<float>& wf,
  const size_t n,
  const int isign) const
{
  /*
  The Cohen-Daubechies-Vial Wavelet 4 Filter on the interval. 
  */
  const float C0=0.4829629131445341, C1=0.8365163037378077,
  C2=0.2241438680420134, C3=-0.1294095225512603;
  const float R00=0.603332511928053,R01=0.690895531839104,
  R02=-0.398312997698228,R10=-0.796543516912183,R11=0.546392713959015,
  R12=-0.258792248333818,R20=0.0375174604524466,R21=0.457327659851769,
  R22=0.850088102549165,R23=0.223820356983114,R24=-0.129222743354319,
  R30=0.0100372245644139,R31=0.122351043116799,R32=0.227428111655837,
  R33=-0.836602921223654,R34=0.483012921773304,R43=0.443149049637559,
  R44=0.767556669298114,R45=0.374955331645687,R46=0.190151418429955,
  R47=-0.194233407427412,R53=0.231557595006790,R54=0.401069519430217,
  R55=-0.717579999353722,R56=-0.363906959570891,R57=0.371718966535296,
  R65=0.230389043796969,R66=0.434896997965703,R67=0.870508753349866,
  R75=-0.539822500731772,R76=0.801422961990337,R77=-0.257512919478482;

  size_t nh, i, j;
  if (n < 8) return;
  nh = n >> 1;
  std::vector<float> wksp(n);
  if (isign >= 0) {
    wksp[0] = R00*wf[0]+R01*wf[1]+R02*wf[2];
    wksp[nh] = R10*wf[0]+R11*wf[1]+R12*wf[2];
    wksp[1] = R20*wf[0]+R21*wf[1]+R22*wf[2]+R23*wf[3]+R24*wf[4];
    wksp[nh+1] = R30*wf[0]+R31*wf[1]+R32*wf[2]+R33*wf[3]+R34*wf[4];
      for (i=2,j=3;j<n-4;j+=2,i++) {
        wksp[i] = C0*wf[j]+C1*wf[j+1]+C2*wf[j+2]+C3*wf[j+3];
        wksp[i+nh] = C3*wf[j]-C2*wf[j+1]+C1*wf[j+2]-C0*wf[j+3];
      }
    wksp[nh-2] = R43*wf[n-5]+R44*wf[n-4]+R45*wf[n-3]+R46*wf[n-2]+R47*wf[n-1];
    wksp[n-2] = R53*wf[n-5]+R54*wf[n-4]+R55*wf[n-3]+R56*wf[n-2]+R57*wf[n-1];
    wksp[nh-1] = R65*wf[n-3]+R66*wf[n-2]+R67*wf[n-1];
    wksp[n-1] = R75*wf[n-3]+R76*wf[n-2]+R77*wf[n-1];
  } else {
    wksp[0] = R00*wf[0]+R10*wf[nh]+R20*wf[1]+R30*wf[nh+1];
    wksp[1] = R01*wf[0]+R11*wf[nh]+R21*wf[1]+R31*wf[nh+1];
    wksp[2] = R02*wf[0]+R12*wf[nh]+R22*wf[1]+R32*wf[nh+1];
    if (n == 8) {
      wksp[3] = R23*wf[1]+R33*wf[5]+R43*wf[2]+R53*wf[6];
      wksp[4] = R24*wf[1]+R34*wf[5]+R44*wf[2]+R54*wf[6];
    } else {
      wksp[3] = R23*wf[1]+R33*wf[nh+1]+C0*wf[2]+C3*wf[nh+2];
      wksp[4] = R24*wf[1]+R34*wf[nh+1]+C1*wf[2]-C2*wf[nh+2];
      wksp[n-5] = C2*wf[nh-3]+C1*wf[n-3]+R43*wf[nh-2]+R53*wf[n-2];
      wksp[n-4] = C3*wf[nh-3]-C0*wf[n-3]+R44*wf[nh-2]+R54*wf[n-2];
    }
    for (i=2,j=5;i<nh-3;i++) {
      wksp[j++] = C2*wf[i]+C1*wf[i+nh]+C0*wf[i+1]+C3*wf[i+nh+1];
      wksp[j++] = C3*wf[i]-C0*wf[i+nh]+C1*wf[i+1]-C2*wf[i+nh+1];
    }
    wksp[n-3] = R45*wf[nh-2]+R55*wf[n-2]+R65*wf[nh-1]+R75*wf[n-1];
    wksp[n-2] = R46*wf[nh-2]+R56*wf[n-2]+R66*wf[nh-1]+R76*wf[n-1];
    wksp[n-1] = R47*wf[nh-2]+R57*wf[n-2]+R67*wf[nh-1]+R77*wf[n-1];
  }
  for (i=0; i<n; ++i) wf[i] = wksp[i];
  return;
}

void sigproc_tools::Daubechies4::condition(
  std::vector<float>& wf,
  const size_t n,
  const int isign) const
{
  float t0, t1, t2, t3;
  if (n < 4) return;
  if (isign >= 0) {
    t0 = 0.324894048898962*wf[0]+0.0371580151158803*wf[1];
    t1 = 1.00144540498130*wf[1];
    t2 = 1.08984305289504*wf[n-2];
    t3 = -0.800813234246437*wf[n-2]+2.09629288435324*wf[n-1];
    wf[0]=t0; wf[1]=t1; wf[n-2]=t2; wf[n-1]=t3;
  } else {
    t0 = 3.07792649138669*wf[0]-0.114204567242137*wf[1];
    t1 = 0.998556681198888*wf[1];
    t2 = 0.917563310922261*wf[n-2];
    t3 = 0.350522032550918*wf[n-2]+0.477032578540915*wf[n-1];
    wf[0]=t0; wf[1]=t1; wf[n-2]=t2; wf[n-1]=t3;
  }
}


#endif
