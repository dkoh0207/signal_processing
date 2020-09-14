//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class sigproc_tools::Morph1D+;
#pragma link C++ class sigproc_tools::Morph2D+;
#pragma link C++ class sigproc_tools::MiscUtils+;
#pragma link C++ class sigproc_tools::Denoising+;
#pragma link C++ class sigproc_tools::AdaptiveWiener+;
#pragma link C++ class sigproc_tools::Deconvolution+;
#pragma link C++ class sigproc_tools::Morph1DFast+;
#pragma link C++ class sigproc_tools::Morph2DFast+;
#pragma link C++ class sigproc_tools::MorphologicalCNC+;
#pragma link C++ class sigproc_tools::FourierCNC+;
#pragma link C++ class FrequencyFilters::sigproc_tools+;
//ADD_NEW_CLASS ... do not change this line
#endif

