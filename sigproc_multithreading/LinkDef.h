//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class sigproc_multithreading::CoherentNoiseCorrection+;
#pragma link C++ class sigproc_multithreading::ROIFinder+;
#pragma link C++ class sigproc_multithreading::Morph1DFast+;
#pragma link C++ class sigproc_multithreading::SigprocParallelDefs+;
#pragma link C++ class sigproc_multithreading::Morph2DFast+;
#pragma link C++ class sigproc_multithreading::MiscUtils+;
//ADD_NEW_CLASS ... do not change this line
#endif
