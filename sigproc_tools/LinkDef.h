//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class sigproc_tools::sample+;
#pragma link C++ class Morph1D::sigproc_tools+;
#pragma link C++ class Morph2D::sigproc_tools+;
#pragma link C++ class MiscUtils::sigproc_tools+;
#pragma link C++ class Denoising::sigproc_tools+;
//ADD_NEW_CLASS ... do not change this line
#endif




