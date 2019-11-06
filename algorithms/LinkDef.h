//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class algorithms::sample+;
#pragma link C++ class algorithms::NoiseRemoval+;
#pragma link C++ class algorithms::MorphologicalFilter+;
#pragma link C++ class std::vector<std::vector<double> >+;
//ADD_NEW_CLASS ... do not change this line
#endif



