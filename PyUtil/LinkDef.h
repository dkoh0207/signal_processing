//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ namespace sproc+;
#pragma link C++ namespace sproc::pyutil+;
#pragma link C++ class sproc::pyutil::load_pyutil+;
#ifndef __CINT__
#pragma link C++ function pyutil::as_float32_vector(PyObject*)+;
#pragma link C++ function pyutil::as_float32_vector_2d(PyObject*)+;
#endif
//ADD_NEW_CLASS ... do not change this line
#endif
