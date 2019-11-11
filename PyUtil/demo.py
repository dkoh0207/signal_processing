from sproc import sproc
import numpy as np
import time,sys

t0 = time.time()
verbose=False
ctr=0
if 'verbose' in sys.argv:
    verbose=True

# This loop test 1D and 2D array conversions.
while 1:
    
    # Test 1D array
    np_array = np.random.random(size=(1000000)).astype(np.float32)
    # geneate C++ vector
    cpp_vec = sproc.pyutil.as_float32_vector(np_array)
    # generte np array back from C++ vector
    ap_array = sproc.pyutil.as_ndarray(cpp_vec)
    # if verbose, print out the size and mean
    if verbose:
        print(len(np_array),len(cpp_vec),len(ap_array))
        print(np.mean(np_array),np.mean(cpp_vec),np.mean(ap_array))
    # assert: since this is only conversion, neither size nor mean value should change
    assert len(np_array) == len(cpp_vec) and len(np_array) == len(ap_array)
    assert np.mean(np_array) == np.mean(ap_array)
    
    # Test 2D array
    np_array_2d = np.random.random(size=(1000,1000)).astype(np.float32)
    # geneate C++ vector
    cpp_vec_2d = sproc.pyutil.as_float32_vector_2d(np_array_2d)
    # generte np array back from C++ vector    
    ap_array_2d = sproc.pyutil.as_ndarray(cpp_vec_2d)
    # if verbose, print out the size and mean    
    if verbose:
        print(len(np_array_2d),len(cpp_vec),len(ap_array_2d))
        print(np.mean(np_array_2d),np.mean(cpp_vec_2d),np.mean(ap_array_2d))
    # assert: since this is only conversion, neither size nor mean value should change
    assert len(np_array_2d) == len(cpp_vec_2d) and len(np_array_2d) == len(ap_array_2d)
    assert np_array_2d.mean() == ap_array_2d.mean()
    t1 = time.time()
    ctr+=1
    if t1 - t0 > 20: break

print('Tested',ctr,'times during',t1-t0,'seconds!')

