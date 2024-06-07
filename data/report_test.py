# Test the netcdf
from netCDF4 import Dataset
import pandas as pd
from numpy import linalg as la
from numpy import sum as nps
import numpy as np

## 1  Test of the deformation calculations
print(' ')
print(" --- 1  Test of the deformation calculations ---")
print(' ')


file_ref  = 'test_ref/SIDRR_20170901.nc'
file_test = 'test_out/SIDRR_20170901.nc'

data_ref = Dataset(file_ref)
data_test = Dataset(file_test)

for var in data_ref.variables:
    if var in data_test.variables:
        test = np.all(data_ref[var][:] == data_test[var][:])
        if test:
            print(var, ' PASS ')
        elif not test:
            #Check the number of triangles
            if (len(data_ref[var][:]) != len(data_test[var][:])):
                print("Differences in variable %s length : ctrl = %s, test = %s" % (var, len(data_ref[var][:]), len(data_test[var][:])))
            else:
                print("variable %s has length : ctrl = %s, test = %s" % (var, len(data_ref[var][:]), len(data_ref[var][:])))
                dist = la.norm( (data_ref[var][:] - data_test[var][:] )/data_ref[var][:] )
                if var == 'idpair':
                    print(var, data_test[var][:]-data_ref[var][:])

                if abs(dist) < 1e-7 :
                    print(var, ' PASS BUT NOT EQUAL (after all, perfection is perfection, and this isnt)')
                elif abs(dist) < 1e-4 :
                    print(var, ' PASS BUT REALLY NOT EQUAL (but if you squeeze your eyes...)')
                else:
                    print(var, ' FAIL: alternative truths are not accepted.')
                    print('ERROR distance = ', dist)
        # print('data length difference', len(data_ref[var][:]) - len(data_test[var][:]))
        # print('REF dist- test dist: ', ( la.norm( data_ref[var][:]) - la.norm( data_test[var][:] ) ) / len(data_ref[var][:]  ) )
        # print('REF sum - test sum: ', ( nps( data_ref[var][:]) - nps( data_test[var][:]) ) / len(data_ref[var][:] ) )
    else:
        print(var, ' is not present in the tested dataset')

