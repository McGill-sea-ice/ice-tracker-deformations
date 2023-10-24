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


file_ref  = 'test_ref/RCMS1SID_20220101_20220102_dt24_tol24_dx.nc'
file_test = 'test_out/test/05_output/RCMS1SID_20220101_20220102_dt24_tol24_dx.nc'

data_ref = Dataset(file_ref)
data_test = Dataset(file_test)

for var in data_ref.variables:
    if var in data_test.variables:
        test = ( data_ref[var][:] == data_test[var][:] ).all()
        if test:
            print(var, ' PASS ')
        elif not test:
            dist = ( la.norm( data_ref[var][:]) - la.norm( data_test[var][:] ) ) / len( data_ref[var][:]  )
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

## 2 test of the coverage frequency
print(' ')
print(" --- 2 test of the coverage frequency ---")
print(' ')

file_ref  = 'test_ref/RCMS1_20220101_20220102_dt24_tol24_res20_int1_reflat70_coverage_area_timeseries.pkl'
file_test = 'test_out/test/figs/RCMS1_20220101_20220102_dt24_tol24_res20_int1_reflat70_coverage_area_timeseries.pkl'

data_ref = pd.read_pickle(file_ref)
data_test = pd.read_pickle(file_test)

for var in data_ref:
    if var in data_test:
        test = ( data_ref[var][:] == data_test[var][:] ).all()
        if test:
            print(var, ' PASS ')
        elif not test:
            dist = ( la.norm( data_ref[var][:]) - la.norm( data_test[var][:] ) ) / len( data_ref[var][:]  )
            if abs(dist) < 1e-5:
                print(var, ' PASS ')
            else:
                print(var, ' FAIL')
                print('ERROR distance = ', dist)
            print(var, ' FAIL')
        # print('data length difference', len(data_ref[var][:]) - len(data_test[var][:]))
        # print('REF dist- test dist: ', ( la.norm( data_ref[var][:]) - la.norm( data_test[var][:] ) ) / len(data_ref[var][:]  ) )
        #print('REF sum - test sum: ', ( nps(data_ref[var][:].astype('float'))  - nps( np.array(data_test[var][:].astype('float')))) / len(data_ref[var][:] ) )
    else:
        print(var, ' is not present in the tested dataset')
