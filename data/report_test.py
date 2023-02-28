# Test the netcdf
from netCDF4 import Dataset

file_ref  ='test_ref/RCMS1SID_20220101_20220102_dt24_tol24_dx.nc'
file_test ='test_out/test/05_output/RCMS1SID_20220101_20220102_dt24_tol24_dx.nc'

data_ref = Dataset(file_ref)
data_test = Dataset(file_test)

print(data_ref)

for var in data_ref.variables:
    test = ( data_ref[var][:] == data_test[var][:] ).all()
    if test:
        print(var, ' PASS ')
    elif not test:
        print(var, ' FAILED')