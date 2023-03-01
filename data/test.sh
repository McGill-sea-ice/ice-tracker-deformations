# test.sh

TEST_DIR="test_out"
SID_DIR="../src/SeaIceDeformation"
SAC_DIR="../SatelliteCoverage"
DATA_DIR="../../data/"

# 0 Remove test_out if it exists

if [ -d "$TEST_DIR" ]; then
	rm -r $TEST_DIR
fi

# 1 run the SID default script

cd $SID_DIR

if [ -e "namelist.ini" ]; then
	mv "namelist.ini" "namelist.ini_bak"
fi

python main.py

if [ -e "namelist.ini_bak" ]; then
	mv "namelist.ini_bak" "namelist.ini"
fi

# 2 run the SAC coverage_frequency_map.py

cd $SAC_DIR

if [ -e "options.ini" ]; then
	mv "options.ini" "options.ini_bak"
fi

python coverage_frequency_map.py

# 3 run netcdf_tools.py

python netcdf_tools.py

if [ -e "options.ini_bak" ]; then
	mv "options.ini_bak" "options.ini"
fi

# 4 run the report_test.py

cd $DATA_DIR

python report_test.py

# end of script
