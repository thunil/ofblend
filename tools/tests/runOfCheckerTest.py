#!/usr/bin/python
#
# Run optical flow checkerboard tests , similar to unit tests in tools/tests. Use as 
# "python runOfCheckerTest.py | grep Result " to check outcome
# (note - needs ../testdata dir! a bit ugly..)
# optional: set evn var MANTA_GEN_TEST_DATA=1 for generating test data
#
import os, sys, shutil

doRealRun = 1;
manta  = "/Users/sinithue/devel/manta/buildMaster/manta"
#infile = "../../scenes/ofCheckerboardTest.py"
infile = "./ofCheckerboardTest.py"

call1 = ('"%s"  "%s"  1  ' % (manta, infile) ); 
print("CALL  %s\n"%(call1));
if 1 and doRealRun:
	os.system(call1);
else: 
	exit(1);

call1 = ('"%s"  "%s"  2  ' % (manta, infile) ); 
print("CALL  %s\n"%(call1));
os.system(call1);

call1 = ('"%s"  "%s"  3  ' % (manta, infile) ); 
print("CALL  %s\n"%(call1));
os.system(call1);

call1 = ('"%s"  "%s"  4  ' % (manta, infile) ); 
print("CALL  %s\n"%(call1));
os.system(call1);


call1 = ('"%s"  "%s"  5  ' % (manta, infile) ); 
print("CALL  %s\n"%(call1));
os.system(call1);

call1 = ('"%s"  "%s"  6  ' % (manta, infile) ); 
print("CALL  %s\n"%(call1));
os.system(call1);

	
