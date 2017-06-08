############################################################################################################################
#
# Data load setup
#
############################################################################################################################

# img out sample calls
#
# gen blended sequence for 2 inputs:
#./manta ../mantaflowgit/scenes/pcgTestOfLoad.py doOnlyLoad 1 doSingleLoad 1 dataSet 201 showUi 1 blendIters 1 blendSingleAlpha 0.5 regen2dSlices 1 dispWidth 500 dispHeight 500 camposz 0.8 singleSel   0 0 0   0 1 0
#
# gen blended sequence for 3 inputs:
#./manta ../mantaflowgit/scenes/pcgTestOfLoad.py doOnlyLoad 1 doSingleLoad 1 dataSet 202 showUi 1 blendIters 1 regen2dSlices 1 dispWidth 500 dispHeight 500 camposz 0.8 singleSel3   0 0 0   0 1 0   1 0 0      blendSingleAlpha 0.5   blendThirdAlpha 0.5
#
# run single grid data
#./manta ../mantaflowgit/scenes/pcgTestOfLoad.py dataSet 202  dosingleselect 1  singleSel 0 1 0  0 2 0
#

import sys
from manta import *
from math  import *

setDebugLevel(1); # set to 2 for kernel msgs, 3 for ranges


dataSet       = 2
paramSet      = 5
showUi        = True
multiStep     = 3
doLsProject   = False
doFinalProject= True
doStoreResult = False
writePngs     = False

# load only modes
doOnlyLoad     = False
doSingleSelect = False
singleSel      = [0,0,0, 1,1,0,  -1,-1,-1]; 
doLoadThird    = False
doConstOverride = 0; # for testing only, paper ex2, turn of SDF!

# configure blending
blendIters       = -9;
blendSingleAlpha = -1;
blendThirdAlpha = -1;

# generate 2d slices from 3d version
regen2dSlices = False
# window settings, for img output
dispWidth  = 800
dispHeight = 600
camPosZ    = 1.2

# data input 
dim = 2
res = 10
resLoad        = 12
isLevelset     = True
frameIterStart = 0
frameIterEnd   = 12
frameIterMult  = 10
frameIterOffset= 10
orderTime  = 1
orderSpace = 1
projSizeThreshold = 959
doLoadGrid = 0
thirdDimFactor = 1.0

# modify filename for some runs...
dataSetSuffix = ""

# fit data into grid with 10% border, used in OF solve for reset
autoBorder  = 0.1
loadOffset  = Vec3(0,0,0);

gridLoadX  = -1;  gridLoadY = -1;  gridLoadZ = -1;
gridOtherX = -1; gridOtherY = -1; gridOtherZ = -1;
gridThirdX = -1; gridThirdY = -1; gridThirdZ = -1;


# check for parameters (ignore case), note: bool given as 0/1
# eg call with: "./manta scene.py   dataSet 2 showUi 1 doOnlyLoad 0 "
skipNext = 0;
for iter in range(1, len(sys.argv)):
	if skipNext>0:
		skipNext = skipNext-1; continue; 
	if(sys.argv[iter].lower() == "showui"):
		showUi = (int(sys.argv[iter+1])!=0);
		print( "showUi set to %d" % (showUi) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "dataset"):
		dataSet = int(sys.argv[iter+1])
		print( "dataSet set to %d" % (dataSet) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "paramset"):
		paramSet = int(sys.argv[iter+1])
		print( "paramSet set to %d" % (paramSet) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "dostoreresult"):
		doStoreResult = (int(sys.argv[iter+1])!=0);
		print( "doStoreResult set to %d" % (doStoreResult) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "doonlyload"):
		doOnlyLoad = (int(sys.argv[iter+1])!=0);
		print( "doOnlyLoad set to %d" % (doOnlyLoad) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "dosingleload") or (sys.argv[iter].lower() == "dosingleselect"):
		doSingleSelect = (int(sys.argv[iter+1])!=0);
		print( "doSingleSelect set to %d" % (doSingleSelect) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "singlesel"):
		for i in range(0, 6):
			singleSel[i] = int(sys.argv[iter+1+i]);
		print( "singleSel set to %d,%d,%d and %d,%d,%d" % (singleSel[0],singleSel[1],singleSel[2],singleSel[3],singleSel[4],singleSel[5]) ); 
		skipNext = 6;
	elif(sys.argv[iter].lower() == "singlesel3"): # double interpol, 9 values
		for i in range(0, 9):
			singleSel[i] = int(sys.argv[iter+1+i]);
		print( "singleSel3 set to " + str(singleSel) ); 
		skipNext = 9;
	elif(sys.argv[iter].lower() == "blenditers"):
		blendIters = int(sys.argv[iter+1])
		print( "blendIters set to %d" % (blendIters) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "blendsinglealpha"):
		blendSingleAlpha = float(sys.argv[iter+1])
		print( "blendSingleAlpha set to %f" % (blendSingleAlpha) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "blendthirdalpha"):
		blendThirdAlpha = float(sys.argv[iter+1])
		print( "blendThirdAlpha set to %f" % (blendThirdAlpha) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "regen2dslices"):
		regen2dSlices = (int(sys.argv[iter+1])!=0);
		print( "regen2dSlices set to %d" % (regen2dSlices) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "dispwidth"):
		dispWidth = int(sys.argv[iter+1])
		print( "dispWidth set to %d" % (dispWidth) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "dispheight"):
		dispHeight = int(sys.argv[iter+1])
		print( "dispHeight set to %d" % (dispHeight) ); skipNext = 1;
	elif(sys.argv[iter].lower() == "camposz"):
		camPosZ = float(sys.argv[iter+1])
		print( "camPosZ set to %d" % (camPosZ) ); skipNext = 1;
	else:
		print( "Error: unknown parameter '" + sys.argv[iter] +"' " );
		exit(1);



############################################################################################################################

# liquids

# 1-3 , 64^2 in pcgData64
if dataSet==1:
	datPath = "pcgData64/phi"
	frameIterEnd = 8
	res = 64
	dat0 = 3
	dat1 = 1  
	resLoad = 64
	#res = 100
	#frameIterOffset= 60; # more interesting start frame
	#dat0 = 1; dat1 = 2; # NT_DEBUG

# 4-6 , 140^2 in pcgData140
elif dataSet==2:
	datPath = "pcgData140/phi"
	res = 140
	resLoad = 140
	res = 168
	frameIterStart = 4 # debug , tough start frame

	# variants
	#dat0 = 4; dat1 = 5; # rel. easy
	dat0 = 4; dat1 = 6;  # tough!
	#dat0 = 6;	dat1 = 4;  # tougher , reverse
	#dat0 = 5; dat1 = 4; # NT_DEBUG new


# border (autoBord,autoBnd) reduce test data
elif dataSet==3:
	autoBorder  = 0.0
	datPath = "pcgDataBordred/phibv"
	res = 120
	dat0 = 1; 
	dat1 = 3; 
	resLoad = 100
	frameIterStart = 2 
	frameIterMult  = 1
	frameIterOffset= 1

# smoke

# 2d plume 1-4 (note 1&2 are too low...?)
elif dataSet==5:
	datPath = "pcgDataSmoke64/dens"
	#frameIterOffset = 5;
	frameIterOffset = 30;
	frameIterMult   = 5;
	frameIterEnd    = 8;
	res = 64;
	isLevelset = False;
	dat0 = 3;
	dat1 = 4;
	resLoad = 64;
	res = 80;
	frameIterEnd    = 18; frameIterOffset = 5 # longer , all frames
	#frameIterOffset = 60 # start later

# 2d paper example
elif dataSet==6:
	datPath = "pcgData2dPaper/phi"
	res = 140
	dat0 = 1; 
	dat1 = 2;
	resLoad = 128
	res = 280;
	frameIterStart = 2 

# 3d

# 7-9 , 64^3 in pcgData64 (3d versions of 1-3)
elif dataSet==10:
	datPath = "pcgData64threed/phi"
	dim = 3; res = 64
	dat0 = 7
	dat1 = 9
	resLoad = 64
	res = 16; # NT_DBEUG 4d
	# shorter
	frameIterStart = 1; frameIterEnd   = 7
	frameIterMult  = 20; frameIterOffset= 10
	#res = 124;	frameIterMult   = 5;	frameIterEnd    = 80; 


# 3d plume 5,6
elif dataSet==15:
	datPath = "pcgDataSmoke64threed/dens"
	frameIterOffset = 50
	frameIterMult   = 5
	dim = 3; res = 64
	isLevelset = False
	dat0 = 5
	dat1 = 6
	resLoad = 64
	# full run
	frameIterOffset = 5; frameIterEnd = 18

# 3d volumes for dat 1-3 , 64^2 in pcgData64conv3d , frame #0 is interval 1, #1 is interval 2
elif dataSet==101:
	datPath = "pcgData64conv3d/phito"
	dim = 3;
	res = 64;
	dat0 = 1;
	dat1 = 3;
	dat0 = 3; dat1 = 1; # reverse
	resLoad = 64;
	res = 76;
	frameIterOffset= 0;
	frameIterMult  = 1;
	frameIterEnd   = 2;
	#resLoad = 64; res = 40; # smaller
	frameIterOffset= 1; frameIterEnd   = 1;

elif dataSet==102:
	# special test, run a case from 202 single, for loop
	dat0 = 1; dat1 = 2; # fixed for grids!
	datPath = "pcgDataGrid01/phito02" 
	doLoadGrid = 1;
	dim        = 3;
	resLoad    = 80;
	res        = 96;
	frameIterOffset= 0; frameIterMult  = 1;
	frameIterStart = 0; frameIterEnd   = 1; 
	gridLoadX,  gridLoadY,  gridLoadZ  = (1,5,0)
	gridLoadX,  gridLoadY,  gridLoadZ  = (1,1,0) # tougher!
	gridOtherX, gridOtherY, gridOtherZ = (3,5,0)
	res        = 44; # small! 
	res        = 64; # ca. 10s

elif dataSet==103:
	# special test, run a case from 202 single, for loop
	dat0 = 1; dat1 = 2; # fixed for grids!
	datPath = "pcgDataFlip2d/phifl2d" 
	doLoadGrid = 0;
	dim        = 3;
	resLoad    = 100;
	res        = 120;
	frameIterOffset= 0; frameIterMult  = 1;
	frameIterStart = 0; frameIterEnd   = 1; 
	res        = 64; # small! 

elif dataSet==104:
	#NT_DEBUG print("ERR fix time below"); exit(1);
	# two d case with longer time
	dat0 = 1; dat1 = 2; 
	datPath = "pcgDataFlip2d/phifl2dLonger" 
	doLoadGrid = 0;
	dim        = 3;
	resLoad    = 50;
	res        = 60;
	frameIterOffset= 0; frameIterMult  = 1;
	frameIterStart = 0; frameIterEnd   = 1; 
	#res        = 34; # small! 
	thirdDimFactor = 1.5


# 3d volumes for 2d flip sims, grid 10x10;   to02: 80 , to03: 46  (pcgTestDatagenArrFlip)
elif dataSet==200:
	dat0 = 1; dat1 = 2; # fixed for grids!
	datPath = "pcgDataGrid01/phito03" # note, suffix added!
	doLoadGrid = 1;
	dim        = 3;
	resLoad    = 46;
	res        = 56;
	#res        = 32; # smaller NT DEBUG
	frameIterOffset= 0; frameIterMult  = 1;
	frameIterStart = 0; frameIterEnd   = 1;
elif dataSet==201 or dataSet==202 or dataSet==203:
	# 201 = 4 drops low
	# 202 = 6 drops higher
	# 203 = 6 drops higher , 2stepOf
	dat0 = 1; dat1 = 2; # fixed for grids!
	datPath = "pcgDataGrid01/phito02" 
	doLoadGrid = 1;
	dim        = 3;
	resLoad    = 80;
	res        = 96;
	frameIterOffset= 0; frameIterMult  = 1;
	frameIterStart = 0; frameIterEnd   = 1;

	#res        = 34; # small! 
	res        = 44; 
	#res        = 54; 
	#res        = 96; 

	# modify load
	if dataSet==202:
		dataSetSuffix = "B_"
	# spcecial = 2 step OF!
	if dataSet==203:
		#dataSetSuffix = "C_" # org, fastMarch = False
		dataSetSuffix = "D_"  # t2, fastMarch = True


else:
	print( "Error , no valid dataSet chosen!"); exit(1);

# automatically compute border offset
loadOffset  = Vec3(0.)
if autoBorder>0:
	resOffset   = int(res * autoBorder);
	loadOffset  = Vec3(resOffset,resOffset,resOffset);

# repeat the first frame?
doRepeat = False; 
repeatLen = int(res*autoBorder*1.0); 
if 1 and (dataSet>=200 or dataSet==102 or dataSet==103 or dataSet==104): #  NT_DEBUG , check repeat for long
	doRepeat = True;
	loadOffset.z = loadOffset.z + repeatLen


# ====== other params ... ====

# maximal region around surface to include from levelset
# modified below...
maxDist = 5

# solving accuracy, default
cgAccuracy = 1e-04

# cfl condition for advection (old 10)
cfl = 999; # NT_DBEUG off 

# dont correct outer layer in levelset 
corrOuter = False;
# use fast-marching , or simple-reinit?
#doFastMarch = True;
doFastMarch = False; # default for now

# misc params, to be initialized
lsFactor = 0
wSmooth = 0
wEnergy = 0
wPostBlur = 0
maxDist = 0
blurSelection = 0

# for projection step (used with doLsProject)
projMaxDist = 4.0; # 4 & 40 - good defaults...
projMaxIter = 40;

gsStd  = Vec3(res,res,res*thirdDimFactor)
if (dim==2):
	gsStd.z=1

solv           = Solver(name='main', gridSize = gsStd, dim=dim )
solv.timestep  = 1.0
times          = Timings()

# setup relics  ...
flags = solv.create(FlagGrid)
flags.initDomain()
flags.setConst(TypeFluid) # set all to TypeFluid

# testing , preview!
vDiff = solv.create(VecGrid)
i0blend  = solv.create(LevelsetGrid)
i0reinit = solv.create(LevelsetGrid)

# allocate grid data
i0      = solv.create(LevelsetGrid)
i1      = solv.create(LevelsetGrid)
i0org   = solv.create(LevelsetGrid)
i1org   = solv.create(LevelsetGrid)
vel     = solv.create(VecGrid)
velTmp  = solv.create(VecGrid);
# debug multi level solve
debvel5 = solv.create(VecGrid)
debvel4 = solv.create(VecGrid)
debvel3 = solv.create(VecGrid)
debvel2 = solv.create(VecGrid)
debvel1 = solv.create(VecGrid)
debls5 = solv.create(LevelsetGrid)
debls4 = solv.create(LevelsetGrid)
debls3 = solv.create(LevelsetGrid)
debls2 = solv.create(LevelsetGrid)
debls1 = solv.create(LevelsetGrid)
#rhs   = solv.create(RealGrid)
deberr1 = Vec3(0,0,0)
deberr2 = Vec3(0,0,0)
deberr3 = Vec3(0,0,0)
deberr4 = Vec3(0,0,0)
deberr5 = Vec3(0,0,0)
deberr6 = Vec3(0,0,0)
deberr7 = Vec3(0,0,0)
deberr8 = Vec3(0,0,0)
# optional for 2-param interpol
velThird   = 0
iDiffThird = 0

# only for velocity extrapol test , shouldnt be needed in general
velMac = ""

# test motion
i0adv   = solv.create(LevelsetGrid) # advected

# multi step OF
velTmp2  = solv.create(VecGrid)
iTmp1    = solv.create(LevelsetGrid)

iDiff     = solv.create(RealGrid)
iDiffBack = solv.create(RealGrid)

# other solver data
meshDest = solv.create(Mesh)
meshRef1 = solv.create(Mesh)
meshRef0 = solv.create(Mesh)

# solver object for load & size convert
gsload = vec3(resLoad,resLoad,resLoad*thirdDimFactor)
if (dim==2):
	gsload.z=1
sload  = Solver(name='loaddata', gridSize = gsload, dim=dim)
iload  = sload.create(RealGrid)

# 3d mesh env for 2d data
doMeshFrom2d = True
if doMeshFrom2d:
	t2gs = vec3(res,res,10)
	#if (dim==2): t2gs.z=1
	t2dsolv  = Solver(name='mesh3d', gridSize = t2gs, dim=3)
	t2dphi   = t2dsolv.create(LevelsetGrid)

# create preview mesh for 2d data
def previewMesh(phi, meshOut):
	for i in range(0,10):
		placeGrid2d(phi,t2dphi,dstz=i)
	t2dphi.setBound(0.1,0);
	t2dphi.createMesh(meshOut);

#///

if showUi and (GUI):
	gui = Gui()
	gui.windowSize(dispWidth,dispHeight);
	gui.setCamPos (0,0,-camPosZ);
	gui.show()
	if not regen2dSlices:
		gui.pause();
	else:
		gui.toggleHideGrids();

avgErrBest  = 1e10
screenNum   = 0
blendScrNum = 0


# helper functions
def extrapolLs(target, maxt=-1):
	# extrapolation distance, default & override
	d = (maxDist+1);
	if maxt>=0:
		d = maxt;

	if 1: # for testing, disable extrapol
		if doFastMarch:
			target.reinitMarching(flags=flags, maxTime=d , correctOuterLayer=corrOuter ); 
		else:
			extrapolateLsSimple(phi=target, distance=d)  
			extrapolateLsSimple(phi=target, distance=d, inside=true)

	if doConstOverride: # for testing, turn into 0-1 surface
		print("WARNING - const override test active!");
		constOverride(phi=target, val=-(1./lsFactor) );

	# extrapolLs done

# process newly loaded data for of
def doPostProcess(iTarget):
	if(isLevelset):
		extrapolLs(iTarget);
		remapLsValues( grid=iTarget, max=maxDist ) 
		iTarget.multConst(lsFactor)
	else:
		# convert to LS:
		iTarget.multConst(-1.0)
		iTarget.addConst(0.5) # org
		iTarget.addConst(0.3) # larger off

		# workaround , convert to LS...
		if 1:
			extrapolLs(iTarget);
			remapLsValues( grid=iTarget, max=maxDist )
			iTarget.multConst(lsFactor)
		#unused: simpleBlur(iTarget, preImgBlur)

	#iTarget.multConst(-1.0); 
	#iTarget.setBound( int(maxDist*-0.25), int(res*autoBorder*0.5) ); # old, hard reset
	iTarget.setBoundNeumann(int(res*autoBorder*0.5));
	# end post process
	

# load and process data file
def getDataName(dataId, num, attachSuffix=True ):
	suffix = "";

	# hard coded offset for paper 2d example
	# try1
	if 0:
		if dataSet==6 and dataId==1:
			num = 100 # img 7
		if dataSet==6 and dataId==2:
			num = 70 + 10
	# try2
	if 1:
		if dataSet==6 and dataId==1: 
			#num = 60 
			num = num + 20 
		if dataSet==6 and dataId==2:
			#num = 40 
			num = num + 0 
	# end - hard coded offset for paper 2d example

	if attachSuffix:
		suffix = ".uni";
	if(doLoadGrid<=0):
		return("./%s%02d_%04d%s" % (datPath, dataId, num, suffix ) );
	else:
		# grid data
		if  (dataId==1):
			return("./%s_x%dy%dz%d_%04d%s" % (datPath, gridLoadX,gridLoadY,gridLoadZ ,num, suffix ) );
		elif(dataId==2):
			return("./%s_x%dy%dz%d_%04d%s" % (datPath, gridOtherX,gridOtherY,gridOtherZ ,num, suffix ) );
		elif(dataId==3):
			return("./%s_x%dy%dz%d_%04d%s" % (datPath, gridThirdX,gridThirdY,gridThirdZ ,num, suffix ) );
		else:
			print( "Error - invalid grid dataId %d \n"%(dataId) ); 
			exit(1);
	return("-error-");

def loadData(iTarget, dataId, num ):
	global iload;
	if 1:
		# regular load
		iload.load( getDataName(dataId, num) );

		if dataSet==6: # 2d paper example
			iload.setBound(0.1, 3 ); # paper example - cut off flat side

	else:
		# fixed init of sphere for debugging
		drop = solv.create(Sphere, center=vec3(0.4,0.5,0.0)*vec3(res,res,res), radius=res*0.25); 
		if dataId==2:
			drop = solv.create(Sphere, center=vec3(0.4+0.05,0.5,0.0)*vec3(res,res,res), radius=res*0.25); 
		iload = drop.computeLevelset()
		# fixed init of sphere, end

	interpolateGrid(iTarget, iload, offset=loadOffset , scale=vec3(1.-2*autoBorder) );

	if doRepeat: 
		repeatFrame( iTarget, loadOffset.z + 1./float(res), repeatLen, bnd=0 ); 

	iTarget.setBound(0.1, int(res*autoBorder) ); # reset border of LS! 


def getResultName(dataId0, dataId1, num, suffix ):
	name = "";
	if(doLoadGrid<=0):
		name = "%s_res%02d_%s.uni" % (getDataName(dataId0, num, False), dataId1, suffix);
	else:
		name = "%s_to_x%dy%dz%d_%s.uni" % (getDataName(dataId0, num, False), gridOtherX,gridOtherY,gridOtherZ, suffix);
	return name;


def storeResult(iTarget, dataId0, dataId1, num , suffix):
	mySuffix = dataSetSuffix+suffix
	if doStoreResult:
		iTarget.save( getResultName(dataId0, dataId1, num, mySuffix) );

def loadResult(iTarget, dataId0, dataId1, num , suffix):
	mySuffix = dataSetSuffix+suffix
	iTarget.load( getResultName(dataId0, dataId1, num, mySuffix) );


# small helper to prepare final LS
def createFinal(iTarget, iSrc, copyFromSrc=True ):
	if(copyFromSrc):
		iTarget.copyFrom(iSrc);
	iTarget.setBound(-0.1, int(res*autoBorder) ); # reset border of LS!
	if(isLevelset):
		dummy = 1.;
		#extrapolLs(iTarget); # only slight gain in accuracy...

# simple helper to select 1 run from grid
def checkSingleSelect(testRunX,testRunY,testRunZ, otherRunX,otherRunY,otherRunZ):
	if(testRunX == singleSel[0] and
	   testRunY == singleSel[1] and
	   testRunZ == singleSel[2] and
	   otherRunX== singleSel[3] and
	   otherRunY== singleSel[4] and
	   otherRunZ== singleSel[5]): 
		return True;
	return False;


############################################################################################################################

accuAvgErr = 0.; # global results of runseq
def runSeq():
	global screenNum, blendScrNum, avgErrBest
	global accuAvgErr
	global gridOtherX, gridOtherY, gridOtherZ, gridThirdX, gridThirdY, gridThirdZ;

	accuAvgErr = 0.;
	accuAvgCnt = 0.;

	#main loop, check different tests
	for frameIter in range(frameIterStart,frameIterEnd): 

		# =========================================== setup

		i0.setConst(0.);
		i1.setConst(0.);
		print( "\nframeIter: %dd "%(frameIter) )

		# load & backup input data
		fileNum = frameIterOffset + frameIter*frameIterMult
		loadData( i0, dat0, fileNum );
		loadData( i1, dat1, fileNum );
		i0org.copyFrom(i0); 
		i1org.copyFrom(i1); 
		if 1 and isLevelset: # why extrapolate original?  -> data should match
			extrapolLs(i0org);
			extrapolLs(i1org);
		doPostProcess(i0);
		doPostProcess(i1);
		
		if not doOnlyLoad:
			# =========================================== step

			vel.setConst(vec3(0,0,0));

			opticalFlowMultiscale3d(vel=vel, i0=i0, i1=i1,   wSmooth=wSmooth, wEnergy=wEnergy  , postVelBlur=wPostBlur 
				, cgAccuracy=cgAccuracy , blurType=blurSelection , cfl=cfl, orderTime=orderTime,orderSpace=orderSpace  
				, resetBndWidth=autoBorder , multiStep=multiStep , projSizeThresh=projSizeThreshold, minGridSize=6, doFinalProject=doFinalProject 
				, dv1=debvel1,dv2=debvel2,dv3=debvel3,dv4=debvel4,dv5=debvel5 , dr1=debls1,dr2=debls2,dr3=debls3,dr4=debls4,dr5=debls5 );

			# use level set projection
			if doLsProject: 
				maxLSPNorm = 0.50; # zero defo where normals disagree (old: skip vary normals & extrapolate if <1 )
				maxLSPDist = 20;   # skip when defo further than this
				lspBlur    = wPostBlur * 0.5;
				#lspBlur    = 0.01; maxLSPDist = 5; # settings, test blur fill in

				# control by loop
				maxLSPNorm = -1. # doesnt matter
				maxLSPDist = projMaxDist;
				#lspBlur    = wPostBlur * 0.5; 
				lspBlur    = 2.
				maxIter    = projMaxIter;

				print( "\nAdding level-set projection, params %d %f , blur %d, iter %d " % (maxLSPDist, maxLSPNorm, lspBlur , maxIter) );

				corrVelsOf3d( velTmp , vel, i0, i0, i1 , maxLSPDist, maxLSPNorm , lspBlur, autoBorder , maxIter ); 
				#vel.addScaled(velTmp, Vec3(-1.0) );

			#if 0:
				#maxVel = vel.getMaxValue(); print( "Max vel "+str(maxVel) )# check... 


			# =========================================== advect forward 
			i0adv.copyFrom(i0);  # calcErrOrg
			advectSemiLagrangeCfl(cfl=cfl, flags=flags, vel=vel, grid=i0adv, order=orderTime)
			i0adv.setBoundNeumann(0); 

			# calcErrOrg 
			createFinal( i0reinit, i0adv );

			# =========================================== calc error
			# load reference surface again
	
			# just for comparison
			if (dim==3) and showUi:
				i1org.createMesh(meshRef1)
				i0org.createMesh(meshRef0)

			avgErr1 = 0.

			# create full SDF for advected i0
			if(isLevelset): 
				# mark disagreements in levelsets (also compute error
				avgErr1 = calcLsDiff3d(i0=i0adv, i1=i1, correction=20. ) * 1.; # calcErrOrg
				#? i0reinit.multConst(-1.); createFinal(i0reinit, i0reinit, False); # calcErrOrg

				if 1 and (dim==3) and showUi:
					iTmp1.copyFrom(i0reinit); 
					iTmp1.multConst( -( maxDist/(lsFactor) ));
					iTmp1.createMesh(meshDest)

			else:
				# non levelsets...
				avgErr1 = calcDensDiff(i0=i0reinit, i1=iTmp1); 

			print( "Frame %d %d->%d: avg curr error %f " % (fileNum, dat0,dat1, avgErr1 ) ) # per frame error

			# =========================================== 


			# recons 3d mesh surface for 2d
			if doMeshFrom2d and dim==2:
				iTmp1.copyFrom(i0reinit); 
				iTmp1.multConst( -( maxDist/(lsFactor) ));
				#simpleBlur(iTmp1, 1.);
				previewMesh(iTmp1, meshDest);
				previewMesh(i1org, meshRef1);
				previewMesh(i0org, meshRef0);


			# store difference grid
			iDiff.copyFrom(i1org); 
			iDiff.sub(i0reinit);
			#iDiff.setBound(0., int(res*autoBorder) ); # reset border 

			# store result
			storeResult( iDiff, dat0,dat1, fileNum, "dif" );
			storeResult( vel,   dat0,dat1, fileNum, "vel" );

		else:
			# load only
			loadResult( iDiff, dat0,dat1, fileNum, "dif" );
			loadResult( vel  , dat0,dat1, fileNum, "vel" );
			avgErr1 = 0.;


		# load data for third 
		if doLoadThird:
			# a bit ugly - bend "other" selector towards third
			(tmpx, tmpy, tmpz) = (gridOtherX, gridOtherY, gridOtherZ)
			(gridOtherX, gridOtherY, gridOtherZ) = getGridDataId(gridThirdX, gridThirdY, gridThirdZ);

			loadResult( iDiffThird, dat0,dat1, fileNum, "dif" );
			loadResult( velThird  , dat0,dat1, fileNum, "vel" );

			(gridOtherX, gridOtherY, gridOtherZ) = (tmpx, tmpy, tmpz)

		accuAvgErr += avgErr1 
		accuAvgCnt += 1.0

		solv.step()
		if writePngs and showUi: # output
			gui.screenshot( 'out02ms3_%05d.png' % screenNum ); 
			screenNum += 1

		# =========================================== 
		# try blending...
		maxBl = blendIters
		for blendIt in range(0,maxBl): 
			if blendSingleAlpha>=0.:
				alpha = blendSingleAlpha;
			else:
				alpha = (1.0*blendIt) / (1.0*(maxBl-1.0));
			beta = blendThirdAlpha; # optional

			i0blend.copyFrom(i0org); 

			# blending: 0 off, 1 SL, 2 wDiff , 3 2param
			blendMode = 1
			if not doLoadThird:
				print( "Blending %f "%(alpha) )
			else:
				print( "Blending sec %f, third %f "%(alpha, beta) )
				blendMode = 3;

			# show org i1 last
			if alpha>(1.0+1e-04): 
				blendMode = 99;

			# ... blend!
			if blendMode==0:
				# naive LS blend
				i0blend.multConst( 1.0 - alpha )
				iTmp1.copyFrom(i1org);
				iTmp1.multConst( alpha )
				i0blend.add( iTmp1 );
			elif blendMode==1:
				# just forward
				advectSemiLagrangeCfl(cfl=cfl, flags=flags, vel=vel, grid=i0blend, order=orderTime, velFactor=alpha)
			elif blendMode==2:
				# forward & backw
				advectSemiLagrangeCfl(cfl=cfl, flags=flags, vel=vel, grid=i0blend, order=orderTime, velFactor=alpha )
				if 1 and (isLevelset):
					createFinal( i0blend, 0 , copyFrom_i0adv=False ); # damn, needed?
				diffDuration = 0.3
				#wdiff = 1.0 * alpha # org identical
				wdiff = (1.0/diffDuration) * alpha - (1.0/diffDuration-1.) # diffDuration% shift
				# only add if >0
				if wdiff>0.:
					iTmp1.copyFrom( iDiff );   
					iTmp1.multConst( wdiff )
					advectSemiLagrangeCfl(cfl=cfl, flags=flags, vel=vel, grid=iTmp1, order=orderTime, velFactor=-(1.-alpha) )
					i0blend.add( iTmp1 );
					print( "Adding diff vector, weight = %f "%(wdiff) )

			elif blendMode==3:
				# 2-param, just forward
				advectSemiLagrangeCfl(cfl=cfl, flags=flags, vel=vel, grid=i0blend, order=orderTime, velFactor=alpha)
				advectSemiLagrangeCfl(cfl=cfl, flags=flags, vel=velThird, grid=i0blend, order=orderTime, velFactor=beta)

			elif blendMode==99:
				print( "No blend! Showing org1");
				i0blend.copyFrom(i1org); 

			else:
				print( "Invalid blend mode! %d"%blendMode ); exit(1);



			# reinit...
			if(isLevelset):
				extrapolLs(i0blend, res/4);
				i0blend.setBound(0.1, int(res*autoBorder) ); # reset border 
				if (dim==3) and showUi:
					i0blend.createMesh(meshDest)
				if (dim==2) and showUi:
					previewMesh(i0blend, meshDest);
				# optionally save mesh...?
				if 0: # output
					meshDest.save( 'mesh_%d_%04d.obj' % (blendIt,blendScrNum) )
			else:
				i0blend.copyFrom(i0blend, False)
				if 0: # output
					i0blend.save( 'outblxy_%d_%04d.uni' % (blendIt,blendScrNum) )

			# output


			# re-generate 2d slices for whole 2d anim, output - similar to preview Mesh , only for solves with time dim
			if regen2dSlices:
				# note - skip region that is only safety border (ie, res* autoBorder*(1+1) )
				# once for blank region, and once in the beginning for the first-frame-repeat
				for z in range( int(res* autoBorder*2.0)-1, int(res*thirdDimFactor-(res* autoBorder)) ):
					getSliceFrom3d(src=i0blend, srcz=z, dst=t2dphi, dstz=0);
					for dstz in range(1,10):
						getSliceFrom3d(src=t2dphi, srcz=0, dst=t2dphi, dstz=dstz);
					t2dphi.setBound(0.1,0);
					t2dphi.createMesh(meshDest);
					if writePngs and showUi: # output
						#gui.screenshot( 'out_%d_%04d.png' % (blendIt,z) ); 
						if not doLoadThird:
							gui.screenshot( 'outbl%02d_%f_%04d.png' % (dataSet,alpha,z) ); 
						else:
							gui.screenshot( 'outbl%02d_%f_%f_%04d.png' % (dataSet,alpha,beta,z) ); 
					solv.step();
			else:
				# normal display / screenshot
				solv.step()
				if 0 and writePngs and showUi: # output
					gui.screenshot( 'outbl01_%d_%04d.png' % (blendIt,blendScrNum) ); 


		if maxBl>0:
			print("Blending done");
			solv.step(); # stop after blending is done?
		# done blending 

		blendScrNum = blendScrNum + 1


	accuAvgErr /= accuAvgCnt
	print( "Accu avg error %f  (params lsf %f, maxd %f, smo %f, ene %f, blr %f, cgacc %f, rI %d, rO %d)   d/pSet:%d,%d|%d,%d " % (accuAvgErr, lsFactor, maxDist, wSmooth, wEnergy, wPostBlur, cgAccuracy, res, resLoad, dataSet, paramSet, dat0,dat1 ) )
	if(accuAvgErr < avgErrBest):
		avgErrBest = accuAvgErr
	times.display()
	exit(1); # quit after one cycle 
		

############################################################################################################################

def loopParamSearch() :
	global lsFactor, wSmooth, wEnergy, wPostBlur;
	global maxDist, cgAccuracy, blurSelection;
	global projMaxDist , projMaxIter;

	# prepare loop ranges
	loopLenMd  = 1
	loopLenLfs = 1
	loopLenBlr = 1

	if paramSet==6: # search
		loopLenMd  = 3
		loopLenLfs = 5
		loopLenBlr = 5

	if paramSet==7: # search
		loopLenMd  = 5
		loopLenLfs = 5
		loopLenBlr = 5


	# main loop, vary parameters
	for testRunMd in range(0,loopLenMd):
		for testRunLsf in range(0,loopLenLfs): # 5
			for testRunBlr in range(0,loopLenBlr): # 10

				# run1, main params (range 4,6,6)
				if paramSet==1: # search 64
					lsFactor  = pow(10, (1.-testRunLsf) )
					wSmooth   = 0.00002
					wEnergy   = 0.001
					wPostBlur = 1.0 * testRunBlr
					maxDist   = 0.5 * pow(5, (testRunMd) ) 

				# good run 102 accu avg error 3.562538  (params lsf 0.100000, maxd 12.500000, smo 0.000020, ene 0.001000, blr 3.000000)
				elif paramSet==2: # search 64 p2
					lsFactor  = 0.1
					wSmooth   = 0.00002
					wEnergy   = 0.001
					wPostBlur = 3.0 
					maxDist   = 12.0
					# run2
					wSmooth   = 0.000001 * pow(5, testRunBlr)
					wEnergy   = 0.0001   * pow(5, testRunLsf)

				elif paramSet==3: # search 140
					lsFactor  = 0.1
					wEnergy   = 0.000125
					wPostBlur = 3.0 * testRunBlr
					maxDist   = 0.5 * pow(5, (testRunMd) ) 
					wSmooth   = 0.000001 * pow(5, testRunLsf)

				# 1/25 of smoothing amount, 3x more blur!
				elif paramSet==4:
					lsFactor  = 0.1
					wSmooth   = 0.000001
					wEnergy   = 0.000125
					wPostBlur = 9.0 
					maxDist   = 12.0

				# good from 64 , smoke
				elif paramSet==5:
					wPostBlur  = 3.0 
					wSmooth    = 0.01
					wEnergy    = 0.0005
					# 1 for median, 2 for simple
					#wPostBlur  = 2.0  

					wPostBlur  = 1.0 * testRunBlr
					blurSelection  = testRunLsf

					#wPostBlur  = 1.0 
					blurSelection  = 1

					maxDist    = 1. + 10. *testRunLsf
					lsFactor   = 0.1 * (testRunBlr+1)
					maxDist    = 20; 
					lsFactor   = 0.1
					wPostBlur  = 6. # compromise?

					#cgAccuracy = 0.003
					cgAccuracy = 0.03 

				# exhaustive search
				elif paramSet==6:
					maxDist    = 20; 
					lsFactor   = 0.1;
					cgAccuracy = 0.01;
					blurSelection = 1;
					wPostBlur = 2.0       *         testRunBlr;
					wSmooth   = 0.00001   * pow(10, testRunLsf);
					wEnergy   = 0.0001    * pow(10, testRunMd);

				# search project
				elif paramSet==7:
					wSmooth    = 0.01
					wEnergy    = 0.0005 
					blurSelection  = 1 
					maxDist    = 20; 
					lsFactor   = 0.1
					wPostBlur  = 6. 
					cgAccuracy = 0.001 
					wPostBlur   = 10. / (testRunBlr+1);
					projMaxDist = 20  / (testRunLsf+1);
					projMaxIter = 200 / (testRunMd+1);
					# search 1: maxdist 20 bad!
					wPostBlur   = 2.* (testRunBlr+0.);
					projMaxDist = 5-(testRunLsf);
					projMaxIter = 5 * (testRunMd+1);
					# good:
					wPostBlur   = 4;  # times 0.5 later -> 2
					projMaxIter = 40;
					projMaxDist = 4;
					# test fill in , NT_DEBUG , doFillIn
					# minimal difference?
					projMaxIter = 3 * (testRunBlr+0);

				# liquid 64 3d
				elif paramset==10:
					wSmooth    = 0.000025
					wEnergy    = 0.000125
					wPostBlur  = 3.0 
					maxDist    = 12.0
					lsFactor   = 0.1
					#cgAccuracy = pow(10, -testRunBlr)
					#cgAccuracy = 0.05
					cgAccuracy = 0.5   # low, but still ok...
					#cgAccuracy = 10.0 # solve effectively off

				else:
					print( "Error , no valid paramset chosen!"); exit(1);

				print( "\nparamSearch loop %d %d %d "%(testRunBlr,testRunMd,testRunLsf) )
				runSeq();

			print( "Accu avg error === " ) # loop sep

	print( "Accu avg error best so far %f " % avgErrBest )


############################################################################################################################

def idString(testRunX,testRunY,testRunZ, otherRunX,otherRunY,otherRunZ):
	ret = "";
	# uni directional graph
	#if   otherRunZ < testRunZ:
		#ret = ("%d_%d_%d vs %d_%d_%d" % (otherRunX,otherRunY,otherRunZ, testRunX,testRunY,testRunZ   ));
	#elif otherRunY < testRunY:
		#ret = ("%d_%d_%d vs %d_%d_%d" % (otherRunX,otherRunY,otherRunZ, testRunX,testRunY,testRunZ   ));
	#elif otherRunX < testRunX:
		#ret = ("%d_%d_%d vs %d_%d_%d" % (otherRunX,otherRunY,otherRunZ, testRunX,testRunY,testRunZ   ));
	#else:
		#ret = ("%d_%d_%d vs %d_%d_%d" % (testRunX,testRunY,testRunZ,    otherRunX,otherRunY,otherRunZ));

	# bi-directional
	ret = ("%d_%d_%d vs %d_%d_%d" % (testRunX,testRunY,testRunZ,    otherRunX,otherRunY,otherRunZ));

	return ret;

def getGridDataId(ix,iy,iz):
	outx,outy,outz = (-1,-1,-1)
	if dataSet==200 or dataSet==201: 
		outx = ix * 2;
		outy = iy * 2;
		outz = iz * 2;
	elif dataSet==202:
		outx = 1 + ix * 2;
		outy = 1 + iy * 2;
		outz = 0; # 2 params
	elif dataSet==203:
		outx = 1 + ix * 2;
		outy = 1 + iy * 2;
		outz = 0; # 2 params
	return (outx,outy,outz)

def loopGrid() :
	global lsFactor, wSmooth, wEnergy, wPostBlur;
	global maxDist, cgAccuracy, blurSelection;
	global gridLoadX, gridLoadY, gridLoadZ;
	global gridOtherX, gridOtherY, gridOtherZ;
	global gridThirdX, gridThirdY, gridThirdZ;
	global doLoadThird, velThird, iDiffThird;

	# sanity check for grid single load
	#if doSingleSelect and not doOnlyLoad:
		#print ("Error, use single load only with load mode"); exit(1);

	errAllRuns = 0.
	# prepare loop ranges
	loopLenZ = 1;
	loopLenY = 1;
	loopLenX = 1;

	if dataSet==200 or dataSet==201: 
		loopLenZ = 1
		loopLenY = 2
		loopLenX = 2
	if dataSet==202:
		loopLenZ = 1
		loopLenY = 3
		loopLenX = 2
	if dataSet==203:
		loopLenZ = 1
		loopLenY = 3
		loopLenX = 2

	gridThirdX, gridThirdY, gridThirdZ = (singleSel[6],singleSel[7],singleSel[8]) 
	if doSingleSelect and gridThirdX>=0 and gridThirdY>=0 and gridThirdZ>=0:
		print("Loading third data point for interpolation")
		doLoadThird = True;
		# allocate grids
		velThird   = solv.create(VecGrid)
		iDiffThird = solv.create(LevelsetGrid)

	# fix params
	maxDist   = 20
	lsFactor  = 0.1
	wSmooth   = 0.01
	wEnergy   = 0.0005
	#wPostBlur = 2
	wPostBlur = 6
	#cgAccuracy = 0.003 # doesnt pay off
	cgAccuracy = 0.03
	blurSelection = 1
	#cgAccuracy = 11.11 # NT DEBUG off

	gv = {};

	# main loop, vary parameters
	for testRunZ in range(0,loopLenZ):
		for testRunY in range(0,loopLenY): 
			for testRunX in range(0,loopLenX): 
				# calculate offsets into data
				gridLoadX, gridLoadY, gridLoadZ = getGridDataId(testRunX,testRunY,testRunZ );

				for otherRunZ in range(0,loopLenZ):
					for otherRunY in range(0,loopLenY): 
						for otherRunX in range(0,loopLenX): 
							if 1 and testRunX==otherRunX and testRunY==otherRunY and testRunZ==otherRunZ:
								continue;
							#gridOtherX = otherRunX * 2; gridOtherY = otherRunY * 2; gridOtherZ = otherRunZ * 2;
							gridOtherX, gridOtherY, gridOtherZ = getGridDataId(otherRunX,otherRunY,otherRunZ );

							#print( "\ngridSearch (%d %d %d) -> (%d %d %d)  "%(testRunX,testRunY,testRunZ, otherRunX,otherRunY,otherRunZ) )
							#print( "    sel    (%d %d %d) -> (%d %d %d)  "%(singleSel[0],singleSel[1],singleSel[2],singleSel[3],singleSel[4],singleSel[5]) );
							if (not doSingleSelect) or (doSingleSelect and checkSingleSelect(testRunX,testRunY,testRunZ, otherRunX,otherRunY,otherRunZ) ):
								#print( "\ngridSearch (%d %d %d) -> (%d %d %d)  "%(testRunX,testRunY,testRunZ, otherRunX,otherRunY,otherRunZ) )
								print( "\ngridSearch (%d %d %d) -> (%d %d %d)  "%(gridLoadX,gridLoadY,gridLoadZ, gridOtherX,gridOtherY,gridOtherZ) )

								runSeq();

								# store error
								print( idString(testRunX,testRunY,testRunZ, otherRunX,otherRunY,otherRunZ) );
								gv[ idString(testRunX,testRunY,testRunZ, otherRunX,otherRunY,otherRunZ) ] = accuAvgErr;
								errAllRuns = errAllRuns + accuAvgErr;
				# sub loop done
	# main loop done

	if not doOnlyLoad:
		if 0:
			print( "\nResults:");
			for item in gv:
				print("%s: %f" % (item, gv[item]));

		if 1:
			sorted_list = [x for x in gv.iteritems()] 
			sorted_list.sort(key=lambda x: x[0]) # sort by id string
			#sorted_list.sort(key=lambda x: x[1]) # sort by value 
			#sorted_list.reverse() # to reverse the sort
			print( "\nResults sorted:");
			for item in sorted_list:
				print("%s: %f" % (item[0], item[1]));

	print( "\nTotal added err: %f"%(errAllRuns));
	# grid done



############################################################################################################################


# select loop type
if dataSet<200:
	loopParamSearch();
else:
	loopGrid();


