#
# Checkerboard test init 
# set activateChecker to modes 1-6 below , usually run with runOfCheckerTest.py script
# (note - very high values around 40k to quickly show errors)
#
import sys;
import os
import shutil
from   manta import *
from helperInclude import *

# hard coded for now...
sys.path.append("../../scenes/") 
from   ofHelpers import *

activateChecker = 3; # default loading
showGui         = 1;

# run as test?
testMode = False;

# override prefix from cmd line
if(len(sys.argv)>=2):
	testMode = True;
	showGui  = 0;
	activateChecker = int(sys.argv[1]);
print("Mode: %d " % activateChecker);

# main modes
res = 64; # default
#res = 44; # smaller test!
if activateChecker == 1: res = 80; # gen volumes
if activateChecker == 2: res = 44; # load & gen defo vel
if activateChecker == 3: res = 94; # load both
if activateChecker == 4: res = 64; # slice mode , like mode 5, but alpha 0
if activateChecker == 5: res = 64; # slice mode , with advect for alpha 0.5
if activateChecker == 6: res = 64; # slice mode opt, otherwise same as mode 5, note: compares with data from 5!


# note - first gen of input volumes
#        are loaded with borders at sides
# then - calc defo from this
# then - load both for output 


# win/mac paths
dataPath="."
dataPath="../testdata/"
if int(os.getenv('OFBLEND_WIN', 0))==1:
	dataPath="I:/dataOfBlend/"
elif int(os.getenv('OFBLEND_WIN', 0))==2:
	dataPath="/Users/sinithue/devel/manta/buildtbb/"
elif int(os.getenv('OFBLEND_WIN', 0))==3:
	# default imac
	dataPath="/Users/sinithue/devel/manta/buildImac/"

screenPrefix = "checkbt03"; 

doRunOf     = 0;
doLoad      = 1;
doSave      = 0;
doLoadAndProject = 0; # note, only makes sense when doLsProject=true, and doRunOf=False
doCacheInputs = 0;
loadSmoke = 0; # extrapolate and gen levelset, or load/recons org smoke data

# output settings, note sync with checkerTest!
# pass together as:  writePngs  , triMesh    , projPpm    , writeUni   , sdfSmooth  , meshPrefix 
writePngs  =0 
triMesh    =0 
projPpm    =0 
writeUni   =0 
sdfSmooth  =0
meshPrefix = "fluidsurface_final";


#doMeshOut = False; 
#meshPrefix = "fluidsurface_final";

doSingleFrameOut = 0; # special output mode, for single frame morph

#doDiffAdd = 0; # different modes!
yFactor = 1.0

scrPre2 = "cbt";
yFactor = 1.2; # a bit more difficult

debugSkipLoad = 999; # off

fourdFactor = 1.5;

# separate time scaling (actual speed, instead of grid size!?)
loadTimeScale = 1.0;

# interval for outpu
outInt = 1;
#outInt = 0.50;  print("Note - increased outInt to %f \n" % outInt);

fourdFactor     = 1.6;

blendAlpha      = 1.0;
blendAlpha      = 0.0;  print("\nWARN! advect effectively off (i0)\n"); # off
#blendAlpha      = 1.0;  # notable difference to i0
# simply amplifies diff in phi blendAlpha      = 1.2;  print("\nWARN! advect with %f \n" % blendAlpha ); # off

if 1 and activateChecker>2:
	loadTimeScale = 1.;  # phiout SL, 22 , 40, 59, 79, 98 gone
	#loadTimeScale = 0.5;  # phi and vel out of sync , check , used?
	print("\nWARN! outInterval to %f , time scale load %f \n"%(outInt,loadTimeScale) );

print("Blend alpha advect %f"%blendAlpha); # off


# note - defo "could" have different time dim (normally not)
fourdFactorDefo = fourdFactor;
# for data load later on
fourdFactorLoad = fourdFactor;

# data file names
fileLoadPrefix = "xxx%03d"; 
fileOutPrefix  = "xxx%03d";  
cacheprefix = "prep";

inpath         = "%spcgDataCheckertest/"%dataPath
defoRes        = 44;
fileLoadPrefix = "%sdefoCheckerC_%s" % ( inpath,"%03d"); cacheprefix="prep";  
fileOutPrefix  = fileLoadPrefix;  

# for special sliced load of large data sets
fileSliceLoadPrefix = "%spcgDataCheckertest/checker100_slice%s.uni" % (dataPath, "%04d" ) # "final" xl data
fileSliceIdxStart   =   0
fileSliceIdxEnd     = 160



# main params
multiStep         = 3;
doLsProject       = 0; # proj on/off
doFinalProject    = True;
projSizeThreshold = 74;
# projection params
#maxLSPDist = 20; 
maxLSPDist = 10; 
lspBlur    = 2.;

doSlicedLoad = 0; # 0=off, 1=load single, 2=load uni sequence

# ===================
# ===================
# ===================

if activateChecker==1: # data gen
	doLoad = 2; doRunOf = 0; doSave = 0;
if activateChecker==2: # load , gen defo
	doLoad = 1; doRunOf = 0; doSave = 0; loadSmoke = 2;
if activateChecker==3: # load both
	doLoad = 1; doRunOf = 0; doSave = 0; loadSmoke = 2;
if activateChecker==4: # load both
	doSlicedLoad = 2; doLoad = 0; doRunOf = 0; doSave = 0; loadSmoke = 2;
	#debugSkipLoad = 20; print("\n\nWARN debug skip!\n\n");	# for quick tests
if activateChecker==5: # load both
	# uses loadSlice_dumpOut and loadAdvectTimeSlice 
	doSlicedLoad = 2; doLoad = 0; doRunOf = 0; doSave = 0; loadSmoke = 2; blendAlpha = 0.5;
if activateChecker==6: # load both
	# uses loadAdvectTimeSlice_OptRun and loadSlice_dumpOutOpt
	doSlicedLoad = 3; doLoad = 0; doRunOf = 0; doSave = 0; loadSmoke = 2; blendAlpha = 0.5;



# ===================

# output data
fileId = res;
if (not doRunOf):
	fileId = defoRes;
fileLoadPrefix= fileLoadPrefix%(fileId);
fileOutPrefix = fileOutPrefix%(fileId);



# ========================
# registration params

autoBorder    = 0.1;
autoBrdRes    = int(res * autoBorder);
postVelBlur   = 4;
repeatStartDist = float(res) * autoBorder * loadTimeScale; # keep float
print("Repeat start "+str(repeatStartDist) );

# resolution of data
fourdFactorLoad = fourdFactor; # assume all are the same!
resLoad         = 80;
dataname_i0     = "%spcgDataCheckertest/checker80b.uni" % dataPath;
dataname_i1     = "%spcgDataCheckertest/checker80b.uni" % dataPath;

# dont thicken for smoke!
#thickenSurface = -0.5 * (resLoad/(res*1.0)); # shift by 1 cell?

lsFactor    = -1.0/(20.0) * 0.1;
maxDist     = 20;
cgAccuracy  = 0.01;
#cgAccuracy  = 10.01; print("\n\nWARN! cg effectively off...\n\n"); # off

maxDist     = 40; print("Note - increased maxDist to %f " % maxDist);
orgErr      = 0.;

# solver setup
dim = 3;
(gs,gs4d) = initRes( res, yFactor, fourdFactorDefo, loadTimeScale);

print("Main resolution "+str(gs4d) );
solv           = Solver(name='main', gridSize = gs, dim=dim, fourthDim=int(res*fourdFactor) );
solv.timestep  = 1.0;

times  = Timings();

# setup relics  ...
flags = solv.create(FlagGrid)
flags.initDomain()
flags.fillGrid()

# allocate grid data
i0   = solv.create(Grid4Real)
if 1 or not doSlicedLoad:
	i1   = solv.create(Grid4Real)
	i0adv= solv.create(Grid4Real)

# preview intermediate step
i0adv_org=""
doIadvOrg = False;
if(doIadvOrg):
	i0adv_org= solv.create(Grid4Real)
if 1 or not doSlicedLoad:
	vel  = solv.create(Grid4Vec4)
if (doRunOf or doLoadAndProject) and doLsProject:
	velTmp = solv.create(Grid4Vec4)

if doCacheInputs==1 and loadSmoke==1:
	print("Turn off caching for smoke!?!"); exit(1);

# rescale time dim for load
if (loadTimeScale != 1.) and (doSlicedLoad<2):
	itimescale = solv.create(Grid4Real);
	# dont mod loadRescale.t!

# loader
(gsload,gs4dload) = initRes( resLoad, yFactor, fourdFactorDefo, 1.); # note - no loadTimeScale here
sload  = Solver(name='loader', gridSize = gsload, dim=dim, fourthDim=gs4dload.t);
print("Data load dimension "+str(gs4dload)+" " );

# loading memory only required if preprocessede files are not available
loadMemRequired = True;
if doSlicedLoad==1: 
	prepname = "%s_%s%04d.uni" %(dataname_i0, cacheprefix, res);
	if doCacheInputs and (os.path.isfile( prepname )):
		loadMemRequired = False;
if loadMemRequired:
	iload  = sload.create(Grid4Real)

# effective time range should not be finer than org data
timerangeStart = (int(autoBrdRes*1.5)+3) ;
timerangeStart = (int(autoBrdRes*1.0)+3) ; # for liqsm
timerangeStop  = int( (res*fourdFactor-2.0*autoBrdRes-1)*loadTimeScale ); # shorten end
#print("Core time range %d to %d; res3d %d , durationT %d " % ( timerangeStart,timerangeStop, res, resLoad*fourdFactor ) );

timerangeStart = 0; timerangeStop  = int(res*fourdFactor*loadTimeScale)-1; print("DEBUG whole time range out \n"); # debug, show whole
# output for retry render of liqsm
#timerangeStart = 60; timerangeStop  = 400; print("Retry 2 render range \n"); # render retry 2, res320
if testMode:
	timerangeStart = 40; timerangeStop  = 45;

if 0 and activateChecker==3:
	timerangeStart = (int(autoBrdRes*1.5)+1) ; # autobrd + repeStart

# init UI
gui = 0;
if showGui and (GUI):
	gui = Gui()
	gui.show()
	gui.setCamPos (0,0,-1.2);

	# debug
	#gui.pause()


# output

#phiout = solv.create(LevelsetGrid)
phiout = solv.create(RealGrid)
if(doIadvOrg):
	phiout_org = solv.create(RealGrid)
if not doSlicedLoad:
	phiout1 = solv.create(RealGrid)
	phiout0 = solv.create(RealGrid)

velout  = solv.create(Vec3Grid)
veltout = solv.create(RealGrid)

# grids for debugging, see below
#phioutDeb = solv.create(RealGrid) # for debugging only
#phioutCmp = solv.create(RealGrid)

# checker debug
phioutSL = solv.create(RealGrid)
velTmp_out  = solv.create(Vec3Grid)
velTmp_tout = solv.create(RealGrid)

phiDiff = solv.create(RealGrid)
velDiff = solv.create(RealGrid)

meshout = solv.create(Mesh); # just as a helper

# not always needed, but allocate solver
(gsDefoLoad,gs4dDefoLoad) = initRes( defoRes, yFactor, fourdFactorDefo, loadTimeScale);
print("Defo (vel) load dimension "+str(gsDefoLoad)+" "+str(defoRes*fourdFactorDefo) );
sDefoload  = Solver(name='defo_loader', gridSize = gsDefoLoad, dim=dim, fourthDim=int(gs4dDefoLoad.t) ); 



# ======================================================================================================================


def my_range(start, end, step):
	while start <= end:
		yield start
		start += step

def initOthers(dst0, dst1, dstvel, dst4d, dstvelTmp, dst4dTmp, t):
	getSliceFrom4d(src=i0, dst=dst0, srct=t); 
	getSliceFrom4d(src=i1, dst=dst1, srct=t); 
	getSliceFrom4dVec(src=vel, dst=dstvel, dstt=dst4d, srct=t); 
	#getSliceFrom4dVec(src=velTmp, dst=dstvelTmp, dstt=dst4dTmp, srct=t); 

# load data sets & registrations

def postProcChecker(iTarget, isSlice=False ):
	# start frame , note - make sure it hits full data!
	repeatFrame4d( iTarget, loadOffset.t + 1./float(res), repeatStartDist, bnd=0 ); 

	if 1:
		print("Resetting full border %d"%autoBrdRes);
		# reset full border
		if loadSmoke==1:
			iTarget.setBound(0.125,autoBrdRes); # neutralize offset
		elif loadSmoke==2:
			print("Load checker handling");
		else:
			iTarget.setBound(0.,1*autoBrdRes); 

	# convert from "LS" to smoke
	if loadSmoke==1:
		iTarget.addConst(-0.125) 
		iTarget.multConst(-1.0 * lsFactor); # undo mult later on!
		#print("Min max %f   %f" % (iTarget.getMax(), iTarget.getMin()) ); 
	elif loadSmoke==2:
		if not isSlice:
			iTarget.addConst(1.0);
			iTarget.multConst(-0.5);
			print("Load checker handling");
	else: 
		#iTarget.addConst( thickenSurface );
		iTarget.addConst(-0.1); print("NOTE - special smoke data iso surface offset active!");

		# extrapol 
		iTarget.setBound(0.1, autoBrdRes ); # reset border of LS! 

		# extrapolate
		extrap4dLsSimple(phi=iTarget, distance=maxDist)  
		extrap4dLsSimple(phi=iTarget, distance=maxDist, inside=true)

		iTarget.multConst( lsFactor );


def loadPhiData(iTarget, fname ):
	# load & resize
	global iload;

	prepname = "%s_%s%04d.uni" %(fname,cacheprefix,res);
	if doCacheInputs and (os.path.isfile( prepname )):
		if (loadTimeScale != 1.): 
			itimescale.load(prepname);
			#sizeTimescale = Vec4(res, res*yFactor, res, (res*fourdFactor*loadTimeScale) );
			sizeTimescale = Vec4(1.,1.,1., (loadTimeScale) );
			print("Load time scale %f " % ( sizeTimescale.t ) );
			interpolateGrid4d(iTarget, itimescale, offset=Vec4(0) , scale=sizeTimescale );
		else:
			# easy, direct load...
			iTarget.load(prepname);
		return;

	iload.load( fname ); # getDataName(dataId, num) ); 
	if loadSmoke==1:
		iload.setBound(0.125,1); # neutralize offset
	if loadSmoke==2:
		iload.setBound(-1.,0); # reset outer bord
	interpolateGrid4d(iTarget, iload, offset=loadOffset , scale=loadScale );

	postProcChecker(iTarget);

	if doCacheInputs:
		iTarget.save(prepname);
	#print("Min max %f   %f phi final" % (iTarget.getMax(), iTarget.getMin()) ); 


def loadPhiDataSliced(phiDst, phiSliceName ):
	# improved but slow version of load! old one had NNB in time...
	loadPlaceGrid4d(fname=phiSliceName , phi=phiDst, offset=loadOffset , scale=loadScale , \
		fileIdxStart=fileSliceIdxStart, fileIdxEnd=fileSliceIdxEnd , debugSkipLoad=debugSkipLoad , spread=1. , \
		overrideSize=gs4d , overrideTimeOff=0. , overrideGoodRegion=0, loadTimeScale=loadTimeScale ); 

	postProcChecker(phiDst, True); 
	# #print("Min max after load   %f   %f" % (iTarget.getMax(), iTarget.getMin()) ); 




# ======================================================================================================================

def dumpOut(iadv, name, interval=1):
	if writePngs or writeUni:
		backupSources( "%s_0000"%(name) );
	if triMesh:
		backupSources( "%s_0000"%(meshPrefix) );

	meshCnt = 0;
	for t in my_range( timerangeStart,timerangeStop, interval):
		print("Step %f , range %f to %f " %(t, timerangeStart,timerangeStop) );
		getSliceFrom4d(src=iadv, dst=phiout, srct=t); 
		if(doIadvOrg):
			getSliceFrom4d(src=i0adv_org, dst=phiout_org, srct=t); 

		# show orgs , debugging!
		if not triMesh:
			getSliceFrom4d(src=i0   , dst=phiout0, srct=t); 
			getSliceFrom4d(src=i1   , dst=phiout1, srct=t); 
			#phiout0.createMesh(mesh0); phiout1.createMesh(mesh1); 

		initOthers(phiout0, phiout1, velout,veltout, velTmp_out,velTmp_tout, t);

		frameNum = int(t/interval);
		if (frameNum != t):
			print("NYI- only integer time interval supported here !"); exit(1);
	
		frameNum = int(t/interval); 
		writeResults(name, frameNum,meshCnt, phiout,meshout,gui,  writePngs  , triMesh    , projPpm    , writeUni   , sdfSmooth  , meshPrefix  );
		solv.step();


# similar to dumpOut, but reduced functionality
def loadSlice_dumpOut( dst, phiSrc, defoName, outName, interval=1):
	if writePngs or writeUni:
		backupSources( "%s_0000"%(outName) );
	if triMesh:
		backupSources( "%s_0000"%(meshPrefix) );

	meshCnt = 0;
	#loadAdvectTimeSlice_OptInit( defoName, dst, phiSrc, timerangeStart+0.5*interval, blendAlpha , loadTimeScale );
	for t in my_range( timerangeStart,timerangeStop, interval):
		loadAdvectTimeSlice       ( 99, defoName, dst, phiSrc, t+0.0*interval, blendAlpha , loadTimeScale ,  \
			defoOffset, defoScale, defoFactor,  overrideSize=gs4d,  debugVel=velout, debugVelT=veltout );

		# manually reset bound for 10 cells below.
		#
		# this is a hack, and necessary because:
		# optrun assumes at least a ten cell border -> cut off here for sanity check
		# this is only necessary to compare loadSlice_dumpOut and loadSlice_dumpOutOpt
		#
		# this is necessary for modes 5 and 6, as mode 6 should recreate the mode 5 data...
		# otherwise, we dont need this & shouldnt do it
		if activateChecker==5 or activateChecker==6:
			phiout.setBound( 0, 9 )

		# just show loaded data for debugging:
		#getSliceFrom4d(src=phiSrc, dst=phiout, srct=t); 

		# get regular data for comparison
		getSliceFrom4d   (src=i1 , dst=phioutSL  , srct=t ); # compare , note - only shows up later!
		getSliceFrom4dVec(src=vel, dst=velTmp_out, dstt=velTmp_tout, srct=t); 

		# calc "error" between the two versions
		calcObfDiff( phiout, phioutSL, phiDiff , velout, velTmp_out, veltout, velTmp_tout, velDiff, autoBrdRes );

		frameNum = int(t/interval);
		writeResults(outName, frameNum,meshCnt, dst,meshout,gui,  writePngs  , triMesh    , projPpm    , writeUni   , sdfSmooth  , meshPrefix  );
		solv.step();


# reduced version, using optimized functions
def loadSlice_dumpOutOpt( dst, phiSrc, defoName, outName, interval=1):
	meshCnt = 0;
	lats = 99; 
	loadAdvectTimeSlice_OptInit( lats, defoName, False,False, -1. ); # partial load and third load off
	for t in my_range( timerangeStart,timerangeStop, interval):
		loadAdvectTimeSlice_OptRun ( lats, defoName, dst, phiSrc, t+0.0*interval, blendAlpha , loadTimeScale ,  \
			defoOffset, defoScale, defoFactor,  overrideSize=gs4d,  debugVel=velout, debugVelT=veltout );

		frameNum = int(t/interval);
		writeResults(outName, frameNum, meshCnt, dst,meshout,gui,  writePngs  , triMesh    , projPpm    , writeUni   , sdfSmooth  , meshPrefix  );
		solv.step();
	loadAdvectTimeSlice_Finish( lats );



# ======================================================================================================================
# main


# first init offsets
(loadOffset,loadScale,defoOffset,defoScale,defoFactor) = initOffsets( gs4d,gs4dDefoLoad, autoBorder, repeatStartDist, loadTimeScale);
print("Defo (vel) scaling factor "+str(defoFactor)+", size scale " +str(defoScale)+" " ); 

# load input data

if doLoad==2:
	initTestCheckerboard(i0, vel);
	initTestCheckerboard(i1, vel);
	if 1 and (res==80):
		if testMode: 
			if getGenRefFileSetting():
				i0.save("checker80b.uni");  

			getSliceFrom4d(src=i0, dst=phiout, srct=55); 
			doTestGrid( sys.argv[0], "A_phi"  , solv, phiout   , threshold=1e-08 , thresholdStrict=1e-14  )
			solv.step(); solv.step(); solv.step();
		else:
			i0.save("checker80b.uni");  

		exit(1); # note - manually set right res above!

else:
	# i0
	if not (doSlicedLoad==2 or doSlicedLoad==3):
		loadPhiData(i0, dataname_i0);
	else:
		loadPhiDataSliced(i0, fileSliceLoadPrefix );
		# regular data, for comparison...
		loadPhiData(i1, dataname_i0);
		loadVelRescaled(vel, fileLoadPrefix, sDefoload, defoOffset, defoScale, defoFactor);

		if blendAlpha>0.:
			# compare...
			#i1adv.copyFrom(i1);
			advect4d(vel=vel, grid=i1, dtFac=blendAlpha);
		print("Max defo loaded %f " % (vel.getMaxValue()) ); 

	orgErr = 1.;

	# i1
	if not doSlicedLoad:
		loadPhiData(i1, dataname_i1);

		erri0i1 = calcLsDiff4d(i0=i0, i1=i1 ); 
		orgErr = erri0i1;
		print("Error between inputs %f " % (erri0i1) ); 

	# load previous defo data
	if doLoad==1:
		loadVelRescaled(vel, fileLoadPrefix, sDefoload, defoOffset, defoScale, defoFactor); 

# checkboard specials
if activateChecker>0 and activateChecker<4:
	lsFactor = 1.;
if activateChecker==4:
	lsFactor = -1./lsFactor; # saved data is pre-multiplied, neutralize...

if activateChecker==2: # load , just for display, offset vels!
	#vel.setConst(Vec4(1,1,1,1));
	initVecFromScalar(source=i0, target=vel);
	if 1 and (res==44):
		vel.save("defoCheckerX_044_vel.uni"); 

if 0 and activateChecker>=2: # load , just for display, offset vels!
	# for display
	vel.addConst(Vec4(1,1,1,1));
	vel.multConst(Vec4(0.5,0.5,0.5,0.5));

if loadMemRequired:
	del iload;

# run of!

if doRunOf or doLoadAndProject:
	if loadSmoke == 1:
		print("ERROR - loadSmoke was active! Needs to be disabled for actual OF runs"); exit(1);

	# run!
	if doRunOf:
		opticalFlowMultiscale4d(vel=vel, i0=i0, i1=i1,  wSmooth=0.001, wEnergy=0.0001 , cgAccuracy=cgAccuracy , \
			postVelBlur=postVelBlur, cfl=999, multiStep=multiStep , minGridSize=20, projSizeThresh=projSizeThreshold,
			doFinalProject=doFinalProject, resetBndWidth=autoBorder);
	else:
		loadVelRescaled(vel, fileLoadPrefix, sDefoload, defoOffset, defoScale, defoFactor);

	if 1 and doLsProject:
		i0adv.copyFrom(i0);
		advect4d(vel=vel, grid=i0adv);
		orgErr = calcLsDiff4d(i0=i0adv, i1=i1 ); 

		corrVelsOf4d( velTmp , vel, i0, i0, i1 , threshPhi=maxLSPDist, threshNorm=-1. , postVelBlur=lspBlur, resetBndWidth=autoBorder); 
		print("Max defo proj %f " % (velTmp.getMaxValue()) ); 

		# test advect 2
		if(doIadvOrg):
			i0adv_org.copyFrom(i0adv); # backup , compare w,w/o proj

	# main step done
	times.display();
	solv.step() 

	if doSave:
		vel.save(  "%s_vel.uni" % fileOutPrefix );
		#i0adv.save("%s_adv.uni" % fileOutPrefix );
		backupSources(fileLoadPrefix);


if not doSlicedLoad:
	print("Max vel %f " % (vel.getMaxValue()) ); 
	
	# do advect
	i0adv.copyFrom(i0);
	if 1: 
		if 1 or blendAlpha<(1.+1e-03):
			#vel.multConst(Vec4(blendAlpha)); print("PARTIAL ADVECT ACTIVE %f !\n" % blendAlpha); # debug, 50% interpol
			print("PARTIAL ADVECT ACTIVE %f !\n" % blendAlpha); 
			advect4d(vel=vel, grid=i0adv, dtFac=blendAlpha);
		else:
			i0adv.copyFrom(i1);  print("SHOW i1 ACTIVE!\n");  # debug, show input 1 , target
			#i0adv.copyFrom(i0);  print("SHOW i0 ACTIVE!\n"); # debug, show input 0 , source

	i0.multConst(    1.0/lsFactor );
	i1.multConst(    1.0/lsFactor );
	i0adv.multConst( 1.0/lsFactor );
else:
	# slices - only prepare i0
	i0.multConst(    1.0/lsFactor );



# ======================================================================================================================
# get 3d slices and output
# "normal" preview mode
if 1 and ( showGui or triMesh or projPpm or writeUni or testMode): 
	if blendAlpha<1.:
		print("PARTIAL ADVECT ACTIVEb %f !\n" % blendAlpha); # debug, 50% interpol

	if not doSlicedLoad:
		dumpOut              ( i0adv, screenPrefix, outInt);
	elif doSlicedLoad==2:
		loadSlice_dumpOut    ( phiout, i0, "%s_vel.uni" % fileLoadPrefix, screenPrefix, outInt);
	elif doSlicedLoad==3:
		loadSlice_dumpOutOpt ( phiout, i0, "%s_vel.uni" % fileLoadPrefix, screenPrefix, outInt);


# done!
solv.step()

if testMode:
	#if activateChecker==1: # done above!

	if activateChecker==2: 
		initOthers(phiout0, phiout1, velout,veltout, velTmp_out,velTmp_tout, 30);
		doTestGrid( sys.argv[0], "B_vel"  , solv, velout   , threshold=1e-08 , thresholdStrict=1e-14  )
		doTestGrid( sys.argv[0], "B_velt" , solv, veltout  , threshold=1e-08 , thresholdStrict=1e-14  )

	# alpha = 0. by default

	if activateChecker==3: 
		doTestGrid( sys.argv[0], "C_phi"  , solv, phiout   , threshold=1e-08 , thresholdStrict=1e-14  )
		doTestGrid( sys.argv[0], "C_phi0" , solv, phiout0  , threshold=1e-08 , thresholdStrict=1e-14  )
		doTestGrid( sys.argv[0], "C_phi1" , solv, phiout1  , threshold=1e-08 , thresholdStrict=1e-14  )
		doTestGrid( sys.argv[0], "C_vel"  , solv, velout   , threshold=1e-08 , thresholdStrict=1e-14  )
		doTestGrid( sys.argv[0], "C_velt" , solv, veltout  , threshold=1e-08 , thresholdStrict=1e-14  )

	if activateChecker==4: 
		doTestGrid( sys.argv[0], "D_phi"  , solv, phiout   , threshold=1e-08 , thresholdStrict=1e-14  )
		doTestGrid( sys.argv[0], "D_vel"  , solv, velout   , threshold=1e-08 , thresholdStrict=1e-14  )
		doTestGrid( sys.argv[0], "D_velt" , solv, veltout  , threshold=1e-08 , thresholdStrict=1e-14  )

	# alpha != 0. for the following (only phi matters then)

	if activateChecker==5: 
		doTestGrid( sys.argv[0], "E_phi"  , solv, phiout   , threshold=1e-08 , thresholdStrict=1e-14  )

	if activateChecker==6: 
		if getGenRefFileSetting():
			print("Skipping file gen for mode 5, use mode 4 to generate."); exit(1);

		# debug code only - load manually , and display
		#if 1:
			#print("Debug E_phi init/test")
			#phioutDeb.load("../testdata/checkerTest_E_phi.uni")
			#phioutCmp.copyFrom(phiout)
			#phioutCmp.addScaled(phioutDeb,-1.) 
			#solv.step()

		# reuse data from 4:
		doTestGrid( sys.argv[0], "E_phi"  , solv, phiout   , threshold=1e-08 , thresholdStrict=1e-14  )
		# note - new optimized version can give different results at outer border


