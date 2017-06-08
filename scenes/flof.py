#
# FlOF example scene
#
import sys;
import os
from manta   import *
from ofHelpers import *
paramUsed = [];

#
# quickstart guide:
#
# first gen, e.g.:      "manta .../dataGen2Drop.py  px 0  save4d 1  savehi 1 ; manta dataGen2Drop.py  px 1  save4d 1  savehi 1 "
#
# then gen defo         "manta .../flof.py dataid0 0  dataid1 1 mode 1 "
# check visually        "manta .../flof.py dataid0 0  dataid1 1 mode 2  showgui 1 "
# test full ouput       "manta .../flof.py dataid0 0  dataid1 1 mode 3  showgui 1 alpha 100 "
#
# gen other defo        "manta .../flof.py dataid0 1  dataid1 0 mode 1 "
# test 2way defo        "manta .../flof.py dataid0 1  dataid1 0 mode 2  showgui 1 alpha 100 "
#
# final 2way interpol   "manta .../flof.py  mode 3  showgui 1  twoway 1   alpha 50 writepngs 1 "
# (with png output)
# 

setDebugLevel(1); # set to 2 for kernel msgs

# main prefix for outputs
screenPrefixOrg = "out"; 
screenPrefixOrg = getParam("prefix",screenPrefixOrg, paramUsed);


# read cmd line & os params

dataPath="./"
dataPathId = int(os.getenv('OFBLEND_PATH', 0));
dataPathId = int( getParam("datapathid",dataPathId, paramUsed) );

# win/mac paths
if dataPathId==1:
	# default mac
	dataPath="/Users/sinithue/dataOfBlend/"
elif dataPathId==2:
	dataPath="I:/dataOfBlend/"
print("Using data path %s " % dataPath);

# optionally separate path for xl slice data
dataPathXl = int( getParam("datapathxl", 0, paramUsed) );
if dataPathXl==1:
	# default mac
	dataPathXl="/Users/sinithue/dataOfBlend/"
	print("Using data path xl %s " % dataPathXl);
elif dataPathXl==2:
	dataPathXl="I:/dataOfBlend/"
	print("Using data path xl %s " % dataPathXl);
else:
	dataPathXl=dataPath

# read params
# 1: source data id (can be interpreted differently)
# 2: target data id  
# 3: output mode  1=run OF, 2=show defo, 3=load and gen hires, 4=hires 2way
dataId0 = int( getParam("dataid0", 0, paramUsed) );
dataId1 = int( getParam("dataid1", 1, paramUsed) );
mainMode = int( getParam("mode",   1, paramUsed) );
	

# main settings container for optical flow & data processing	
ofsd = { }
ofsd["postProcMode"]    = 0
ofsd["isosurfaceOff"]   = 0
ofsd["isosurfaceMult"]  = 1.
ofsd["fourdFactor"]     = 0
ofsd["yFactor"]         = 1.
ofsd["autoBorder"]      = 0
ofsd["autoBrdRes"]      = 0
ofsd["lsFactor"]        = 1.
ofsd["maxDist"]         = 0
ofsd["cgAccuracy"]      = 0
ofsd["resBackup"]       = 0
ofsd["repeatStartDist"] = 0
ofsd["postVelBlur"]     = 0
ofsd["maxLSPDist"]      = 0
ofsd["lspBlur"]         = 0


# when loading hi-res SDF data, rescale values for new grid size
rescaleSdfSlices = False
sdfIsoOff        = 0.

# output settings:
# writePngs   png screen shot
# writeUni    write volumes
# projPpm     debugging, projections as ppm
# triMesh     generate triangle triMesh
# meshPrefix  filenames for blender
# sdfSmooth   post processing for output, smooth sdf values
# output settings, note sync with checkerTest!
writePngs = 0
triMesh   = 0
projPpm   = 0
writeUni  = 0
sdfSmooth = 0
sdfSmoothBord = 0
meshPrefix = "fluidsurface_final";

# remember whether we're dealing with SDFs or smoke data (from post-proc mode, which can be reset!)
haveLevelsetData = False;
# slight temporal filtering , good for liquids
doTimeSdfMerge = False
timeSdfIsFirst = True

setup       = 1;
resHiresOut = 100; 

# optional second interpolation dimension
doThirdLoad = dataId2 = 0;
doThirdLoad = int( getParam("thirdload", 0, paramUsed) );

# current data point positions (indexed by dataidX)
currp = [ vec3(0) , vec3(0) , vec3(0) , vec3(0) ];

if setup==1: # new two drop example

	ofsd["fourdFactor"] = 1.5 
	resLoad             = 40 
	resLoadXl           = 2*resLoad 
	resDefo             = 32 # could be anything...
	ofsd["postProcMode"]= 1 

	haveLevelsetData 	= True
	rescaleSdfSlices    = True
	sdfIsoOff           = -0.5 
	sdfSmooth           = 2
	doTimeSdfMerge      = True

	dataPrefix          = "out"   # main data prefix
	dataPrefixXl        = "outxl" # hires slice data prefix

	fileLoadPrefix = "%sdefo01_%03d_%03d_%03d" % ( dataPath, dataId0,dataId1,resDefo); 
	fileLoadInv    = "%sdefo01_%03d_%03d_%03d" % ( dataPath, dataId1,dataId0,resDefo); 
	fileOutPrefix  = fileLoadPrefix;  

	dataname_i0     = "%s/%s_r%03d_x%03d.uni" % (dataPathXl,dataPrefix, resLoad, dataId0 )
	dataname_i1     = "%s/%s_r%03d_x%03d.uni" % (dataPathXl,dataPrefix, resLoad, dataId1 )

	fileSliceLoadPrefix0= "%s/%s_r%03d_x%03d_%s.uni" % (dataPathXl, dataPrefixXl, resLoadXl, dataId0, "%04d" ) 
	fileSliceLoadPrefix1= "%s/%s_r%03d_x%03d_%s.uni" % (dataPathXl, dataPrefixXl, resLoadXl, dataId1, "%04d" )
	fileSliceIdxStart   = 0;
	fileSliceIdxEnd     = resLoad*4.5;

else:
	print("Error unkown setup %d " % (setup) ); exit(1);


# set other params more or less automatically

# main run mode
doRunOf     = False;
showGui     = 0;
doLoad      = 1;
doSave      = 0;
doCacheInputs    = 0;

# main resolution , set later on depending on mode
res = -1; 

# sliced loading: 0=off, 1=load single, 2=load uni sequence
doSlicedLoad    = 0; 
doOptLoadadv    = 1; # use optimized slice load version?
#doOptLoadadv    = 0; # NT_DEBUG off

# different modes for 2way interpol  (0=off, 1=on, 2=naive blend)
doTwoWayBlend   = 0; 

# main blending param
blendAlpha      = 1.0; 

# faster load for debugging
debugSkipLoad   = 9999; # = off

# for two dim interpol (optional, not for all setups)
thirdAlpha  = 0.;

# interval for output 
outputInterval = 1;
outputInterval = float( getParam("interval",outputInterval, paramUsed) );
if(outputInterval!=1.0): print("Note - changed outputInterval to %f " % outputInterval);

# get "partial load" fraction
partialLoad = -1.; # off
partialLoad = float( getParam("partialload",partialLoad, paramUsed) );

overrideToff   = 0.;  # accumulate offset

if mainMode==1:
	# run OF and save OF data
	res = resDefo;
	doLoad  = False; 
	doRunOf = True; 
	doSave  = True;   
	# optional - enable caching

elif mainMode==2:
	res = resLoad;
	# allow main resolution override
	res = int( getParam("resload",res, paramUsed) );

	# apply defo from OF run to 4d data for visual check
	doLoad = 1; 
	doSlicedLoad = 0; 
	doCacheInputs = 1;  
	doRunOf = False; 

elif mainMode==3:
	# load sliced hi res data and output deformed density files
	doLoad = 0; 
	doSlicedLoad = 2; 
	ofsd["postProcMode"] = 0; # no post-processing! leave original data
	if partialLoad>0.:
		ofsd["postProcMode"] = -1; # also no rep. start frame
	doRunOf = False;  
	#writeUni = 1; # uncomment to always write

	# allow resHiresOut override
	resHiresOut = int( getParam("hiresout",resHiresOut, paramUsed) );
	res = resHiresOut;


else:
	print("Unknown run mode %d " % mainMode); exit(1);

if mainMode>1:
	blendAlpha = int( getParam("alpha", int(100.*blendAlpha), paramUsed) ) * 0.01;
	if doThirdLoad:
		thirdAlpha = int( getParam("thirdalpha", int(100.*thirdAlpha), paramUsed) ) * 0.01;
print("Blend alpha debug: %f %f %d "%(blendAlpha,thirdAlpha,doThirdLoad));


# calculate blend weights for 2dim 3point interpol using barycentric coords
def recalcUvw():
	global uvw, blendAlpha, thirdAlpha;

	# problem specific geometry
	if setup==4:
		r1 = Vec3(dataId0, depth0, 0);
		r2 = Vec3(dataId1, depth1, 0);
		r3 = Vec3(dataId2, depth2, 0);
	elif setup==5:
		r1 = currp[0]; r2 = currp[1]; r3 = currp[2]
	else:
		print("Missing coords"); exit(1);

	pos = r1 + blendAlpha* (r2-r1) + thirdAlpha * (r3-r1);

	# uvw = Vec3(triu,triv,triw); # "return" value
	uvw = world2uvw( pos, r1,r2,r3 );
	print("UVW for setup %f %f %f from alphas %f,%f" %(uvw.x,uvw.y,uvw.z, blendAlpha, thirdAlpha) );
	#exit(1);

if doThirdLoad:
	recalcUvw()


# copy to settings struct
ofsd["resBackup"] = res;

# override defaults
writePngs    = int( getParam("writepngs",writePngs, paramUsed) ); 
writeUni     = int( getParam("writeuni" ,writeUni , paramUsed) ); 

doTwoWayBlend = int( getParam("twoway"  ,0             , paramUsed) );
doTWUnion     = int( getParam("twunion" ,0             , paramUsed) ); 
if doTWUnion and not doTwoWayBlend:
	doTwoWayBlend = 1

# output meshes? or ppms? (for debugging)
projPpm      = int( getParam("projppm"  ,projPpm  , paramUsed) ); 
triMesh      = int( getParam("trimesh"  ,triMesh  , paramUsed) ); 

debugSkipLoad = int( getParam("debugskip",debugSkipLoad, paramUsed) ); 
if(debugSkipLoad<9999): print("\nWARNING debug skip enabled %d!\n"%debugSkipLoad); # for quick tests

# note - defo "could" have different time dim (normally not)
ofsd["fourdFactorDefo"] = ofsd["fourdFactor"];
# for data load later on
ofsd["fourdFactorLoad"] = ofsd["fourdFactor"];

# data file names for caching
cacheprefix = "prep";

suffixThirdAlpha = ""
if doThirdLoad:
	suffixThirdAlpha = ("b%03d" % (thirdAlpha*100.));

# adjust output path
screenPrefix = "%s/%s_f%01dt%01d_a%03d%s" % (dataPath,screenPrefixOrg,dataId0,dataId1, int(blendAlpha*100.), suffixThirdAlpha );
meshPrefix = "%s/%s" % (dataPath, meshPrefix);


# ========================
# main OF params

multiStep           = 3;
ofsd["cgAccuracy"]  = 0.01;
dim                 = 3; # not really, four in practice...

# projection params
ofsd["maxLSPDist"]     = 10; 
ofsd["lspBlur"]        = 2.;

# misc parameters
ofsd["autoBorder"]     = 0.1;
ofsd["autoBrdRes"]     = int(res * ofsd["autoBorder"]);
ofsd["postVelBlur"]    = 4;
sdfSmoothBord = ofsd["autoBrdRes"]; # just to pass value along...

ofsd["repeatStartDist"] = float(res) * ofsd["autoBorder"] * 1.; # keep float
ofsd["lsFactor"]    = -0.1/20.0;
ofsd["maxDist"]     =  40; 
orgErr      = 0.;


# solver setup
(gs,gs4d)     = initRes( res, ofsd["yFactor"], ofsd["fourdFactorDefo"], 1.)
alloc4thDim   = gs4d.t
if partialLoad>0.:
	alloc4thDim   = gs4d.t * partialLoad / 1.;
	print("Partial sliced load active, alloc: %d  " % (alloc4thDim) );
alloc4thDim = int(alloc4thDim)

solv          = Solver(name='main', gridSize = gs, dim=dim, fourthDim=alloc4thDim );
solv.timestep = 1.0;
print("Main dimension "+vec4intString( Vec4(gs4d.x,gs4d.y,gs4d.z,alloc4thDim) ) +" , org length " +str(gs4d.t) );

timing  = Timings();

# setup relics  ...
flags = solv.create(FlagGrid)
flags.initDomain()
flags.fillGrid()

# allocate grid data
i0   = solv.create(Grid4Real)
if not doSlicedLoad:
	i1   = solv.create(Grid4Real)
	i0adv= solv.create(Grid4Real)
	if doTwoWayBlend>0 or doThirdLoad>0:
		i1adv    = solv.create(Grid4Real)
		velInv   = solv.create(Grid4Vec4)
	if doThirdLoad>0:
		i2       = solv.create(Grid4Real)
		thirdVel = solv.create(Grid4Vec4); 
else:
	if doTwoWayBlend>0 or doThirdLoad>0:
		i1       = solv.create(Grid4Real)
		i2       = solv.create(Grid4Real)

# preview intermediate step
i0adv_org=""
doIadvOrg = False;
if(doIadvOrg):
	i0adv_org = solv.create(Grid4Real)
if not doSlicedLoad:
	vel  = solv.create(Grid4Vec4)
if doRunOf:
	velTmp = solv.create(Grid4Vec4)

if doCacheInputs==1 and ofsd["postProcMode"]==3:
	print("Turn off caching for smoke post-proc!"); exit(1);

# loader
(gsload,gs4dload) = initRes( resLoad, ofsd["yFactor"], ofsd["fourdFactorDefo"], 1.); 
sload  = Solver(name='loader', gridSize = gsload, dim=dim, fourthDim=resLoad*ofsd["fourdFactorLoad"]);
print("Data load dimension "+vec4intString(gs4dload)+" " );

# loading memory only required if preprocessede files are not available
loadMemRequired = True;
if doSlicedLoad==1: 
	prepname = "%s_%s%04d.uni" %(dataname_i0, cacheprefix, res);
	if doCacheInputs and (os.path.isfile( prepname )):
		loadMemRequired = False;
if loadMemRequired:
	iload  = sload.create(Grid4Real)


# effective time range should not be finer than org data (optionally full time range)
timerangeStart = int(ofsd["autoBrdRes"]*1.0) + ofsd["repeatStartDist"] + 1
timerangeStop  = int( (res*ofsd["fourdFactor"] -  ofsd["autoBrdRes"]-1)*1. );

timerangeStart = int( getParam("start",timerangeStart, paramUsed) );
timerangeStop  = int( getParam("stop" ,timerangeStop , paramUsed) );

# optionally, allow as 0..1 factor
startfac = float( getParam("startfac", -1., paramUsed) );
stopfac  = float( getParam("stopfac" , -1., paramUsed) );
if startfac>=0.:
	timerangeStart = int(res*ofsd["fourdFactor"] * startfac);
if stopfac>=0.:
	timerangeStop = int(res*ofsd["fourdFactor"] * stopfac);

# init offset
if timerangeStop>0. and partialLoad>0.:
	alloc4thDim   = gs4d.t * partialLoad / 1.;
	# init first offset , use 1/2 of partial load window
	overrideToff = -( timerangeStart - int( alloc4thDim * (1./2.) ) );

showGui      = int( getParam("showgui",showGui, paramUsed) );
guiHideGrids = int( getParam("hidegrids", 1   , paramUsed) ); # on by default...

# ensure all cmd line params were used
checkUnusedParam(paramUsed);

# init UI
gui = 0;
if showGui and (GUI):
	gui = Gui()
	gui.show()
	gui.setCamPos (0,0,-1.2);
	if guiHideGrids:
		gui.toggleHideGrids();
	#gui.pause()
else:
	writePngs = 0; # make sure it's off...


# output data
phiout  = solv.create(LevelsetGrid)
meshout = solv.create(Mesh)
if(doIadvOrg):
	phiout_org = solv.create(LevelsetGrid)

# time merge test for SDFs
if doTimeSdfMerge:
	dstPrev1 = solv.create(LevelsetGrid)
	dstPrev2 = solv.create(LevelsetGrid)

# allocate debugging grids
if not doSlicedLoad:
	phiout1 = solv.create(LevelsetGrid)
	phiout0 = solv.create(LevelsetGrid)

	mesh1   = solv.create(Mesh)
	mesh0   = solv.create(Mesh)

	velout  = solv.create(Vec3Grid)
	veltout = solv.create(RealGrid)

	velTmp_out  = solv.create(Vec3Grid)
	velTmp_tout = solv.create(RealGrid)
else:
	# velocity debug
	velout  = solv.create(VecGrid)
	veltout = solv.create(RealGrid)

	if doTwoWayBlend==1:
		phiout1 = solv.create(LevelsetGrid)
	if doThirdLoad:
		phiout1 = solv.create(LevelsetGrid)
		phiout2 = solv.create(LevelsetGrid)

	if doTWUnion==1:
		phioutTmp = solv.create(LevelsetGrid); # tmp grid for union interpolation
	if doTWUnion==1 and (doThirdLoad):
		dstTmp = solv.create(LevelsetGrid); # another tmp grid needed...



# not always needed, but allocate solver for loading defo data
(gsDefoLoad,gs4dDefoLoad) = initRes( resDefo, ofsd["yFactor"], ofsd["fourdFactorDefo"], 1.) 
sDefoload  = Solver(name='defo_loader', gridSize = gsDefoLoad, dim=dim, fourthDim=int(gs4dDefoLoad.t) )
print("Defo (vel) load dimension "+vec4intString(gs4dDefoLoad)+" " )


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


# most parameters hidden in ofsd dictionary
def postProcessPhiData(iTarget, isSlice=False , loadOffset=vec4(0) , ofsd={} ):  
	printMinMax = True;
	if(printMinMax): print("Min max before post proc: %f to %f" % (iTarget.getMin(), iTarget.getMax()) ); 

	if ofsd["postProcMode"]==-1:
		print("Post processing off, no rep. start. fr.");
		return; # special for streaming load

	# start frame , note - make sure it hits full data!
	repeatFrame4d( iTarget, loadOffset.t + 1./float(ofsd["resBackup"]), ofsd["repeatStartDist"], bnd=0 ); 

	if   ofsd["postProcMode"]==0:
		print("Post processing off");

	elif ofsd["postProcMode"]==1:
		# extrapolate levelset
		print("Post processing as level set, offset %f, mult %f " %(ofsd["isosurfaceOff"],ofsd["isosurfaceMult"]) );

		if ofsd["isosurfaceOff"] != 0.:
			iTarget.addConst(ofsd["isosurfaceOff"]); 
		if ofsd["isosurfaceMult"] != 1.:
			iTarget.multConst(ofsd["isosurfaceMult"]); 

		iTarget.setBound(0.1, ofsd["autoBrdRes"] ); # reset border of LS to outside

		# extrapolate
		extrap4dLsSimple(phi=iTarget, distance=ofsd["maxDist"])  
		extrap4dLsSimple(phi=iTarget, distance=ofsd["maxDist"], inside=True)

		iTarget.multConst( ofsd["lsFactor"] );

	elif ofsd["postProcMode"]==3:
		print("Post processing with smokeWLT special!");

		# convert from "LS" to smoke
		iTarget.addConst  (-0.125) 
		iTarget.multConst (-1.0 * ofsd["lsFactor"]); # undo mult later on!
		iTarget.setBound  (0. , ofsd["autoBrdRes"]); 
		# note - does special set bound later on!

	if(printMinMax): print("Min max after post proc: %f to %f" % (iTarget.getMin(), iTarget.getMax()) ); 

# load data sets & registrations

# regular 4d data load (also handles caching)
def loadPhiData(iTarget, fname ):
	# load & resize
	global iload, loadOffset;

	prepname = "%s_%s%04d.uni" %(fname,cacheprefix,res);
	if doCacheInputs and (os.path.isfile( prepname )):
		iTarget.load(prepname);
		return;

	iload.load( fname ); 

	# special - reset sides before rescaling!
	if ofsd["postProcMode"]==3:
		iload.setBound(0.125,1); # neutralize offset

	interpolateGrid4d(iTarget, iload, offset=loadOffset , scale=loadScale );

	postProcessPhiData( iTarget, False, loadOffset, ofsd )

	if doCacheInputs:
		iTarget.save(prepname);

# sliced load from 3d sequence (for hi-res data)
def loadPhiDataSliced( phiDst, phiSliceName ):
	# loop over given range and load
	loadPlaceGrid4d(fname=phiSliceName , phi=phiDst, offset=loadOffset , scale=loadScale , \
		fileIdxStart=fileSliceIdxStart, fileIdxEnd=fileSliceIdxEnd , debugSkipLoad=debugSkipLoad , spread=1. , \
		overrideSize=gs4d , overrideTimeOff=overrideToff , overrideGoodRegion=0 , loadTimeScale=1. , rescaleSdfValues=rescaleSdfSlices, sdfIsoOff=sdfIsoOff, repeatStartFrame=ofsd["autoBorder"] ); 

	postProcessPhiData( phiDst, True, loadOffset, ofsd )

# simple quadratic ease in/out , t = time, b = startvalue, c = change in value, d = duration
def easeInOut (t, b, c, d):
	t = t / (d*0.5);
	if t < 1.: 
		return( c/2.*t*t + b );
	t = t-1;
	return( -c/2. * (t*(t-2.) - 1.) + b);

# blend weight
def blendW (t, l):
	w0 = ( 1.0+l - 2.0*t) / (1.0+l);
	w1 = (-1.0+l + 2.0*t) / (1.0+l);
	if(w0<0.): w0=0.;
	if(w1<0.): w1=0.;
	if(w0>1.): w0=1.;
	if(w1>1.): w1=1.;
	wu = 1. - w0 - w1;
	print("Blend for %f:     w0 %f, wu %f, w1 %f " % (t,w0,wu,w1) );
	return (w0,wu,w1);


# ======================================================================================================================

# apply deformation to whole 4d volume (using global params)
def deform4d():
	i0adv.copyFrom(i0);

	if blendAlpha<(0.+1e-03) and not doThirdLoad: 
		print("Showing i0!\n");  # debug info, show input 1 , target

	elif blendAlpha<(1.+1e-03) or doThirdLoad: # for third-load, always apply (dont skip to i1)
		print("Applying defo with alpha=%f !\n" % blendAlpha); 

		if blendAlpha>(1.0-1e-03) and ofsd["postProcMode"]==0:
			massBefore = debugGridAvg4d(i0adv,ofsd["autoBrdRes"]);

		advect4d(vel=vel, grid=i0adv, dtFac=blendAlpha);

		if doTwoWayBlend==1:
			# reverse advect
			i1adv.copyFrom(i1);
			advect4d(vel=velInv, grid=i1adv, dtFac= (1.-blendAlpha) );
			i0adv.multConst(1.-blendAlpha);
			i0adv.addScaled(i1adv, blendAlpha);

		elif doTwoWayBlend==2: # naive blend, for comparison only
			print("\nShowing naive blend result!\n" );
			i0adv.multConst(1.-blendAlpha);
			i0adv.addScaled(i1, blendAlpha);

		# only for full advect, recalc diff
		if blendAlpha>(1.0-1e-03) and (orgErr>1e-03):
			if ofsd["postProcMode"]==0:
				massAfter = debugGridAvg4d(i0adv,ofsd["autoBrdRes"]);
				print("Avg mass before/after adv %f,%f; factor %f" % (massBefore,massAfter,(massBefore/massAfter)) );
			avgErr = calcLsDiff4d(i0=i0adv, i1=i1 ); 
			print("Final ls diff = %f , org fac %f\n" % (avgErr, (avgErr/orgErr) ) );

	else:
		# debug, show input 1 , target
		i0adv.copyFrom(i1);  print("Showing i1!\n");  

	# second defo for two-dim data (note - copied from first defo)
	if doThirdLoad:
		if thirdAlpha<(0.+1e-03): 
			print("Not active: vel2!\n"); 

		print("Applying defo2 with third-alpha=%f !\n" % thirdAlpha); 
		advect4d(vel=vel,      grid=thirdVel, dtFac=blendAlpha ); # align
		advect4d(vel=thirdVel, grid=i0adv,    dtFac=thirdAlpha );



def dumpOut(iadv, outName, interval=1, nrOverride=-1):
	meshCnt = 0;
	for t in my_range( timerangeStart,timerangeStop, interval):
		getSliceFrom4d(src=iadv, dst=phiout, srct=t); 
		phiout.createMesh(meshout);
		if(doIadvOrg):
			getSliceFrom4d(src=i0adv_org, dst=phiout_org, srct=t); 

		# show orgs , debugging!
		if not triMesh:
			getSliceFrom4d(src=i0   , dst=phiout0, srct=t); 
			getSliceFrom4d(src=i1   , dst=phiout1, srct=t); 
			phiout0.createMesh(mesh0); phiout1.createMesh(mesh1); 

		initOthers(phiout0, phiout1, velout,veltout, velTmp_out,velTmp_tout, t);

		if nrOverride<0:
			frameNum = int(t/interval);
			if (frameNum != t):
				print("NYI- only integer time interval supported here !"); exit(1); 
			writeResults(outName, frameNum, meshCnt, phiout,meshout,gui, writePngs  , triMesh    , projPpm    , writeUni   , sdfSmooth,sdfSmoothBord  , meshPrefix  ); 
			meshCnt = meshCnt+1; 
		else:
			print("Nr override active!");
			writeResults(outName, nrOverride,nrOverride, phiout, meshout,gui, writePngs  , triMesh    , projPpm    , writeUni   , sdfSmooth,sdfSmoothBord  , meshPrefix  );
		solv.step();

# similar to dumpOut, but more functionality (optimized loads, multiple defos)
isInitialized = False;
def loadSlice_dumpOut( dst, phiSrc, defoName0, outName, interval=1, phiSliceName="", \
					   defoName1="", dst1=0, phiSrc1=0, phiSliceName1="", \
					   defoName2="",   dst2=0, phiSrc2=0, phiSliceName2="", \
					   doCleanup=True , nrOverride=-1 , \
					   defoName3="",   dst3=0, phiSrc3=0, phiSliceName3="" , defoAniFac = 1.):
	global overrideToff, isInitialized, timeSdfIsFirst;
	meshCnt =  0;
	lats1   =  7; # arbitrary IDs
	lats2   = 13; # actual value doesnt matter
	lats    = [7,            13,            73,                11,               17   ]; # last one unused for now...
	latName = [defoName0,    defoName1,     defoName2,         defoName3,        "-"  ];
	latSrc  = [phiSrc,       phiSrc1,       phiSrc2,           phiSrc3,          0    ];
	latDst  = [dst,          dst1,          dst2,              dst3,             0    ];
	latSlic = [phiSliceName, phiSliceName1, phiSliceName2,     phiSliceName3,    0    ];
	latNum = 0;

	# offset by fraction of single frame
	frameOffset = 0.5;  

	loadFunc = loadAdvectTimeSlice;
	if doOptLoadadv: loadFunc = loadAdvectTimeSlice_OptRun; # optimized version

	latNum = 1;
	if doTwoWayBlend: latNum = 2;
	if doThirdLoad  : latNum = 3;
	useDefoVols = (latNum>2)

	if not isInitialized:
		for i in range(0,latNum):
			print("Advect init %d %s " % (i,latName[i]) );
			loadAdvectTimeSlice_OptInit(    lats[i], latName[i] , useDefoVols, False, partialLoad ); 
			if useDefoVols:
				loadAdvectTimeSlice_OptAdd( lats[i], latName[(i+1)%latNum] );
		isInitialized = True

	for t in my_range( timerangeStart,timerangeStop, interval):
		print("Time %f, in range %f to %f, inter %f " % (t,timerangeStart,timerangeStop,interval) );

		# activate partial loading , shift & load whenever necessary
		print("Pos in slice %d , dim %d , toff %d " % ( (t+overrideToff) , alloc4thDim, overrideToff) );

		keep      = int(alloc4thDim-2); 
		keepCheck = int(alloc4thDim*0.5)+1	

		if partialLoad>0. and (t + overrideToff) > keepCheck: 
			shift = (alloc4thDim - keep); # shift amount (remainder)
			overrideToff = overrideToff - shift;
			print("Pos shift %d " % (shift) );

			for i in range(0,latNum):
				shiftForwGrid4d( phi=latSrc[i] , overrideGoodRegion=keep );

			for i in range(0,latNum): loadPlaceGrid4d(fname=latSlic[i] , phi=latSrc[i], offset=loadOffset , scale=loadScale , \
				fileIdxStart=fileSliceIdxStart, fileIdxEnd=fileSliceIdxEnd , debugSkipLoad=debugSkipLoad , spread=1. , \
				overrideSize=gs4d , overrideTimeOff=overrideToff , overrideGoodRegion=keep, loadTimeScale=1. , rescaleSdfValues=rescaleSdfSlices ,sdfIsoOff=sdfIsoOff,  repeatStartFrame=ofsd["autoBorder"] ); 

		# interpolation go...

		if not doThirdLoad:
			# 2 data sets , 1d line interpolations

			loadFunc ( lats[0], latName[0], latDst[0], latSrc[0], t+frameOffset*interval, blendAlpha , 1. ,  \
				defoOffset, defoScale, defoFactor, overrideSize=gs4d, overrideTimeOff=overrideToff, zeroVel=False , thirdAlpha = thirdAlpha, bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );

			alphaInv = 1. - blendAlpha;
			if doTwoWayBlend==2: # recreate naive blend , alpha = 0 for second data
				alphaInv = 0.;

			if doTwoWayBlend==1 and doTWUnion:
				loadFunc ( lats[1], latName[1], latDst[1], latSrc[1], t+frameOffset*interval, alphaInv , 1. ,  defoOffset, defoScale, defoFactor, gs4d, overrideToff, zeroVel=False, bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );
				(w0,wu,w1) = blendW(blendAlpha, 0.1);

				# blend union
				phioutTmp.copyFrom(dst); 
				phioutTmp.join(latDst[1]); 

				dst.multConst(           w0);
				dst.addScaled(phioutTmp, wu);
				dst.addScaled(latDst[1], w1);
				print("Blended union with weights %f,%f,%f " % (w0,wu,w1) );

			elif doTwoWayBlend==1 or doTwoWayBlend==2:
				loadFunc ( lats[1], latName[1], latDst[1], latSrc[1], t+frameOffset*interval, alphaInv , 1. ,  defoOffset, defoScale, defoFactor, gs4d, overrideToff, zeroVel=False, bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );

				# include additional blending factors
				dst.multConst(           (1.-blendAlpha) );
				dst.addScaled(latDst[1], (   blendAlpha) );

		elif doThirdLoad: # doThirdLoad , 3 data sets

			if doTwoWayBlend==0: # simple interpol , 1 src data point only
				loadFunc ( lats[0], latName[0], latDst[0], latSrc[0], t+frameOffset*interval, (uvw.y + uvw.z) , 1. ,defoOffset, defoScale, defoFactor, overrideSize=gs4d, overrideTimeOff=overrideToff, zeroVel=False , thirdAlpha = uvw.z , bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );

			if doTwoWayBlend==1 and doTWUnion:
				print("Advect stanton 2dim with uvws %f %f %f " % (uvw.x,uvw.y,uvw.z) );
				loadFunc ( lats[0], latName[0], latDst[0], latSrc[0], t+frameOffset*interval, (uvw.y + uvw.z) , 1. ,defoOffset, defoScale, defoFactor, overrideSize=gs4d, overrideTimeOff=overrideToff, zeroVel=False , thirdAlpha = uvw.z , bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );
				loadFunc ( lats[1], latName[1], latDst[1], latSrc[1], t+frameOffset*interval, (uvw.z + uvw.x) , 1. ,defoOffset, defoScale, defoFactor, overrideSize=gs4d, overrideTimeOff=overrideToff, zeroVel=False , thirdAlpha = uvw.x , bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );
				loadFunc ( lats[2], latName[2], latDst[2], latSrc[2], t+frameOffset*interval, (uvw.x + uvw.y) , 1. ,defoOffset, defoScale, defoFactor, overrideSize=gs4d, overrideTimeOff=overrideToff, zeroVel=False , thirdAlpha = uvw.y , bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );

				# calculate blendW 
				(w0,w1,w2,w_uv,w_vw,w_uw,w_uvw) = weights2dim(uvw)
				sum  =w0+w1+w2+w_uv+w_vw+w_uw+w_uvw;
				print( "uvw-summ %4.3f,%4.3f,%4.3f     w0=%4.3f, w1=%4.3f, w2=%4.3f, w_uv=%4.3f, w_vw=%4.3f, w_uw=%4.3f, w_uvw=%4.3f , sum %f" % (uvw.x,uvw.y,uvw.z,  w0,w1,w2,w_uv,w_vw,w_uw,w_uvw, sum) );
				# calculate unions w weight >0
				dstTmp.setConst(0.);

				if w_uv>1e-03: 
					phioutTmp.copyFrom (latDst[0]);
					phioutTmp.join     (latDst[1]);
					dstTmp.addScaled(phioutTmp, w_uv );
				if w_vw>1e-03: 
					phioutTmp.copyFrom (latDst[1]);
					phioutTmp.join     (latDst[2]);
					dstTmp.addScaled(phioutTmp, w_vw );
				if w_uw>1e-03: 
					phioutTmp.copyFrom (latDst[0]);
					phioutTmp.join     (latDst[2]);
					dstTmp.addScaled(phioutTmp, w_uw );

				if w_uvw>1e-03: 
					phioutTmp.copyFrom (latDst[0]);
					phioutTmp.join     (latDst[1]);
					phioutTmp.join     (latDst[2]);
					dstTmp.addScaled(phioutTmp, w_uvw ); 

				# last, add orgs
				if w0>1e-03: dstTmp.addScaled(latDst[0], w0 );
				if w1>1e-03: dstTmp.addScaled(latDst[1], w1 );
				if w2>1e-03: dstTmp.addScaled(latDst[2], w2 );

				dst.copyFrom(dstTmp); # finally copy over, overwrite deformed i0

			elif doTwoWayBlend==1: # simple interpol
				print("Advect 2dim with uvws %f %f %f " % (uvw.x,uvw.y,uvw.z) );
				loadFunc ( lats[0], latName[0], latDst[0], latSrc[0], t+frameOffset*interval, (uvw.y + uvw.z) , 1. ,defoOffset, defoScale, defoFactor, overrideSize=gs4d, overrideTimeOff=overrideToff, zeroVel=False , thirdAlpha = uvw.z , bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );
				loadFunc ( lats[1], latName[1], latDst[1], latSrc[1], t+frameOffset*interval, (uvw.z + uvw.x) , 1. ,defoOffset, defoScale, defoFactor, overrideSize=gs4d, overrideTimeOff=overrideToff, zeroVel=False , thirdAlpha = uvw.x , bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );
				loadFunc ( lats[2], latName[2], latDst[2], latSrc[2], t+frameOffset*interval, (uvw.x + uvw.y) , 1. ,defoOffset, defoScale, defoFactor, overrideSize=gs4d, overrideTimeOff=overrideToff, zeroVel=False , thirdAlpha = uvw.y , bordSkip=ofsd["autoBrdRes"] , defoAniFac=defoAniFac );

				if 0:
					# NT DEBUG, only dst , leave unmodified
					print("Note: No blend, only dst active!");
				else:
					dst.multConst(           uvw.x );
					dst.addScaled(latDst[1], uvw.y );
					dst.addScaled(latDst[2], uvw.z );

			if doTwoWayBlend==2: # not for naive blend!
				print("No blend, naive on");

		if doTimeSdfMerge: # "motion blur", temporal blur effect
			if timeSdfIsFirst:
				dstPrev2.copyFrom(dst)
				dstPrev1.copyFrom(dst)
				timeSdfIsFirst = False 
			dstPrev2.copyFrom(dstPrev1)
			dstPrev1.copyFrom(dst)
			dst.join(dstPrev2)
			simpleBlurSpecial( dst, 1 , 0.); 
			print("Time merge sdf mode enabled!");

		if haveLevelsetData == 1: # for levelset data reset border for LS
			#print("Set bound %d "%(ofsd["autoBrdRes"]+50));
			dst.setBound(0.5, ofsd["autoBrdRes"]+1); 

		if nrOverride<0:
			frameNum = int(t/interval); 
			writeResults(outName, frameNum,meshCnt,      dst, meshout,gui,  writePngs  , triMesh    , projPpm    , writeUni   , sdfSmooth,sdfSmoothBord  , meshPrefix );
			meshCnt = meshCnt+1; 
		else:
			print("Nr override active!");
			writeResults(outName, nrOverride,nrOverride, dst, meshout,gui,  writePngs  , triMesh    , projPpm    , writeUni   , sdfSmooth,sdfSmoothBord  , meshPrefix );

		#timing.display()
		solv.step(); 

	# cleanup globals
	if doOptLoadadv and doCleanup:
		for i in range(0,latNum):
			loadAdvectTimeSlice_Finish( lats[i] );


# ======================================================================================================================
# main


# load input data

# first init offsets
(loadOffset,loadScale,defoOffset,defoScale,defoFactor) = initOffsets( gs4d,gs4dDefoLoad, ofsd["autoBorder"], ofsd["repeatStartDist"], 1.);

print("Data scaling factor "+str(loadScale)+", offset " +str(loadOffset)+" " ); 

# i0
if not doSlicedLoad==2:
	loadPhiData(i0, dataname_i0);
else:
	loadPhiDataSliced(i0, fileSliceLoadPrefix0 );
	if doTwoWayBlend==1:
		loadPhiDataSliced(i1, fileSliceLoadPrefix1 );
orgErr = 1.;

# i1
if not doSlicedLoad:
	loadPhiData(i1, dataname_i1);

	erri0i1 = calcLsDiff4d(i0=i0, i1=i1 ); 
	orgErr = erri0i1;
	print("Error between inputs %f " % (erri0i1) ); 

# i2
if doThirdLoad:
	if not doSlicedLoad:
		loadPhiData(i2, dataname_i2);
	else:
		loadPhiDataSliced(i2, fileSliceLoadPrefix2 );

# load previous defo data
if doLoad:
	loadVelRescaled(vel, fileLoadPrefix, sDefoload, defoOffset, defoScale, defoFactor);

	if doTwoWayBlend==1:
		# load reverse defo  
		loadVelRescaled(velInv, fileLoadInv , sDefoload, defoOffset, defoScale, defoFactor);
		print("Two-way (1dim) interpolation active" );

	if doThirdLoad==1:
		loadVelRescaled(thirdVel, fileLoadThird1 , sDefoload, defoOffset, defoScale, defoFactor);
		print("Two-dim interpolation active" );

	vel.setBoundNeumann(ofsd["autoBrdRes"]+0); 
	if doTwoWayBlend==1: velInv.setBoundNeumann(ofsd["autoBrdRes"]+0); 
	if doThirdLoad>0:    thirdVel.setBoundNeumann(ofsd["autoBrdRes"]+0); 

if loadMemRequired:
	del iload;

# run of!

if doRunOf:
	if ofsd["postProcMode"]>=2:
		print("ERROR - smoke post-proc was active! Needs to be disabled for actual OF runs"); exit(1);

	# run!
	if doRunOf and 1:
		opticalFlowMultiscale4d(vel=vel, i0=i0, i1=i1,  wSmooth=0.001, wEnergy=0.0001 , cgAccuracy=ofsd["cgAccuracy"] , \
			postVelBlur=ofsd["postVelBlur"], cfl=999, multiStep=multiStep , minGridSize=20, 
			doFinalProject=True, resetBndWidth=ofsd["autoBorder"] );

	# main OF step done
	solv.step() 

	if doSave:
		vel.save(  "%s_vel.uni" % fileOutPrefix );


# before applying defo, undo scaling for OF
if not doSlicedLoad:
	i0.multConst(    1.0/ofsd["lsFactor"] );
	i1.multConst(    1.0/ofsd["lsFactor"] );

	if doTWUnion or doThirdLoad:
		print("These modes are only supported with slices!"); exit(1);


# ======================================================================================================================
# get 3d slices and output result in some way

if ( showGui or triMesh or projPpm or writeUni or writePngs ): 

	doClean = True;
	defoAniFac = 1.0 
	steps = 1;

	# iterate for single step mode
	nrOverride = -1;
	for step in range(steps):
		# apply deformation once
		if step==0 and not doSlicedLoad:
			deform4d()

		# actual output modes
		if not doSlicedLoad:
			dumpOut(i0adv, screenPrefix, outputInterval, nrOverride);
		else:
			if doThirdLoad>0:
				loadSlice_dumpOut( dst=phiout, phiSrc=i0, defoName0=("%s_vel.uni" % fileLoadPrefix), outName=screenPrefix, \
					interval=outputInterval, phiSliceName=fileSliceLoadPrefix0 ,
					defoName1=("%s_vel.uni" % fileLoadThird1), dst1=phiout1, phiSrc1=i1, phiSliceName1=fileSliceLoadPrefix1 ,
					defoName2=("%s_vel.uni" % fileLoadThird2), dst2=phiout2, phiSrc2=i2, phiSliceName2=fileSliceLoadPrefix2 , nrOverride=nrOverride, doCleanup=doClean, defoAniFac=defoAniFac  );
			elif doTwoWayBlend>0:
				loadSlice_dumpOut( dst=phiout, phiSrc=i0, defoName0=("%s_vel.uni" % fileLoadPrefix), outName=screenPrefix, \
					interval=outputInterval, phiSliceName=fileSliceLoadPrefix0,   \
					defoName1=("%s_vel.uni" % fileLoadInv), dst1=phiout1, phiSrc1=i1, phiSliceName1=fileSliceLoadPrefix1 , nrOverride=nrOverride, doCleanup=doClean , defoAniFac=defoAniFac  );
			else:
				loadSlice_dumpOut( dst=phiout, phiSrc=i0, defoName0=("%s_vel.uni" % fileLoadPrefix), outName=screenPrefix, \
					interval=outputInterval, phiSliceName=fileSliceLoadPrefix0 , nrOverride=nrOverride, doCleanup=doClean , defoAniFac=defoAniFac  );
else:
	print("No frame output.");


# done!
solv.step()


