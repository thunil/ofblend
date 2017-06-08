#
# 2 Drop data generation setup
# run with: ./manta .../dataGen2Drop.py   px 0    save4d 1  savehi 1
#           ./manta .../dataGen2Drop.py   px 1    save4d 1  savehi 1
#
# add "res 80" for high res data
# 
import sys;
import os
import shutil
from manta import *
from ofHelpers import *
paramUsed = [];

# output prefix
outPrefix = "out"; 
outPrefix = getParam("prefix",outPrefix, paramUsed);

# UI on/off
showGui           = int( getParam("showgui", 0, paramUsed) ); 

# modes
do4dOutput        = 1;
timeinterval      = 3; # only save every n'th frame for 4d volume
timeintervalHires = 1; # only save every n'th frame for hi res slices

doHires  =   int( getParam("hires"   , 1, paramUsed) );  # up res version?
upresFac = float( getParam("upresfac", 2, paramUsed) );  # increase resolution by which factor?

doOutput4dVol    = int( getParam("save4d"     , 0, paramUsed) ); # save full 4d volume
doSaveHires      = int( getParam("savehi"     , 0, paramUsed) ); # save hi res slices?
doOutput4dSlices = int( getParam("save4dhi"   , 0, paramUsed) ); # in combi with doHires, and do4dOutput only 
doWritePngs      = int( getParam("writepngs"  , 0, paramUsed) ); # screenshots?

# for load, only for hires
scrPrefixIn = "<<empty>>";

# ===

# read params

res        = int( getParam("res", 40, paramUsed) );
dropPosPx  = int( getParam("px",  0, paramUsed) );
if res<1:
	print("Give at least 4 params: manta gen.py  res R   px X py Y pz Z  " ); exit(1);

checkUnusedParam(paramUsed);

# solver params
dim = 3;
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

s.timestepMin = 0.1   # time step range
s.timestepMax = 2.0
s.cfl         = 2.0   
timesteps     = res * 1.5 
s.frameLength = 30./timesteps

gravVec           = vec3(0,-0.0025,0); 

s.timestep   = (s.timestepMax+s.timestepMin)*0.5 
partDisc     = 2
minParticles = pow(partDisc,dim)
timings = Timings()

# size of particles 
radiusFactor      = 1.4
radiusFactorHires = 1.0

# prepare grids and particles
flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)

vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
phiObs   = s.create(LevelsetGrid)
phiTmp   = s.create(LevelsetGrid)

pp       = s.create(BasicParticleSystem) 
pVel     = pp.create(PdataVec3) 

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

fractions = s.create(MACGrid)


if do4dOutput:
	resf = res;
	if doOutput4dSlices:
		resf = int(upresFac*res)
		print("Note, 4d out & do hires slices! Extrap & slices output at full resolution..."); 
	fourthDim = timesteps
		
	datout = Solver(name='main', gridSize = vec3(resf,resf,resf), dim=3, fourthDim=fourthDim )
	phiout = datout.create(Grid4Real)
	print("4d res: %d,%d,%d,%d " % (resf,resf,resf,fourthDim) );


# init geo ==============================================================================================

flags.initDomain(boundaryWidth=1)

phi.setConst(999.);
phiObs.setConst(999.);

if 1:
	# falling drops 
	fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.15,1.0)) # basin

	# drop heights
	height1 = 0.5
	height2 = 0.7 

	height1 = 0.5
	height2 = 0.25 

	if dropPosPx==0:
		# v01
		dropCenter1 = vec3(0.33, height1, 0.33)
		dropCenter2 = vec3(0.66, height2, 0.66)
	elif dropPosPx==1:
		# v02 , swapped 
		dropCenter1 = vec3(0.66, height1, 0.66)
		dropCenter2 = vec3(0.33, height2, 0.33)
	else:
		print("Unknown 2 drop version! %d " % dropPosPx); exit(1);

	rad = 0.1
	fluidDrop1  = s.create(Sphere, center=gs*dropCenter1, radius=res*rad)
	fluidDrop2  = s.create(Sphere, center=gs*dropCenter2, radius=res*rad)

	# new

	rad = 0.08
	if dropPosPx==0:
		# v01
		dropCenter1 = vec3(0.25, 0.7 , 0.33)
		dropCenter2 = vec3(0.75, 0.4 , 0.66)

		fluidDrop1  = s.create(Sphere, center=gs*dropCenter1, radius=res*(rad*1.5))
		fluidDrop2  = s.create(Box,    p0=gs*(dropCenter2-Vec3(rad)), p1=gs*(dropCenter2+Vec3(rad)) )
	elif dropPosPx==1:
		# v02 , swapped 
		dropCenter1 = vec3(0.45, 0.7 , 0.33)
		dropCenter2 = vec3(0.55, 0.4 , 0.66)

		fluidDrop1  = s.create(Box,    p0=gs*(dropCenter1-Vec3(rad)), p1=gs*(dropCenter1+Vec3(rad)) )
		fluidDrop2  = s.create(Sphere, center=gs*dropCenter2, radius=res*(rad*1.5))
	else:
		print("Unknown 2 drop version! %d " % dropPosPx); exit(1);

	phi.join( fluidBasin.computeLevelset() );
	phi.join( fluidDrop1.computeLevelset() );
	phi.join( fluidDrop2.computeLevelset() );
	
	radiusFactorHires = 1.5


if 1: # sample pool and initial region
	pVel.setSource( vel, isMAC=True )
	sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=partDisc, randomness=0.06 )
	mapGridToPartsVec3(source=vel, parts=pp, target=pVel )
	flags.updateFromLevelset(phi)


# surface
sres = int(upresFac*res)
srs = Solver(name='surf', gridSize = vec3(sres,sres,sres) , dim=3)
if doHires>0:
	sflags    = srs.create(FlagGrid)
	sphi      = srs.create(LevelsetGrid)
	spp       = srs.create(BasicParticleSystem) 
	spindex = srs.create(ParticleIndexSystem) 
	sgpi    = srs.create(IntGrid)
	sflags.initDomain(boundaryWidth=0)

	# meshes, show fine
	smesh = Mesh(parent=srs); # fine
	mesh  = Mesh(parent=s  ); # low res

else:
	mesh  = Mesh(parent=s  ); # low res, show
	smesh = Mesh(parent=srs); # fine


if showGui and (GUI):	
	gui = Gui()
	gui.setCamPos (0,0,-1.2);
	gui.toggleHideGrids();
	gui.show( dim==2 )
	#gui.pause()

lastFrame = -1
updateFractions(flags=flags, phiObs=phiObs, fractions=fractions)

# backup...
filenameSuffix   =  'r%03d_x%03d' % ( res, dropPosPx); 
filenameSuffixXl =  'r%03d_x%03d' % (sres, dropPosPx); 

#main loop
while s.frame < 9999:
	
	if 1:
		maxVel = vel.getMaxValue()
		s.adaptTimestep( maxVel )
		print("Time %f, maxvel %f, timestep %f " %(s.timeTotal, maxVel, s.timestep) );

		pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=True ) 
		deleteTopParts( parts=pp, maxHeight = res*0.95 );
		limitPvel( parts=pp, pvel=pVel, max=7.);

		# make sure we have velocities throughout liquid region
		mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
		extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) 
		markFluidCells( parts=pp, flags=flags )

		# create approximate surface level set, resample particles
		gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
		averagedParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor , 1, 1 );
		extrapolateLsSimple(phi, 4);

		if doHires>0 and upresFac>0.:
			# hi res surface!
			scalePartPos(pp, upresFac);
			sphi.clear()
			gridParticleIndex( parts=pp , flags=sflags, indexSys=spindex, index=sgpi )
			if upresFac>=3.0:
				print("Radius fac 2! for upres3...");
				averagedParticleLevelset( pp, spindex, sflags, sgpi, sphi , 2.0*radiusFactor*radiusFactorHires , 1, 1 ) 
			elif upresFac>=2.0:
				averagedParticleLevelset( pp, spindex, sflags, sgpi, sphi , 1.5*radiusFactor*radiusFactorHires , 1, 1 ) 
			else: # upresFac 1.0:
				averagedParticleLevelset( pp, spindex, sflags, sgpi, sphi , 1.0*radiusFactor*radiusFactorHires , 1, 1 ) 

			extrapolateLsSimple(phi=sphi, distance=2, inside=False)
			extrapolateLsSimple(phi=sphi, distance=2, inside=True)

			sphi.setBound(value=0., boundaryWidth=1)
			scalePartPos(pp, 1./upresFac);

		addGravity(flags=flags, vel=vel, gravity=gravVec )

		setBoundMAC(vel, vec3(0.), 2, True); 
		solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)

		# make sure we have proper velocities
		extrapolateMACSimple( flags=flags, vel=vel, distance=(int(maxVel)+25) ); 
		setBoundMAC(vel, vec3(0.), 2, True);

		flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )

		pVel.setSource( vel, isMAC=True ) 
		adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor ) 

	else:
		sphi.load( '%sxl_%s_%04d.uni' % (scrPrefixIn, filenameSuffixXl,s.frame) ); # v04



	# --- sim done


	# store slice
	t   = s.frame 
	tlo = t/timeinterval;
	thi = t/timeintervalHires;
	if t>=0 and (lastFrame != s.frame) and (t%timeinterval==0):
		print("Timestep %d/%d "%(tlo, timesteps) );
		if 1 and do4dOutput:
			print("placing grid to t=%d " % (tlo) );
			if doHires>0:
				if doOutput4dSlices:
					placeGrid3d(sphi,phiout,dstt=tlo) 
				else:
					interpolateGrid( target=phiTmp, source=sphi , orderSpace=2 )
					placeGrid3d(phiTmp,phiout,dstt=tlo) 
			else:
				placeGrid3d(phi,phiout,dstt=tlo) 

			if(tlo >= timesteps):
				# extrap 4d 
				if 1:
					maxDist = 10;
					print("4d extrapol with distance %f " % maxDist);
					extrap4dLsSimple(phi=phiout, distance=maxDist)  
					extrap4dLsSimple(phi=phiout, distance=maxDist, inside=True)

				if 1 and doOutput4dSlices:
					# output slices again
					for i in range(int(timesteps)):
						getSliceFrom4d(src=phiout, dst=sphi, srct=i); 
						sphi.save( '%sxle_%s_%04d.uni' % (outPrefix, filenameSuffixXl,i) );

				# full 4d volume
				if doOutput4dVol:
					phiout.save( '%s_%s.uni' % (outPrefix, filenameSuffix) );

	# output 3d slices
	if t>=0 and (lastFrame != s.frame) and (t%timeintervalHires==0):
		# save hi res slices (only for do4dOutput=0)
		if doSaveHires:
			quantizeGrid( grid = sphi, step = (1.0/256.0) )
			sphi.save( '%sxl_%s_%04d.uni' % (outPrefix, filenameSuffixXl, thi) )

	if(tlo >= timesteps):
		exit(1);

	# mesh gen, only for UI
	if 1 and showGui and (dim==3):
		if doHires>0:
			sphi.createMesh(smesh)
		else:
			phi.createMesh(mesh)

	lastFrame = s.frame;
	
	#timings.display()
	#s.printMemInfo()
	s.step()

	if doWritePngs and showGui and (GUI):
		gui.screenshot( '%s_%s_%04d.png' % (outPrefix, filenameSuffix, s.frame) );


