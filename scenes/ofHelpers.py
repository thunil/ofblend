#
# Helper functions for optical flow
#
import sys
import os
import shutil
from   manta import *

def vec4ToVec3(vec4d):
	ret   = vec3(0.)
	ret.x = vec4d.x;
	ret.y = vec4d.y;
	ret.z = vec4d.z;
	return ret;

def vec4ToVec3i(vec4d):
	ret   = vec3(0)
	ret.x = int(vec4d.x);
	ret.y = int(vec4d.y);
	ret.z = int(vec4d.z);
	return ret;

def vec4intString(vec4d):
	return ("[%d,%d,%d,%d]" % ( int(vec4d.x), int(vec4d.y), int(vec4d.z), int(vec4d.t) ) );

def initRes(res, yFactor, fourdFactor, loadTimeScale):
	ret = vec4(0.)
	ret.x   = (res * 1       );
	ret.y   = (res * yFactor );
	ret.z   = (res * 1       );
	ret.t   = (res * fourdFactor * loadTimeScale);
	ret3d   = vec4ToVec3i(ret);
	return (ret3d,ret);

# ======================================================================================================================
# read parameters

def getParam(name, default, paramUsed):
	while( len(paramUsed)<len(sys.argv) ):
		paramUsed.append(0); 
	for iter in range(1, len(sys.argv)):
		if(sys.argv[iter].lower() == name) and (iter+1<len(paramUsed)):
			paramUsed[iter] = paramUsed[iter+1] = 1;
			return sys.argv[iter+1];
	return default;

def checkUnusedParam(paramUsed):
	err = False;
	for iter in range(1, len(sys.argv)):
		if(paramUsed[iter]==0):
			print("Error: param %d '%s' not used!" % (iter,sys.argv[iter]) ); 
			err = True;
	if err:
		exit(1);

# ======================================================================================================================
# calculate offset for phi load
# used for regular & sliced versions
def initOffsets(gs4d, gs4dDefo, autoBorder, repeatStartDist, loadTimeScale):

	loadOffset   = vec4(0.); # res * autoBorder);
	loadOffset.x = gs4d.x * autoBorder;
	loadOffset.y = gs4d.y * autoBorder;
	loadOffset.z = gs4d.z * autoBorder;
	loadOffset.t = gs4d.t * autoBorder;
	loadOffset.t = loadOffset.t + repeatStartDist 

	loadScale    = vec4( (1. - 2.0*autoBorder) );
	loadScale.t  = loadScale.t * loadTimeScale;

	# and again for deformation, plus re-scaling of magnitude (defoFactor)

	defoOffset = vec4(0.); # off

	defoScale  = vec4(1., 1., 1., loadTimeScale); 

	defoFactor = vec4(0.)
	defoFactor.x   = gs4d.x / gs4dDefo.x;
	defoFactor.y   = gs4d.y / gs4dDefo.y;
	defoFactor.z   = gs4d.z / gs4dDefo.z;
	defoFactor.t   = gs4d.t / gs4dDefo.t;

	return (loadOffset,loadScale,defoOffset,defoScale,defoFactor);



# ======================================================================================================================
# loading functions

# load and adapt 4d vel/defo data
def loadVelRescaled(velgrid, prefix, defoSolver, defoOffset, defoScale, defoFactor):
	# rescale load (only vel, readvect phi)
	# size for surface data & defo could be different
	ipol_vel   = defoSolver.create(Grid4Vec4)

	ipol_vel.load( "%s_vel.uni" % prefix );
	interpolateGrid4dVec(velgrid,   ipol_vel, offset=defoOffset , scale=defoScale );
	del(ipol_vel);

	velgrid.multConst(defoFactor); 



# ======================================================================================================================
# outputs


def writeResults(name, t, meshCnt, grid, meshout, gui,  writePngs  , triMesh    , projPpm    , writeUni   , sdfSmooth,sdfSmoothBord  , meshPrefix ):
	#print( "Output filename: %s_%04d" % (name,t) );
	if sdfSmooth>0:
		print("Output sdf smoothing active %d, bord %d " % (sdfSmooth, sdfSmoothBord) );
		simpleBlurSpecial( grid, 1, 0.   , sdfSmoothBord-1 ); 
		simpleBlurSpecial( grid, 1,-999. , sdfSmoothBord-1 ); # all

	if writePngs or triMesh or gui: # show mesh...
		grid.createMesh(meshout, notiming=True ); 

	if writePngs:
		gui.screenshot( '%s_%04d.png' % (name,t), notiming=True  );

	if triMesh: 
		meshout.save( "%s_%04d.bobj.gz" % (meshPrefix, meshCnt), notiming=True  );

	if projPpm:
		#projectPpmFull( grid, '%s_%04d.ppm' % (name,t) , 0, 10.1);  # smoke
		projectPpmFull( grid, '%s_%04d.ppm' % (name,t) , 1, 1.); # surface test

	if writeUni:
		grid.save( "%s_%04d.uni" % (name, t), notiming=True  );


# ======================================================================================================================
# helper functions for UVW conversions

def world2uvw(p, r1,r2,r3):
	#print("w2uvw " + str(p) + " tri " + str(r1)+" "+ str(r2)+" "+ str(r3)); 

	# calc ((r1-rr3) (r2-r3))^-1 (p-r3)
	r13 = r1-r3;
	r23 = r2-r3;
	ma = r13.x; mc = r13.y;
	mb = r23.x; md = r23.y;
	if abs(ma*md-mb*mc)<1e-05:
		print(("Invalid coords? %f ") % (ma*md-mb*mc)); exit(1);
	mdet = 1./(ma*md-mb*mc);
	ina = mdet * md;
	inb = mdet * -mb;
	inc = mdet * -mc;
	ind = mdet * ma;
	#print("InvM %f %f   %f  %f " %(ina,inb,inc,ind) );
	pr3 = p-r3;
	triu = ina*pr3.x + inb*pr3.y;
	triv = inc*pr3.x + ind*pr3.y;
	triw = (1.-triu-triv) 
	#print("UVW %f %f %f" %(triu,triv,triw) );
	eps = 1e-03; # accept round off...
	if triu<(0.-eps) or triv<(0.-eps) or triw<(0.-eps) or triu>(1.+eps) or triv>(1.+eps) or triw>(1.+eps):
		print("w2uvw " + str(p) + " tri " + str(r1)+" "+ str(r2)+" "+ str(r3) + " , uvw %f,%f,%f " % (triu,triv,triw) ); 
		print("Invalid UVW!"); exit(1);
	# make sure...
	clamp01(triu);
	clamp01(triv);
	triw = (1.-triu-triv) 

	#print("Pos " + str(triu*r1+ triv*r2+ triw*r3) );  # double check
	uvw = vec3(triu,triv,triw); # "return" value
	return uvw;

def clamp01(v):
	if v<0.: return 0.;
	if v>1.: return 1.;
	return v;

def weights2dim(uvw):
	u = uvw.x;
	v = uvw.y;
	w = uvw.z;
	w0    = w1    = w2    = w_uv  = w_vw  = w_uw  = w_uvw = 0.;
	
	p=vec3(u,v,w);
	if u>=0.5:
		r1=vec3(1.0,0.0,0.0);
		r2=vec3(0.5,0.5,0.0);
		r3=vec3(0.5,0.0,0.5);
		uvw = world2uvw(p, r1,r2,r3);
		w0    = uvw.x;
		w_uv  = uvw.y;
		w_uw  = uvw.z;
	elif v>=0.5:
		r1=vec3(0.0,1.0,0.0);
		r2=vec3(0.5,0.5,0.0);
		r3=vec3(0.0,0.5,0.5);
		uvw = world2uvw(p, r1,r2,r3);
		w1    = uvw.x;
		w_uv  = uvw.y;
		w_vw  = uvw.z;
	elif w>=0.5:
		r1=vec3(0.0,0.0,1.0);
		r2=vec3(0.5,0.0,0.5);
		r3=vec3(0.0,0.5,0.5);
		uvw = world2uvw(p, r1,r2,r3);
		w2    = uvw.x;
		w_uw  = uvw.y;
		w_vw  = uvw.z;
	else:
		# inner 3 triangles
		if u >= v and u >= w:
			#print "tri D";
			r1=vec3(0.5,0.5,0.0);
			r2=vec3(0.5,0.0,0.5);
			r3=vec3(1./3.,1./3.,1./3.);
			uvw = world2uvw(p, r1,r2,r3);
			w_uv  = uvw.x;
			w_uw  = uvw.y;
			w_uvw = uvw.z;
		elif v >= u and v >= w:
			#print "tri E";
			r1=vec3(0.5,0.5,0.0);
			r2=vec3(0.0,0.5,0.5);
			r3=vec3(1./3.,1./3.,1./3.);
			uvw = world2uvw(p, r1,r2,r3);
			w_uv  = uvw.x;
			w_vw  = uvw.y;
			w_uvw = uvw.z;
		elif w >= u and w >= v:
			#print "tri F";
			r1=vec3(0.5,0.0,0.5);
			r2=vec3(0.0,0.5,0.5);
			r3=vec3(1./3.,1./3.,1./3.);
			uvw = world2uvw(p, r1,r2,r3);
			w_uw  = uvw.x;
			w_vw  = uvw.y;
			w_uvw = uvw.z;
		else:
			print "This should never happen...???";

	return (w0,w1,w2,w_uv,w_vw,w_uw,w_uvw);


