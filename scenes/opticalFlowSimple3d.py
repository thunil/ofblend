#
# simple 3d optical flow example
#
from manta import *

# solver params
res = 60
#res = 100 
dim = 3
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)


flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)
i0       = s.create(LevelsetGrid)
i1       = s.create(LevelsetGrid) 
vel      = s.create(MACGrid)

# show surface?
mesh0     = s.create(Mesh)
mesh1     = s.create(Mesh)
meshOut   = s.create(Mesh)

# setup source & target 
box1 = s.create(Box, p0=gs*vec3(0.3), p1=gs*vec3(0.5)) 
phi.setConst( 999. )
phi.join( box1.computeLevelset() ) 
i0.copyFrom(phi)

box2 = s.create(Box, p0=gs*vec3(0.4), p1=gs*vec3(0.7)) 
phi.setConst( 999. )
phi.join( box2.computeLevelset() ) 
i1.copyFrom(phi)

# prepare, control SDF value range
dist = int(res/4.)

extrapolateLsSimple(phi=i0, distance=dist)
extrapolateLsSimple(phi=i0, distance=dist, inside=True)

extrapolateLsSimple(phi=i1, distance=dist)
extrapolateLsSimple(phi=i1, distance=dist, inside=True)

i0.multConst(0.0005*float(dist))
i1.multConst(0.0005*float(dist))

# go...

if 1 and (GUI):	
	gui = Gui()
	gui.show()
	gui.pause()

while s.frame < 9999:

	# compute deformation
	opticalFlowMultiscale3d( i0=i0, i1=i1, vel=vel, wSmooth=0.5, wEnergy=0.0001, multiStep=4 )

	# apply
	phi.copyFrom(i0)
	advectSemiLagrangeCfl(flags=flags, vel=vel, grid=phi, order=1, velFactor=(float)(1.0), cfl=999) # cfl off

	if 1: # show meshes
		i0.createMesh(mesh0)
		i1.createMesh(mesh1)
		phi.createMesh(meshOut)

	s.step()

