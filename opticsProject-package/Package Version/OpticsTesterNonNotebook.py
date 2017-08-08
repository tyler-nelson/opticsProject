from Driver.Driver import driver
from Lens.Lens import Lens
from Mirror.Mirror import Mirror
from Ray.Ray import Ray, RayGenerator
from functools import partial
from Components.Presets import circle, line

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import warnings

plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 20})
mpl.rcParams['figure.figsize'] = 8, 8
warnings.filterwarnings("ignore")


# fix problem with updating objIndex
# bug with plano concave r1 = 0, r2 < 0
# bug at -1.5 for l2, it's because our rays were generated in the l2 lens
# should improve efficiency of ray propagating algorithm by using a flag to
# indicate whether a ray is worth working on


Ray.rayCount = 0

myRayList = RayGenerator.beam_source(4, -3, 6, -3, 11)

f1 = lambda x: -((x-5)**2+1)**(1/2)+8
f2 = lambda x: 1

l = [Lens.custom(surfaceList=[[[3], [f2(3), f1(3)]],
                              [[3, 7], [f2, f1]], [[7], [f2(7), f1(7)]]])]
m = []
theDriver = driver(lambda x: 1, myRayList.theRayList,
                   listOfLenses=l, listOfMirrors=m)
theDriver.rayFunc(threshold=0)
fig = plt.figure()

theDriver.plotter(fig)
plt.ylim(-3, 14)
plt.xlim(1, 9)
plt.legend()
plt.tight_layout()
plt.show()

# -----------------------------------------------------------------------------

# fix problem with updating objIndex
# bug with plano concave r1 = 0, r2 < 0
# bug at -1.5 for l2, it's because our rays were generated in the l2 lens
# should improve efficiency of ray propagating algorithm by using a flag to
# indicate whether a ray is worth working on

# TODO: BROKEN

Ray.rayCount = 0

myRayList = RayGenerator.beam_source(4, -3, 6, -3, 11)
myLens = Lens.gaussianOptics(r1=5, r2=-5, diameter=3, thickness=0.5, xOffSet=0,
                             yOffSet=-1, nFunction=lambda x: 1.5)
myLens2 = Lens.gaussianOptics(r1=-3, r2=3, diameter=3, thickness=0.5,
                              xOffSet=2, yOffSet=-4.75, nFunction=lambda x: 1.5)
myMirror = Mirror.custom(surfaceList=[[[0, 10], [lambda x: -(x-5)**2+8]]])
l = [myLens]#, myLens2]
m = [myMirror]

theDriver = driver(lambda x: 1, myRayList.theRayList,
                   listOfLenses=l, listOfMirrors=m)
theDriver.rayFunc(threshold=1e-1)
print("driver ray")
fig = plt.figure()

theDriver.plotter(fig)
plt.ylim(-3, 14)
plt.xlim(1, 9)
plt.legend()
plt.tight_layout()
plt.show()

# END TODO

# -----------------------------------------------------------------------------

Ray.rayCount = 0

myRayList = RayGenerator.point_source(5, 8+2**0.5, 11, 75*pi/180+pi,
                                      105*pi/180+pi)

f1 = partial(circle, False, 8, 1, 5)
f2 = lambda x: 1

l = [Lens.custom(surfaceList=[[[1], [f2(1), f1(1)]], [[1, 9], [f2, f1]],
                              [[9], [f2(9), f1(9)]]],
                 nFunction=lambda x: 2**0.5)]
m = []
theDriver = driver(lambda x: 1, myRayList.theRayList,
                   listOfLenses=l, listOfMirrors=m)
fig = plt.figure()

theDriver.plotter(fig)
plt.ylim(-3, 20)
plt.xlim(1, 9)
plt.legend()
plt.tight_layout()
plt.show()
'''
# -----------------------------------------------------------------------------
'''
Ray.rayCount = 0

myRayList = RayGenerator.point_source(5, 3, 11, 75*pi/180+pi, 105*pi/180+pi)

f1 = lambda x: -(4/5*(x-5)**2+4)**(1/2)+0
f2 = lambda x: (4/5*(x-5)**2+4)**(1/2)-9

l = [Lens.custom(surfaceList=[[[1], [f2(1), f1(1)]], [[1, 9], [f2, f1]],
                              [[9], [f2(9), f1(9)]]], nFunction=lambda x: 1.5)]
m = []
theDriver = driver(lambda x: 1, myRayList.theRayList,
                   listOfLenses=l, listOfMirrors=m)
theDriver.rayFunc(threshold=0)
fig = plt.figure()

theDriver.plotter(fig)
plt.ylim(-13, 4)
plt.xlim(0, 10)
plt.legend()
plt.tight_layout()
plt.show()


Ray.rayCount = 0

myRayList = RayGenerator.point_source(2+5, 0, 11, -30*pi/180+pi,30*pi/180+pi)

f1 = lambda x: -(-5/9*(x-5)**2+5)**(1/2)
f2 = lambda x: (-5/9*(x-5)**2+5)**(1/2)

l = [Lens.custom(surfaceList=[[[2], [f2(2), f1(2)]], [[2, 8], [f2, f1]],
                              [[8], [f2(8), f1(8)]]], nFunction=lambda x: 1.5)]
theDriver = driver(lambda x: 1, myRayList.theRayList,
                   listOfLenses=l, listOfMirrors=[])
theDriver.rayFunc(threshold=1e-3, segments=32)
fig = plt.figure()

theDriver.plotter(fig)
plt.ylim(-5, 5)
plt.xlim(0, 10)
plt.legend()
plt.tight_layout()
plt.show()

# ------------------------------------------------------------------------------

Ray.rayCount = 0

myRayList = RayGenerator.point_source(5, -2, 11, -30*pi/180+pi/2,
                                      30*pi/180+pi/2)

f1 = lambda x: -(-9/5*(x-5)**2+9)**(1/2)
f2 = lambda x: (-9/5*(x-5)**2+9)**(1/2)

l = [Lens.custom(surfaceList=[[[5-5**(1/2)], [f2(5-5**(1/2)), f1(5-5**(1/2))]],
                              [[5-5**(1/2), 5+5**(1/2)], [f2, f1]],
                              [[5+5**(1/2)], [f2(5+5**(1/2)), f1(5+5**(1/2))]]],
                 nFunction=lambda x: 1.5)]
theDriver = driver(lambda x: 1, myRayList.theRayList,
                   listOfLenses=l, listOfMirrors=[])
theDriver.rayFunc(threshold=1e-3, segments=32)
fig = plt.figure()

theDriver.plotter(fig)
plt.ylim(-5, 5)
plt.xlim(0, 10)
plt.legend()
plt.tight_layout()
plt.show()

# ------------------------------------------------------------------------------

Ray.rayCount = 0

myRayList = RayGenerator.point_source(2+5, 0, 11, -30*pi/180+pi, 30*pi/180+pi)

f1 = lambda x: (-(x-5-xOffSet)**2 + 9)**(1/2) + 1 +yOffSet
f2 = lambda x: -(-(x-5-xOffSet)**2 + 9)**(1/2) + 8 + yOffSet
l = [myLens]
m = []
theDriver = driver(lambda x: 1, myRayList.theRayList,
                   listOfLenses=l, listOfMirrors=m)

fig = plt.figure()

theDriver.plotterJustLens(fig)
plt.ylim(-5, 5)
plt.xlim(0, 10)
plt.legend()
plt.tight_layout()
plt.show()

myRayList = RayGenerator.beam_source(3, 0, 7, 0, 10)
#myLens = Lens.gaussianOptics(r1=5, r2=-5, diameter=3, thickness=0.5,
#                             xOffSet=0, yOffSet=0, nFunction=lambda x: 1.5)
myLens = Lens.predefined("meniscus concave")
theDriver = driver(lambda x: 1, myRayList.theRayList,
                   listOfLenses=[myLens], listOfMirrors=[])

theDriver.rayFunc(threshold=1e-3, segments=8)
fig = plt.figure()

theDriver.plotter(fig)
plt.ylim(-2, 8)
plt.xlim(0, 10)
plt.legend()
plt.tight_layout()
plt.show()
