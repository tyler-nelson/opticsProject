from numpy import pi
import numpy as np
import mpmath

'''
The ray class is esentially a container for a variety of attributes associated
with a light ray

We have broken down which category each property best fits into
this typically will indicate which part of the program will use the field
although there may be overlap with other categories.

Summary of fields:
The light ray is defined in a mathematical sense as a normal ray which has an
end point (source location) and an angle that provides the ray direction
(denoted angle) measured CCW from the horizontal.

Due to our method of determining where a ray will strike a surface.
'''


class Ray:
    rayCount = 0

    def __init__(self, x, y, angle, wavelength=0, identity=42, intensity=1,
                 objIndex=-1, sPolarization=0.5, pPolarization=0.5):

        # vector quantities
        self.sourceLocation = (x, y)
        self.angle = angle
        self.unitVector = [np.cos(angle), np.sin(angle)]

        # line quantities
        self.slope = (self.unitVector[1]/self.unitVector[0])
        self.yIntercept = y - x*self.slope  # y = mx + b --> b = y - mx
        self.xIntercept = - self.yIntercept/self.slope  # x = -b/m

        # plot quantities
        # keep track of the prevous interface points
        self.locationHistory = [self.sourceLocation]

        # light quantities
        self.wavelength = wavelength
        self.intensity = intensity
        self.intensityList = [intensity]
        self.parallel = pPolarization
        self.perpendicular = sPolarization

        # program quantities
        self.objIndex = -1
        # self.identity = identity #this can be set up in the ray generator
        self.identity = Ray.rayCount
        Ray.rayCount += 1
        self.reflected = 0
        self.worthUsing = 1

    def update_location_history(self, x, y, angle, objIndex):

        # update vector quantities:
        self.sourceLocation = (x, y)
        self.angle = angle

        ang = mpmath.mpmathify(angle)
        x0, y0 = mpmath.cos(ang), mpmath.sin(ang)
        self.unitVector = [float(x0), float(y0)]

        self.slope = float(mpmath.fdiv(y0, x0))
        self.yIntercept = y - x*self.slope  # y = mx + b --> b = y - mx
        self.xIntercept = - self.yIntercept/self.slope  # x = -b/m

        # update plot quantities
        self.locationHistory.append((x, y))

        # update program quantity
        self.objIndex = objIndex

    def rPerp(angleI, angleT, nI, nT):
        return (nI*np.cos(angleI)-nT*np.cos(angleT)) /\
                (nI*np.cos(angleI)+nT*np.cos(angleT))

    def rPara(angleI, angleT, nI, nT):
        return (nT*np.cos(angleI)-nI*np.cos(angleT)) /\
                (nT*np.cos(angleI)+nI*np.cos(angleT))

    def tPerp(angleI, angleT, nI, nT):
        return 2*nI*np.cos(angleI)/(nI*np.cos(angleI)+nT*np.cos(angleT))

    def tPara(angleI, angleT, nI, nT):
        return 2*nI*np.cos(angleI)/(nT*np.cos(angleI)+nI*np.cos(angleT))

    # should be irradience
    def update_intensity_and_polarization(self, theta1, theta2, n_1, n_2):
        t_perp_sqr = Ray.tPerp(theta1, theta2, n_1, n_2)**2
        t_para_sqr = Ray.tPara(theta1, theta2, n_1, n_2)**2
        self.perpendiular = t_perp_sqr/(t_perp_sqr+t_para_sqr)
        self.parallel = 1 - self.perpendicular

        if(theta1 == 0):
            self.intensity *= (self.perpendicular + self.parallel)*n_2 /\
                            n_1*(2*n_1/(n_1+n_2))**2
        else:
            self.intensity *= n_2*np.cos(theta2)/(n_1*np.cos(theta1)) * \
                                         (self.perpendicular*t_perp_sqr +
                                          self.parallel*t_para_sqr)
        self.intensityList.append(self.intensity)

    def printString(self):
        return "{0}\n{1}\n".format(self.slope, self.xIntercept)


class RayGenerator:
    """
    This class acts as a container for a ray list to be given to a driver class
    along some optical elements for those rays to interact with.

    There are three predefined methods of construction that provide frequently
    seen ray sources. Those are explained in turn.
    """
    def __init__(self, rayList):
        self.theRayList = rayList

    @classmethod
    def point_source(cls, x0, y0, numberOfRays, startAngle, stopAngle,
                     startLambda=1, endLambda=1, sPolarization=0.5,
                     pPolarization=0.5):
        """
        Overview :
        This sets up a diverging source eminating from a point (x0,y0) sweeping
        from start angle to stop angle in evenly spaced intervals based on how
        many rays are given.

        Inputs:
        x0,y0 -                 (x,y) coordinates of the source
        startAngle, stopAngle - angles to begin and end the generation process
                                (in radians)
        startLambda, endLambda - the wavelengths to begin and end at for the
                                 rays (they're divided in the same way the
                                 angles are)
        sPolarization, pPolarization - the fraction of perpendicular or
                                       parallel light, these must add to 1

        Outputs:
        a ray list which fits these criteria
        """
        if startAngle == 0 and stopAngle == 2*np.pi:
            theAngleList = np.linspace(startAngle, stopAngle, numberOfRays+1)
        else:
            theAngleList = np.linspace(startAngle, stopAngle, numberOfRays)

        return cls([Ray(x0, y0, (angle + 2*pi if angle < 0 else angle), 0,
                        sPolarization=sPolarization,
                        pPolarization=pPolarization)
                    for angle in theAngleList])

    @classmethod
    def beam_source(cls, x0, y0, x1, y1, numberOfRays, startLambda=1,
                    endLambda=1, sPolarization=0.5, pPolarization=0.5):
        """
        Overview - this creates a line segment determined by points (x0,y0) and
        (x1,y1) that the rays will emerge perpendicularly from at evenly space
        intervals based on the number of rays. If only one ray is given, the
        ray will emminate from the (x0,y0) point. Since this simulator was
        originally based around having light move from top to bottom, hence the
        light always is oriented to move upward.

        Inputs:
        x0,y0 - start point of the line segment
        x1,y1 - end point of the line segment
        numberOfRays - self evident
        OTHER PARAMETERS - explained above

        Outputs:
        a list of rays which satify these conditions

        """
        if x1 - x0 == 0 and y1 - y0 == 0:
            raise(Exception("beam width must be non-zero"))
        rayAngle = np.arctan2((x1-x0), -(y1-y0))
        xList = np.linspace(x0, x1, numberOfRays)
        yList = np.linspace(y0, y1, numberOfRays)

        return cls([Ray(p0, p1, rayAngle, 0, identity=index,
                    sPolarization=sPolarization,
                    pPolarization=pPolarization)
                    for index, (p0, p1) in enumerate(zip(xList, yList))])

    @classmethod
    def converging_source(cls, x0, y0, x1, y1, numberOfRays, startAngle,
                          stopAngle, startLambda=1, endLambda=1,
                          sPolarization=0.5, pPolarization=0.5):
        """
        Overview - this constructor seeks to combine the beam source with the
        point source to create a converging source from a line.
        Inputs:
        x0,y0      - define the start of the line
        x1,y1      - define the end of the line
        startAngle - angle of the leftmost ray
        stopAngle  - angle of the rightmost ray (with rays inbetween having
                     angles that evenly divide the range)
        OTHER PARAMETERS - see point source constructor

        Outputs:
        a ray list that meets these specifications

        """
        if x1-x0 == 0 and y1-y0 == 0:
            raise(Exception("converging beam width must be non-zero"))
        xList = np.linspace(x0, x1, numberOfRays)
        yList = np.linspace(y0, y1, numberOfRays)
        angleList = np.linspace(startAngle, stopAngle, numberOfRays)

        return cls([Ray(p0, p1, angle, 0, sPolarization=sPolarization,
                        pPolarization=pPolarization)
                    for p0, p1, angle in zip(xList, yList, angleList)])
