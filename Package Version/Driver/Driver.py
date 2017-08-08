import mpmath
import numpy as np
from numpy import pi
import scipy.optimize as opt

from Ray.Ray import Ray


class driver:

    def __init__(self, n, listOfRays, listOfLenses=[], listOfMirrors=[]):
        self.nMedium = n
        self.listOfRays = listOfRays
        self.listOfMirrors = listOfMirrors
        self.listOfLenses = listOfLenses
        self.surfaceList = []
        for i in listOfLenses:
            self.surfaceList += i.surfaceList
        for i in listOfMirrors:
            self.surfaceList += i.surfaceList
        for i in listOfRays:
            i.objIndex = self.locatingIndexFunction(i.sourceLocation)

    def locatingIndexFunction(self, sourceLocation):
        """
        Overview: this function seeks to determine which index this location
                  belongs to

        inputs:  sourceLocation, (implicit list of lenses, list of mirrors)
        outputs: the label of the lenses in the map that contains this location
                 or -1 to indicate that it's in the surrounding medium

        Procedure:
        1. Iterate through the map of lenses and labels, choose a lens l_i
        2. Iterate through the list of surfaces associated with this lens,
           choose s_j
        ### ASSUMING NOT VERTICAL LINE ###
        3. if check s_j's domain bounds the x coordinate, if not proceed to
           next surface
        4. check if pairs of surfaces defined over the domain bound the y
           position by evaluating them
        ### VERTICAL LINE ####
        3. check if the line is equal within tolerance to the x coordinate
        4. check if the y position is bounded by the pairs of y coordinates
           that decribe the vertical line
        #### END OF CASES ####
        5. return if successful otherwise continue
        """
        x,y = sourceLocation
        listOfLenses = self.listOfLenses
        for i in range(len(self.listOfLenses)):
            for orderedPair in listOfLenses[i].surfaceList:
                domain, surfaces = orderedPair
                if(len(domain) == 1):
                    if(len(surfaces) == 2):
                        if(x > domain[0]-0.001 and x < domain[0]+0.001 and y > surfaces[0] and y < surfaces[1]):
                            return i
                    else:
                        print("this is the pathelogical broken lens case")
                        raise(Exception("I haven't implemented this yet"))
                else:
                    if(x > domain[0] and x < domain[1]):#this could still be pathelogical
                        for m in range(0, len(surfaces), 2):
                            f1, f2 = surfaces[m], surfaces[m+1] #I assume that odd numbers of surfaces shouldn't be possible
                            lB, uB = sorted([complex(f1(x)).real,complex(f2(x)).real])
                            if(lB < y and uB > y ):
                                return i
        return -1

    def rayFunc(self,threshold=0,segments=8,extraRays=1,passLimit=12,digits=100):
        """
        NOTE:
        The beginning of this function is comprised of definitions of helper
        functions that the rayFunc needs to work
        each of these functions will be described under their definition.

        Inputs: threshold(Optional), segments(optional) (an instance of the
        class must call this method)
        Outputs: none

        ***
        These input/output lists are not really accurate as to what the
        function is doing since an instance of the object driver is modifying
        its own data members they don't need to be parameters for the function.
        A more accurate description is, the driver is modifying the list of
        rays that was provided on construction of the instance by having those
        rays interact with the various optical elements provided in the list of
        lenses and list of mirrors given to the class initially. Hence the
        input and output are a list of rays.
        ***

        Description of procedure for rayFunc:
        1. Iterate through the list of rays, suppose we are on the ith step and
           have a ray r_i that we want to simulate
        2. Iterate through all surfaces provided to determine which the ray
           will impinge upon(if any).
           These surfaces are those which have a point that lies on the
           extension of the ray. The ray will interact with the surface that
           minimizes the Opitcal Path Length(OPL). Since we only work with
           homogenous media (that is, n is not a function of spatial
           coordinates within an element), the shortest path length will be a
           straight line.
        3. Determine if the next medium is a lens, mirror, or the background by
           checking the position of the interaction point and searching for
           which element contains this point
        4. determine angle of intersection by using
           arccosine(r_i.unitVector dotproduct surfaceNormalVector)
        #### MIRROR ONLY ####
        5. Use law of reflection and a rotation matrix to specify the
           omponents of the reflected ray
        ->continue reading at the end of lens portion
        #### LENS ONLY ####
        5. Determine critical angle if n1 > n2
        6.a. If not total internal reflection:
            i.   Use Snell's Law to calculate the angle relative to the normal
                 of the refracted ray
            ii.  Use a rotation matrix to determine the components of the
                 refracted ray
            iii. Determine the intensity of the transmited (refracted ray) and,
                 by conservation of energy, the intensity of
                 the internally reflected ray as well
            iv. Calculate the s and p components of the reflected and
                transmitted ray
        6.b.  Use law of reflection and a rotation matrix to specify the
              components of the reflected ray
        #### END OF CASES ####
        7. update the ray object and record the normal vector

        """
        def normalVectorFunc(function, value):
            """
            This function returns the normal vector of a function at a
            specified point using a derivative based approach

            inputs:  function, location
            outputs: an ordered pair corresponding to the xy components of the
                     normal vector

            NOTE: The direction of the normal vector is ambiguous for the ray
            which requires a check later to ensure the correct orientation
            """

            # Value is the x value where the given ray intersects the lens
            # surface
            try:
                fprime = float(mpmath.diff(function, value))
                return [-fprime/(1+fprime**2)**(0.5), 1/(1+fprime**2)**(0.5)]
            except TypeError:
                raise(Exception("Unexpectedly recieved complex number."))

        # expected inputs for line: [slope, x intercept], surface is a
        # function, value is a reasonable estimate for where the zero occurs
        def rayIntersectSurface(line, surface, value):
            """
            Given two functions, f1 and f2, f1 intersected f2 when
            f1(x) - f2(x) = 0
            This function is a wrapper for the fsolve routine in the scipy
            optimize suite

            inputs:
            line : a list containing the slope and x intersecpt of our ray
            surface : the surface that the ray might hit
            value: a guess to help fsolve find a root

            outputs:
            according to fsolve: a list of intersection points or the last
            value used if the search was unsuccessful
            """
            return [float(mpmath.re(mpmath.findroot(
                lambda x: surface(x)-x*line[0]+line[0]*line[1],
                value, tol=1e-10, verify=False)))]

        def refractedRayFunction(normalVector, angle, isNegative):
            """
            A function that calculates the components of the possible reflected
            or refracted ray depending on whether isNegative is true or false.
            If false, i.e. equal to 0 then the reflected ray is calculated.
            Otherwise the refracted ray is calculated. Either of these are
            achieved by rotating the normal vector.

            NOTE:
            The rotation matrix method is ambigious about whether the angle
            that we are rotating about is positive or negative.
            Without a sure way to calculate the correct one only, we calculate
            both possibilities and take the one that has the greater
            dot product to the ray's unit vector

            input:
            normalVector, angle to rotate by, isNegative - a flag to indicate
            whether we are using the normal vector or
            minus the normal vector

            output:
            two lists containing the components of the rotated rays
            """
            xN, yN = normalVector
            if(isNegative):
                xN, yN = -xN, -yN
            refractedRay1 = [xN*np.cos(-theta2) - yN*np.sin(-theta2),
                             xN*np.sin(-theta2) + yN*np.cos(-theta2)]
            refractedRay2 = [xN*np.cos(theta2) - yN*np.sin(theta2),
                             xN*np.sin(theta2) + yN*np.cos(theta2)]
            return refractedRay1, refractedRay2
        # lets make it more modular :D

        def surfaceFinderPart(ray, listOfElements, segments):
            """
            Inputs:
            ray : a ray that we are seeking to have intersect a surface
            listOfElements: a list containing either lenses or mirrors that we
            are going to search their list of surfaces

            Outputs:
            xMin, yMin represent the (x,y) of the nearest location of
            intersection. minDistSurface is the function that minimized the OPL
            for this set of optical elements. minDistance is the length of the
            OPL.

            Procedure: (for not vertical lines)
            1. iterate through the list of elements, select one say e_i
            2. iterate through its list of surfaces, s_j
            ### FOR NOT VERTICAL LINES ###
            3. compare the domain of the surface to the point's x location
            4.a. if the x location is not bounded by the domain iterate to next
                 surface
            4.b  if the x location is bounded, break the domain into 8 (or more
                 if desired) evenly spaced intervals
            5. for each interval check if the rayIntersectSurface function
               returns a valid point (that is a point which is along
               the extension of the ray)
            ### FOR VERTICAL LINES ### !If there are multiple vertical lines
            which share a x location, this should iterate through each line
            segment
            3. check if x location of the point is within tolerance equal to
               the vertical line, if not iterate to next lens
            4. check if y location is bounded by the line segment, if not
               iterate to next lens
            ### END OF CASES ###
            6. if a valid point is found, compare its OPL to the best OPL so
               far, if better, update the quantities
            """
            x0, y0 = ray.sourceLocation
            minDistance = 1e10
            minDistSurface = -1
            xMin, yMin = -1, -1
            for element in listOfElements:
                listOfSurfaces = element.surfaceList
                for s in range(len(listOfSurfaces)):
                    domain, surfaces = listOfSurfaces[s]
                    # should add the shifted case, where there are two line
                    # segments
                    if(len(domain) == 1):  # this is a vertical line
                        pointOfIntersection = ray.slope*domain[0] +\
                                              ray.yIntercept
                        realPart = [mpmath.re(i) for i in surfaces]
                        surfacesSorted = sorted(realPart)
                        if(pointOfIntersection >= surfacesSorted[0]
                           and pointOfIntersection <= surfacesSorted[1]):
                            distance = ((domain[0]-x0)**2 +
                                        (pointOfIntersection-y0)**2)**(1/2)
                            checkAngle = mpmath.atan2(pointOfIntersection -
                                                      y0, domain[0]-x0)
                            rayAngle = ray.angle
                            if(checkAngle < 0):
                                checkAngle += 2*pi
                            if(checkAngle > rayAngle+0.001
                               or checkAngle < rayAngle-0.001):
                                # This ray intersection is clearly bad so
                                # disregard it
                                continue
                            if(distance < minDistance and distance > 0.001):
                                minDistance = distance
                                xMin, yMin = domain, pointOfIntersection
                                # I use this as the label to differentiate
                                # between the vertical line case and the rest
                                minDistSurface = [domain[0], surfacesSorted[:2]]
                        elif(len(surfaces) == 4):
                            print("update this part with the pieces from above")
                            if(pointOfIntersection >= surfaces[2]
                               and pointOfIntersection <= surfaces[3]):
                                distance = ((domain[0]-x0)**2 +
                                            (pointOfIntersection-y0)**2)**(1/2)

                                if(distance < minDistance
                                   and distance > 0.001):
                                    minDistance = distance
                                    xMin, yMin = domain, pointOfIntersection
                                    minDistSurface = [domain[0], surfaces[2:4]]
                    # should probably include error checking to ensure that the
                    # length of any of these domain lists is either 1 or 2
                    else:  # we're not working with a vertical line
                        # given a domain (a,b), my strategy is to use fsolve
                        # with an intial guess of (a+b)/2, if the ray happens
                        # to already be on the surface, we refine the search by
                        # splitting the surface into eighths
                        a, b = domain

                        # somewhat arbitarily chose to split the domain into
                        # eighths
                        initial_guess = [i*(a+b)/segments
                                         for i in range(segments)]
                        for guess in initial_guess:
                            # we have to loop over all the surfaces defined
                            # on this domain
                            for i in range(len(surfaces)):
                                surface_i = surfaces[i]
                                listOfPointsOfIntersection = rayIntersectSurface([ray.slope, ray.xIntercept], surface_i, guess)
                                for pointOfIntersection in listOfPointsOfIntersection:

                                    # we've got a valid intersection
                                    if(pointOfIntersection > a
                                       and pointOfIntersection < b):
                                        y1 = surface_i(pointOfIntersection)
                                        distance = ((pointOfIntersection-x0)**2
                                                    + (y1-y0)**2)**(1/2)
                                        checkAngle = mpmath.atan2(y1-y0,
                                                                  pointOfIntersection-x0)
                                        rayAngle = ray.angle
                                        if(checkAngle < 0):
                                            checkAngle += 2*pi

                                        if(checkAngle > rayAngle+0.001
                                           or checkAngle < rayAngle-0.001):
                                            # This ray intersection is clearly
                                            # bad so disregard it
                                            continue
                                        if(distance < minDistance
                                           and distance > 0.001):
                                            # Put in a min distance, this might
                                            # fix the ray not moving problem
                                            minDistance = distance
                                            xMin, yMin = pointOfIntersection, y1
                                            minDistSurface = [surface_i]
            return xMin, yMin, minDistSurface, minDistance

        mp.prec = digits
        segments = int(segments)
        if(segments <= 1):
            segments = 2
        passes = 0
        listOfRays = self.listOfRays
        while(passes < passLimit):
            listOfThingsToAdd = []
            for i in range(len(listOfRays)):
                ray = listOfRays[i]
                if(ray.intensity < threshold or not ray.worthUsing):
                    continue
                isMirror = 0
                xMin1, yMin1, minDistSurface1, distance1 = surfaceFinderPart(ray, self.listOfLenses,segments=segments)
                xMin2, yMin2, minDistSurface2, distance2 = surfaceFinderPart(ray, self.listOfMirrors,segments=segments)
                if(distance1 > 0 and distance2 > 0):
                    if(distance1 < distance2):
                        xMin, yMin, minDistSurface = xMin1, yMin1, minDistSurface1
                    else:
                        isMirror = 1
                        xMin, yMin, minDistSurface = xMin2, yMin2, minDistSurface2
                elif(distance1 > 0):
                    xMin, yMin, minDistSurface = xMin1, yMin1, minDistSurface1
                elif(distance2 > 0):
                    isMirror = 1
                    xMin, yMin, minDistSurface = xMin2, yMin2, minDistSurface2
                else:
                    # no points of intersection, may as well move on to the
                    # next ray
                    ray.worthUsing = 0
                    continue
                # we have now found the next point of intersection or exhausted
                # all possibilities so we can move on to the rest of this
                # before we proceed we again need to have the special case
                # handled for vertical lines
                if(minDistSurface == -1):
                    ray.worthUsing = 0
                    continue
                if(len(minDistSurface) == 2):  # the special case
                    xMin, normalVector = xMin[0], [-1, 0]
                else:
                    normalVector = normalVectorFunc(minDistSurface[0], xMin)
                if(np.dot(normalVector, ray.unitVector) > 0):
                    normalVector = [-i for i in normalVector]
                theta1 = np.arccos(np.dot(ray.unitVector, normalVector))

                if(not isMirror):
                    if(len(minDistSurface) == 2):
                        nextObj1 = self.locatingIndexFunction((xMin+0.001,
                                                               yMin))
                        nextObj2 = self.locatingIndexFunction((xMin-0.001,
                                                               yMin))
                    else:
                        nextObj1 = self.locatingIndexFunction((xMin, yMin +
                                                               0.0001))
                        nextObj2 = self.locatingIndexFunction((xMin, yMin -
                                                               0.0001))
                    if(nextObj1 != ray.objIndex):
                        nextObj = nextObj1
                    else:
                        nextObj = nextObj2

                    if(ray.objIndex == -1 and nextObj == -1):
                        continue
                    elif(ray.objIndex == -1):
                        n1 = self.nMedium
                        n2 = self.mapOfLenses.get(nextObj).refractiveIndex
                    elif(nextObj == -1):
                        n1 = self.mapOfLenses.get(ray.objIndex).refractiveIndex
                        n2 = self.nMedium
                    else:
                        n1 = self.mapOfLenses.get(ray.objIndex).refractiveIndex
                        n2 = self.mapOfLenses.get(nextObj).refractiveIndex
                    wavelength = ray.wavelength
                    n_1, n_2 = n1(wavelength), n2(wavelength)
                    if(theta1 > pi/2):
                        theta1 = pi - theta1
                    # internal reflection case need to do this still
                    if(abs(n_1) > abs(n_2) and theta1 >= np.arcsin(n_2/n_1)):#internal reflection case need to do this still
                        theta2 = theta1
                        refractedRay1, refractedRay2 = refractedRayFunction(normalVector, theta1, 0)
                        if(np.dot(refractedRay1, unitVector) > np.dot(refractedRay2, unitVector)):#don't use absolute values here
                            refractedRay = refractedRay1
                        else:
                            refractedRay = refractedRay2
                        ray.intensityList.append(ray.intensity)
                        nextObj = ray.objIndex
                    else:
                        i0 = ray.intensity
                        theta2 = 0 if(theta1 == 0) else np.arcsin(n_1/n_2*np.sin(theta1))
                        refractedRay1, refractedRay2 = refractedRayFunction(normalVector, theta2, 1)
                        if(abs(np.dot(refractedRay1, unitVector)) > abs(np.dot(refractedRay2, unitVector))):
                            refractedRay = refractedRay1
                        else:
                            refractedRay = refractedRay2
                        reflect1, reflect2 = refractedRayFunction(normalVector, theta1, 0)
                        if(np.dot(reflect1, unitVector) > np.dot(reflect2, unitVector)):
                            reflect = reflect1
                        else:
                            reflect = reflect2
                        ray.update_intensity_and_polarization(theta1, theta2, n_1, n_2)
                        reflectedXaxisTheta = np.arctan2(reflect[1],reflect[0])
                        if(reflectedXaxisTheta < 0):
                            reflectedXaxisTheta += 2*np.pi
                        Rperp = Ray.rPerp(theta1,theta2,n_1,n_2)**2
                        Rpara = Ray.rPara(theta1,theta2,n_1,n_2)**2
                        s = Rperp/(Rperp + Rpara)
                        newRay = Ray(xMin, yMin, reflectedXaxisTheta,
                                     intensity = abs(i0 - ray.intensity),
                                     objIndex = ray.objIndex,sPolarization=s, pPolarization=(1-s))
                        #comment this part out to look at the other behavior
                        listOfThingsToAdd.append(newRay)
                else:
                    # We can add the loss function bit later
                    theta2 = theta1
                    refractedRay1, refractedRay2 = refractedRayFunction(
                        normalVector, theta2, 1)

                    # Don't use absolute values here
                    if(np.dot(refractedRay1, ray.unitVector) >
                       np.dot(refractedRay2, ray.unitVector)):
                        refractedRay = refractedRay1
                    else:
                        refractedRay = refractedRay2
                    ray.intensityList.append(ray.intensity)
                    nextObj = ray.objIndex

                ray.normalVectorHistory.append([[xMin, normalVector[0]+xMin],
                                                [yMin, normalVector[1]+yMin]])
                refractedXaxisTheta = float(mpmath.atan2(
                    mpmath.mpmathify(refractedRay[1]),
                    mpmath.mpmathify(refractedRay[0])))

                if(refractedXaxisTheta < 0):
                    refractedXaxisTheta += 2*np.pi
                #ray.angleHistory.append(theta1*180/pi)
                #ray.angleHistory.append(theta2*180/pi)
                ray.update_location_history(xMin, yMin,
                                            refractedXaxisTheta, nextObj)
            if(extraRays):
                listOfRays+=listOfThingsToAdd
            passes += 1

    # This method will probably keep growing however it's not too bad right now
    def plotter(self, fig):
        sub = fig.add_subplot(111)
        surfaceList = self.surfaceList
        # First we have to handle graphing the lens or lenses
        for s in range(len(surfaceList)):
            domain, surfaces = surfaceList[s]

            # This will need some refining for the 4 point case
            if(len(domain) == 1):
                x = domain[0]
                xList = [x]*len(surfaces)
                sub.plot(xList, surfaces)
            else:
                x = np.linspace(domain[0], domain[1], 1000)

                # This will probably need special handling for constant
                # functions
                for function in surfaces:
                    tester_vec = np.frompyfunc(function, 1, 1)
                    tester = tester_vec(x[:20])

                    if(type(tester) is int or type(tester) is float
                       or type(tester) is complex):
                        yList = [function(domain[0])]*len(x)
                    else:
                        yList = tester_vec(x)
                    realPiece = np.frompyfunc(mpmath.re, 1, 1)
                    sub.plot(x, realPiece(yList))
        rayList = self.listOfRays
        for i in range(len(rayList)):
            xAx, yAx = zip(*rayList[i].locationHistory)
            xAxList, yAxList = list(xAx), list(yAx)
            ray = rayList[i]
            angle = ray.angle
            slope = ray.slope
            ray1, ray2 = [xAxList[-1]-5, yAxList[-1]-5*slope], [xAxList[-1]+5,
                                                                yAxList[-1]+5
                                                                * slope]
            checkAngle = mpmath.atan2(ray1[1]-yAxList[-1], ray1[0]-xAxList[-1])
            if(checkAngle < 0):
                checkAngle += 2*pi
            if(checkAngle > angle + 0.001 or checkAngle < angle - 0.001):
                xAxList.append(ray2[0])
                yAxList.append(ray2[1])
            else:
                xAxList.append(ray1[0])
                yAxList.append(ray1[1])

            # Change this piece for different plot configurations
            for j in range(len(ray.intensityList)):
                sub.plot(xAxList[j:j+2], yAxList[j:j+2], 'k',
                         alpha=(ray.intensityList[j]))

    def plotterJustLens(self, fig):  # for debugging
        sub = fig.add_subplot(111)
        surfaceList = self.surfaceList
        for s in range(len(surfaceList)):
            domain, surfaces = surfaceList[s]
            # this will need some refining for the 4 point case
            if(len(domain) == 1):
                x = domain[0]
                xList = [x]*len(surfaces)
                sub.plot(xList, surfaces)
            else:
                x = np.linspace(domain[0], domain[1], 1000)
                # This will probably need special handling for constant
                # functions
                # This seems to be a constant function, should probably put in
                # more
                for function in surfaces:
                    if(function(1) == function(2)
                       and function(1) == function(1.5)):
                        # robust handling when the domain for the contain
                        # function isn't so nice
                        yList = [function(1)]*len(x)
                    else:
                        yList = function(x)
                    sub.plot(x, yList)
