import scipy.optimize as opt
from functools import partial
from Components.Presets import hyperbola, circle, circle_circle, line


class Lens:
    def __init__(self, surfaceList, n=lambda x: 1.5):

        self.refractiveIndex = n
        for domain, surfaces in surfaceList:
            for i in range(len(domain)):
                if(isinstance(domain[i], complex)):
                    domain[i] = domain[i].real
            if(len(domain) == 1):
                for j in range(len(surfaces)):
                    if(isinstance(surfaces[j], complex)):
                        surfaces[j] = surfaces[j].real
        self.surfaceList = surfaceList

    @classmethod
    def predefined(cls, name, nFunction=lambda x: 1.5, xOffSet=0, yOffSet=0,
                   r1=0, r2=0):
        """
        I should probably include the r1, r2, height, and width values in these
        examples
        ***NOTE***
        Conventions:
        1. the rays typically move vertically from - to +
        When we say one thing is before another we mean that the former has a
        lesser "y" value than the latter
        Inputs:
        name - the type of predefined lens desired
        nFunction - an index of refraction function based on wavelength
        xOffSet, yOffSet - horizontal and vertical displacements
        Output:
        An instance of the class that the specified lens
        Description:
        The primary purpose of this construtor is to illustrate the formate of
        building a lens. The portions that aren't labeled hyperbolic or
        elliptic follow the conventions of gaussian optics which are detailed
        in that constructor
        """

        # should probably change the convex parts so that they reflect the
        # updated methodology
        if(name == "biconvex"):
            f1 = partial(circle, False, 1 + yOffSet, 3, 5-xOffSet)
            f2 = partial(circle, True, 5 + yOffSet, 3, 5-xOffSet)
            x0, x1 = circle_circle(r1, r2, 2, xOffSet)
            return cls([[[x0, x1], [f1, f2]]], n=nFunction)
        elif(name == "biconcave"):
            f1 = partial(circle, False, 1 + yOffSet, 3, 5-xOffSet)
            f2 = partial(circle, True, 8 + yOffSet, 3, 5-xOffSet)
            mySurfaceList = [[[3+xOffSet], [f1(3+xOffSet), f2(3+xOffSet)]],
                             [[3+xOffSet, 7+xOffSet], [f1, f2]],
                             [[7+xOffSet], [f1(7+xOffSet), f2(7+xOffSet)]]]
            return cls(mySurfaceList, n=nFunction)
        elif(name == "planar convex"):
            f1 = partial(circle, True, 5 + yOffSet, 3, 5-xOffSet)
            f2 = partial(line, 3 + yOffSet)
            x0 = opt.fsolve(lambda x: f1(x)-f2(x), 2+xOffSet,
                            xtol=1e-10, maxfev=200)
            x1 = opt.fsolve(lambda x: f1(x)-f2(x), 7+xOffSet,
                            xtol=1e-10, maxfev=200)
            return cls([[[x0[0], x1[0]], [f1, f2]]], n=nFunction)
        elif(name == "planar concave"):
            f1 = partial(circle, False, 1 + yOffSet, 3, 5-xOffSet)
            f2 = partial(line, 6 + yOffSet)
            mySurfaceList = [[[3+xOffSet], [f1(3+xOffSet), f2(3+xOffSet)]],
                             [[3+xOffSet, 7+xOffSet], [f1, f2]],
                             [[7+xOffSet], [f1(7+xOffSet), f2(7+xOffSet)]]]
            return cls(mySurfaceList, n=nFunction)
        elif(name == "meniscus convex"):
            f1 = partial(circle, True, 6 + yOffSet, 4, 5-xOffSet)
            f2 = partial(circle, True, 9 + yOffSet, 6, 5-xOffSet)
            x0 = opt.fsolve(lambda x: f1(x)-f2(x), 2+xOffSet,
                            xtol=1e-10, maxfev=200)
            x1 = opt.fsolve(lambda x: f1(x)-f2(x), 8+xOffSet,
                            xtol=1e-10, maxfev=200)
            return cls([[[x0[0], x1[0]], [f1, f2]]], n=nFunction)
        elif(name == "meniscus concave"):
            # I could save using the append opperation and memory
            # initialization by just expliciting spelling out what I wanted my
            # list to be
            f1 = partial(circle, True, 4.5 + yOffSet, 3, 5-xOffSet)
            f2 = partial(circle, True, 8 + yOffSet, 7, 5-xOffSet)
            mySurfaceList = [[[3+xOffSet], [f1(3+xOffSet), f2(3 + xOffSet)]],
                             [[3+xOffSet, 7+xOffSet], [f1, f2]],
                             [[7+xOffSet], [f1(7+xOffSet), f2(7 + xOffSet)]]]
            return cls(mySurfaceList, n=nFunction)
        elif(name == "hyberbolic biconvex"):
            f1 = partial(hyperbola, True, yOffSet, 5, 5-xOffSet)
            f2 = partial(hyperbola, False, yOffSet, 5, 5-xOffSet)
            mySurfaceList = [[[2+xOffSet], [f2(2+xOffSet), f1(2+xOffSet)]],
                             [[2+xOffSet, 8+xOffSet], [f2, f1]],
                             [[8+xOffSet], [f2(8+xOffSet), f1(8+xOffSet)]]]
            return cls(mySurfaceList, n=nFunction)
        elif(name == "hyberbolic plano convex"):
            f1 = partial(hyperbola, True, 8 + yOffSet, 1, 5-xOffSet)
            f2 = partial(line, 1 + yOffSet)
            mySurfaceList = [[[3+xOffSet], [f2(3+xOffSet), f1(3+xOffSet)]],
                             [[3+xOffSet, 7+xOffSet], [f2, f1]],
                             [[7+xOffSet], [f2(7+xOffSet), f1(7+xOffSet)]]],
            return cls(mySurfaceList, n=nFunction)
        elif(name == "ellipse biconvex"):
            f1 = partial(hyperbola, True, yOffSet, 5, 5-xOffSet)
            f2 = partial(hyperbola, False, yOffSet, 5, 5-xOffSet)
            mySurfaceList = [[[2+xOffSet], [f2(2+xOffSet), f1(2+xOffSet)]],
                             [[2+xOffSet, 8+xOffSet], [f2, f1]],
                             [[8+xOffSet], [f2(8+xOffSet), f1(8+xOffSet)]]]
            return cls(mySurfaceList, n=nFunction)
        else:
            raise(Exception("Unknown lens type, please enter one of the \
                            following:\nbiconvex, biconcave, planar convex, \
                            planar concave, meniscus convex, meniscus \
                            concave"))

    @classmethod
    def gaussianOptics(cls, r1=0, r2=0, thickness=0, diameter=0, xOffSet=0,
                       yOffSet=0, nFunction=lambda x: 1.5):
        """
        ***NOTE***
        Conventions:
        1. the rays typically move vertically from - to +. When we say one
           thing is before another we mean that the former has a lesser "y"
           value than the latter
        2. if the arc that defines a surface of a lens comes before the center,
           that radius is positive, and negative otherwise
        3. in place of using infinity for describing planar surfaces, we let
           the radius go to zero because this is an otherwise unused value
        Inputs:
        r1 - radius of curvature of the first surface (in the vertical
             direction)
        r2 - radius of curvature of the second surface
        thickness - the vertical distance between the vertices of the two
                    surfaces
        diameter -  the length of the domain of the lens (their horizontal
                    extent)
        xOffSet, yOffSet - horiztonal and vertical displacements respectively
        nFunction - a index of refraction function dependent on wavelength
        Output:
        An instance on the lens class with a surface built by the inputs.
        ***NOTE***
        Due to some quirk in the program, the lens doesn't always behave
        correctly when it is centered on the y axis, therefore an offset is
        built in to ensure that this isn't normally encountered
        Description:
        Using the conventions of gaussian optics and our surface representation
        scheme we build the lens specified by the parameters given. There are
        six lenses that may be produced. Those are as follows:
        biconvex, biconcave, meniscus convex, meniscus concave, plano concave,
        plano convex.
        Below are some helper functions to condense the code and prevent having
        to rewrite the same code.
        After those are a variety of cases for producing different lenses based
        on the inputs.
        """

        def abrvFunc(r1, r2, f1, f2, diameter, x0=0, x1=0, xOffSet=0,
                     thickness=0):
            # This is a helper function to make the code more readable
            # (or less eye watering if you prefer)
            rmax = max(r1, r2)
            if((x0 or x1) and abs(x1-x0) <= diameter):
                raise(Exception("this diameter exceeds the limitations imposed\
                                by the specified diameter"))
            xmin = rmax - diameter/2
            xmax = rmax + diameter/2
            mySurfaceList = [[[xmin+xOffSet], sorted(
                [f1(xmin+xOffSet), f2(xmax+xOffSet)])]]
            mySurfaceList.append([[xmin+xOffSet, xmax+xOffSet], [f1, f2]])
            mySurfaceList.append([[xmax+xOffSet], sorted(
                [f1(xmax+xOffSet), f2(xmax+xOffSet)])])
            return mySurfaceList

        def widthHandler(f1, f2, x0, x1, surfaceList):
            f1y0, f1y1, f2y0, f2y1 = f1(x0), f1(x1), f2(x0), f2(x1)
            if(f1y0 < f2y0 - 0.001 or f1y0 > f2y0 + 0.001):
                surfaceList.insert(0, [[x0], sorted([f1y0, f2y0])])
            if(f1y1 < f2y1 - 0.001 or f1y1 > f2y1 + 0.001):
                surfaceList.append([[x1], sorted([f1y1, f2y1])])

        if(((r1 or r2 or diameter) and thickness)
           or ((r1 or r2) and diameter) or (r1 and r2)):
            # This should say that at least two of the following
            # (r1, r2, width, height) have to be given
            if(thickness < 0 or diameter < 0):  # deal with negative width
                raise(Exception("width and height are constrained to be \
                                greater than or equal to zero"))

            elif(thickness == 0 and diameter == 0):
                raise(Exception("at minimum a width or height must be \
                                specified for a lens along with r1 or r2"))

            elif(r1 > 0 and r2 > 0 and r1 < r2):
                # meniscus convex case 1
                if(not thickness):
                    thickness = abs(r1-(r1**2 - diameter**2/4)**(1/2) -
                                    r2 + (r2**2-diameter**2/4)**(1/2))

                f1 = partial(circle, True,
                             r2+1 + yOffSet - r2 + r1 - thickness, r1,
                             r2-xOffSet)
                f2 = partial(circle, True, r2 + 1 + yOffSet, r2, r2-xOffSet)
                x0, x1 = circle_circle(r1, r2, thickness, xOffSet)
                if(diameter):
                    return cls(abrvFunc(r1, r2, f1, f2, diameter, x0, x1,
                                        xOffSet=xOffSet),
                               n=nFunction)
                surfaceList = [[[x0, x1], [f1, f2]]]
                widthHandler(f1, f2, x0, x1, surfaceList)
                return cls(surfaceList, n=nFunction)

            elif(r1 < 0 and r2 < 0 and r2 > r1):
                # meniscus convex case 2
                if(not thickness):
                    thickness = abs(r1-(r1**2 - diameter**2/4)**(1/2) - r2 +
                                    (r2**2-diameter**2/4)**(1/2))
                f1 = partial(circle, False, 1 + yOffSet, r1, r1 + xOffSet)
                r1, r2 = -r1, -r2
                f2 = partial(circle, False, 1 + yOffSet + r1 + thickness - r2,
                             r2, r1 + xOffSet)
                x0, x1 = circle_circle(r1, r2, thickness, xOffSet)
                if(diameter):
                    return cls(abrvFunc(r1, r2, f1, f2, diameter, x0, x1,
                                        xOffSet=xOffSet), n=nFunction)
                surfaceList = [[[x0, x1], [f1, f2]]]
                widthHandler(f1, f2, x0, x1, surfaceList)
                return cls(surfaceList, n=nFunction)

            elif(r1 > 0 and r2 == 0):
                # plano-convex case 1
                if(not thickness):
                    thickness = abs(r1-(r1**2-diameter**2/4)**(1/2))
                f1 = partial(circle, 1 + yOffSet + r1, r1, r1 + xOffSet)
                f2 = lambda x: 1 + yOffSet + thickness
                d = r1 - thickness
                b = -2*(r1+xOffSet)
                c = d**2 + (r1+xOffSet)**2 - r1**2
                x0 = (-b-(b**2 - 4*c)**(1/2))/2
                x1 = (-b+(b**2 - 4*c)**(1/2))/2

                if(diameter):
                    return cls(abrvFunc(r1, r2, f1, f2, diameter, x0, x1,
                                        xOffSet=xOffSet), n=nFunction)
                surfaceList = [[[x0, x1], [f1, f2]]]
                widthHandler(f1, f2, x0, x1, surfaceList)
                return cls(surfaceList, n=nFunction)

            elif(r1 > 0 and r2 < 0):
                # biconvex
                r2 = -r2
                if(not thickness):
                    thickness = abs(r1-(r1**2-diameter**2/4)**(1/2)) +\
                        abs(r2-(r2**2-diameter**2/4)**(1/2))
                rmax = max(r1, r2)
                f1 = partial(circle, False,
                             rmax + 1 + thickness - r1 - r2 + yOffSet, r2,
                             rmax - xOffSet)
                f2 = partial(circle, True,
                             rmax + 1 + yOffSet, r1, rmax - xOffSet)
                d = r1 + r2 - thickness
                x0, x1 = circle_circle(r1, r2, thickness, xOffSet)
                if(diameter):
                    return cls(abrvFunc(r1, r2, f1, f2, diameter,
                                        x0, x1, xOffSet=xOffSet),
                               n=nFunction)
                surfaceList = [[[x0, x1], [f1, f2]]]
                widthHandler(f1, f2, x0, x1, surfaceList)
                return cls(surfaceList, n=nFunction)

            elif(r1 == 0 and r2 > 0):
                # plano-concave case 1
                if(diameter == 0 or thickness == 0):
                    raise(Exception("Insufficient information provided for \
                                    concave class object, both diameter and \
                                    thickness must be postive"))
                f1 = lambda x: 1 + yOffSet - thickness
                f2 = partial(circle, True, 1 + r2 + yOffSet, r2, r2-xOffSet)
                return cls(abrvFunc(r1, r2, f1, f2, diameter, xOffSet=xOffSet),
                           n=nFunction)

            elif(r1 == 0 and r2 == 0):
                # rectangle
                mySurfaceList = [[[1 + xOffSet],
                                  [1 + yOffSet, thickness + 1 + yOffSet]]]
                mySurfaceList.append([[1 + xOffSet, diameter + 1 + xOffSet],
                                      [lambda x: 1 + yOffSet,
                                       lambda x: thickness + 1 + yOffSet]])
                mySurfaceList.append([[diameter + 1 + xOffSet],
                                      [1 + yOffSet, thickness + 1 + yOffSet]])
                return cls(mySurfaceList, n=nFunction)

            elif(r1 == 0 and r2 < 0):
                # plano-convex case 2
                if(not thickness):
                    thickness = abs(r2-(r2**2-diameter**2/4)**(1/2))
                r2 = -r2
                f1 = partial(circle, True, 1 + yOffSet, r2, r2-xOffSet)
                f2 = lambda x: 1 + yOffSet + r2 - thickness
                b = -2*(r2+xOffSet)
                c = (r2 - thickness)**2 + (r2+xOffSet)**2 - r2**2
                x0 = (-b-(b**2 - 4*c)**(1/2))/2
                x1 = (-b+(b**2 - 4*c)**(1/2))/2
                if(diameter):
                    return cls(abrvFunc(r1, r2, f1, f2, diameter, x0, x1,
                                        xOffSet=xOffSet), n=nFunction)
                surfaceList = [[[x0, x1], [f1, f2]]]
                widthHandler(f1, f2, x0, x1, surfaceList)
                return cls(surfaceList, n=nFunction)

            elif(r1 < 0 and r2 > 0):
                # biconcave
                if(diameter == 0 or thickness == 0):
                    raise(Exception("Insufficient information provided for \
                                    concave class object, both diameter and \
                                    thickness must be positive"))
                r1 = -r1
                rmax = max(r1, r2)
                f1 = partial(circle, False, rmax + 1 + yOffSet, r1,
                             rmax + xOffSet)
                f2 = partial(circle, True,
                             rmax+1 + thickness + r1 + r2 + yOffSet,
                             r2, rmax + xOffSet)
                return cls(abrvFunc(r1, r2, f1, f2, diameter, xOffSet=xOffSet),
                           n=nFunction)

            elif(r1 < 0 and r2 == 0):
                # plano-concave case 2
                if(diameter == 0 or thickness == 0):
                    raise(Exception("Insufficient information provided for \
                                    concave class object, both diameter and \
                                    thickness must be positive"))
                r1 = -r1
                f1 = partial(circle, False, 1 + yOffSet, r1, r1 + xOffSet)
                f2 = lambda x: 1 + r1 + thickness + yOffSet
                return cls(abrvFunc(r1, r2, f1, f2, diameter, xOffSet=xOffSet),
                           n=nFunction)

            elif(r1 > 0 and r2 > 0):
                # meniscus concave case 1
                if(diameter == 0 or thickness == 0):
                    raise(Exception("Insufficient information provided for \
                                    concave class object, both diameter and \
                                    thickness must be positive"))
                f1 = partial(circle, True, r1 + 1 + yOffSet, r1, r1 + xOffSet)
                f2 = partial(circle, True,
                             r1 + 1 - thickness + r2 - r1 + yOffSet, r2,
                             r1 + xOffSet)
                if(diameter < 2*r2):  # reasonable diameter
                    return cls(abrvFunc(r1, r2, f1, f2, diameter,
                                        xOffSet=xOffSet), n=nFunction)
                raise(Exception("This is an unreasonable diameter for a \
                                meniscus concave lens"))
            elif(r1 < 0 and r2 < 0):  # meniscus concave case 2
                r1, r2 = -r1, -r2
                if(diameter == 0 or thickness == 0):
                    raise(Exception("Insufficient information provided for \
                                    concave class object, both diameter and \
                                    thickness must be positive"))
                f1 = partial(circle, False, 1 + yOffSet, r1, r2 + xOffSet)
                f2 = partial(circle, False, 1 - thickness - r2 + r1 + yOffSet,
                             r2, r2 + xOffSet)
                if(diameter < 2*r1):
                    # reasonable diameter
                    return cls(abrvFunc(r1, r2, f1, f2, diameter,
                                        xOffSet=xOffSet), n=nFunction)
                raise(Exception("This is an unreasonable diameter for a \
                                meniscus concave lens"))
        raise(Exception("Not enough information supplied or diameter and \
                        thickness are 0 which is impossible"))

    @classmethod
    def custom(cls, surfaceList, nFunction=lambda x: 1.5):
        return cls(surfaceList, n=nFunction)
