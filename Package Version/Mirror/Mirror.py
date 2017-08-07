from functools import partial
from Components.Presets import circle, line


class Mirror:
    def __init__(self, surfaceList=[], lossFunction=lambda i: 0):
        self.surfaceList = surfaceList
        self.lossFunction = lossFunction

    @classmethod
    def predefined(cls, name, lossFunction=lambda i: 0, xOffSet=0, yOffSet=0):
        """
        Inputs:
        name - name of the desired mirror
        lossFunction - a function to reduce the intensity of incident light
        #### NOT IMPLEMENTED YET
        xOffSet, yOffSet - horizontal and vertical displacements

        Outputs:
        A mirror with the specified shape

        Description:
        This constructor provides several example mirrors that could be used on
        their own or studied to learn how to build another mirror
        """
        if(name == "concave up"):
            f = partial(circle, True, 5+yOffSet, 9, xOffSet)
            return cls(surfaceList=[[[xOffSet, 3+xOffSet], [f]]],
                       lossFunction=lossFunction)
        elif(name == "concave down"):
            f = partial(circle, False, 5+yOffSet, 16, xOffSet)
            return cls(surfaceList=[[[xOffSet, 4+xOffSet], [f]]],
                       lossFunction=lossFunction)
        elif(name == "horizontal"):
            f = partial(line, 4 + yOffSet)
            return cls(surfaceList=[[[xOffSet, 5+xOffSet], [f]]],
                       lossFunction=lossFunction)
        elif(name == "vertical"):
            return cls(surfaceList=[[[3+xOffSet], [1+yOffSet, 6+yOffSet]]],
                       lossFunction=lossFunction)
        elif(name == "parabolic down"):
            return cls(surfaceList=[[[xOffSet, 10+xOffSet],
                                     [lambda x: -(x-5-xOffSet)**2+8+yOffSet]]])
        elif(name == "hyperbolic down"):
            return cls(surfaceList=[[[xOffSet, 10+xOffSet],
                                     [lambda x: -((x-5-xOffSet)**2+1)**(1/2)
                                      + 8 + yOffSet]]])
        else:
            raise(Exception("Unknown predefined case\nThe acceptable cases\
                            are: concave up, concave down, vertical, and\
                             horizontal"))

    @classmethod
    def gaussianOptics(cls, R=0, vertical=0, diameter=0, xOffSet=0, yOffSet=0,
                       lossFunction=lambda i: 0, x0=0, y0=0, x1=0, y1=0):
        """
        Inputs:
        R - radius of curvature of the mirror, negative indicates that the
                   center is below the arc of the mirror
        vertical - if R = 0 then we have a planar mirror, vertical is a flag to
                   indicate the orientation of the mirror
        diameter - specifies the horizontal extent of the mirror (or vertical
                   extent if R = 0 and vertical = 1)
        xOffSet, yOffSet - self explanatory
        lossFunction - a function to describe how much energy is lost upon a
                       reflection ### not implemented

        Outputs:
        An instance of mirror with the traits required by the parameters.

        Description:
        Specifies a mirror by the radius of curvature, if the radius is
        negative, that indicates that the center is below the arc of the mirror
        (or that the center has a lesser y coordinate compared to the arc)

        """
        print("Gauss")
        if(diameter or R):
            if(R > 0):
                f = partial(circle, True, yOffSet, R, xOffSet)
                if(diameter):
                    if(diameter < 2*R):
                        return cls(surfaceList=[[[xOffSet-diameter/2,
                                                  xOffSet+diameter/2],
                                                 [f]]],
                                   lossFunction=lossFunction)
                    else:
                        raise(Exception("This diameter exceeds the domain of\
                                        the mirror surface"))
                return cls([[[xOffSet-R, xOffSet+R], [f]]],
                           lossFunction=lossFunction)
            elif(R < 0):
                R = -R

                f = partial(circle, False, yOffSet, R, xOffSet)
                if(diameter):
                    if(diameter < 2*R):
                        return cls(surfaceList=[[[xOffSet-diameter/2,
                                                  xOffSet+diameter/2],
                                                 [f]]],
                                   lossFunction=lossFunction)
                    else:
                        raise(Exception("This diameter exceeds the domain of\
                                        the mirror surface"))
                return cls(surfaceList=[[[xOffSet-R, xOffSet+R], [f]]],
                           lossFunction=lossFunction)
            elif(diameter):
                if(vertical):
                    return cls(surfaceList=[[[xOffSet], [yOffSet-diameter/2,
                                                         yOffSet+diameter/2]]])
                return cls(surfaceList=[[[xOffSet-diameter/2,
                                          xOffSet+diameter/2],
                                         [lambda x: 3+yOffSet]]])
            # else....
            return cls(surfaceList=[[[x0, x1], lambda x: y1-y0/(x1-x0)]],
                       lossFunction=lossFunction)

    @classmethod
    def custom(cls, surfaceList, lossFunction=lambda i: 0):
        return cls(surfaceList=surfaceList, lossFunction=lossFunction)
