def circle(isNeg, yOffSet, rSquared, x0, x1):
    sign = -1 if isNeg else 1
    return sign*(-(x1-x0)**2 + rSquared*rSquared)**(1/2) + yOffSet


def line(offset, x):
    return offset + x


def hyperbola(isNeg, yOffSet, r, x0, x1):
    sign = -1 if isNeg else 1
    return sign*(-5/9 * (x1 - x0)**2 + r)**(1/2) + yOffSet


def circle_circle(r1, r2, thickness, xOffSet):
    d = -r1 + r2 + thickness
    if(d < 0):
        d = -d
    elif(d == 0):
        d = 0.001
    a = (1/d)*(4*d**2*r1**2 - (d**2 - r2**2 + r1**2)**2)**(1/2)
    return -a/2 + r2 + xOffSet, a/2 + r2 + xOffSet
