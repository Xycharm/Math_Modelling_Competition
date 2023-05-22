"""algorithm for optimization"""
from geomentry_structure import *
import math
import cmath

def single_quartic(a0, b0, c0, d0, e0):
    ''' Reduce the quartic equation to to form:
        x**4 + a*x**3 + b*x**2 + c*x + d = 0'''
    a, b, c, d = b0 / a0, c0 / a0, d0 / a0, e0 / a0

    # Some repeating variables
    a0 = 0.25 * a
    a02 = a0 * a0

    # Coefficients of subsidiary cubic equation
    p = 3 * a02 - 0.5 * b
    q = a * a02 - b * a0 + 0.5 * c
    r = 3 * a02 * a02 - b * a02 + c * a0 - d

    # One root of the cubic equation
    z0, _, _ = single_cubic(1, p, r, p * r - 0.5 * q * q)

    # Additional variables
    s = cmath.sqrt(2 * p + 2 * z0.real + 0j)
    if s == 0:
        t = z0 * z0 + r
    else:
        t = -q / s

    # Compute roots by quadratic equations
    r0, r1 = single_quadratic(1, s, z0 + t)
    r2, r3 = single_quadratic(1, -s, z0 - t)

    return r0 - a0, r1 - a0, r2 - a0, r3 - a0

def single_quadratic(a0, b0, c0):
    ''' Reduce the quadratic equation to to form:
        x**2 + a*x + b = 0 '''
    a, b = b0 / a0, c0 / a0

    # Some repeating variables
    a0 = -0.5*a
    delta = a0*a0 - b
    sqrt_delta = cmath.sqrt(delta)

    # Roots
    r1 = a0 - sqrt_delta
    r2 = a0 + sqrt_delta

    return r1, r2

def single_cubic(a0, b0, c0, d0):
    ''' Reduce the cubic equation to to form:
        x**3 + a*x**2 + b*x + c = 0 '''
    a, b, c = b0 / a0, c0 / a0, d0 / a0

    # Some repeating constants and variables
    third = 1./3.
    a13 = a*third
    a2 = a13*a13
    sqr3 = math.sqrt(3)

    # Additional intermediate variables
    f = third*b - a2
    g = a13 * (2*a2 - b) + c
    h = 0.25*g*g + f*f*f

    def cubic_root(x):
        ''' Compute cubic root of a number while maintaining its sign'''
        if x.real >= 0:
            return x**third
        else:
            return -(-x)**third

    if f == g == h == 0:
        r1 = -cubic_root(c)
        return r1, r1, r1

    elif h <= 0:
        j = math.sqrt(-f)
        k = math.acos(-0.5*g / (j*j*j))
        m = math.cos(third*k)
        n = sqr3 * math.sin(third*k)
        r1 = 2*j*m - a13
        r2 = -j * (m + n) - a13
        r3 = -j * (m - n) - a13
        return r1, r2, r3

    else:
        sqrt_h = cmath.sqrt(h)
        S = cubic_root(-0.5*g + sqrt_h)
        U = cubic_root(-0.5*g - sqrt_h)
        S_plus_U = S + U
        S_minus_U = S - U
        r1 = S_plus_U - a13
        r2 = -0.5*S_plus_U - a13 + S_minus_U*sqr3*0.5j
        r3 = -0.5*S_plus_U - a13 - S_minus_U*sqr3*0.5j
        return r1, r2, r3

def PCP(a, b, O):
    """
    using PCP method,return the optimal reflexive point of point a,b on circle O
    """
    x_O = O.o.x
    y_O = O.o.y
    x_a = a.x - x_O
    x_b = b.x - x_O
    y_a = a.y - y_O
    y_b = b.y - y_O
    r = O.r
    A = 4*(x_a*x_a*x_b*x_b+x_b*x_b*y_a*y_a+x_a*x_a*y_b*y_b+y_a*y_a*y_b*y_b)
    B = (-4)*r*(x_a*x_a*x_b+x_a*x_b*x_b+x_b*y_a*y_a+x_a*y_b*y_b)
    C = r*r*x_a*x_a+2*r*r*x_a*x_b+r*r*x_b*x_b-4*x_a*x_a*x_b*x_b+r*r*y_a*y_a-4*x_b*x_b*y_a*y_a+2*r*r*y_a*y_b+r*r*y_b*y_b-4*x_a*x_a*y_b*y_b-4*y_a*y_a*y_b*y_b
    D = 2*r*(2*x_a*x_a*x_b+2*x_a*x_b*x_b+x_b*y_a*y_a+x_a*y_b*y_b-x_a*y_a*y_b-x_b*y_a*y_b)
    E = x_a*x_a*y_b*y_b+2*x_a*x_b*y_a*y_b+x_b*x_b*y_a*y_a-r*r*x_b*x_b-2*r*r*x_a*x_b-r*r*x_a*x_a
    result= single_quartic(A,B,C,D,E)
    #function to calculate distance
    def distm(xp, yp):
        return math.sqrt((x_a - xp) ** 2 + (y_a - yp) ** 2) + math.sqrt((x_b - xp) ** 2 + (y_b - yp) ** 2)
    mdist = math.inf
    for i in range(len(result)):
        if result[i].imag==0:
            item = result[i].real
            x_p1=item*r
            x_p2=item*r
            y_p1=r*math.sqrt(1-item*item)
            y_p2=-r*math.sqrt(1 - item * item)

            if distm(x_p1, y_p1) < mdist:
                minx = x_p1
                miny = y_p1
                mdist=distm(x_p1, y_p1)

            if distm(x_p2, y_p2) < mdist:
                minx = x_p2
                miny = y_p2
                mdist = distm(x_p2, y_p2)

    point_m = Point(minx+x_O, miny+y_O)
    return point_m


def intersection(a, O, b):
    """
    calculate the intersection between a circle and a line segment
    """
    x_O = O.o.x
    y_O = O.o.y
    x_a = a.x - x_O
    x_b = b.x - x_O
    y_a = a.y - y_O
    y_b = b.y - y_O
    r = O.r
    if x_a != x_b:
        k = (y_a - y_b) / (x_a - x_b)
        b = y_a - k * x_a
        if k ** 2 * b ** 2 - (1 + k ** 2) * (b ** 2 - r ** 2) >= 0:
            x1 = (-k * b - math.sqrt(k ** 2 * b ** 2 - (1 + k ** 2) * (b ** 2 - r ** 2))) / (1 + k ** 2)
            x2 = (-k * b + math.sqrt(k ** 2 * b ** 2 - (1 + k ** 2) * (b ** 2 - r ** 2))) / (1 + k ** 2)
            xi = x1 if (x_a <= x1 <= x_b or x_b <= x1 <= x_a) else x2
            yi = k * xi + b
            return Point(xi + x_O, yi + y_O)
        else:
            return
    else:
        if -r < x_a < r:
            y1 = math.sqrt(r ** 2 - x_a ** 2)
            y2 = -math.sqrt(r ** 2 - x_a ** 2)
            yi = y1 if (y_b <= y1 <= y_a or y_a <= y2 <= y_b) else y2
            return Point(x_a + x_O, yi + y_O)
        else:
            return


def area(a, O, b):
    """
    judge the C/R area
    """
    x_O = O.o.x
    y_O = O.o.y
    x_a = a.x - x_O
    x_b = b.x - x_O
    y_a = a.y - y_O
    y_b = b.y - y_O
    r = O.r
    if x_b ** 2 + y_b ** 2 < r ** 2:
        return 1
    elif x_b * x_a + y_b * y_a < r ** 2:
        if intersection(a, O, b) != None:
            return 1
        else:
            return 0
    return 0


def GetoptPi(x, y, O):
    """
    optimize the point between other 2-adjacent two points
    """
    if O.isOnCircle(x, O) == 1:
        if area(x, O, y) == 1:
            return intersection(x, O, y)
        else:
            return PCP(x, y, O)
    else:
        if O.isOnCircle(y, O) == -1:
            return y
        else:
            return intersection(x, O, y)


def Length(l):
    """
    calculate the length of a traverse chain
    """
    sum_l = 0
    for m in range(len(l) - 1):
        sum_l += math.sqrt((l[m].x - l[m + 1].x) ** 2 + (l[m].y - l[m + 1].y) ** 2)
    return sum_l
