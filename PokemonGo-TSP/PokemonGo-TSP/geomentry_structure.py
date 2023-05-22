"""
data structure of geometric elements
"""


class Point(object):
    def __init__(self, x_value, y_value):
        super(Point, self).__init__()
        self.x = x_value
        self.y = y_value


class Line(object):
    def __init__(self, a, b, c):
        super(Line, self).__init__()
        self.a = a
        self.b = b
        self.c = c


class Circle(object):
    def __init__(self, point_o, r):
        super(Circle, self).__init__()
        self.o = point_o
        self.r = r

    def isOnCircle(self, point, circle):  # to judge the relation between circle and point
        px = point.x
        py = point.y
        cx = circle.o.x
        cy = circle.o.y
        r = circle.r
        if (px - cx) ** 2 + (py - cy) ** 2 > r * r:
            return 1
        elif (px - cx) ** 2 + (py - cy) ** 2 == r * r:
            return 0
        else:
            return -1
