
# -*- coding: utf-8 -*-

u'''Exact ellipsoidal geodesy, intended I{for testing purposes only}.

Ellipsoidal geodetic (lat-/longitude) L{LatLon} and geocentric
(ECEF) L{Cartesian} classes and functions L{areaOf}, L{intersections2},
L{isclockwise}, L{nearestOn} and L{perimeterOf} based on module
L{geodsolve}, a wrapper invoking I{Karney}'s U{GeodSolve
<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>} utility.
'''

# from pygeodesy.datums import _WGS84  # from .ellipsoidalBase
from pygeodesy.ellipsoidalBase import CartesianEllipsoidalBase, \
                                     _nearestOn, _WGS84
from pygeodesy.ellipsoidalBaseDI import LatLonEllipsoidalBaseDI, _TOL_M, \
                                       _intersection3, _intersections2
# from pygeodesy.errors import _xkwds  # from .karney
from pygeodesy.karney import fabs, _polygon,  Property_RO, _xkwds
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, _ALL_OTHER
from pygeodesy.points import _areaError, ispolar  # PYCHOK exported
# from pygeodesy.props import Property_RO  # from .karney

# from math import fabs  # from .karney

__all__ = _ALL_LAZY.ellipsoidalGeodSolve
__version__ = '25.05.28'


class Cartesian(CartesianEllipsoidalBase):
    '''Extended to convert exact L{Cartesian} to exact L{LatLon} points.
    '''

    def toLatLon(self, **LatLon_and_kwds):  # PYCHOK LatLon=LatLon, datum=None
        '''Convert this cartesian point to an exact geodetic point.

           @kwarg LatLon_and_kwds: Optional L{LatLon} and L{LatLon} keyword
                                   arguments as C{datum}.  Use C{B{LatLon}=...,
                                   B{datum}=...} to override this L{LatLon}
                                   class or specify C{B{LatLon}=None}.

           @return: The geodetic point (L{LatLon}) or if C{B{LatLon} is None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon_and_kwds}} argument.
        '''
        kwds = _xkwds(LatLon_and_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianEllipsoidalBase.toLatLon(self, **kwds)


class LatLon(LatLonEllipsoidalBaseDI):
    '''An ellipsoidal L{LatLon} like L{ellipsoidalKarney.LatLon} but using (exact)
       geodesic I{wrapper} L{GeodesicSolve} to compute the geodesic distance,
       initial and final bearing (azimuths) between two given points or the
       destination point given a start point and an (initial) bearing.
    '''

    @Property_RO
    def Equidistant(self):
        '''Get the prefered azimuthal equidistant projection I{class} (L{EquidistantGeodSolve}).
        '''
        return _MODS.azimuthal.EquidistantGeodSolve

    @Property_RO
    def geodesicx(self):
        '''Get this C{LatLon}'s (exact) geodesic (L{GeodesicSolve}).
        '''
        return self.datum.ellipsoid.geodsolve

    geodesic = geodesicx  # for C{._Direct} and C{._Inverse}

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to exact cartesian (ECEF) coordinates.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}} and other keyword
                  arguments, ignored if C{B{Cartesian} is None}.  Use C{B{Cartesian}=Class}
                  to override this L{Cartesian} class or set C{B{Cartesian}=None}.

           @return: The cartesian (ECEF) coordinates (L{Cartesian}) or if C{B{Cartesian} is
                    None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with
                    C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian}}, B{C{datum}} or other B{C{Cartesian_datum_kwds}}.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian, datum=self.datum)
        return LatLonEllipsoidalBaseDI.toCartesian(self, **kwds)


def areaOf(points, datum=_WGS84, wrap=True, polar=False):
    '''Compute the area of an (ellipsoidal) polygon or composite.

       @arg points: The polygon points (L{LatLon}[], L{BooleanFHP} or L{BooleanGH}).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the B{C{points}} (C{bool}).
       @kwarg polar: Use C{B{polar}=True} if the polygon encloses a pole (C{bool}), see
                     function L{ispolar<pygeodesy.points.ispolar>} and U{area of a polygon
                     enclosing a pole<https://GeographicLib.SourceForge.io/C++/doc/
                     classGeographicLib_1_1GeodesicExact.html#a3d7a9155e838a09a48dc14d0c3fac525>}.

       @return: Area (C{meter}, same as units of the B{C{datum}}'s ellipsoid axes, I{squared}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid C{B{wrap}=False}, unwrapped, unrolled longitudes not supported.

       @see: Functions L{pygeodesy.areaOf}, L{ellipsoidalExact.areaOf}, L{ellipsoidalKarney.areaOf},
             L{sphericalNvector.areaOf} and L{sphericalTrigonometry.areaOf}.
    '''
    return fabs(_polygon(datum.ellipsoid.geodsolve, points, True, False, wrap, polar))


def intersection3(start1, end1, start2, end2, height=None, wrap=False,  # was=True
                  equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''I{Iteratively} compute the intersection point of two lines, each defined
       by two (ellipsoidal) points or by an (ellipsoidal) start point and an
       (initial) bearing from North.

       @arg start1: Start point of the first line (L{LatLon}).
       @arg end1: End point of the first line (L{LatLon}) or the initial bearing
                  at the first point (compass C{degrees360}).
       @arg start2: Start point of the second line (L{LatLon}).
       @arg end2: End point of the second line (L{LatLon}) or the initial bearing
                  at the second point (compass C{degrees360}).
       @kwarg height: Optional height at the intersection (C{meter}, conventionally)
                      or C{None} for the mean height.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the B{C{start2}}
                    and B{C{end*}} points (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (I{class} or function
                           L{pygeodesy.equidistant}) or C{None} for the preferred
                           C{B{start1}.Equidistant}.
       @kwarg tol: Tolerance for convergence and for skew line distance and length
                   (C{meter}, conventionally).
       @kwarg LatLon: Optional class to return the intersection points (L{LatLon})
                      or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword arguments,
                           ignored if C{B{LatLon} is None}.

       @return: An L{Intersection3Tuple}C{(point, outside1, outside2)} with C{point}
                a B{C{LatLon}} or if C{B{LatLon} is None}, a L{LatLon4Tuple}C{(lat,
                lon, height, datum)}.

       @raise IntersectionError: Skew, colinear, parallel or otherwise
                                 non-intersecting lines or no convergence
                                 for the given B{C{tol}}.

       @raise TypeError: Invalid or non-ellipsoidal B{C{start1}}, B{C{end1}},
                         B{C{start2}} or B{C{end2}} or invalid B{C{equidistant}}.

       @note: For each line specified with an initial bearing, a pseudo-end point
              is computed as the C{destination} along that bearing at about 1.5
              times the distance from the start point to an initial gu-/estimate
              of the intersection point (and between 1/8 and 3/8 of the authalic
              earth perimeter).

       @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
             calculating-intersection-of-two-circles>} and U{Karney's paper
             <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME
             BOUNDARIES} for more details about the iteration algorithm.
    '''
    return _intersection3(start1, end1, start2, end2, height=height, wrap=wrap,
                          equidistant=equidistant, tol=tol, LatLon=LatLon, **LatLon_kwds)


def intersections2(center1, radius1, center2, radius2, height=None, wrap=False,  # was=True
                   equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''I{Iteratively} compute the intersection points of two circles, each defined
       by an (ellipsoidal) center point and a radius.

       @arg center1: Center of the first circle (L{LatLon}).
       @arg radius1: Radius of the first circle (C{meter}, conventionally).
       @arg center2: Center of the second circle (L{LatLon}).
       @arg radius2: Radius of the second circle (C{meter}, same units as
                     B{C{radius1}}).
       @kwarg height: Optional height for the intersection points (C{meter},
                      conventionally) or C{None} for the I{"radical height"}
                      at the I{radical line} between both centers.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{center2}}
                    (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (I{class} or
                           function L{pygeodesy.equidistant}) or C{None} for
                           the preferred C{B{center1}.Equidistant}.
       @kwarg tol: Convergence tolerance (C{meter}, same units as B{C{radius1}}
                   and B{C{radius2}}).
       @kwarg LatLon: Optional class to return the intersection points (L{LatLon})
                      or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword arguments,
                           ignored if C{B{LatLon} is None}.

       @return: 2-Tuple of the intersection points, each a B{C{LatLon}} instance
                or L{LatLon4Tuple}C{(lat, lon, height, datum)} if C{B{LatLon} is
                None}.  For abutting circles, both points are the same instance,
                aka the I{radical center}.

       @raise IntersectionError: Concentric, antipodal, invalid or non-intersecting
                                 circles or no convergence for the B{C{tol}}.

       @raise TypeError: Invalid or non-ellipsoidal B{C{center1}} or B{C{center2}}
                         or invalid B{C{equidistant}}.

       @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{height}}.

       @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
             calculating-intersection-of-two-circles>}, U{Karney's paper
             <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME BOUNDARIES},
             U{Circle-Circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>} and
             U{Sphere-Sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>}
             intersections.
    '''
    return _intersections2(center1, radius1, center2, radius2, height=height, wrap=wrap,
                           equidistant=equidistant, tol=tol, LatLon=LatLon, **LatLon_kwds)


def isclockwise(points, datum=_WGS84, wrap=True, polar=False):
    '''Determine the direction of a path or polygon.

       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the B{C{points}} (C{bool}).
       @kwarg polar: Use C{B{polar}=True} if the C{B{points}} enclose a pole (C{bool}),
                     see function U{ispolar<pygeodeys.points.ispolar>}.

       @return: C{True} if B{C{points}} are clockwise, C{False} otherwise.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: The B{C{points}} enclose a pole or zero area.

       @see: L{pygeodesy.isclockwise}.
    '''
    a = _polygon(datum.ellipsoid.geodsolve, points, True, False, wrap, polar)
    if a < 0:
        return True
    elif a > 0:
        return False
    raise _areaError(points)


def nearestOn(point, point1, point2, within=True, height=None, wrap=False,
              equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''I{Iteratively} locate the closest point on the geodesic between
       two other (ellipsoidal) points.

       @arg point: Reference point (C{LatLon}).
       @arg point1: Start point of the geodesic (C{LatLon}).
       @arg point2: End point of the geodesic (C{LatLon}).
       @kwarg within: If C{True}, return the closest point I{between}
                      B{C{point1}} and B{C{point2}}, otherwise the
                      closest point elsewhere on the geodesic (C{bool}).
       @kwarg height: Optional height for the closest point (C{meter},
                      conventionally) or C{None} or C{False} for the
                      interpolated height.  If C{False}, the closest
                      takes the heights of the points into account.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll both
                    B{C{point1}} and B{C{point2}} (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (I{class}
                           or function L{pygeodesy.equidistant}) or C{None}
                           for the preferred C{B{point}.Equidistant}.
       @kwarg tol: Convergence tolerance (C{meter}).
       @kwarg LatLon: Optional class to return the closest point
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: Closest point, a B{C{LatLon}} instance or if C{B{LatLon}
                is None}, a L{LatLon4Tuple}C{(lat, lon, height, datum)}.

       @raise TypeError: Invalid or non-ellipsoidal B{C{point}}, B{C{point1}}
                         or B{C{point2}} or invalid B{C{equidistant}}.

       @raise ValueError: No convergence for the B{C{tol}}.

       @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
             calculating-intersection-of-two-circles>} and U{Karney's paper
             <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME
             BOUNDARIES} for more details about the iteration algorithm.
    '''
    return _nearestOn(point, point1, point2, within=within, height=height, wrap=wrap,
                      equidistant=equidistant, tol=tol, LatLon=LatLon, **LatLon_kwds)


def perimeterOf(points, closed=False, datum=_WGS84, wrap=True):
    '''Compute the perimeter of an (ellipsoidal) polygon or composite.

       @arg points: The polygon points (L{LatLon}[], L{BooleanFHP} or L{BooleanGH}).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the B{C{points}} (C{bool}).

       @return: Perimeter (C{meter}, same as units of the B{C{datum}}'s ellipsoid axes).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid C{B{wrap}=False}, unwrapped, unrolled longitudes not
                          supported or C{B{closed}=False} with C{B{points}} a composite.

       @see: Functions L{pygeodesy.perimeterOf}, L{ellipsoidalExact.perimeterOf},
             L{ellipsoidalKarney.perimeterOf}, L{sphericalNvector.perimeterOf}
             and L{sphericalTrigonometry.perimeterOf}.
    '''
    return _polygon(datum.ellipsoid.geodsolve, points, closed, True, wrap, False)


__all__ += _ALL_OTHER(Cartesian, LatLon,  # classes
                      areaOf,  # functions
                      intersection3, intersections2, isclockwise, ispolar,
                      nearestOn, perimeterOf)

# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
