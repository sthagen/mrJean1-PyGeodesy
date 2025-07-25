
# -*- coding: utf-8 -*-

u'''(INTERNAL) Private ellipsoidal base classes C{CartesianEllipsoidalBase}
and C{LatLonEllipsoidalBase}.

A pure Python implementation of geodesy tools for ellipsoidal earth models,
transcoded in part from JavaScript originals by I{(C) Chris Veness 2005-2024}
and published under the same MIT Licence**, see for example U{latlon-ellipsoidal
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

# from pygeodesy.azimuthal import EquidistantExact, EquidistantKarney  # _MODS
from pygeodesy.basics import _isin, _xinstanceof
from pygeodesy.constants import EPS, EPS0, EPS1, _0_0, _0_5
from pygeodesy.cartesianBase import CartesianBase  # PYCHOK used!
# from pygeodesy.css import toCss  # _MODS
from pygeodesy.datums import Datum, Datums, _earth_ellipsoid, _ellipsoidal_datum, \
                             Transform, _WGS84,  _EWGS84  # _spherical_datum
# from pygeodesy.dms import parse3llh  # _MODS
# from pygeodesy.elevations import elevation2, geoidHeight2  # _MODS
# from pygeodesy.ellipsoidalBaseDI import _intersect3, _intersections2, _nearestOn2  # _MODS
# from pygeodesy.ellipsoids import _EWGS84  # from .datums
from pygeodesy.errors import _IsnotError, RangeError, _TypeError, _xattr, _xellipsoidal, \
                             _xellipsoids, _xError, _xkwds, _xkwds_not
# from pygeodesy.etm import etm, toEtm8  # _MODS
# from pygeodesy.fmath import favg  # _MODS
from pygeodesy.interns import NN, _COMMA_, _ellipsoidal_
from pygeodesy.latlonBase import LatLonBase, _trilaterate5,  fabs, _Wrap
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
# from pygeodesy.lcc import toLcc  # _MODS
# from pygeodesy.namedTuples import Vector3Tuple  # _MODS
# from pygeodesy.osgr import toOsgr  # _MODS
# from pygeodesy.points import isenclosedBy  # _MODS
from pygeodesy.props import deprecated_method, deprecated_property_RO, \
                            Property_RO, property_doc_, property_RO, _update_all
# from pygeodesy.trf import RefFrame, _toRefFrame  # _MODS
from pygeodesy.units import Epoch, _isDegrees, Radius_, _1mm as _TOL_M
# from pygeodesy import ups, utm, utmups  # MODS
# from pygeodesy.utmupsBase import _lowerleft  # MODS
# from pygeodesy.utily import _Wrap  # from .latlonBase
# from pygeodesy.vector3d import _intersects2  # _MODS

# from math import fabs  # from .latlonBase

__all__ = _ALL_LAZY.ellipsoidalBase
__version__ = '25.07.21'


class CartesianEllipsoidalBase(CartesianBase):
    '''(INTERNAL) Base class for ellipsoidal C{Cartesian}s.
    '''
    _datum   = _WGS84  # L{Datum}
    _epoch   =  None   # overriding .reframe.epoch (C{float})
    _reframe =  None   # reference frame (L{RefFrame})

    def __init__(self, x_xyz, y=None, z=None, reframe=None, epoch=None,
                                                       **datum_ll_name):
        '''New ellispoidal C{Cartesian...}.

           @kwarg reframe: Optional reference frame (L{RefFrame}).
           @kwarg epoch: Optional epoch to observe for B{C{reframe}} (C{scalar}),
                         a non-zero, fractional calendar year; silently ignored
                         if C{B{reframe}=None}.

           @raise TypeError: Non-scalar B{C{x_xyz}}, B{C{y}} or B{C{z}} coordinate
                             or B{C{x_xyz}} not a C{Cartesian} L{Ecef9Tuple},
                             L{Vector3Tuple} or L{Vector4Tuple} or B{C{datum}} is
                             not a L{Datum}, B{C{reframe}} is not a L{RefFrame} or
                             B{C{epoch}} is not C{scalar} non-zero.

           @see: Class L{CartesianBase<CartesianBase.__init__>} for more details.
        '''
        CartesianBase.__init__(self, x_xyz, y=y, z=z, **datum_ll_name)
        if reframe:
            self.reframe = reframe
            self.epoch = epoch

#   def __matmul__(self, other):  # PYCHOK Python 3.5+
#       '''Return C{NotImplemented} for C{c_ = c @ datum}, C{c_ = c @ reframe} and C{c_ = c @ Transform}.
#       '''
#       RefFrame = _MODS.trf.RefFrame
#       return NotImplemented if isinstance(other, (Datum, RefFrame, Transform)) else \
#             _NotImplemented(self, other)

    @deprecated_method
    def convertRefFrame(self, reframe2, reframe, epoch=None):
        '''DEPRECATED, use method L{toRefFrame}.'''
        return self.toRefFrame(reframe2, reframe=reframe, epoch=epoch)  # PYCHOK no cover

    @property_RO
    def ellipsoidalCartesian(self):
        '''Get this C{Cartesian}'s ellipsoidal class.
        '''
        return type(self)

    @property_doc_(''' this cartesian's observed or C{reframe} epoch (C{float}).''')
    def epoch(self):
        '''Get this cartesian's observed or C{reframe} epoch (C{Epoch}) or C{None}.
        '''
        return self._epoch or (self.reframe.epoch if self.reframe else None)

    @epoch.setter  # PYCHOK setter!
    def epoch(self, epoch):
        '''Set or clear this cartesian's observed epoch, a fractional
           calendar year (L{Epoch}, C{scalar} or C{str}) or C{None}.

           @raise TRFError: Invalid B{C{epoch}}.
        '''
        self._epoch = None if epoch is None else Epoch(epoch)

    def intersections2(self, radius, center2, radius2, sphere=True,
                                                       Vector=None, **Vector_kwds):
        '''Compute the intersection of two spheres or circles, each defined by a
           cartesian center point and a radius.

           @arg radius: Radius of this sphere or circle (same units as this point's
                        coordinates).
           @arg center2: Center of the second sphere or circle (C{Cartesian}, L{Vector3d},
                         C{Vector3Tuple} or C{Vector4Tuple}).
           @arg radius2: Radius of the second sphere or circle (same units as this and
                         the B{C{other}} point's coordinates).
           @kwarg sphere: If C{True}, compute the center and radius of the intersection
                          of two I{spheres}.  If C{False}, ignore the C{z}-component and
                          compute the intersection of two I{circles} (C{bool}).
           @kwarg Vector: Class to return intersections (C{Cartesian}, L{Vector3d} or
                          C{Vector3Tuple}) or C{None} for an instance of this (sub-)class.
           @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword arguments,
                               ignored if C{B{Vector} is None}.

           @return: If C{B{sphere} is True}, a 2-tuple of the C{center} and C{radius} of
                    the intersection of the I{spheres}.  The C{radius} is C{0.0} for
                    abutting spheres (and the C{center} is aka the I{radical center}).

                    If C{B{sphere} is False}, a 2-tuple with the two intersection points
                    of the I{circles}.  For abutting circles, both points are the same
                    instance (aka the I{radical center}).

           @raise IntersectionError: Concentric, invalid or non-intersecting spheres or circles.

           @raise TypeError: Invalid B{C{center2}}.

           @raise UnitError: Invalid B{C{radius}} or B{C{radius2}}.

           @see: U{Sphere-Sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>},
                 U{Circle-Circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>}
                 Intersection and function L{pygeodesy.radical2}.
        '''
        try:
            return _MODS.vector3d._intersects2(self, Radius_(radius=radius),
                                            center2, Radius_(radius2=radius2),
                                            sphere=sphere, clas=self.classof,
                                            Vector=Vector, **Vector_kwds)
        except (TypeError, ValueError) as x:
            raise _xError(x, center=self, radius=radius, center2=center2, radius2=radius2)

    @property_doc_(''' this cartesian's reference frame (L{RefFrame}).''')
    def reframe(self):
        '''Get this cartesian's reference frame (L{RefFrame}) or C{None}.
        '''
        return self._reframe

    @reframe.setter  # PYCHOK setter!
    def reframe(self, reframe):
        '''Set or clear this cartesian's reference frame (L{RefFrame}) or C{None}.

           @raise TypeError: The B{C{reframe}} is not a L{RefFrame}.
        '''
        _set_reframe(self, reframe)

    def toLatLon(self, datum=None, height=None, **LatLon_and_kwds):  # PYCHOK signature
        '''Convert this cartesian to a I{geodetic} (lat-/longitude) point.

           @see: Method L{toLatLon<cartesianBase.CartesianBase.toLatLon>}
                 for further details.
        '''
        kwds = LatLon_and_kwds
        if kwds:
            kwds = _xkwds(kwds, reframe=self.reframe, epoch=self.epoch)
        return CartesianBase.toLatLon(self, datum=datum, height=height, **kwds)

    def toRefFrame(self, reframe2, reframe=None, epoch=None, epoch2=None, **name):
        '''Convert this point to an other reference frame and epoch.

           @arg reframe2: Reference frame to convert I{to} (L{RefFrame}).
           @kwarg reframe: Optional reference frame to convert I{from} (L{RefFrame}),
                           overriding this point's reference frame.
           @kwarg epoch: Optional epoch (L{Epoch}, C{scalar} or C{str}), overriding
                         this point's C{epoch or B{reframe}.epoch}.
           @kwarg epoch2: Optional epoch to observe for the converted point (L{Epoch},
                          C{scalar} or C{str}), otherwise B{C{epoch}}.
           @kwarg name: Optional C{B{name}=NN} (C{str}), overriding C{B{reframe2}.name}.

           @return: The converted point (ellipsoidal C{Cartesian}) or if conversion
                    C{isunity}, this point or a copy of this point if the B{C{name}}
                    is non-empty.

           @raise TRFError: This point's C{reframe} is not defined, invalid B{C{epoch}}
                            or B{C{epoch2}} or conversion from this point's C{reframe}
                            to B{C{reframe2}} is not available.

           @raise TypeError: B{C{reframe2}} or B{C{reframe}} not a L{RefFrame}.
        '''
        return _MODS.trf._toRefFrame(self, reframe2, reframe=reframe, epoch=epoch,
                                                     epoch2=epoch2, **name)

    @deprecated_method
    def toTransforms_(self, *transforms, **datum):  # PYCHOK no cover
        '''DEPRECATED on 2024.02.14, use method C{toTransform}.'''
        r = self
        for t in transforms:
            r = r.toTransform(t)
        return r.dup(**datum) if datum else r


class LatLonEllipsoidalBase(LatLonBase):
    '''(INTERNAL) Base class for ellipsoidal C{LatLon}s.
    '''
    _datum          = _WGS84  # L{Datum}
    _elevation2to   =  None   # _elevation2 timeout (C{secs})
    _epoch          =  None   # overriding .reframe.epoch (C{float})
    _gamma          =  None   # UTM/UPS meridian convergence (C{degrees})
    _geoidHeight2to =  None   # _geoidHeight2 timeout (C{secs})
    _reframe        =  None   # reference frame (L{RefFrame})
    _scale          =  None   # UTM/UPS scale factor (C{float})
    _toLLEB_args    = ()      # Etm/Utm/Ups._toLLEB arguments

    def __init__(self, latlonh, lon=None, height=0, datum=_WGS84, reframe=None,
                                          epoch=None, wrap=False, **name):
        '''Create an ellipsoidal C{LatLon} point from the given lat-, longitude
           and height on the given datum, reference frame and epoch.

           @arg latlonh: Latitude (C{degrees} or DMS C{str} with N or S suffix) or
                         a previous C{LatLon} instance provided C{B{lon}=None}.
           @kwarg lon: Longitude (C{degrees} or DMS C{str} with E or W suffix) or C(None),
                       indicating B{C{latlonh}} is a C{LatLon}.
           @kwarg height: Optional height above (or below) the earth surface (C{meter},
                          same units as the datum's ellipsoid axes).
           @kwarg datum: Optional, ellipsoidal datum to use (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg reframe: Optional reference frame (L{RefFrame}).
           @kwarg epoch: Optional epoch to observe for B{C{reframe}} (C{scalar}), a
                         non-zero, fractional calendar year, but silently ignored if
                         C{B{reframe}=None}.
           @kwarg wrap: If C{True}, wrap or I{normalize} B{C{lat}} and B{C{lon}} (C{bool}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise RangeError: Value of C{lat} or B{C{lon}} outside the valid range and
                              L{rangerrors} set to C{True}.

           @raise TypeError: If B{C{latlonh}} is not a C{LatLon}, B{C{datum}} is not a
                             L{Datum}, B{C{reframe}} is not a L{RefFrame} or B{C{epoch}}
                             is not C{scalar} non-zero.

           @raise UnitError: Invalid B{C{lat}}, B{C{lon}} or B{C{height}}.
        '''
        LatLonBase.__init__(self, latlonh, lon=lon, height=height, wrap=wrap, **name)
        if not _isin(datum, None, self._datum, _EWGS84):
            self.datum = _ellipsoidal_datum(datum, name=self.name)
        if reframe:
            self.reframe = reframe
            self.epoch = epoch

#   def __matmul__(self, other):  # PYCHOK Python 3.5+
#       '''Return C{NotImplemented} for C{ll_ = ll @ datum} and C{ll_ = ll @ reframe}.
#       '''
#       RefFrame = _MODS.trf.RefFrame
#       return NotImplemented if isinstance(other, (Datum, RefFrame)) else \
#             _NotImplemented(self, other)

    def antipode(self, height=None):
        '''Return the antipode, the point diametrically opposite
           to this point.

           @kwarg height: Optional height of the antipode, height
                          of this point otherwise (C{meter}).

           @return: The antipodal point (C{LatLon}).
        '''
        lla = LatLonBase.antipode(self, height=height)
        if lla.datum != self.datum:
            lla.datum = self.datum
        return lla

    @deprecated_property_RO
    def convergence(self):
        '''DEPRECATED, use property C{gamma}.'''
        return self.gamma  # PYCHOK no cover

    @deprecated_method
    def convertDatum(self, datum2):
        '''DEPRECATED, use method L{toDatum}.'''
        return self.toDatum(datum2)

    @deprecated_method
    def convertRefFrame(self, reframe2):
        '''DEPRECATED, use method L{toRefFrame}.'''
        return self.toRefFrame(reframe2)

    @property_doc_(''' this points's datum (L{Datum}).''')
    def datum(self):
        '''Get this point's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this point's datum I{without conversion} (L{Datum}).

           @raise TypeError: The B{C{datum}} is not a L{Datum} or not ellipsoidal.
        '''
        _xinstanceof(Datum, datum=datum)
        if not datum.isEllipsoidal:
            raise _IsnotError(_ellipsoidal_, datum=datum)
        if self._datum != datum:
            _update_all(self)
            self._datum = datum

    def distanceTo2(self, other, wrap=False):
        '''I{Approximate} the distance and (initial) bearing between this
           and an other (ellipsoidal) point based on the radii of curvature.

           I{Suitable only for short distances up to a few hundred Km
           or Miles and only between points not near-polar}.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} the B{C{other}}
                        point (C{bool}).

           @return: An L{Distance2Tuple}C{(distance, initial)}.

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: Incompatible datum ellipsoids.

           @see: Method L{Ellipsoid.distance2} and U{Local, flat earth
                 approximation<https://www.EdWilliams.org/avform.htm#flat>}
                 aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>}
                 formula.
        '''
        p = self.others(other)
        if wrap:  # PYCHOK no cover
            p = _Wrap.point(p)
        E = self.ellipsoids(other)
        return E.distance2(*(self.latlon + p.latlon))

    @Property_RO
    def _elevation2(self):
        '''(INTERNAL) Get elevation and data source.
        '''
        return _MODS.elevations.elevation2(self.lat, self.lon,
                                           timeout=self._elevation2to)

    def elevation2(self, adjust=True, datum=None, timeout=2):
        '''Return elevation of this point for its or the given datum, ellipsoid
           or sphere.

           @kwarg adjust: Adjust the elevation for a B{C{datum}} other than
                          C{NAD83} (C{bool}).
           @kwarg datum: Optional datum overriding this point's datum (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2}, L{a_f2Tuple} or C{scalar}
                         radius).
           @kwarg timeout: Optional query timeout (C{seconds}).

           @return: An L{Elevation2Tuple}C{(elevation, data_source)} or
                    C{(None, error)} in case of errors.

           @note: The adjustment applied is the difference in geocentric earth
                  radius between the B{C{datum}} and C{NAV83} upon which the
                  L{elevations.elevation2} is based.

           @note: NED elevation is only available for locations within the U{Conterminous
                  US (CONUS)<https://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.elevation2} and method C{Ellipsoid.Rgeocentric}
                 for further details and possible C{error}s.
        '''
        if self._elevation2to != timeout:
            self._elevation2to = timeout
            LatLonEllipsoidalBase._elevation2._update(self)
        return self._Radjust2(adjust, datum, self._elevation2)

    def ellipsoid(self, datum=_WGS84):
        '''Return the ellipsoid of this point's datum or the given datum.

           @kwarg datum: Default datum (L{Datum}).

           @return: The ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).
        '''
        return _xattr(self, datum=datum).ellipsoid

    @property_RO
    def ellipsoidalLatLon(self):
        '''Get this C{LatLon}'s ellipsoidal class.
        '''
        return type(self)

    def ellipsoids(self, other):
        '''Check the type and ellipsoid of this and an other point's datum.

           @arg other: The other point (C{LatLon}).

           @return: This point's datum ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: Incompatible datum ellipsoids.
        '''
        self.others(other, up=2)  # ellipsoids' caller

        E = self.ellipsoid()
        try:  # other may be Sphere, etc.
            e = other.ellipsoid()
        except AttributeError:
            try:  # no ellipsoid method, try datum
                e = other.datum.ellipsoid
            except AttributeError:
                e = E  # no datum, XXX assume equivalent?
        return _xellipsoids(E, e)

    @property_doc_(''' this point's observed or C{reframe} epoch (C{float}).''')
    def epoch(self):
        '''Get this point's observed or C{reframe} epoch (L{Epoch}) or C{None}.
        '''
        return self._epoch or (self.reframe.epoch if self.reframe else None)

    @epoch.setter  # PYCHOK setter!
    def epoch(self, epoch):
        '''Set or clear this point's observed epoch, a fractional
           calendar year (L{Epoch}, C{scalar} or C{str}) or C{None}.

           @raise TRFError: Invalid B{C{epoch}}.
        '''
        self._epoch = None if epoch is None else Epoch(epoch)

    @Property_RO
    def Equidistant(self):
        '''Get the prefered azimuthal equidistant projection I{class} (L{EquidistantKarney} or L{EquidistantExact}).
        '''
        try:
            _ = self.datum.ellipsoid.geodesic
            return _MODS.azimuthal.EquidistantKarney
        except ImportError:  # no geographiclib
            return _MODS.azimuthal.EquidistantExact  # XXX no longer L{azimuthal.Equidistant}

    @Property_RO
    def _etm(self):
        '''(INTERNAL) Get this C{LatLon} point as an ETM coordinate (L{pygeodesy.toEtm8}).
        '''
        return self._toX8(_MODS.etm.toEtm8)

    @property_RO
    def gamma(self):
        '''Get this point's UTM or UPS meridian convergence (C{degrees}) or
           C{None} if not available or not converted from L{Utm} or L{Ups}.
        '''
        return self._gamma

    @Property_RO
    def _geoidHeight2(self):
        '''(INTERNAL) Get geoid height and model.
        '''
        return _MODS.elevations.geoidHeight2(self.lat, self.lon, model=0,
                                             timeout=self._geoidHeight2to)

    def geoidHeight2(self, adjust=False, datum=None, timeout=2):
        '''Return geoid height of this point for its or the given datum, ellipsoid
           or sphere.

           @kwarg adjust: Adjust the geoid height for a B{C{datum}} other than
                          C{NAD83/NADV88} (C{bool}).
           @kwarg datum: Optional datum overriding this point's datum (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2}, L{a_f2Tuple} or C{scalar}
                         radius).
           @kwarg timeout: Optional query timeout (C{seconds}).

           @return: A L{GeoidHeight2Tuple}C{(height, model_name)} or
                    C{(None, error)} in case of errors.

           @note: The adjustment applied is the difference in geocentric earth
                  radius between the B{C{datum}} and C{NAV83/NADV88} upon which
                  the L{elevations.geoidHeight2} is based.

           @note: The geoid height is only available for locations within the U{Conterminous
                  US (CONUS)<https://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.geoidHeight2} and method C{Ellipsoid.Rgeocentric}
                 for further details and possible C{error}s.
        '''
        if self._geoidHeight2to != timeout:
            self._geoidHeight2to = timeout
            LatLonEllipsoidalBase._geoidHeight2._update(self)
        return self._Radjust2(adjust, datum, self._geoidHeight2)

    def intermediateTo(self, other, fraction, height=None, wrap=False):  # PYCHOK no cover
        '''I{Must be overloaded}.'''
        self._notOverloaded(other, fraction, height=height, wrap=wrap)

    def intersection3(self, end1, start2, end2, height=None, wrap=False,  # was=True
                                          equidistant=None, tol=_TOL_M):
        '''I{Iteratively} compute the intersection point of two geodesic lines, each
           given as two points or as a start point and a bearing from North.

           @arg end1: End point of this line (C{LatLon}) or the initial bearing at
                      this point (compass C{degrees360}).
           @arg start2: Start point of the second line (this C{LatLon}).
           @arg end2: End point of the second line (this C{LatLon}) or the initial
                      bearing at B{C{start2}} (compass C{degrees360}).
           @kwarg height: Optional height at the intersection (C{meter}, conventionally)
                          or C{None} for the mean height.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{start2}} and
                        both B{C{end*}} points (C{bool}).
           @kwarg equidistant: An azimuthal equidistant projection (I{class} or function
                               L{pygeodesy.equidistant}), or C{None} for this point's
                               preferred C{.Equidistant}.
           @kwarg tol: Tolerance for convergence and skew line distance and length
                       (C{meter}, conventionally).

           @return: An L{Intersection3Tuple}C{(point, outside1, outside2)} with C{point}
                    a C{LatLon} instance.

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found, but only if
                               C{B{equidistant}=}L{EquidistantKarney}.

           @raise IntersectionError: Skew, colinear, parallel or otherwise non-intersecting
                                     lines or no convergence for the given B{C{tol}}.

           @raise TypeError: Invalid B{C{end1}}, B{C{start2}} or B{C{end2}}.

           @note: For each line specified with an initial bearing, a pseudo-end point is
                  computed as the C{destination} along that bearing at about 1.5 times the
                  distance from the start point to an initial gu-/estimate of the intersection
                  point (and between 1/8 and 3/8 of the C{authalic} earth perimeter).

           @see: I{Karney's} U{intersect.cpp<https://SourceForge.net/p/geographiclib/
                 discussion/1026621/thread/21aaff9f/>}, U{The B{ellipsoidal} case<https://
                 GIS.StackExchange.com/questions/48937/calculating-intersection-of-two-circles>}
                 and U{Karney's paper<https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section
                 B{14. MARITIME BOUNDARIES} for more details about the iteration algorithm.
        '''
        try:
            s2 = self.others(start2=start2)
            return _MODS.ellipsoidalBaseDI._intersect3(self, end1,
                                                       s2,   end2,
                                                       height=height, wrap=wrap,
                                                       equidistant=equidistant, tol=tol,
                                                       LatLon=self.classof, datum=self.datum)
        except (TypeError, ValueError) as x:
            raise _xError(x, start1=self, end1=end1, start2=start2, end2=end2,
                                          height=height, wrap=wrap, tol=tol)

    def intersections2(self, radius1, center2, radius2, height=None, wrap=False,  # was=True
                                                   equidistant=None, tol=_TOL_M):
        '''I{Iteratively} compute the intersection points of two circles, each
           defined by a center point and a radius.

           @arg radius1: Radius of this circle (C{meter}, conventionally).
           @arg center2: Center of the other circle (this C{LatLon}).
           @arg radius2: Radius of the other circle (C{meter}, same units as
                         B{C{radius1}}).
           @kwarg height: Optional height for the intersection points (C{meter},
                          conventionally) or C{None} for the I{"radical height"}
                          at the I{radical line} between both centers.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{center2}}
                        (C{bool}).
           @kwarg equidistant: An azimuthal equidistant projection (I{class} or
                               function L{pygeodesy.equidistant}) or C{None}
                               for this point's preferred C{.Equidistant}.
           @kwarg tol: Convergence tolerance (C{meter}, same units as
                       B{C{radius1}} and B{C{radius2}}).

           @return: 2-Tuple of the intersection points, each a C{LatLon}
                    instance.  For abutting circles, both intersection
                    points are the same instance, aka the I{radical center}.

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found, but only if
                               C{B{equidistant}=}L{EquidistantKarney}.

           @raise IntersectionError: Concentric, antipodal, invalid or
                                     non-intersecting circles or no
                                     convergence for the given B{C{tol}}.

           @raise TypeError: Invalid B{C{center2}} or B{C{equidistant}}.

           @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{height}}.

           @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
                 calculating-intersection-of-two-circles>}, U{Karney's paper
                 <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME BOUNDARIES},
                 U{circle-circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>} and
                 U{sphere-sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>}
                 intersections.
        '''
        try:
            c2 = self.others(center2=center2)
            return _MODS.ellipsoidalBaseDI._intersections2(self, radius1,
                                                           c2,   radius2,
                                                           height=height, wrap=wrap,
                                                           equidistant=equidistant, tol=tol,
                                                           LatLon=self.classof, datum=self.datum)
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, center=self, radius1=radius1, center2=center2, radius2=radius2,
                                                           height=height, wrap=wrap, tol=tol)

    def isenclosedBy(self, points, wrap=False):
        '''Check whether a polygon or composite encloses this point.

           @arg points: The polygon points or clips (C{LatLon}[],
                        L{BooleanFHP} or L{BooleanGH}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{points}} (C{bool}).

           @return: C{True} if this point is inside the polygon or composite,
                    C{False} otherwise.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not C{LatLon}.

           @raise ValueError: Invalid B{C{point}}, lat- or longitude.

           @see: Functions L{pygeodesy.isconvex}, L{pygeodesy.isenclosedBy}
                 and L{pygeodesy.ispolar} especially if the B{C{points}} may
                 enclose a pole or wrap around the earth I{longitudinally}.
        '''
        return _MODS.points.isenclosedBy(self, points, wrap=wrap)

    @property_RO
    def iteration(self):
        '''Get the most recent C{intersections2} or C{nearestOn} iteration
           number (C{int}) or C{None} if not available/applicable.
        '''
        return self._iteration

    def midpointTo(self, other, height=None, fraction=_0_5, wrap=False):
        '''Find the midpoint on a geodesic between this and an other point.

           @arg other: The other point (C{LatLon}).
           @kwarg height: Optional height for midpoint, overriding the
                          mean height (C{meter}).
           @kwarg fraction: Midpoint location from this point (C{scalar}),
                            may be negative or greater than 1.0.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: Midpoint (C{LatLon}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: Invalid B{C{height}}.

           @see: Methods C{intermediateTo} and C{rhumbMidpointTo}.
        '''
        return self.intermediateTo(other, fraction, height=height, wrap=wrap)

    def nearestOn(self, point1, point2, within=True, height=None, wrap=False,  # was=True
                                        equidistant=None, tol=_TOL_M):
        '''I{Iteratively} locate the closest point on the geodesic (line)
           between two other (ellipsoidal) points.

           @arg point1: Start point of the geodesic (C{LatLon}).
           @arg point2: End point of the geodesic (C{LatLon}).
           @kwarg within: If C{True}, return the closest point I{between} B{C{point1}} and
                          B{C{point2}}, otherwise the closest point elsewhere on the geodesic
                          (C{bool}).
           @kwarg height: Optional height for the closest point (C{meter}, conventionally)
                          or C{None} or C{False} for the interpolated height.  If C{False},
                          the closest takes the heights of the points into account.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll both B{C{point1}} and
                        B{C{point2}} (C{bool}).
           @kwarg equidistant: An azimuthal equidistant projection (I{class} or function
                               L{pygeodesy.equidistant}) or C{None} for this point's
                               preferred C{Equidistant}, like L{Equidistant}.
           @kwarg tol: Convergence tolerance (C{meter}, conventionally).

           @return: Closest point (C{LatLon}).

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found, but only if
                               C{B{equidistant}=}L{EquidistantKarney}.

           @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{equidistant}}.

           @raise ValueError: Datum or ellipsoid of B{C{point1}} or B{C{point2}} is incompatible
                              or no convergence for the given B{C{tol}}.

           @see: I{Karney}'s U{intercept.cpp<https://SourceForge.net/p/geographiclib/
                 discussion/1026621/thread/21aaff9f/>}, U{The B{ellipsoidal} case<https://
                 GIS.StackExchange.com/questions/48937/calculating-intersection-of-two-circles>}
                 and U{Karney's paper<https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section
                 B{14. MARITIME BOUNDARIES} for details about the iteration algorithm.
        '''
        try:
            t = _MODS.ellipsoidalBaseDI._nearestOn2(self, point1, point2, within=within,
                                                          height=height, wrap=wrap,
                                                          equidistant=equidistant,
                                                          tol=tol, LatLon=self.classof)
        except (TypeError, ValueError) as x:
            raise _xError(x, point=self, point1=point1, point2=point2, within=within,
                                         height=height, wrap=wrap, tol=tol)
        return t.closest

    def parse(self, strllh, height=0, datum=None, epoch=None, reframe=None,
                                        sep=_COMMA_, wrap=False, **name):
        '''Parse a string consisting of C{"lat, lon[, height]"},
           representing a similar, ellipsoidal C{LatLon} point.

           @arg strllh: Lat, lon and optional height (C{str}), see function
                        L{pygeodesy.parse3llh}.
           @kwarg height: Optional, default height (C{meter} or C{None}).
           @kwarg datum: Optional datum (L{Datum}), overriding this datum
                         I{without conversion}.
           @kwarg epoch: Optional datum (L{Epoch}), overriding this epoch
                         I{without conversion}.
           @kwarg reframe: Optional reference frame (L{RefFrame}), overriding
                           this reframe I{without conversion}.
           @kwarg sep: Optional separator (C{str}).
           @kwarg wrap: If C{True}, wrap or I{normalize} the lat- and
                        longitude (C{bool}).
           @kwarg name: Optional C{B{name}=NN} (C{str}), overriding this name.

           @return: The similar point (ellipsoidal C{LatLon}).

           @raise ParseError: Invalid B{C{strllh}}.
        '''
        a, b, h = _MODS.dms.parse3llh(strllh, height=height, sep=sep, wrap=wrap)
        return self.classof(a, b, height=h, datum=datum   or self.datum,
                                            epoch=epoch   or self.epoch,
                                          reframe=reframe or self.reframe, **name)

    def _Radjust2(self, adjust, datum, meter_text2):
        '''(INTERNAL) Adjust an C{elevation} or C{geoidHeight} with
           difference in Gaussian radii of curvature of the given
           datum and NAD83 ellipsoids at this point's latitude.

           @note: This is an arbitrary, possibly incorrect adjustment.
        '''
        if adjust:  # Elevation2Tuple or GeoidHeight2Tuple
            m, t = meter_text2
            if isinstance(m, float) and fabs(m) > EPS:  # PYCHOK no cover
                n = Datums.NAD83.ellipsoid.rocGauss(self.lat)
                if n > EPS0:
                    # use ratio, datum and NAD83 units may differ
                    E = self.ellipsoid() if _isin(datum, None, self.datum) else \
                      _earth_ellipsoid(datum)
                    r = E.rocGauss(self.lat)
                    if r > EPS0 and fabs(r - n) > EPS:  # EPS1
                        m *= r / n
                        meter_text2 = meter_text2.classof(m, t)
        return self._xnamed(meter_text2)

    @property_doc_(''' this point's reference frame (L{RefFrame}).''')
    def reframe(self):
        '''Get this point's reference frame (L{RefFrame}) or C{None}.
        '''
        return self._reframe

    @reframe.setter  # PYCHOK setter!
    def reframe(self, reframe):
        '''Set or clear this point's reference frame (L{RefFrame}) or C{None}.

           @raise TypeError: The B{C{reframe}} is not a L{RefFrame}.
        '''
        _set_reframe(self, reframe)

    @Property_RO
    def scale(self):
        '''Get this point's UTM grid or UPS point scale factor (C{float})
           or C{None} if not converted from L{Utm} or L{Ups}.
        '''
        return self._scale

    def _toX8(self, toX8, **kwds):
        '''(INTERNAL) Return toX8(self, ...).
        '''
        return toX8(self, **_xkwds(kwds, datum=self.datum, name=self.name))

    def toCartesian(self, height=None, **Cartesian_and_kwds):  # PYCHOK signature
        '''Convert this point to cartesian, I{geocentric} coordinates,
           also known as I{Earth-Centered, Earth-Fixed} (ECEF).

           @see: Method L{toCartesian<latlonBase.LatLonBase.toCartesian>}
                 for further details.
        '''
        kwds = Cartesian_and_kwds
        if kwds:
            kwds = _xkwds(kwds, reframe=self.reframe, epoch=self.epoch)
        return LatLonBase.toCartesian(self, height=height, **kwds)

    def toCss(self, **toCss_kwds):
        '''Convert this C{LatLon} point to a Cassini-Soldner location.

           @kwarg toCss_kwds: Optional keyword arguments for function
                              L{pygeodesy.toCss}.

           @return: The Cassini-Soldner location (L{Css}).
        '''
        return _MODS.css.toCss(self, **self._name1__(toCss_kwds))

    def toDatum(self, datum2, height=None, **name):
        '''Convert this point to an other datum.

           @arg datum2: Datum to convert I{to} (L{Datum}).
           @kwarg height: Optional height, overriding the
                          converted height (C{meter}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: The converted point (this C{LatLon}) or a copy
                    of this point if B{C{datum2}} matches this
                    point's C{datum}.

           @raise TypeError: Invalid B{C{datum2}}.
        '''
        n  =  self._name__(name)
        d2 = _ellipsoidal_datum(datum2, name=n)
        if self.datum == d2:
            r = self.copy(name=n)
        else:
            kwds = _xkwds_not(None, LatLon=self.classof, name=n,
                                    epoch=self.epoch, reframe=self.reframe)
            c = self.toCartesian().toDatum(d2)
            r = c.toLatLon(datum=d2, height=height, **kwds)
        return r

    def toEtm(self, **toEtm8_kwds):
        '''Convert this C{LatLon} point to an ETM coordinate.

           @kwarg toEtm8_kwds: Optional keyword arguments for
                               function L{pygeodesy.toEtm8}.

           @return: The ETM coordinate (L{Etm}).
        '''
        return self._etm if not toEtm8_kwds else \
               self._toX8(_MODS.etm.toEtm8, **toEtm8_kwds)

    def toLcc(self, **toLcc_kwds):
        '''Convert this C{LatLon} point to a Lambert location.

           @kwarg toLcc_kwds: Optional keyword arguments for
                              function L{pygeodesy.toLcc}.

           @return: The Lambert location (L{Lcc}).
        '''
        return _MODS.lcc.toLcc(self, **self._name1__(toLcc_kwds))

    def toMgrs(self, center=False, **toUtmUps_kwds):
        '''Convert this C{LatLon} point to an MGRS coordinate.

           @kwarg center: If C{True}, try to I{un}-center MGRS
                          to its C{lowerleft} (C{bool}) or by
                          C{B{center} meter} (C{scalar}).
           @kwarg toUtmUps_kwds: Optional keyword arguments for
                                 method L{toUtmUps}.

           @return: The MGRS coordinate (L{Mgrs}).

           @see: Methods L{toUtmUps} and L{toMgrs<pygeodesy.utmupsBase.UtmUpsBase.toMgrs>}.
        '''
        return self.toUtmUps(center=center, **toUtmUps_kwds).toMgrs(center=False)

    def toOsgr(self, kTM=False, **toOsgr_kwds):
        '''Convert this C{LatLon} point to an OSGR coordinate.

           @kwarg kTM: If C{True}, use I{Karney}'s Krüger method from module
                       L{ktm}, otherwise I{Ordinance Survery}'s recommended
                       formulation (C{bool}).
           @kwarg toOsgr_kwds: Optional keyword arguments for function
                               L{pygeodesy.toOsgr}.

           @return: The OSGR coordinate (L{Osgr}).
        '''
        return _MODS.osgr.toOsgr(self, kTM=kTM, **self._name1__(toOsgr_kwds))

    def toRefFrame(self, reframe2, reframe=None, epoch=None, epoch2=None, height=None, **name):
        '''Convert this point to an other reference frame and epoch.

           @arg reframe2: Reference frame to convert I{to} (L{RefFrame}).
           @kwarg reframe: Optional reference frame to convert I{from} (L{RefFrame}),
                           overriding this point's reference frame.
           @kwarg epoch: Optional epoch (L{Epoch}, C{scalar} or C{str}), overriding
                         this point's C{epoch or B{reframe}.epoch}.
           @kwarg epoch2: Optional epoch to observe for the converted point (L{Epoch},
                          C{scalar} or C{str}), otherwise B{C{epoch}}.
           @kwarg height: Optional height, overriding the converted height (C{meter}).
           @kwarg name: Optional C{B{name}=NN} (C{str}), overriding C{B{reframe2}.name}.

           @return: The converted point (ellipsoidal C{LatLon}) or if conversion
                    C{isunity}, this point or a copy of this point if the B{C{name}}
                    is non-empty.

           @raise TRFError: This point's C{reframe} is not defined, invalid B{C{epoch}}
                            or B{C{epoch2}} or conversion from this point's C{reframe}
                            to B{C{reframe2}} is not available.

           @raise TypeError: B{C{reframe2}} or B{C{reframe}} not a L{RefFrame}.
        '''
        return _MODS.trf._toRefFrame(self, reframe2, reframe=reframe, epoch=epoch,
                                           epoch2=epoch2, height=height, **name)

    def toTransform(self, transform, inverse=False, datum=None, **LatLon_kwds):
        '''Apply a Helmert transform to this geodetic point.

           @arg transform: Transform to apply (L{Transform} or L{TransformXform}).
           @kwarg inverse: Apply the inverse of the Helmert transform (C{bool}).
           @kwarg datum: Datum for the transformed point (L{Datum}), overriding
                         this point's datum but I{not} taken it into account.
           @kwarg LatLon_kwds: Optional keyword arguments for the transformed
                               point, like C{B{height}=...}.

           @return: A transformed point (C{LatLon}) or a copy of this point if
                    C{B{transform}.isunity}.

           @raise TypeError: Invalid B{C{transform}}.
        '''
        _xinstanceof(Transform, transform=transform)
        d = datum or self.datum
        if transform.isunity:
            r = self.dup(datum=d, **LatLon_kwds)
        else:
            c = self.toCartesian()
            c = c.toTransform(transform, inverse=inverse, datum=d)
            r = c.toLatLon(LatLon=self.classof, **_xkwds(LatLon_kwds, height=self.height))
        return r

    def toUps(self, center=False, **toUps8_kwds):
        '''Convert this C{LatLon} point to a UPS coordinate.

           @kwarg center: If C{True}, I{un}-center the UPS to its
                          C{lowerleft} (C{bool}) or by C{B{center}
                          meter} (C{scalar}).
           @kwarg toUps8_kwds: Optional keyword arguments for
                               function L{pygeodesy.toUps8}.

           @return: The UPS coordinate (L{Ups}).
        '''
        u = self._ups if (not toUps8_kwds) and self._upsOK() else \
            self._toX8(_MODS.ups.toUps8, **toUps8_kwds)
        return _lowerleft(u, center)

    def toUtm(self, center=False, **toUtm8_kwds):
        '''Convert this C{LatLon} point to a UTM coordinate.

           @kwarg center: If C{True}, I{un}-center the UTM to its
                          C{lowerleft} (C{bool}) or by C{B{center}
                          meter} (C{scalar}).
           @kwarg toUtm8_kwds: Optional keyword arguments for function
                               L{pygeodesy.toUtm8}.

           @return: The UTM coordinate (L{Utm}).

           @note: For the highest accuracy, use method L{toEtm} and
                  class L{pygeodesy.Etm} instead of L{pygeodesy.Utm}.
        '''
        u = self._utm if not toUtm8_kwds else \
            self._toX8(_MODS.utm.toUtm8, **toUtm8_kwds)
        return _lowerleft(u, center)

    def toUtmUps(self, pole=NN, center=False, **toUtmUps8_kwds):
        '''Convert this C{LatLon} point to a UTM or UPS coordinate.

           @kwarg pole: Optional top/center of UPS (stereographic)
                        projection (C{str}, 'N[orth]' or 'S[outh]').
           @kwarg center: If C{True}, I{un}-center the UTM or UPS to
                          its C{lowerleft} (C{bool}) or by C{B{center}
                          meter} (C{scalar}).
           @kwarg toUtmUps8_kwds: Optional keyword arguments for
                                  function L{pygeodesy.toUtmUps8}.

           @return: The UTM or UPS coordinate (L{Utm} or L{Ups}).
        '''
        x = not toUtmUps8_kwds
        if x and self._utmOK():
            u = self._utm
        elif x and self._upsOK(pole):
            u = self._ups
        else:  # no cover
            utmups = _MODS.utmups
            u = self._toX8(utmups.toUtmUps8, pole=pole, **toUtmUps8_kwds)
            if isinstance(u, utmups.Utm):
                self._update(False, _utm=u)  # PYCHOK kwds
            elif isinstance(u, utmups.Ups):
                self._update(False, _ups=u)  # PYCHOK kwds
            else:
                _xinstanceof(utmups.Utm, utmups.Ups, toUtmUps8=u)
        return _lowerleft(u, center)

    @deprecated_method
    def to3xyz(self):  # PYCHOK no cover
        '''DEPRECATED, use method C{toEcef}.

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @note: Overloads C{LatLonBase.to3xyz}
        '''
        r = self.toEcef()
        return _MODS.namedTuples.Vector3Tuple(r.x, r.y, r.z, name=self.name)

    def triangulate(self, bearing1, other, bearing2, **height_wrap_tol):
        '''I{Iteratively} locate a point given this, an other point and a bearing
           from North at each point.

           @arg bearing1: Bearing at this point (compass C{degrees360}).
           @arg other: The other point (C{LatLon}).
           @arg bearing2: Bearing at the B{C{other}} point (compass C{degrees360}).
           @kwarg height_wrap_tol: Optional keyword arguments C{B{height}=None},
                         C{B{wrap}=False} and C{B{tol}}, see method L{intersection3}.

           @return: Triangulated point (C{LatLon}).

           @see: Method L{intersection3} for further details.
        '''
        if _isDegrees(bearing1) and _isDegrees(bearing2):
            r = self.intersection3(bearing1, other, bearing2, **height_wrap_tol)
            return r.point
        raise _TypeError(bearing1=bearing1, bearing2=bearing2 **height_wrap_tol)

    def trilaterate5(self, distance1, point2, distance2, point3, distance3,
                           area=True, eps=EPS1, wrap=False):
        '''Trilaterate three points by I{area overlap} or I{perimeter intersection}
           of three intersecting circles.

           @arg distance1: Distance to this point (C{meter}), same units as B{C{eps}}).
           @arg point2: Second center point (C{LatLon}).
           @arg distance2: Distance to point2 (C{meter}, same units as B{C{eps}}).
           @arg point3: Third center point (C{LatLon}).
           @arg distance3: Distance to point3 (C{meter}, same units as B{C{eps}}).
           @kwarg area: If C{True}, compute the area overlap, otherwise the perimeter
                        intersection of the circles (C{bool}).
           @kwarg eps: The required I{minimal overlap} for C{B{area}=True} or the
                       I{intersection margin} for C{B{area}=False} (C{meter},
                       conventionally).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{point2}}
                        and B{C{point3}} (C{bool}).

           @return: A L{Trilaterate5Tuple}C{(min, minPoint, max, maxPoint, n)} with
                    C{min} and C{max} in C{meter}, same units as B{C{eps}}, the
                    corresponding trilaterated points C{minPoint} and C{maxPoint}
                    as I{ellipsoidal} C{LatLon} and C{n}, the number of trilatered
                    points found for the given B{C{eps}}.

                    If only a single trilaterated point is found, C{min I{is} max},
                    C{minPoint I{is} maxPoint} and C{n=1}.

                    If C{B{area}=False}, C{min} and C{max} represent the nearest
                    respectively farthest intersection margin.

                    If C{B{area}=True}, C{min} and C{max} are the smallest respectively
                    largest I{radial} overlap found.

                    If C{B{area}=True} and all 3 circles are concentric, C{n=0} and
                    C{minPoint} and C{maxPoint} are the B{C{point#}} with the smallest
                    B{C{distance#}} C{min} respectively largest B{C{distance#}} C{max}.

           @raise IntersectionError: Trilateration failed for the given B{C{eps}},
                                     insufficient overlap for C{B{area}=True}, no
                                     circle intersections for C{B{area}=False} or
                                     all circles are (near-)concentric.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @raise ValueError: Coincident B{C{points}} or invalid B{C{distance1}},
                              B{C{distance2}} or B{C{distance3}}.

           @note: Ellipsoidal trilateration invokes methods C{LatLon.intersections2}
                  and C{LatLon.nearestOn} based on I{Karney}'s Python U{geographiclib
                  <https://PyPI.org/project/geographiclib>} if installed, otherwise
                  the accurate (but slower) C{ellipsoidalExact.LatLon} methods.
        '''
        return _trilaterate5(self,                       distance1,
                             self.others(point2=point2), distance2,
                             self.others(point3=point3), distance3,
                             area=area, eps=eps, wrap=wrap)

    @Property_RO
    def _ups(self):  # __dict__ value overwritten by method C{toUtmUps}
        '''(INTERNAL) Get this C{LatLon} point as UPS coordinate (L{Ups}),
           see L{pygeodesy.toUps8}.
        '''
        return self._toX8(_MODS.ups.toUps8)  # pole=NN, falsed=True

    def _upsOK(self, pole=NN, falsed=True, **unused):
        '''(INTERNAL) Check matching C{Ups}.
        '''
        try:
            u = self._ups
        except RangeError:
            return False
        return falsed and (u.pole == pole[:1].upper() or not pole)

    @Property_RO
    def _utm(self):  # __dict__ value overwritten by method C{toUtmUps}
        '''(INTERNAL) Get this C{LatLon} point as UTM coordinate (L{Utm}),
           see L{pygeodesy.toUtm8}.
        '''
        return self._toX8(_MODS.utm.toUtm8)

    def _utmOK(self):
        '''(INTERNAL) Check C{Utm}.
        '''
        try:
            _ = self._utm
        except RangeError:
            return False
        return True


def _lowerleft(utmups, center):
    '''(INTERNAL) Optionally I{un}-center C{utmups}.
    '''
    if _isin(center, False, 0, _0_0):
        u = utmups
    elif _isin(center, True):
        u = utmups._lowerleft
    else:
        u = _MODS.utmupsBase._lowerleft(utmups, center)
    return u


def _nearestOn(point, point1, point2, within=True, height=None, wrap=False,  # was=True
                      equidistant=None, tol=_TOL_M, **LatLon_and_kwds):
    '''(INTERNAL) Get closest point, imported by .ellipsoidalExact,
       -GeodSolve, -Karney and -Vincenty to embellish exceptions.
    '''
    try:
        p = _xellipsoidal(point=point)
        t = _MODS.ellipsoidalBaseDI._nearestOn2(p, point1, point2, within=within,
                                                   height=height, wrap=wrap,
                                                   equidistant=equidistant,
                                                   tol=tol, **LatLon_and_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, point=point, point1=point1, point2=point2)
    return t.closest


def _set_reframe(inst, reframe):
    '''(INTERNAL) Set or clear an instance's reference frame.
    '''
    if reframe is not None:
        _xinstanceof(_MODS.trf.RefFrame, reframe=reframe)
        inst._reframe = reframe
    elif inst.reframe is not None:
        inst._reframe = None


__all__ += _ALL_DOCS(CartesianEllipsoidalBase, LatLonEllipsoidalBase)

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
