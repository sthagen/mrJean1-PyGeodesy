
# -*- coding: utf-8 -*-

u'''Formulary of basic geodesy functions and approximations.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.basics import _copysign, _isin  # _args_kwds_count2
# from pygeodesy.cartesianBase import CartesianBase  # _MODS
from pygeodesy.constants import EPS, EPS0, EPS1, PI, PI2, PI3, PI_2, R_M, \
                               _0_0s, float0_, isnon0, remainder, _umod_PI2, \
                               _0_0, _0_125, _0_25, _0_5, _1_0, _2_0, _4_0, \
                               _90_0, _180_0, _360_0
from pygeodesy.datums import Datum, Ellipsoid, _ellipsoidal_datum, \
                            _mean_radius, _spherical_datum, _WGS84,  _EWGS84
# from pygeodesy.ellipsoids import Ellipsoid, _EWGS84  # from .datums
from pygeodesy.errors import IntersectionError, LimitError, limiterrors, \
                            _TypeError, _ValueError, _xattr, _xError, \
                            _xcallable, _xkwds, _xkwds_pop2
from pygeodesy.fmath import euclid, fdot_, fprod, hypot, hypot2, sqrt0
from pygeodesy.fsums import fsumf_,  Fmt, unstr
# from pygeodesy.internals import typename  # from .named
from pygeodesy.interns import _delta_, _distant_, _inside_, _SPACE_, _too_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _name__, _name2__, _NamedTuple, _xnamed,  typename
from pygeodesy.namedTuples import Bearing2Tuple, Distance4Tuple, LatLon2Tuple, \
                                  Intersection3Tuple, PhiLam2Tuple
# from pygeodesy.streprs import Fmt, unstr  # from .fsums
# from pygeodesy.triaxials import _hartzell3  # _MODS
from pygeodesy.units import _isDegrees, _isHeight, _isRadius, Bearing, Degrees_, \
                             Distance, Distance_, Height, Lamd, Lat, Lon, Meter_, \
                             Phid, Radians, Radians_, Radius, Radius_, Scalar, _100km
from pygeodesy.utily import acos1, asin1, atan2, atan2b, degrees2m, hav, _loneg, \
                            m2degrees, tan_2, sincos2, sincos2_, _Wrap
# from pygeodesy.vector3d import _otherV3d  # _MODS
# from pygeodesy.vector3dBase import _xyz_y_z3  # _MODS
# from pygeodesy import ellipsoidalExact, ellipsoidalKarney, vector3d, \
#                       sphericalNvector, sphericalTrigonometry  # _MODS

from contextlib import contextmanager
from math import atan, cos, degrees, fabs, radians, sin, sqrt  # pow

__all__ = _ALL_LAZY.formy
__version__ = '25.05.12'

_RADIANS2 =  radians(_1_0)**2  # degree to radians-squared
_ratio_   = 'ratio'
_xline_   = 'xline'


def angle2chord(rad, radius=R_M):
    '''Get the chord length of a (central) angle or I{angular} distance.

       @arg rad: Central angle (C{radians}).
       @kwarg radius: Mean earth radius (C{meter}, conventionally), datum (L{Datum}) or ellipsoid
                      (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) to use or C{None}.

       @return: Chord length (C{meter}, same units as B{C{radius}} or if C{B{radius} is None}, C{radians}).

       @see: Function L{chord2angle}, method L{intermediateChordTo<sphericalNvector.LatLon.intermediateChordTo>} and
             U{great-circle-distance<https://WikiPedia.org/wiki/Great-circle_distance#Relation_between_central_angle_and_chord_length>}.
    '''
    d = _isDegrees(rad, iscalar=False)
    r =  sin((radians(rad) if d else rad) / _2_0) * _2_0
    return (degrees(r) if d else r) if radius is None else (_mean_radius(radius) * r)


def _anti2(a, b, n_2, n, n2):
    '''(INTERNAL) Helper for C{antipode} and C{antipode_}.
    '''
    r = remainder(a, n) if fabs(a) > n_2 else a
    if r == a:
        r = -r
        b += n
    if fabs(b) > n:
        b = remainder(b, n2)
    return float0_(r, b)


def antipode(lat, lon, **name):
    '''Return the antipode, the point diametrically opposite to a given
       point in C{degrees}.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg name: Optional C{B{name}=NN} (C{str}).

       @return: A L{LatLon2Tuple}C{(lat, lon)}.

       @see: Functions L{antipode_} and L{normal} and U{Geosphere
             <https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return LatLon2Tuple(*_anti2(lat, lon, _90_0, _180_0, _360_0), **name)


def antipode_(phi, lam, **name):
    '''Return the antipode, the point diametrically opposite to a given
       point in C{radians}.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).
       @kwarg name: Optional C{B{name}=NN} (C{str}).

       @return: A L{PhiLam2Tuple}C{(phi, lam)}.

       @see: Functions L{antipode} and L{normal_} and U{Geosphere
             <https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return PhiLam2Tuple(*_anti2(phi, lam, PI_2, PI, PI2), **name)


def bearing(lat1, lon1, lat2, lon2, **final_wrap):
    '''Compute the initial or final bearing (forward or reverse azimuth) between two
       (spherical) points.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg final_wrap: Optional keyword arguments for function L{pygeodesy.bearing_}.

       @return: Initial or final bearing (compass C{degrees360}) or zero if both points
                coincide.
    '''
    r = bearing_(Phid(lat1=lat1), Lamd(lon1=lon1),
                 Phid(lat2=lat2), Lamd(lon2=lon2), **final_wrap)
    return degrees(r)


def bearing_(phi1, lam1, phi2, lam2, final=False, wrap=False):
    '''Compute the initial or final bearing (forward or reverse azimuth) between two
       (spherical) points.

       @arg phi1: Start latitude (C{radians}).
       @arg lam1: Start longitude (C{radians}).
       @arg phi2: End latitude (C{radians}).
       @arg lam2: End longitude (C{radians}).
       @kwarg final: If C{True}, return the final, otherwise the initial bearing (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{phi2}} and B{C{lam2}}
                    (C{bool}).

       @return: Initial or final bearing (compass C{radiansPI2}) or zero if both points
                coincide.

       @see: U{Bearing<https://www.Movable-Type.co.UK/scripts/latlong.html>}, U{Course
             between two points<https://www.EdWilliams.org/avform147.htm#Crs>} and
             U{Bearing Between Two Points<https://web.Archive.org/web/20020630205931/
             https://MathForum.org/library/drmath/view/55417.html>}.
    '''
    db, phi2, lam2 = _Wrap.philam3(lam1, phi2, lam2, wrap)
    if final:  # swap plus PI
        phi1, lam1, phi2, lam2, db = phi2, lam2, phi1, lam1, -db
        r = PI3
    else:
        r = PI2
    sa1, ca1, sa2, ca2, sdb, cdb = sincos2_(phi1, phi2, db)

    x = ca1 * sa2 - sa1 * ca2 * cdb
    y = sdb * ca2
    return _umod_PI2(atan2(y, x) + r)  # .utily.wrapPI2


def _bearingTo2(p1, p2, wrap=False):  # for points.ispolar, sphericalTrigonometry.areaOf
    '''(INTERNAL) Compute initial and final bearing.
    '''
    try:  # for LatLon_ and ellipsoidal LatLon
        return p1.bearingTo2(p2, wrap=wrap)
    except AttributeError:
        pass
    # XXX spherical version, OK for ellipsoidal ispolar?
    t = p1.philam + p2.philam
    i = bearing_(*t, final=False, wrap=wrap)
    f = bearing_(*t, final=True,  wrap=wrap)
    return Bearing2Tuple(degrees(i), degrees(f),
                         name__=_bearingTo2)


def chord2angle(chord, radius=R_M):
    '''Get the (central) angle from a chord length or distance.

       @arg chord: Length or distance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}, conventionally), datum (L{Datum}) or
                      ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) to use.

       @return: Angle (C{radians} with sign of B{C{chord}}) or C{0} if C{B{radius}=0}.

       @note: The angle will exceed C{PI} if C{B{chord} > B{radius} * 2}.

       @see: Function L{angle2chord}.
    '''
    m = _mean_radius(radius)
    r =  fabs(chord / (m * _2_0)) if m > 0 else _0_0
    if r:
        i = int(r)
        if i > 0:
            r -= i
            i *= PI
        r = (asin1(r) + i) * _2_0
    return _copysign(r, chord)


def compassAngle(lat1, lon1, lat2, lon2, adjust=True, wrap=False):
    '''Return the angle from North for the direction vector M{(lon2 - lon1,
       lat2 - lat1)} between two points.

       Suitable only for short, not near-polar vectors up to a few hundred
       Km or Miles.  Use function L{pygeodesy.bearing} for longer vectors.

       @arg lat1: From latitude (C{degrees}).
       @arg lon1: From longitude (C{degrees}).
       @arg lat2: To latitude (C{degrees}).
       @arg lon2: To longitude (C{degrees}).
       @kwarg adjust: Adjust the longitudinal delta by the cosine of the mean
                      latitude (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and
                    B{C{lon2}} (C{bool}).

       @return: Compass angle from North (C{degrees360}).

       @note: Courtesy of Martin Schultz.

       @see: U{Local, flat earth approximation
             <https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    d_lon, lat2, lon2 = _Wrap.latlon3(lon1, lat2, lon2, wrap)
    if adjust:  # scale delta lon
        d_lon *= _scale_deg(lat1, lat2)
    return atan2b(d_lon, lat2 - lat1)


def cosineLaw(lat1, lon1, lat2, lon2, corr=0, earth=None, wrap=False,
                                              datum=_WGS84, radius=R_M):
    '''Compute the distance between two points using the U{Law of Cosines
       <https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       formula, optionally corrected.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg corr: Use C{B{corr}=2} to apply the U{Forsythe-Andoyer-Lambert
                    <https://www2.UNB.CA/gge/Pubs/TR77.pdf>}, C{B{corr}=1} for the
                    U{Andoyer-Lambert<https://Books.Google.com/books?id=x2UiAQAAIAAJ>}
                    corrected (ellipsoidal) or keep C{B{corr}=0} for the uncorrected
                    (spherical) C{Law of Cosines} formula (C{int}).
       @kwarg earth: Mean earth radius (C{meter}) or datum (L{Datum}) or ellipsoid
                     (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and B{C{lat2}} and B{C{lon2}}
                    (C{bool}).
       @kwarg datum: Default ellipsiodal B{C{earth}} (and for backward compatibility).
       @kwarg radius: Default spherical B{C{earth}} (and for backward compatibility).

       @return: Distance (C{meter}, same units as B{C{radius}} or the datum's or
                ellipsoid axes).

       @raise TypeError: Invalid B{C{earth}}, B{C{datum}} or B{C{radius}}.

       @raise ValueError: Invalid B{C{corr}}.

       @see: Functions L{cosineLaw_}, L{equirectangular}, L{euclidean}, L{flatLocal} /
             L{hubeny}, L{flatPolar}, L{haversine}, L{thomas} and L{vincentys} and
             method L{Ellipsoid.distance2}.

       @note: See note at function L{vincentys_}.
    '''
    return _dE(cosineLaw_, earth or datum,  wrap, lat1, lon1, lat2, lon2, corr=corr) if corr else \
           _dS(cosineLaw_, earth or radius, wrap, lat1, lon1, lat2, lon2)


def cosineLaw_(phi2, phi1, lam21, corr=0, earth=None, datum=_WGS84):
    '''Compute the I{angular} distance between two points using the U{Law of Cosines
       <https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>} formula,
       optionally corrected.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg corr: Use C{B{corr}=2} to apply the U{Forsythe-Andoyer-Lambert
                    <https://www2.UNB.CA/gge/Pubs/TR77.pdf>}, C{B{corr}=1} for the
                    U{Andoyer-Lambert<https://Books.Google.com/books?id=x2UiAQAAIAAJ>}
                    corrected (ellipsoidal) or keep C{B{corr}=0} for the uncorrected
                    (spherical) C{Law of Cosines} formula (C{int}).
       @kwarg earth: Mean earth radius (C{meter}) or datum (L{Datum}) or ellipsoid
                     (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg datum: Default ellipsoidal B{C{earth}} (and for backward compatibility).

       @return: Angular distance (C{radians}).

       @raise TypeError: Invalid B{C{earth}} or B{C{datum}}.

       @raise ValueError: Invalid B{C{corr}}.

       @see: Functions L{cosineLaw}, L{euclidean_}, L{flatLocal_} / L{hubeny_},
             L{flatPolar_}, L{haversine_}, L{thomas_} and L{vincentys_} and U{Geodesy-PHP
             <https://GitHub.com/jtejido/geodesy-php/blob/master/src/Geodesy/Distance/
             AndoyerLambert.php>}.
    '''
    s2, c2, s1, c1, r, c21 = _sincosa6(phi2, phi1, lam21)
    if corr and isnon0(c1) and isnon0(c2):
        E = _ellipsoidal(earth or datum, cosineLaw_)
        f = _0_25 * E.f
        if f:  # ellipsoidal
            if corr == 1:  # Andoyer-Lambert
                r2 = atan2(E.b_a * s2, c2)
                r1 = atan2(E.b_a * s1, c1)
                s2, c2, s1, c1 = sincos2_(r2, r1)
                r = acos1(s1 * s2 + c1 * c2 * c21)
                if r:
                    sr, _, sr_2, cr_2 = sincos2_(r, r * _0_5)
                    if isnon0(sr_2) and isnon0(cr_2):
                        s  = (sr + r) * ((s1 - s2) / sr_2)**2
                        c  = (sr - r) * ((s1 + s2) / cr_2)**2
                        r += (c  - s) * _0_5 * f

            elif corr == 2:  # Forsythe-Andoyer-Lambert
                sr, cr, s2r, _ = sincos2_(r, r * 2)
                if isnon0(sr) and fabs(cr) < EPS1:
                    s = (s1 + s2)**2 / (_1_0 + cr)
                    t = (s1 - s2)**2 / (_1_0 - cr)
                    x = s + t
                    y = s - t

                    s =  8 * r**2 / sr
                    a = 64 * r  + s * cr * 2  # 16 * r**2 / tan(r)
                    d = 48 * sr + s  # 8 * r**2 / tan(r)
                    b = -2 * d
                    e = 30 * s2r

                    c  = fdot_(30, r,  cr, s,  e, _0_5)  # 8 * r**2 / tan(r)
                    t  = fdot_( a, x,  b,  y,  e,  y**2,  -c, x**2,  d, x * y) * _0_125
                    r += fdot_(-r, x,  sr, y * 3,  t, f) * f
            else:
                raise _ValueError(corr=corr)
    return r


def _d3(wrap, lat1, lon1, lat2, lon2):
    '''(INTERNAL) Helper for _dE, _dS, ....
    '''
    if wrap:
        d_lon, lat2, _ = _Wrap.latlon3(lon1, lat2, lon2, wrap)
        return radians(lat2), Phid(lat1=lat1), radians(d_lon)
    else:  # for backward compaibility
        return Phid(lat2=lat2), Phid(lat1=lat1), radians(lon2 - lon1)


def _dE(fun_, earth, wrap, *lls, **corr):
    '''(INTERNAL) Helper for ellipsoidal distances.
    '''
    E = _ellipsoidal(earth, fun_)
    r =  fun_(*_d3(wrap, *lls), datum=E, **corr)
    return r * E.a


def _dS(fun_, radius, wrap, *lls, **adjust):
    '''(INTERNAL) Helper for spherical distances.
    '''
    r = fun_(*_d3(wrap, *lls), **adjust)
    if radius is not R_M:
        try:  # datum?
            radius = radius.ellipsoid.R1
        except AttributeError:
            pass  # scalar?
        lat1, _, lat2, _ = lls
        radius = _mean_radius(radius, lat1, lat2)
    return r * radius


def _ellipsoidal(earth, where):
    '''(INTERNAL) Helper for distances.
    '''
    return _EWGS84 if _isin(earth, _EWGS84, _WGS84) else (
             earth if  isinstance(earth, Ellipsoid) else
            (earth if  isinstance(earth, Datum)     else  # PYCHOK indent
            _ellipsoidal_datum(earth, name__=where)).ellipsoid)


def equirectangular(lat1, lon1, lat2, lon2, radius=R_M, **adjust_limit_wrap):
    '''Approximate the distance between two points using the U{Equirectangular Approximation
       / Projection<https://www.Movable-Type.co.UK/scripts/latlong.html#equirectangular>}.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum}) or ellipsoid
                      (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg adjust_limit_wrap: Optionally, keyword arguments for function L{equirectangular4}.

       @return: Distance (C{meter}, same units as B{C{radius}} or the datum's
                ellipsoid axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: Function L{equirectangular4} for more details, the available B{C{options}},
             errors, restrictions and other, approximate or accurate distance functions.
    '''
    r = _mean_radius(radius, lat1, lat2)
    t =  equirectangular4(Lat(lat1=lat1), Lon(lon1=lon1),
                          Lat(lat2=lat2), Lon(lon2=lon2),
                        **adjust_limit_wrap)  # PYCHOK 4 vs 2-3
    return degrees2m(sqrt(t.distance2), radius=r)


def _equirectangular(lat1, lon1, lat2, lon2, **adjust_limit_wrap):
    '''(INTERNAL) Helper for classes L{frechet._FrechetMeterRadians} and
       L{hausdorff._HausdorffMeterRedians}.
    '''
    t = equirectangular4(lat1, lon1, lat2, lon2, **adjust_limit_wrap)
    return t.distance2 * _RADIANS2


def equirectangular4(lat1, lon1, lat2, lon2, adjust=True, limit=45, wrap=False):
    '''Approximate the distance between two points using the U{Equirectangular Approximation
       / Projection<https://www.Movable-Type.co.UK/scripts/latlong.html#equirectangular>}.

       This approximation is valid for short distance of several hundred Km or Miles, see
       the B{C{limit}} keyword argument and L{LimitError}.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta by the cosine of the
                      mean latitude (C{bool}).
       @kwarg limit: Optional limit for lat- and longitudinal deltas (C{degrees}) or C{None}
                     or C{0} for unlimited.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and B{C{lon2}}
                    (C{bool}).

       @return: A L{Distance4Tuple}C{(distance2, delta_lat, delta_lon, unroll_lon2)} with
                C{distance2} in C{degrees squared}.

       @raise LimitError: The lat- or longitudinal delta exceeds the B{C{-limit..limit}}
                          range and L{limiterrors<pygeodesy.limiterrors>} is C{True}.

       @see: U{Local, flat earth approximation<https://www.EdWilliams.org/avform.htm#flat>},
             functions L{equirectangular}, L{cosineLaw}, L{euclidean}, L{flatLocal} /
             L{hubeny}, L{flatPolar}, L{haversine}, L{thomas} and L{vincentys} and methods
             L{Ellipsoid.distance2}, C{LatLon.distanceTo*} and C{LatLon.equirectangularTo}.
    '''
    if wrap:
        d_lon, lat2, ulon2 = _Wrap.latlon3(lon1, lat2, lon2, wrap)
    else:
        d_lon,       ulon2 = (lon2 - lon1), lon2
    d_lat = lat2 - lat1

    if limit and limit > 0 and limiterrors():
        d = max(fabs(d_lat), fabs(d_lon))
        if d > limit:
            t = _SPACE_(_delta_, Fmt.PAREN_g(d), Fmt.exceeds_limit(limit))
            s =  unstr(equirectangular4, lat1, lon1, lat2, lon2,
                                         limit=limit, wrap=wrap)
            raise LimitError(s, txt=t)

    if adjust:  # scale delta lon
        d_lon *= _scale_deg(lat1, lat2)

    d2 = hypot2(d_lat, d_lon)  # degrees squared!
    return Distance4Tuple(d2, d_lat, d_lon, ulon2 - lon2)


def euclidean(lat1, lon1, lat2, lon2, radius=R_M, adjust=True, wrap=False):
    '''Approximate the C{Euclidean} distance between two (spherical) points.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum}) or ellipsoid
                      (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg adjust: Adjust the longitudinal delta by the cosine of the mean
                      latitude (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and
                    B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}} or the ellipsoid
                or datum axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: U{Distance between two (spherical) points
             <https://www.EdWilliams.org/avform.htm#Dist>}, functions L{euclid},
             L{euclidean_}, L{cosineLaw}, L{equirectangular}, L{flatLocal} /
             L{hubeny}, L{flatPolar}, L{haversine}, L{thomas} and L{vincentys}
             and methods L{Ellipsoid.distance2}, C{LatLon.distanceTo*} and
             C{LatLon.equirectangularTo}.
    '''
    return _dS(euclidean_, radius, wrap, lat1, lon1, lat2, lon2, adjust=adjust)


def euclidean_(phi2, phi1, lam21, adjust=True):
    '''Approximate the I{angular} C{Euclidean} distance between two (spherical) points.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg adjust: Adjust the longitudinal delta by the cosine of the mean
                      latitude (C{bool}).

       @return: Angular distance (C{radians}).

       @see: Functions L{euclid}, L{euclidean}, L{cosineLaw_}, L{flatLocal_} /
             L{hubeny_}, L{flatPolar_}, L{haversine_}, L{thomas_} and L{vincentys_}.
    '''
    if adjust:
        lam21 *= _scale_rad(phi2, phi1)
    return euclid(phi2 - phi1, lam21)


def excessAbc_(A, b, c):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle from two sides
       and the included (small) angle.

       @arg A: An interior triangle angle (C{radians}).
       @arg b: Frist adjacent triangle side (C{radians}).
       @arg c: Second adjacent triangle side (C{radians}).

       @return: Spherical excess (C{radians}).

       @raise UnitError: Invalid B{C{A}}, B{C{b}} or B{C{c}}.

       @see: Functions L{excessGirard_}, L{excessLHuilier_} and U{Spherical
             trigonometry<https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    A = Radians_(A=A)
    b = Radians_(b=b) * _0_5
    c = Radians_(c=c) * _0_5

    sA, cA, sb, cb, sc, cc = sincos2_(A, b, c)
    s = sA * sb * sc
    c = cA * sb * sc + cc * cb
    return atan2(s, c) * _2_0


def excessCagnoli_(a, b, c):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle using U{Cagnoli's
       <https://Zenodo.org/record/35392>} (D.34) formula.

       @arg a: First triangle side (C{radians}).
       @arg b: Second triangle side (C{radians}).
       @arg c: Third triangle side (C{radians}).

       @return: Spherical excess (C{radians}).

       @raise UnitError: Invalid B{C{a}}, B{C{b}} or B{C{c}}.

       @see: Function L{excessLHuilier_} and U{Spherical trigonometry
             <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    a = Radians_(a=a)
    b = Radians_(b=b)
    c = Radians_(c=c)

    r = _maprod(cos, a * _0_5, b * _0_5, c * _0_5)
    if r:
        s =  fsumf_(a, b, c) * _0_5
        t = _maprod(sin, s, s - a, s - b, s - c)
        r =  asin1(sqrt(t) * _0_5 / r) if t > 0 else _0_0
    return Radians(Cagnoli=r * _2_0)


def excessGirard_(A, B, C):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle using U{Girard's
       <https://MathWorld.Wolfram.com/GirardsSphericalExcessFormula.html>} formula.

       @arg A: First interior triangle angle (C{radians}).
       @arg B: Second interior triangle angle (C{radians}).
       @arg C: Third interior triangle angle (C{radians}).

       @return: Spherical excess (C{radians}).

       @raise UnitError: Invalid B{C{A}}, B{C{B}} or B{C{C}}.

       @see: Function L{excessLHuilier_} and U{Spherical trigonometry
             <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    r = fsumf_(Radians_(A=A),
               Radians_(B=B),
               Radians_(C=C), -PI)
    return Radians(Girard=r)


def excessLHuilier_(a, b, c):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle using U{L'Huilier's
       <https://MathWorld.Wolfram.com/LHuiliersTheorem.html>}'s Theorem.

       @arg a: First triangle side (C{radians}).
       @arg b: Second triangle side (C{radians}).
       @arg c: Third triangle side (C{radians}).

       @return: Spherical excess (C{radians}).

       @raise UnitError: Invalid B{C{a}}, B{C{b}} or B{C{c}}.

       @see: Function L{excessCagnoli_}, L{excessGirard_} and U{Spherical
             trigonometry<https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    a = Radians_(a=a)
    b = Radians_(b=b)
    c = Radians_(c=c)

    s =  fsumf_(a, b, c) * _0_5
    r = _maprod(tan_2, s, s - a, s - b, s - c)
    r =  atan(sqrt(r)) if r > 0 else _0_0
    return Radians(LHuilier=r * _4_0)


def excessKarney(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the surface area of a (spherical) quadrilateral bounded by a
       segment of a great circle, two meridians and the equator using U{Karney's
       <https://MathOverflow.net/questions/97711/the-area-of-spherical-polygons>}
       method.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum}) or ellipsoid
                      (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and
                    B{C{lon2}} (C{bool}).

       @return: Surface area, I{signed} (I{square} C{meter} or the same units as
                B{C{radius}} I{squared}) or the I{spherical excess} (C{radians})
                if C{B{radius}=0} or C{None}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise UnitError: Invalid B{C{lat2}} or B{C{lat1}}.

       @raise ValueError: Semi-circular longitudinal delta.

       @see: Functions L{excessKarney_} and L{excessQuad}.
    '''
    r = excessKarney_(*_d3(wrap, lat1, lon1, lat2, lon2))
    if radius:
        r *= _mean_radius(radius, lat1, lat2)**2
    return r


def excessKarney_(phi2, phi1, lam21):
    '''Compute the I{spherical excess} C{E} of a (spherical) quadrilateral bounded by
       a segment of a great circle, two meridians and the equator using U{Karney's
       <https://MathOverflow.net/questions/97711/the-area-of-spherical-polygons>}
       method.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Spherical excess, I{signed} (C{radians}).

       @raise ValueError: Semi-circular longitudinal delta B{C{lam21}}.

       @see: Function L{excessKarney} and U{Area of a spherical polygon
       <https://MathOverflow.net/questions/97711/the-area-of-spherical-polygons>}.
    '''
    # from: Veness <https://www.Movable-Type.co.UK/scripts/latlong.html>  Area
    # method due to Karney: for each edge of the polygon,
    #
    #                 tan(Δλ / 2) · (tan(φ1 / 2) + tan(φ2 / 2))
    #    tan(E / 2) = -----------------------------------------
    #                           1 +  tan(φ1 / 2) · tan(φ2 / 2)
    #
    # where E is the spherical excess of the trapezium obtained by extending
    # the edge to the equator-circle vector for each edge (see also ***).
    t2 = tan_2(phi2)
    t1 = tan_2(phi1)
    c  = (t1 * t2) + _1_0
    s  = (t1 + t2) *  tan_2(lam21, lam21=None)
    return Radians(Karney=atan2(s, c) * _2_0)


# ***) Original post no longer available, following is a copy of the main part
# <http://OSGeo-org.1560.x6.Nabble.com/Area-of-a-spherical-polygon-td3841625.html>
#
# The area of a polygon on a (unit) sphere is given by the spherical excess
#
#    A = 2 * pi - sum(exterior angles)
#
# However this is badly conditioned if the polygon is small.  In this case, use
#
#    A = sum(S12{i, i+1}) over the edges of the polygon
#
# where S12 is the area of the quadrilateral bounded by an edge of the polygon,
# two meridians and the equator, i.e. with vertices (phi1, lambda1), (phi2,
# lambda2), (0, lambda1) and (0, lambda2).  S12 is given by
#
#    tan(S12 / 2) = tan(lambda21 / 2) * (tan(phi1 / 2) + tan(phi2 / 2)) /
#                                       (tan(phi1 / 2) * tan(phi2 / 2) + 1)
#
#                 = tan(lambda21 / 2) * tanh((Lamb(phi1) + Lamb(phi2)) / 2)
#
# where lambda21 = lambda2 - lambda1 and Lamb(x) is the Lambertian (or the
# inverse Gudermannian) function
#
#    Lambertian(x) = asinh(tan(x)) = atanh(sin(x)) = 2 * atanh(tan(x / 2))
#
# Notes: The formula for S12 is exact, except that...
#      - it is indeterminate if an edge is a semi-circle
#      - the formula for A applies only if the polygon does not include a pole
#        (if it does, then add +/- 2 * pi to the result)
#      - in the limit of small phi and lambda, S12 reduces to the trapezoidal
#        formula, S12 = (lambda2 - lambda1) * (phi1 + phi2) / 2
#      - I derived this result from the equation for the area of a spherical
#        triangle in terms of two edges and the included angle given by, e.g.
#        U{Todhunter, I. - Spherical Trigonometry (1871), Sec. 103, Eq. (2)
#        <http://Books.Google.com/books?id=3uBHAAAAIAAJ&pg=PA71>}
#      - I would be interested to know if this formula for S12 is already known
#      - Charles Karney


def excessQuad(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the surface area of a (spherical) quadrilateral bounded by a segment
       of a great circle, two meridians and the equator.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum}) or ellipsoid
                      (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and
                    B{C{lon2}} (C{bool}).

       @return: Surface area, I{signed} (I{square} C{meter} or the same units as
                B{C{radius}} I{squared}) or the I{spherical excess} (C{radians})
                if C{B{radius}=0} or C{None}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise UnitError: Invalid B{C{lat2}} or B{C{lat1}}.

       @see: Function L{excessQuad_} and L{excessKarney}.
    '''
    r = excessQuad_(*_d3(wrap, lat1, lon1, lat2, lon2))
    if radius:
        r *= _mean_radius(radius, lat1, lat2)**2
    return r


def excessQuad_(phi2, phi1, lam21):
    '''Compute the I{spherical excess} C{E} of a (spherical) quadrilateral bounded
       by a segment of a great circle, two meridians and the equator.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Spherical excess, I{signed} (C{radians}).

       @see: Function L{excessQuad} and U{Spherical trigonometry
             <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    c = cos((phi2 - phi1) * _0_5)
    s = sin((phi2 + phi1) * _0_5) * tan_2(lam21)
    return Radians(Quad=atan2(s, c) * _2_0)


def flatLocal(lat1, lon1, lat2, lon2, datum=_WGS84, scaled=True, wrap=False):
    '''Compute the distance between two (ellipsoidal) points using
       the U{ellipsoidal Earth to plane projection<https://WikiPedia.org/
       wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
       aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid (L{Ellipsoid}, L{Ellipsoid2}
                     or L{a_f2Tuple}) to use.
       @kwarg scaled: Scale prime_vertical by C{cos(B{phi})} (C{bool}), see
                      method L{pygeodesy.Ellipsoid.roc2_}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and
                    B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s or ellipsoid axes).

       @raise TypeError: Invalid B{C{datum}}.

       @note: The meridional and prime_vertical radii of curvature are taken and
              scaled at the mean of both latitude.

       @see: Functions L{flatLocal_} or L{hubeny_}, L{cosineLaw}, L{equirectangular},
             L{euclidean}, L{flatPolar}, L{haversine}, L{thomas} and L{vincentys}, method
             L{Ellipsoid.distance2} and U{local, flat earth approximation
             <https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    t = _d3(wrap, lat1, lon1, lat2, lon2)
    E = _ellipsoidal(datum, flatLocal)
    return E._hubeny_2(*t, scaled=scaled, squared=False) * E.a

hubeny = flatLocal  # PYCHOK for Karl Hubeny


def flatLocal_(phi2, phi1, lam21, datum=_WGS84, scaled=True):
    '''Compute the I{angular} distance between two (ellipsoidal) points using
       the U{ellipsoidal Earth to plane projection<https://WikiPedia.org/
       wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
       aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid (L{Ellipsoid}, L{Ellipsoid2}
                     or L{a_f2Tuple}) to use.
       @kwarg scaled: Scale prime_vertical by C{cos(B{phi})} (C{bool}), see
                      method L{pygeodesy.Ellipsoid.roc2_}.

       @return: Angular distance (C{radians}).

       @raise TypeError: Invalid B{C{datum}}.

       @note: The meridional and prime_vertical radii of curvature are taken and
              scaled I{at the mean of both latitude}.

       @see: Functions L{flatLocal} or L{hubeny}, L{cosineLaw_}, L{flatPolar_},
             L{euclidean_}, L{haversine_}, L{thomas_} and L{vincentys_} and U{local,
             flat earth approximation<https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    E = _ellipsoidal(datum, flatLocal_)
    return E._hubeny_2(phi2, phi1, lam21, scaled=scaled, squared=False)

hubeny_ = flatLocal_  # PYCHOK for Karl Hubeny


def flatPolar(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using the U{polar
       coordinate flat-Earth<https://WikiPedia.org/wiki/Geographical_distance
       #Polar_coordinate_flat-Earth_formula>} formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum}) or ellipsoid
                      (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and B{C{lat2}} and B{C{lon2}}
                    (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}} or the datum's or
                ellipsoid axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: Functions L{flatPolar_}, L{cosineLaw}, L{flatLocal} / L{hubeny},
             L{equirectangular}, L{euclidean}, L{haversine}, L{thomas} and L{vincentys}.
    '''
    return _dS(flatPolar_, radius, wrap, lat1, lon1, lat2, lon2)


def flatPolar_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two (spherical) points using the
       U{polar coordinate flat-Earth<https://WikiPedia.org/wiki/Geographical_distance
       #Polar_coordinate_flat-Earth_formula>} formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{flatPolar}, L{cosineLaw_}, L{euclidean_}, L{flatLocal_} /
             L{hubeny_}, L{haversine_}, L{thomas_} and L{vincentys_}.
    '''
    a = fabs(PI_2 - phi1)  # co-latitude
    b = fabs(PI_2 - phi2)  # co-latitude
    if a < b:
        a, b = b, a
    if a < EPS0:
        a = _0_0
    elif b > 0:
        b  = b / a  # /= chokes PyChecker
        c  = b * cos(lam21) * _2_0
        c  = fsumf_(_1_0, b**2, -fabs(c))
        a *= sqrt0(c)
    return a


def _hartzell(pov, los, earth, **kwds):
    '''(INTERNAL) Helper for C{CartesianBase.hartzell} and C{LatLonBase.hartzell}.
    '''
    if earth is None:
        earth =  pov.datum
    else:
        earth = _spherical_datum(earth, name__=hartzell)
        pov   =  pov.toDatum(earth)
    h = pov.height
    if h < 0:  # EPS0
        t = _SPACE_(Fmt.PARENSPACED(height=h), _inside_)
        raise IntersectionError(pov=pov, earth=earth, txt=t)
    return hartzell(pov, los=los, earth=earth, **kwds) if h > 0 else pov  # EPS0


def hartzell(pov, los=False, earth=_WGS84, **name_LatLon_and_kwds):
    '''Compute the intersection of the earth's surface and a Line-Of-Sight from
       a Point-Of-View in space.

       @arg pov: Point-Of-View outside the earth (C{LatLon}, C{Cartesian},
                 L{Ecef9Tuple} or L{Vector3d}).
       @kwarg los: Line-Of-Sight, I{direction} to earth (L{Los}, L{Vector3d}),
                   C{True} for the I{normal, plumb} onto the surface or C{False}
                   or C{None} to point to the center of the earth.
       @kwarg earth: The earth model (L{Datum}, L{Ellipsoid}, L{Ellipsoid2},
                     L{a_f2Tuple} or a C{scalar} earth radius in C{meter}).
       @kwarg name_LatLon_and_kwds: Optional C{B{name}="hartzell"} (C{str}), class
                   C{B{LatLon}=None} to return the intersection and optionally,
                   additional C{LatLon} keyword arguments, include the B{C{datum}}
                   if different from and to convert from B{C{earth}}.

       @return: The intersection (L{Vector3d}, B{C{pov}}'s C{cartesian type} or
                the given B{C{LatLon}} instance) with attribute C{height} set to
                the distance to the B{C{pov}}.

       @raise IntersectionError: Invalid B{C{pov}} or B{C{pov}} inside the earth or
                                 invalid B{C{los}} or B{C{los}} points outside or
                                 away from the earth.

       @raise TypeError: Invalid B{C{earth}}, C{ellipsoid} or C{datum}.

       @see: Class L{Los}, functions L{tyr3d} and L{hartzell4} and methods
             L{Ellipsoid.hartzell4}, any C{Cartesian.hartzell} and C{LatLon.hartzell}.
    '''
    n, kwds = _name2__(name_LatLon_and_kwds, name__=hartzell)
    try:
        D = _spherical_datum(earth, name__=hartzell)
        r, h, i = _MODS.triaxials._hartzell3(pov, los, D.ellipsoid._triaxial)

        C = _MODS.cartesianBase.CartesianBase
        if kwds:
            c = C(r, datum=D)
            r = c.toLatLon(**_xkwds(kwds, height=h))
        elif isinstance(r, C):
            r.height = h
        if i:
            r._iteration = i
    except Exception as x:
        raise IntersectionError(pov=pov, los=los, earth=earth, cause=x, **kwds)
    return _xnamed(r, n) if n else r


def haversine(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using the U{Haversine
       <https://www.Movable-Type.co.UK/scripts/latlong.html>} formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum}) or ellipsoid
                      (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and
                    B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @raise TypeError: Invalid B{C{radius}}.

       @see: U{Distance between two (spherical) points
             <https://www.EdWilliams.org/avform.htm#Dist>}, functions L{cosineLaw},
             L{equirectangular}, L{euclidean}, L{flatLocal} / L{hubeny}, L{flatPolar},
             L{thomas} and L{vincentys} and methods L{Ellipsoid.distance2},
             C{LatLon.distanceTo*} and C{LatLon.equirectangularTo}.

       @note: See note at function L{vincentys_}.
    '''
    return _dS(haversine_, radius, wrap, lat1, lon1, lat2, lon2)


def haversine_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two (spherical) points using the
       U{Haversine<https://www.Movable-Type.co.UK/scripts/latlong.html>} formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{haversine}, L{cosineLaw_}, L{euclidean_}, L{flatLocal_} /
             L{hubeny_}, L{flatPolar_}, L{thomas_} and L{vincentys_}.

       @note: See note at function L{vincentys_}.
    '''
    h = hav(phi2 - phi1) + cos(phi1) * cos(phi2) * hav(lam21)  # haversine
    return atan2(sqrt0(h), sqrt0(_1_0 - h)) * _2_0  # == asin1(sqrt(h)) * 2


def heightOf(angle, distance, radius=R_M):
    '''Determine the height above the (spherical) earth' surface after
       traveling along a straight line at a given tilt.

       @arg angle: Tilt angle above horizontal (C{degrees}).
       @arg distance: Distance along the line (C{meter} or same units as
                      B{C{radius}}).
       @kwarg radius: Optional mean earth radius (C{meter}).

       @return: Height (C{meter}, same units as B{C{distance}} and B{C{radius}}).

       @raise ValueError: Invalid B{C{angle}}, B{C{distance}} or B{C{radius}}.

       @see: U{MultiDop geog_lib.GeogBeamHt<https://GitHub.com/NASA/MultiDop>}
             (U{Shapiro et al. 2009, JTECH
             <https://Journals.AMetSoc.org/doi/abs/10.1175/2009JTECHA1256.1>}
             and U{Potvin et al. 2012, JTECH
             <https://Journals.AMetSoc.org/doi/abs/10.1175/JTECH-D-11-00019.1>}).
    '''
    r = h = Radius(radius)
    d = fabs(Distance(distance))
    if d > h:
        d, h = h, d

    if d > EPS0:  # and h > EPS0
        d = d / h  # /= h chokes PyChecker
        s = sin(Phid(angle=angle, clip=_180_0))
        s = fsumf_(_1_0, s * d * _2_0, d**2)
        if s > 0:
            return h * sqrt(s) - r

    raise _ValueError(angle=angle, distance=distance, radius=radius)


def heightOrthometric(h_loc, N):
    '''Get the I{orthometric} height B{H}, the height above the geoid, earth surface.

       @arg h_loc: The height above the ellipsoid (C{meter}) or an I{ellipsoidal}
                   location (C{LatLon} or C{Cartesian} with a C{height} or C{h}
                   attribute), otherwise C{0 meter}.
       @arg N: The I{geoid} height (C{meter}), the height of the geoid above the
               ellipsoid at the same B{C{h_loc}} location.

       @return: I{Orthometric} height C{B{H} = B{h} - B{N}} (C{meter}, same units
                as B{C{h}} and B{C{N}}).

       @see: U{Ellipsoid, Geoid, and Orthometric Heights<https://www.NGS.NOAA.gov/
             GEOID/PRESENTATIONS/2007_02_24_CCPS/Roman_A_PLSC2007notes.pdf>}, page
             6 and module L{pygeodesy.geoids}.
    '''
    h = h_loc if _isHeight(h_loc) else _xattr(h_loc, height=_xattr(h_loc, h=0))
    return Height(H=Height(h=h) - Height(N=N))


def horizon(height, radius=R_M, refraction=False):
    '''Determine the distance to the horizon from a given altitude above the
       (spherical) earth.

       @arg height: Altitude (C{meter} or same units as B{C{radius}}).
       @kwarg radius: Optional mean earth radius (C{meter}).
       @kwarg refraction: Consider atmospheric refraction (C{bool}).

       @return: Distance (C{meter}, same units as B{C{height}} and B{C{radius}}).

       @raise ValueError: Invalid B{C{height}} or B{C{radius}}.

       @see: U{Distance to horizon<https://www.EdWilliams.org/avform.htm#Horizon>}.
    '''
    h, r = Height(height), Radius(radius)
    if min(h, r) < 0:
        raise _ValueError(height=height, radius=radius)

    if refraction:
        r *= 2.415750694528  # 2.0 / 0.8279
    else:
        r += r + h
    return sqrt0(r * h)


class _idllmn6(object):  # see also .geodesicw._wargs, .latlonBase._toCartesian3, .vector2d._numpy
    '''(INTERNAL) Helper for C{intersection2} and C{intersections2}.
    '''
    @contextmanager  # <https://www.Python.org/dev/peps/pep-0343/> Examples
    def __call__(self, datum, lat1, lon1, lat2, lon2, small, wrap, s, **kwds):
        try:
            if wrap:
                _, lat2, lon2 = _Wrap.latlon3(lon1, lat2, lon2, wrap)
                kwds = _xkwds(kwds, wrap=wrap)  # for _xError
            m = small if small is _100km else Meter_(small=small)
            n = typename(intersections2 if s else intersection2)
            if datum is None or euclidean(lat1, lon1, lat2, lon2) < m:
                d, m = None, _MODS.vector3d
                _i   = m._intersects2 if s else m._intersect3d3
            elif _isRadius(datum) and datum < 0 and not s:
                d = _spherical_datum(-datum, name=n)
                m = _MODS.sphericalNvector
                _i = m.intersection
            else:
                d = _spherical_datum(datum, name=n)
                if d.isSpherical:
                    m = _MODS.sphericalTrigonometry
                    _i = m._intersects2 if s else m._intersect
                elif d.isEllipsoidal:
                    try:
                        if d.ellipsoid.geodesic:
                            pass
                        m = _MODS.ellipsoidalKarney
                    except ImportError:
                        m = _MODS.ellipsoidalExact
                    _i = m._intersections2 if s else m._intersection3  # ellipsoidalBaseDI
                else:
                    raise _TypeError(datum=datum)
            yield _i, d, lat2, lon2, m, n

        except (TypeError, ValueError) as x:
            raise _xError(x, lat1=lat1, lon1=lon1, datum=datum,
                             lat2=lat2, lon2=lon2, small=small, **kwds)

_idllmn6 = _idllmn6()  # PYCHOK singleton


def intersection2(lat1, lon1, bearing1,
                  lat2, lon2, bearing2, datum=None, wrap=False, small=_100km):  # was=True
    '''I{Conveniently} compute the intersection of two lines each defined by
       a (geodetic) point and a bearing from North, using either ...

       1) L{vector3d.intersection3d3} for B{C{small}} distances (below 100 Km
       or about 0.88 degrees) or if I{no} B{C{datum}} is specified, or ...

       2) L{sphericalTrigonometry.intersection} for a spherical B{C{datum}}
       or a C{scalar B{datum}} representing the earth radius, conventionally
       in C{meter} or ...

       3) L{sphericalNvector.intersection} if B{C{datum}} is a I{negative}
       C{scalar}, (negative) earth radius, conventionally in C{meter} or ...

       4) L{ellipsoidalKarney.intersection3} for an ellipsoidal B{C{datum}}
       and if I{Karney}'s U{geographiclib<https://PyPI.org/project/geographiclib>}
       is installed, otherwise ...

       5) L{ellipsoidalExact.intersection3}, provided B{C{datum}} is ellipsoidal.

       @arg lat1: Latitude of the first point (C{degrees}).
       @arg lon1: Longitude of the first point (C{degrees}).
       @arg bearing1: Bearing at the first point (compass C{degrees360}).
       @arg lat2: Latitude of the second point (C{degrees}).
       @arg lon2: Longitude of the second point (C{degrees}).
       @arg bearing2: Bearing at the second point (compass C{degrees360}).
       @kwarg datum: Optional datum (L{Datum}) or ellipsoid (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) or C{scalar} earth radius
                     (C{meter}, same units as B{C{radius1}} and B{C{radius2}})
                     or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and
                    B{C{lon2}} (C{bool}).
       @kwarg small: Upper limit for small distances (C{meter}).

       @return: Intersection point (L{LatLon2Tuple}C{(lat, lon)}).

       @raise IntersectionError: No or an ambiguous intersection or colinear,
                                 parallel or otherwise non-intersecting lines.

       @raise TypeError: Invalid B{C{datum}}.

       @raise UnitError: Invalid B{C{lat1}}, B{C{lon1}}, B{C{bearing1}}, B{C{lat2}},
                         B{C{lon2}} or B{C{bearing2}}.

       @see: Method L{RhumbLine.intersection2}.
    '''
    b1 = Bearing(bearing1=bearing1)
    b2 = Bearing(bearing2=bearing2)
    with _idllmn6(datum, lat1, lon1, lat2, lon2,
                  small, wrap, False, bearing1=b1, bearing2=b2) as t:
        _i, d, lat2, lon2, m, n = t
        if d is None:
            t, _, _ = _i(m.Vector3d(lon1, lat1, 0), b1,
                         m.Vector3d(lon2, lat2, 0), b2, useZ=False)
            t = LatLon2Tuple(t.y, t.x, name=n)

        else:
            t = _i(m.LatLon(lat1, lon1, datum=d), b1,
                   m.LatLon(lat2, lon2, datum=d), b2,
                     LatLon=None, height=0, wrap=False)
            if isinstance(t, Intersection3Tuple):  # ellipsoidal
                t, _, _ = t
            t = LatLon2Tuple(t.lat, t.lon, name=n)
    return t


def intersections2(lat1, lon1, radius1,
                   lat2, lon2, radius2, datum=None, wrap=False, small=_100km):  # was=True
    '''I{Conveniently} compute the intersections of two circles each defined
       by a (geodetic) center point and a radius, using either ...

       1) L{vector3d.intersections2} for B{C{small}} distances (below 100 Km
       or about 0.88 degrees) or if I{no} B{C{datum}} is specified, or ...

       2) L{sphericalTrigonometry.intersections2} for a spherical B{C{datum}}
       or a C{scalar B{datum}} representing the earth radius, conventionally
       in C{meter} or ...

       3) L{ellipsoidalKarney.intersections2} for an ellipsoidal B{C{datum}}
       and if I{Karney}'s U{geographiclib<https://PyPI.org/project/geographiclib>}
       is installed, otherwise ...

       4) L{ellipsoidalExact.intersections2}, provided B{C{datum}} is ellipsoidal.

       @arg lat1: Latitude of the first circle center (C{degrees}).
       @arg lon1: Longitude of the first circle center (C{degrees}).
       @arg radius1: Radius of the first circle (C{meter}, conventionally).
       @arg lat2: Latitude of the second circle center (C{degrees}).
       @arg lon2: Longitude of the second circle center (C{degrees}).
       @arg radius2: Radius of the second circle (C{meter}, same units as B{C{radius1}}).
       @kwarg datum: Optional datum (L{Datum}) or ellipsoid (L{Ellipsoid}, L{Ellipsoid2}
                     or L{a_f2Tuple}) or C{scalar} earth radius (C{meter}, same units as
                     B{C{radius1}} and B{C{radius2}}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and B{C{lon2}}
                    (C{bool}).
       @kwarg small: Upper limit for small distances (C{meter}).

       @return: 2-Tuple of the intersection points, each a L{LatLon2Tuple}C{(lat, lon)}.
                Both points are the same instance, aka the I{radical center} if the
                circles are abutting

       @raise IntersectionError: Concentric, antipodal, invalid or non-intersecting
                                 circles or no convergence.

       @raise TypeError: Invalid B{C{datum}}.

       @raise UnitError: Invalid B{C{lat1}}, B{C{lon1}}, B{C{radius1}}, B{C{lat2}},
                         B{C{lon2}} or B{C{radius2}}.
    '''
    r1 = Radius_(radius1=radius1)
    r2 = Radius_(radius2=radius2)
    with _idllmn6(datum, lat1, lon1, lat2, lon2,
                  small, wrap, True, radius1=r1, radius2=r2) as t:
        _i, d, lat2, lon2, m, n = t
        if d is None:
            r1 = m2degrees(r1, radius=R_M, lat=lat1)
            r2 = m2degrees(r2, radius=R_M, lat=lat2)

            def _V2T(x, y, _, **unused):  # _ == z unused
                return LatLon2Tuple(y, x, name=n)

            t = _i(m.Vector3d(lon1, lat1, 0), r1,
                   m.Vector3d(lon2, lat2, 0), r2, sphere=False,
                     Vector=_V2T)
        else:
            def _LL2T(lat, lon, **unused):
                return LatLon2Tuple(lat, lon, name=n)

            t = _i(m.LatLon(lat1, lon1, datum=d), r1,
                   m.LatLon(lat2, lon2, datum=d), r2,
                     LatLon=_LL2T, height=0, wrap=False)
    return t


def isantipode(lat1, lon1, lat2, lon2, eps=EPS):
    '''Check whether two points are I{antipodal}, on diametrically
       opposite sides of the earth.

       @arg lat1: Latitude of one point (C{degrees}).
       @arg lon1: Longitude of one point (C{degrees}).
       @arg lat2: Latitude of the other point (C{degrees}).
       @arg lon2: Longitude of the other point (C{degrees}).
       @kwarg eps: Tolerance for near-equality (C{degrees}).

       @return: C{True} if points are antipodal within the
                B{C{eps}} tolerance, C{False} otherwise.

       @see: Functions L{isantipode_} and L{antipode}.
    '''
    return (fabs(lat1 + lat2) <= eps and
            fabs(lon1 + lon2) <= eps) or _isequalTo(
            normal(lat1, lon1), antipode(lat2, lon2), eps)


def isantipode_(phi1, lam1, phi2, lam2, eps=EPS):
    '''Check whether two points are I{antipodal}, on diametrically
       opposite sides of the earth.

       @arg phi1: Latitude of one point (C{radians}).
       @arg lam1: Longitude of one point (C{radians}).
       @arg phi2: Latitude of the other point (C{radians}).
       @arg lam2: Longitude of the other point (C{radians}).
       @kwarg eps: Tolerance for near-equality (C{radians}).

       @return: C{True} if points are antipodal within the
                B{C{eps}} tolerance, C{False} otherwise.

       @see: Functions L{isantipode} and L{antipode_}.
    '''
    return (fabs(phi1 + phi2) <= eps and
            fabs(lam1 + lam2) <= eps) or _isequalTo_(
            normal_(phi1, lam1), antipode_(phi2, lam2), eps)


def _isequalTo(p1, p2, eps=EPS):
    '''Compare 2 point lat-/lons ignoring C{class}.
    '''
    return (fabs(p1.lat - p2.lat) <= eps and
            fabs(p1.lon - p2.lon) <= eps) if eps else (p1.latlon == p2.latlon)


def _isequalTo_(p1, p2, eps=EPS):  # underscore_!
    '''(INTERNAL) Compare 2 point phi-/lams ignoring C{class}.
    '''
    return (fabs(p1.phi - p2.phi) <= eps and
            fabs(p1.lam - p2.lam) <= eps) if eps else (p1.philam == p2.philam)


def isnormal(lat, lon, eps=0):
    '''Check whether B{C{lat}} I{and} B{C{lon}} are within their
       respective I{normal} range in C{degrees}.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg eps: Optional tolerance C{degrees}).

       @return: C{True} if C{(abs(B{lat}) + B{eps}) <= 90} and
                C{(abs(B{lon}) + B{eps}) <= 180}, C{False} otherwise.

       @see: Functions L{isnormal_} and L{normal}.
    '''
    return _loneg(fabs(lon)) >= eps and (_90_0 - fabs(lat)) >= eps  # co-latitude


def isnormal_(phi, lam, eps=0):
    '''Check whether B{C{phi}} I{and} B{C{lam}} are within their
       respective I{normal} range in C{radians}.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).
       @kwarg eps: Optional tolerance C{radians}).

       @return: C{True} if C{(abs(B{phi}) + B{eps}) <= PI/2} and
                C{(abs(B{lam}) + B{eps}) <= PI}, C{False} othwerwise.

       @see: Functions L{isnormal} and L{normal_}.
    '''
    return (PI_2 - fabs(phi)) >= eps and (PI - fabs(lam)) >= eps


def _maprod(fun_, *ts):
    '''(INTERNAL) Helper for C{excessCagnoli_} and C{excessLHuilier_}.
    '''
    return fprod(map(fun_, ts))


def _normal2(a, b, n_2, n, n2):
    '''(INTERNAL) Helper for C{normal} and C{normal_}.
    '''
    if fabs(b) > n:
        b = remainder(b, n2)
    if fabs(a) > n_2:
        r = remainder(a, n)
        if r != a:
            a  = -r
            b -=  n if b > 0 else -n
    return float0_(a, b)


def normal(lat, lon, **name):
    '''Normalize a lat- I{and} longitude pair in C{degrees}.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg name: Optional C{B{name}="normal"} (C{str}).

       @return: L{LatLon2Tuple}C{(lat, lon)} with C{-90 <= lat <= 90}
                and C{-180 <= lon <= 180}.

       @see: Functions L{normal_} and L{isnormal}.
    '''
    return LatLon2Tuple(*_normal2(lat, lon, _90_0, _180_0, _360_0),
                          name=_name__(name, name__=normal))


def normal_(phi, lam, **name):
    '''Normalize a lat- I{and} longitude pair in C{radians}.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).
       @kwarg name: Optional C{B{name}="normal_"} (C{str}).

       @return: L{PhiLam2Tuple}C{(phi, lam)} with C{abs(phi) <= PI/2}
                and C{abs(lam) <= PI}.

       @see: Functions L{normal} and L{isnormal_}.
    '''
    return PhiLam2Tuple(*_normal2(phi, lam, PI_2, PI, PI2),
                          name=_name__(name, name__=normal_))


def _opposes(d, m, n, n2):
    '''(INTERNAL) Helper for C{opposing} and C{opposing_}.
    '''
    d = d % n2  # -20 % 360 == 340, -1 % PI2 == PI2 - 1
    return False if d < m or d > (n2 - m) else (
           True if (n - m) < d < (n  + m) else None)


def opposing(bearing1, bearing2, margin=_90_0):
    '''Compare the direction of two bearings given in C{degrees}.

       @arg bearing1: First bearing (compass C{degrees}).
       @arg bearing2: Second bearing (compass C{degrees}).
       @kwarg margin: Optional, interior angle bracket (C{degrees}).

       @return: C{True} if both bearings point in opposite, C{False} if
                in similar or C{None} if in I{perpendicular} directions.

       @see: Function L{opposing_}.
    '''
    m = Degrees_(margin=margin, low=EPS0, high=_90_0)
    return _opposes(bearing2 - bearing1, m, _180_0, _360_0)


def opposing_(radians1, radians2, margin=PI_2):
    '''Compare the direction of two bearings given in C{radians}.

       @arg radians1: First bearing (C{radians}).
       @arg radians2: Second bearing (C{radians}).
       @kwarg margin: Optional, interior angle bracket (C{radians}).

       @return: C{True} if both bearings point in opposite, C{False} if
                in similar or C{None} if in perpendicular directions.

       @see: Function L{opposing}.
    '''
    m = Radians_(margin=margin, low=EPS0, high=PI_2)
    return _opposes(radians2 - radians1, m, PI, PI2)


def _Propy(func, nargs, kwds):
    '''(INTERNAL) Helper for the C{frechet.[-]Frechet**} and
       C{hausdorff.[-]Hausdorff*} classes.
    '''
    try:
        _xcallable(distance=func)
        # assert _args_kwds_count2(func)[0] == nargs + int(ismethod(func))
        _ = func(*_0_0s(nargs), **kwds)
    except Exception as x:
        t = unstr(func, **kwds)
        raise _TypeError(t, cause=x)
    return func


def _radical2(d, r1, r2, **name):  # in .ellipsoidalBaseDI, .sphericalTrigonometry, .vector3d
    # (INTERNAL) See C{radical2} below
    # assert d > EPS0
    r =  fsumf_(_1_0, (r1 / d)**2, -(r2 / d)**2) * _0_5
    n = _name__(name, name__=radical2)
    return Radical2Tuple(max(_0_0, min(_1_0, r)), r * d, name=n)


def radical2(distance, radius1, radius2, **name):
    '''Compute the I{radical ratio} and I{radical line} of two U{intersecting
       circles<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>}.

       The I{radical line} is perpendicular to the axis thru the centers of
       the circles at C{(0, 0)} and C{(B{distance}, 0)}.

       @arg distance: Distance between the circle centers (C{scalar}).
       @arg radius1: Radius of the first circle (C{scalar}).
       @arg radius2: Radius of the second circle (C{scalar}).
       @kwarg name: Optional C{B{name}=NN} (C{str}).

       @return: A L{Radical2Tuple}C{(ratio, xline)} where C{0.0 <= ratio <=
                1.0} and C{xline} is along the B{C{distance}}.

       @raise IntersectionError: The B{C{distance}} exceeds the sum of
                                 B{C{radius1}} and B{C{radius2}}.

       @raise UnitError: Invalid B{C{distance}}, B{C{radius1}} or B{C{radius2}}.
    '''
    d  = Distance_(distance, low=_0_0)
    r1 = Radius_(radius1=radius1)
    r2 = Radius_(radius2=radius2)
    if d > (r1 + r2):
        raise IntersectionError(distance=d, radius1=r1, radius2=r2,
                                            txt=_too_(_distant_))
    return _radical2(d, r1, r2, **name) if d > EPS0 else \
            Radical2Tuple(_0_5, _0_0, **name)


class Radical2Tuple(_NamedTuple):
    '''2-Tuple C{(ratio, xline)} of the I{radical} C{ratio} and
       I{radical} C{xline}, both C{scalar} and C{0.0 <= ratio <= 1.0}
    '''
    _Names_ = (_ratio_, _xline_)
    _Units_ = ( Scalar,  Scalar)


def _radistance(inst):
    '''(INTERNAL) Helper for the L{frechet._FrechetMeterRadians}
       and L{hausdorff._HausdorffMeterRedians} classes.
    '''
    wrap_, kwds_ = _xkwds_pop2(inst._kwds, wrap=False)
    func_ = inst._func_
    try:  # calling lower-overhead C{func_}
        func_(0, _0_25, _0_5, **kwds_)
        wrap_ = _Wrap._philamop(wrap_)
    except TypeError:
        return inst.distance

    def _philam(p):
        try:
            return p.phi, p.lam
        except AttributeError:  # no .phi or .lam
            return radians(p.lat), radians(p.lon)

    def _func_wrap(point1, point2):
        phi1, lam1 = wrap_(*_philam(point1))
        phi2, lam2 = wrap_(*_philam(point2))
        return func_(phi2, phi1, lam2 - lam1, **kwds_)

    inst._units = inst._units_
    return _func_wrap


def _scale_deg(lat1, lat2):  # degrees
    # scale factor cos(mean of lats) for delta lon
    m = fabs(lat1 + lat2) * _0_5
    return cos(radians(m)) if m < _90_0 else _0_0


def _scale_rad(phi1,  phi2):  # radians, by .frechet, .hausdorff, .heights
    # scale factor cos(mean of phis) for delta lam
    m = fabs(phi1 + phi2) * _0_5
    return cos(m) if m < PI_2 else _0_0


def _sincosa6(phi2, phi1, lam21):  # [4] in cosineLaw
    '''(INTERNAL) C{sin}es, C{cos}ines and C{acos}ine.
    '''
    s2, c2, s1, c1, _, c21 = sincos2_(phi2, phi1, lam21)
    return s2, c2, s1, c1, acos1(s1 * s2 + c1 * c2 * c21), c21


def thomas(lat1, lon1, lat2, lon2, datum=_WGS84, wrap=False):
    '''Compute the distance between two (ellipsoidal) points using U{Thomas'
       <https://apps.DTIC.mil/dtic/tr/fulltext/u2/703541.pdf>} formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid (L{Ellipsoid}, L{Ellipsoid2}
                     or L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and
                    B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s or ellipsoid axes).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{thomas_}, L{cosineLaw}, L{equirectangular}, L{euclidean},
             L{flatLocal} / L{hubeny}, L{flatPolar}, L{haversine}, L{vincentys} and
             method L{Ellipsoid.distance2}.
    '''
    return _dE(thomas_, datum, wrap, lat1, lon1, lat2, lon2)


def thomas_(phi2, phi1, lam21, datum=_WGS84):
    '''Compute the I{angular} distance between two (ellipsoidal) points using
       U{Thomas'<https://apps.DTIC.mil/dtic/tr/fulltext/u2/703541.pdf>} formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Datum (L{Datum}) ?or ellipsoid to use (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).

       @return: Angular distance (C{radians}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{thomas}, L{cosineLaw_}, L{euclidean_}, L{flatLocal_} /
             L{hubeny_}, L{flatPolar_}, L{haversine_} and L{vincentys_} and
             U{Geodesy-PHP<https://GitHub.com/jtejido/geodesy-php/blob/master/
             src/Geodesy/Distance/ThomasFormula.php>}.
    '''
    s2, c2, s1, c1, r, _ = _sincosa6(phi2, phi1, lam21)
    if r and isnon0(c1) and isnon0(c2):
        E = _ellipsoidal(datum, thomas_)
        f =  E.f * _0_25
        if f:  # ellipsoidal
            r1 = atan2(E.b_a * s1, c1)
            r2 = atan2(E.b_a * s2, c2)

            j = (r2 + r1) * _0_5
            k = (r2 - r1) * _0_5
            sj, cj, sk, ck, h, _ = sincos2_(j, k, lam21 * _0_5)

            h =  fsumf_(sk**2, (ck * h)**2, -(sj * h)**2)
            u = _1_0 - h
            if isnon0(u) and isnon0(h):
                r = atan(sqrt0(h / u)) * 2  # == acos(1 - 2 * h)
                sr, cr = sincos2(r)
                if isnon0(sr):
                    u = (sj * ck)**2 * 2 / u
                    h = (sk * cj)**2 * 2 / h
                    x = u + h
                    y = u - h

                    b = r * 2
                    s = r / sr
                    e = 4 * s**2
                    d = 2 * cr
                    a = e * d
                    c = s - (a - d) * _0_5

                    t  = fdot_(a, x,  -b, y,  -d, y**2,  c, x**2,  e, x * y) * _0_25
                    r -= fdot_(s, x,  -1, y,  -t, f) * f * sr
    return r


def vincentys(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using
       U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
       spherical formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum}) or
                      ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple})
                      to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}} and
                    B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @raise UnitError: Invalid B{C{radius}}.

       @see: Functions L{vincentys_}, L{cosineLaw}, L{equirectangular}, L{euclidean},
             L{flatLocal}/L{hubeny}, L{flatPolar}, L{haversine} and L{thomas} and
             methods L{Ellipsoid.distance2}, C{LatLon.distanceTo*} and
             C{LatLon.equirectangularTo}.

       @note: See note at function L{vincentys_}.
    '''
    return _dS(vincentys_, radius, wrap, lat1, lon1, lat2, lon2)


def vincentys_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two (spherical) points using
       U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
       spherical formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{vincentys}, L{cosineLaw_}, L{euclidean_}, L{flatLocal_} /
             L{hubeny_}, L{flatPolar_}, L{haversine_} and L{thomas_}.

       @note: Functions L{vincentys_}, L{haversine_} and L{cosineLaw_} produce
              equivalent results, but L{vincentys_} is suitable for antipodal
              points and slightly more expensive (M{3 cos, 3 sin, 1 hypot, 1 atan2,
              6 mul, 2 add}) than L{haversine_} (M{2 cos, 2 sin, 2 sqrt, 1 atan2, 5
              mul, 1 add}) and L{cosineLaw_} (M{3 cos, 3 sin, 1 acos, 3 mul, 1 add}).
    '''
    s1, c1, s2, c2, s21, c21 = sincos2_(phi1, phi2, lam21)

    c = c2 * c21
    x = s1 * s2 + c1 * c
    y = c1 * s2 - s1 * c
    return atan2(hypot(c2 * s21, y), x)

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
