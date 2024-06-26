
# -*- coding: utf-8 -*-

u'''Various utility functions.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>} and
U{Vector-based geodesy<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _copysign, isinstanceof, isint, isstr, neg
from pygeodesy.constants import EPS, EPS0, INF, NAN, PI, PI2, PI_2, R_M, \
                               _M_KM, _M_NM, _M_SM, _0_0, _1__90, _0_5, _1_0, \
                               _N_1_0, _2__PI, _10_0, _90_0, _N_90_0, _180_0, \
                               _N_180_0, _360_0, _400_0, _copysign_0_0, \
                               _float as _F, _isfinite, isnan, isnear0, \
                               _over, _umod_360, _umod_PI2
from pygeodesy.errors import _ValueError, _xkwds, _xkwds_get1,  _ALL_LAZY, _MODS
from pygeodesy.internals import _passargs  # , _MODS?
from pygeodesy.interns import _edge_, _radians_, _semi_circular_, _SPACE_
# from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS  # from .errors
from pygeodesy.units import Degrees, Degrees_, Feet, Float, Lam, Lamd, \
                            Meter, Meter2, Radians, Radians_

from math import acos, asin, atan2, cos, degrees, fabs, radians, sin, tan  # pow

__all__ = _ALL_LAZY.utily
__version__ = '24.06.15'

# read constant name "_M_Unit" as "meter per Unit"
_M_CHAIN     = _F(  20.1168)     # yard2m(1) * 22
_M_FATHOM    = _F(   1.8288)     # yard2m(1) * 2 or _M_NM * 1e-3
_M_FOOT      = _F(   0.3048)     # Int'l (1 / 3.2808398950131 = 10_000 / (254 * 12))
_M_FOOT_GE   = _F(   0.31608)    # German Fuss (1 / 3.1637560111364)
_M_FOOT_FR   = _F(   0.3248406)  # French Pied-du-Roi or pied (1 / 3.0784329298739)
_M_FOOT_USVY = _F(   0.3048006096012192)  # US Survey (1200 / 3937)
_M_FURLONG   = _F( 201.168)      # 220 * yard2m(1) == 10 * m2chain(1)
# _M_KM      = _F(1000.0)        # kilo meter
# _M_NM      = _F(1852.0)        # nautical mile
# _M_SM      = _F(1609.344)      # statute mile
_M_TOISE     = _F(   1.9490436)  # French toise, 6 pieds (6 / 3.0784329298739)
_M_YARD_UK   = _F(   0.9144)     # 254 * 12 * 3 / 10_000 == 3 * ft2m(1) Int'l


def _abs1nan(x):
    '''(INTERNAL) Bracket C{x}.
    '''
    return _N_1_0 < x < _1_0 or isnan(x)


def acos1(x):
    '''Return C{math.acos(max(-1, min(1, B{x})))}.
    '''
    return acos(x) if _abs1nan(x) else (PI if x < 0 else _0_0)


def acre2ha(acres):
    '''Convert acres to hectare.

       @arg acres: Value in acres (C{scalar}).

       @return: Value in C{hectare} (C{float}).

       @raise ValueError: Invalid B{C{acres}}.
    '''
    # 0.40468564224 == acre2m2(1) / 10_000
    return Float(ha=Float(acres) * 0.40468564224)


def acre2m2(acres):
    '''Convert acres to I{square} meter.

       @arg acres: Value in acres (C{scalar}).

       @return: Value in C{meter^2} (C{float}).

       @raise ValueError: Invalid B{C{acres}}.
    '''
    # 4046.8564224 == chain2m(1) * furlong2m(1)
    return Meter2(Float(acres) * 4046.8564224)


def asin1(x):
    '''Return C{math.asin(max(-1, min(1, B{x})))}.
    '''
    return asin(x) if _abs1nan(x) else _copysign(PI_2, x)


def atan1(y, x=_1_0):
    '''Return C{atan(B{y} / B{x})} angle in C{radians} M{[-PI/2..+PI/2]}
       using C{atan2} for consistency and to avoid C{ZeroDivisionError}.
    '''
    return atan2(-y, -x) if x < 0 else atan2(y, x or _0_0)  # -0. to 0.


def atan1d(y, x=_1_0):
    '''Return C{atan(B{y} / B{x})} angle in C{degrees} M{[-90..+90]}
       using C{atan2d} for consistency and to avoid C{ZeroDivisionError}.

       @see: Function L{pygeodesy.atan2d}.
    '''
    return atan2d(-y, -x) if x < 0 else atan2d(y, x or _0_0)  # -0. to 0.


def atan2b(y, x):
    '''Return C{atan2(B{y}, B{x})} in degrees M{[0..+360]}.

       @see: Function L{pygeodesy.atan2d}.
    '''
    b = atan2d(y, x)
    if b < 0:
        b += _360_0
    return b


def atan2d(y, x, reverse=False):
    '''Return C{atan2(B{y}, B{x})} in degrees M{[-180..+180]},
       optionally I{reversed} (by 180 degrees for C{azimuth}s).

       @see: I{Karney}'s C++ function U{Math.atan2d
             <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>}.
    '''
    if fabs(y) > fabs(x) > 0:
        if y < 0:  # q = 3
            d = degrees(atan2(x, -y)) - _90_0
        else:  # q = 2
            d = _90_0 - degrees(atan2(x, y))
    elif isnan(x) or isnan(y):
        return NAN
    elif y:
        if x > 0:  # q = 0
            d =  degrees(atan2(y, x))
        elif x < 0:  # q = 1
            d = _copysign(_180_0, y) - degrees(atan2(y, -x))
        else:  # x == 0
            d = _copysign(_90_0, y)
    else:
        d = _180_0 if x < 0 else _0_0
    return _azireversed(d) if reverse else d


def _azireversed(azimuth):  # in .rhumbBase
    '''(INTERNAL) Return the I{reverse} B{C{azimuth}} in degrees M{[-180..+180]}.
    '''
    return azimuth + (_N_180_0 if azimuth > 0 else _180_0)


def chain2m(chains):
    '''Convert I{UK} chains to meter.

       @arg chains: Value in chains (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{chains}}.
    '''
    return Meter(Float(chains=chains) * _M_CHAIN)


def circle4(earth, lat):
    '''Get the equatorial or a parallel I{circle of latitude}.

       @arg earth: The earth radius, ellipsoid or datum
                   (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                   L{Datum} or L{a_f2Tuple}).
       @arg lat: Geodetic latitude (C{degrees90}, C{str}).

       @return: A L{Circle4Tuple}C{(radius, height, lat, beta)}
                instance.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{earth}}.

       @raise ValueError: B{C{earth}} or B{C{lat}}.
    '''
    E = _MODS.datums._earth_ellipsoid(earth)
    return E.circle4(lat)


def cot(rad, **error_kwds):
    '''Return the C{cotangent} of an angle in C{radians}.

       @arg rad: Angle (C{radians}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: C{cot(B{rad})}.

       @raise ValueError: L{pygeodesy.isnear0}C{(sin(B{rad})}.
    '''
    s, c = sincos2(rad)
    if c:
        if isnear0(s):
            raise _valueError(cot, rad, **error_kwds)
        c = c / s  # /= chokes PyChecker
    return c


def cot_(*rads, **error_kwds):
    '''Return the C{cotangent} of angle(s) in C{radiansresection}.

       @arg rads: One or more angles (C{radians}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: Yield the C{cot(B{rad})} for each angle.

       @raise ValueError: See L{pygeodesy.cot}.
    '''
    try:
        for r in rads:
            yield cot(r)
    except ValueError:
        raise _valueError(cot_, r, **error_kwds)


def cotd(deg, **error_kwds):
    '''Return the C{cotangent} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: C{cot(B{deg})}.

       @raise ValueError: L{pygeodesy.isnear0}C{(sin(B{deg})}.
    '''
    s, c = sincos2d(deg)
    if c:
        if isnear0(s):
            raise _valueError(cotd, deg, **error_kwds)
        c = c / s  # /= chokes PyChecker
    elif s < 0:
        c = -c  # negate-0
    return c


def cotd_(*degs, **error_kwds):
    '''Return the C{cotangent} of angle(s) in C{degrees}.

       @arg degs: One or more angles (C{degrees}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: Yield the C{cot(B{deg})} for each angle.

       @raise ValueError: See L{pygeodesy.cotd}.
    '''
    try:
        for d in degs:
            yield cotd(d)
    except ValueError:
        raise _valueError(cotd_, d, **error_kwds)


def degrees90(rad):
    '''Convert radians to degrees and wrap M{[-90..+90)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{degrees90}).
    '''
    return wrap90(degrees(rad))


def degrees180(rad):
    '''Convert radians to degrees and wrap M{[-180..+180)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{degrees180}).
    '''
    return wrap180(degrees(rad))


def degrees360(rad):
    '''Convert radians to degrees and wrap M{[0..+360)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{degrees360}).
    '''
    return _umod_360(degrees(rad))


def degrees2grades(deg):
    '''Convert degrees to I{grades} (aka I{gons} or I{gradians}).

       @arg deg: Angle (C{degrees}).

       @return: Angle (C{grades}).
    '''
    return Degrees(deg) * _400_0 / _360_0


def degrees2m(deg, radius=R_M, lat=0):
    '''Convert an angle to a distance along the equator or
       along the parallel at an other (geodetic) latitude.

       @arg deg: The angle (C{degrees}).
       @kwarg radius: Mean earth radius, ellipsoid or datum
                      (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                      L{Datum} or L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Distance (C{meter}, same units as B{C{radius}}
                or equatorial and polar radii) or C{0.0} for
                near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{deg}}, B{C{radius}} or
                          B{C{lat}}.

       @see: Function L{radians2m} and L{m2degrees}.
    '''
    return _Radians2m(Lamd(deg=deg, clip=0), radius, lat)


def fathom2m(fathoms):
    '''Convert I{Imperial} fathom to meter.

       @arg fathoms: Value in fathoms (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{fathoms}}.

       @see: Function L{toise2m}, U{Fathom<https://WikiPedia.org/wiki/Fathom>}
             and U{Klafter<https://WikiPedia.org/wiki/Klafter>}.
    '''
    return Meter(Float(fathoms=fathoms) * _M_FATHOM)


def ft2m(feet, usurvey=False, pied=False, fuss=False):
    '''Convert I{International}, I{US Survey}, I{French} or I{German}
       B{C{feet}} to C{meter}.

       @arg feet: Value in feet (C{scalar}).
       @kwarg usurvey: If C{True}, convert I{US Survey} foot else ...
       @kwarg pied: If C{True}, convert French I{pied-du-Roi} else ...
       @kwarg fuss: If C{True}, convert German I{Fuss}, otherwise
                    I{International} foot to C{meter}.

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{feet}}.
    '''
    return Meter(Feet(feet) * (_M_FOOT_USVY if usurvey else
                              (_M_FOOT_FR   if pied    else
                              (_M_FOOT_GE   if fuss    else _M_FOOT))))


def furlong2m(furlongs):
    '''Convert a furlong to meter.

       @arg furlongs: Value in furlongs (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{furlongs}}.
    '''
    return Meter(Float(furlongs=furlongs) * _M_FURLONG)


def grades(rad):
    '''Convert radians to I{grades} (aka I{gons} or I{gradians}).

       @arg rad: Angle (C{radians}).

       @return: Angle (C{grades}).
    '''
    return Float(grades=Float(rad=rad) * _400_0 / PI2)


def grades400(rad):
    '''Convert radians to I{grades} (aka I{gons} or I{gradians}) and wrap M{[0..+400)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{grades}).
    '''
    return Float(grades400=(grades(rad) % _400_0) or _0_0)  # _umod_400


def grades2degrees(gon):
    '''Convert I{grades} (aka I{gons} or I{gradians}) to C{degrees}.

       @arg gon: Angle (C{grades}).

       @return: Angle (C{degrees}).
    '''
    return Degrees(Float(gon=gon) * _360_0 / _400_0)


def grades2radians(gon):
    '''Convert I{grades} (aka I{gons} or I{gradians}) to C{radians}.

       @arg gon: Angle (C{grades}).

       @return: Angle (C{radians}).
    '''
    return Radians(Float(gon=gon) * PI2 / _400_0)


def km2m(km):
    '''Convert kilo meter to meter (m).

       @arg km: Value in kilo meter (C{scalar}).

       @return: Value in meter (C{float}).

       @raise ValueError: Invalid B{C{km}}.
    '''
    return Meter(Float(km=km) * _M_KM)


def _loneg(lon):
    '''(INTERNAL) "Complement" of C{lon}.
    '''
    return _180_0 - lon


def m2chain(meter):
    '''Convert meter to I{UK} chains.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{chains} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(chain=Meter(meter) / _M_CHAIN)  # * 0.049709695378986715


def m2degrees(distance, radius=R_M, lat=0):
    '''Convert a distance to an angle along the equator or
       along the parallel at an other (geodetic) latitude.

       @arg distance: Distance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      an L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or
                      L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Angle (C{degrees}) or C{INF} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{distance}}, B{C{radius}}
                          or B{C{lat}}.

       @see: Function L{m2radians} and L{degrees2m}.
    '''
    return degrees(m2radians(distance, radius=radius, lat=lat))


def m2fathom(meter):
    '''Convert meter to I{Imperial} fathoms.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{fathoms} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.

       @see: Function L{m2toise}, U{Fathom<https://WikiPedia.org/wiki/Fathom>}
             and U{Klafter<https://WikiPedia.org/wiki/Klafter>}.
    '''
    return Float(fathom=Meter(meter) / _M_FATHOM)  # * 0.546806649


def m2ft(meter, usurvey=False, pied=False, fuss=False):
    '''Convert meter to I{International}, I{US Survey}, I{French} or
       or I{German} feet (C{ft}).

       @arg meter: Value in meter (C{scalar}).
       @kwarg usurvey: If C{True}, convert to I{US Survey} foot else ...
       @kwarg pied: If C{True}, convert to French I{pied-du-Roi} else ...
       @kwarg fuss: If C{True}, convert to German I{Fuss}, otherwise to
                    I{International} foot.

       @return: Value in C{feet} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    # * 3.2808333333333333, US Survey 3937 / 1200
    # * 3.2808398950131235, Int'l 10_000 / (254 * 12)
    return Float(feet=Meter(meter) / (_M_FOOT_USVY if usurvey else
                                     (_M_FOOT_FR   if pied    else
                                     (_M_FOOT_GE   if fuss    else _M_FOOT))))


def m2furlong(meter):
    '''Convert meter to furlongs.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{furlongs} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(furlong=Meter(meter) / _M_FURLONG)  # * 0.00497096954


def m2km(meter):
    '''Convert meter to kilo meter (Km).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in Km (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(km=Meter(meter) / _M_KM)


def m2NM(meter):
    '''Convert meter to nautical miles (NM).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{NM} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(NM=Meter(meter) / _M_NM)  # * 5.39956804e-4


def m2radians(distance, radius=R_M, lat=0):
    '''Convert a distance to an angle along the equator or along the
       parallel at an other (geodetic) latitude.

       @arg distance: Distance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      an L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or
                      L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Angle (C{radians}) or C{INF} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{distance}}, B{C{radius}}
                          or B{C{lat}}.

       @see: Function L{m2degrees} and L{radians2m}.
    '''
    m = circle4(radius, lat).radius
    return INF if m < EPS0 else Radians(Float(distance=distance) / m)


def m2SM(meter):
    '''Convert meter to statute miles (SM).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{SM} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(SM=Meter(meter) / _M_SM)  # * 6.21369949e-4 == 1 / 1609.344


def m2toise(meter):
    '''Convert meter to French U{toises<https://WikiPedia.org/wiki/Toise>}.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{toises} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.

       @see: Function L{m2fathom}.
    '''
    return Float(toise=Meter(meter) / _M_TOISE)  # * 0.513083632632119


def m2yard(meter):
    '''Convert meter to I{UK} yards.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{yards} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(yard=Meter(meter) / _M_YARD_UK)  # * 1.0936132983377078


def NM2m(nm):
    '''Convert nautical miles to meter (m).

       @arg nm: Value in nautical miles (C{scalar}).

       @return: Value in meter (C{float}).

       @raise ValueError: Invalid B{C{nm}}.
    '''
    return Meter(Float(nm=nm) * _M_NM)


def radians2m(rad, radius=R_M, lat=0):
    '''Convert an angle to a distance along the equator or
       along the parallel at an other (geodetic) latitude.

       @arg rad: The angle (C{radians}).
       @kwarg radius: Mean earth radius, ellipsoid or datum
                      (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                      L{Datum} or L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Distance (C{meter}, same units as B{C{radius}}
                or equatorial and polar radii) or C{0.0} for
                near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{rad}}, B{C{radius}} or
                          B{C{lat}}.

       @see: Function L{degrees2m} and L{m2radians}.
    '''
    return _Radians2m(Lam(rad=rad, clip=0), radius, lat)


def _Radians2m(rad, radius, lat):
    '''(INTERNAL) Helper for C{degrees2m} and C{radians2m}.
    '''
    m = circle4(radius, lat).radius
    return _0_0 if m < EPS0 else (rad * m)


def radiansPI(deg):
    '''Convert and wrap degrees to radians M{[-PI..+PI]}.

       @arg deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI})
    '''
    return wrapPI(radians(deg))


def radiansPI2(deg):
    '''Convert and wrap degrees to radians M{[0..+2PI)}.

       @arg deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI2})
    '''
    return _umod_PI2(radians(deg))


def radiansPI_2(deg):
    '''Convert and wrap degrees to radians M{[-3PI/2..+PI/2]}.

       @arg deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI_2})
    '''
    return wrapPI_2(radians(deg))


def _sin0cos2(q, r, sign):
    '''(INTERNAL) 2-tuple (C{sin(r), cos(r)}) in quadrant C{0 <= B{q} <= 3}
       and C{sin} zero I{signed} with B{C{sign}}.
    '''
    if r < PI_2:
        s, c = sin(r), cos(r)
        t = s, c, -s, -c, s
    else:  # r == PI_2
        t = _1_0, _0_0, _N_1_0, _0_0, _1_0
#   else:  # r == 0, testUtility failures
#       t = _0_0, _1_0, _0_0, _N_1_0, _0_0
#   q &= 3
    s = t[q]     or _copysign_0_0(sign)
    c = t[q + 1] or _0_0
    return s, c


def SinCos2(x):
    '''Get C{sin} and C{cos} of I{typed} angle.

       @arg x: Angle (L{Degrees}, L{Radians} or scalar C{radians}).

       @return: 2-Tuple (C{sin(B{x})}, C{cos(B{x})}).
    '''
    return sincos2d(x) if isinstanceof(x, Degrees, Degrees_) else (
           sincos2(x)  if isinstanceof(x, Radians, Radians_) else
           sincos2(float(x)))  # assume C{radians}


def sincos2(rad):
    '''Return the C{sine} and C{cosine} of an angle in C{radians}.

       @arg rad: Angle (C{radians}).

       @return: 2-Tuple (C{sin(B{rad})}, C{cos(B{rad})}).

       @see: U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/
             classGeographicLib_1_1Math.html#sincosd>} function U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             python/geographiclib/geomath.py#l155>} and C++ U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             include/GeographicLib/Math.hpp#l558>}.
    '''
    if _isfinite(rad):
        q = int(rad * _2__PI)  # int(math.floor)
        if q < 0:
            q -= 1
        t = _sin0cos2(q & 3, rad - q * PI_2, rad)
    else:
        t =  NAN, NAN
    return t


def sincos2_(*rads):
    '''Return the C{sine} and C{cosine} of angle(s) in C{radians}.

       @arg rads: One or more angles (C{radians}).

       @return: Yield the C{sin(B{rad})} and C{cos(B{rad})} for each angle.

       @see: function L{sincos2}.
    '''
    for r in rads:
        s, c = sincos2(r)
        yield s
        yield c


def sincos2d(deg, **adeg):
    '''Return the C{sine} and C{cosine} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).
       @kwarg adeg: Optional correction (C{degrees}).

       @return: 2-Tuple (C{sin(B{deg_})}, C{cos(B{deg_})}, C{B{deg_} =
                B{deg} + B{adeg}}).

       @see: U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/
             classGeographicLib_1_1Math.html#sincosd>} function U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             python/geographiclib/geomath.py#l155>} and C++ U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             include/GeographicLib/Math.hpp#l558>}.
    '''
    if _isfinite(deg):
        q = int(deg * _1__90)  # int(math.floor)
        if q < 0:
            q -= 1
        d = deg - q * _90_0
        if adeg:
            t = _xkwds_get1(adeg, adeg=_0_0)
            d = _MODS.karney._around(d + t)
        t = _sin0cos2(q & 3, radians(d), deg)
    else:
        t =  NAN, NAN
    return t


def sincos2d_(*degs):
    '''Return the C{sine} and C{cosine} of angle(s) in C{degrees}.

       @arg degs: One or more angles (C{degrees}).

       @return: Yield the C{sin(B{deg})} and C{cos(B{deg})} for each angle.

       @see: Function L{sincos2d}.
    '''
    for d in degs:
        s, c = sincos2d(d)
        yield s
        yield c


def sincostan3(rad):
    '''Return the C{sine}, C{cosine} and C{tangent} of an angle in C{radians}.

       @arg rad: Angle (C{radians}).

       @return: 3-Tuple (C{sin(B{rad})}, C{cos(B{rad})}, C{tan(B{rad})}).

       @see: Function L{sincos2}.
    '''
    s, c = sincos2(float(rad))
    t = NAN if s is NAN else (_over(s, c) if s else neg(s, neg0=c < 0))
    return s, c, t


def SM2m(sm):
    '''Convert statute miles to meter (m).

       @arg sm: Value in statute miles (C{scalar}).

       @return: Value in meter (C{float}).

       @raise ValueError: Invalid B{C{sm}}.
    '''
    return Meter(Float(sm=sm) * _M_SM)


def tan_2(rad, **semi):  # edge=1
    '''Compute the tangent of half angle.

       @arg rad: Angle (C{radians}).
       @kwarg semi: Angle or edge name and index
                    for semi-circular error.

       @return: M{tan(rad / 2)} (C{float}).

       @raise ValueError: If B{C{rad}} is semi-circular
                          and B{C{semi}} is given.
    '''
    # .formy.excessKarney_, .sphericalTrigonometry.areaOf
    if semi and isnear0(fabs(rad) - PI):
        for n, v in semi.items():
            break
        n = _SPACE_(n, _radians_) if not isint(v) else \
            _SPACE_(_MODS.streprs.Fmt.SQUARE(**semi), _edge_)
        raise _ValueError(n, rad, txt=_semi_circular_)

    return tan(rad * _0_5) if _isfinite(rad) else NAN


def tand(deg, **error_kwds):
    '''Return the C{tangent} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: C{tan(B{deg})}.

       @raise ValueError: If L{pygeodesy.isnear0}C{(cos(B{deg})}.
    '''
    s, c = sincos2d(deg)
    if s:
        if isnear0(c):
            raise _valueError(tand, deg, **error_kwds)
        s = s / c  # /= chokes PyChecker
    elif c < 0:
        s = -s  # negate-0
    return s


def tand_(*degs, **error_kwds):
    '''Return the C{tangent} of angle(s) in C{degrees}.

       @arg degs: One or more angles (C{degrees}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: Yield the C{tan(B{deg})} for each angle.

       @raise ValueError: See L{pygeodesy.tand}.
    '''
    for d in degs:
        yield tand(d, **error_kwds)


def tanPI_2_2(rad):
    '''Compute the tangent of half angle, 90 degrees rotated.

       @arg rad: Angle (C{radians}).

       @return: M{tan((rad + PI/2) / 2)} (C{float}).
    '''
    return tan((rad + PI_2) * _0_5) if _isfinite(rad) else (
           NAN if isnan(rad) else (_N_90_0 if rad < 0 else _90_0))


def toise2m(toises):
    '''Convert French U{toises<https://WikiPedia.org/wiki/Toise>} to meter.

       @arg toises: Value in toises (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{toises}}.

       @see: Function L{fathom2m}.
    '''
    return Meter(Float(toises=toises) * _M_TOISE)


def truncate(x, ndigits=None):
    '''Truncate to the given number of digits.

       @arg x: Value to truncate (C{scalar}).
       @kwarg ndigits: Number of digits (C{int}),
                       aka I{precision}.

       @return: Truncated B{C{x}} (C{float}).

       @see: Python function C{round}.
    '''
    if isint(ndigits):
        p = _10_0**ndigits
        x =  int(x * p) / p
    return x


def unroll180(lon1, lon2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in degrees.

       @arg lon1: Start longitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg wrap: If C{True}, wrap and unroll to the M{(-180..+180]}
                    range (C{bool}).

       @return: 2-Tuple C{(B{lon2}-B{lon1}, B{lon2})} unrolled (C{degrees},
                C{degrees}).

       @see: Capability C{LONG_UNROLL} in U{GeographicLib
             <https://GeographicLib.SourceForge.io/html/python/interface.html#outmask>}.
    '''
    d = lon2 - lon1
    if wrap:
        u = wrap180(d)
        if u != d:
            return u, (lon1 + u)
    return d, lon2


def _unrollon(p1, p2, wrap=False):  # unroll180 == .karney._unroll2
    '''(INTERNAL) Wrap/normalize, unroll and replace longitude.
    '''
    lat, lon = p2.lat, p2.lon
    if wrap and _Wrap.normal:
        lat, lon = _Wrap.latlon(lat, lon)
    _, lon = unroll180(p1.lon, lon, wrap=True)
    if lat != p2.lat or fabs(lon - p2.lon) > EPS:
        p2 = p2.dup(lat=lat, lon=wrap180(lon))
        # p2 = p2.copy(); p2.latlon = lat, wrap180(lon)
    return p2


def _unrollon3(p1, p2, p3, wrap=False):
    '''(INTERNAL) Wrap/normalize, unroll 2 points.
    '''
    w = wrap
    if w:
        w  = _Wrap.normal
        p2 = _unrollon(p1, p2, wrap=w)
        p3 = _unrollon(p1, p3, wrap=w)
        p2 = _unrollon(p2, p3)
    return p2, p3, w  # was wrapped?


def unrollPI(rad1, rad2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in radians.

       @arg rad1: Start longitude (C{radians}).
       @arg rad2: End longitude (C{radians}).
       @kwarg wrap: If C{True}, wrap and unroll to the M{(-PI..+PI]}
                    range (C{bool}).

       @return: 2-Tuple C{(B{rad2}-B{rad1}, B{rad2})} unrolled
                (C{radians}, C{radians}).

       @see: Capability C{LONG_UNROLL} in U{GeographicLib
             <https://GeographicLib.SourceForge.io/html/python/interface.html#outmask>}.
    '''
    r = rad2 - rad1
    if wrap:
        u = wrapPI(r)
        if u != r:
            return u, (rad1 + u)
    return r, rad2


def _valueError(where, x, **kwds):
    '''(INTERNAL) Return a C{_ValueError}.
    '''
    t = _MODS.streprs.Fmt.PAREN(where.__name__, x)
    return _ValueError(t, **kwds)


class _Wrap(object):

    _normal = False  # default

    @property
    def normal(self):
        '''Get the current L{normal} setting (C{True},
           C{False} or C{None}).
        '''
        return self._normal

    @normal.setter  # PYCHOK setter!
    def normal(self, setting):
        '''Set L{normal} to C{True}, C{False} or C{None}.
        '''
        m = _MODS.formy
        t = {True:  (m.normal,        m.normal_),
             False: (self.wraplatlon, self.wraphilam),
             None:  (_passargs,      _passargs)}.get(setting, ())
        if t:
            self.latlon, self.philam = t
            self._normal = setting

    def latlonDMS2(self, lat, lon, **DMS2_kwds):
        if isstr(lat) or isstr(lon):
            kwds = _xkwds(DMS2_kwds, clipLon=0, clipLat=0)
            lat, lon = _MODS.dms.parseDMS2(lat, lon, **kwds)
        return self.latlon(lat, lon)

#   def normalatlon(self, *latlon):
#       return _MODS.formy.normal(*latlon)

#   def normalamphi(self, *philam):
#       return _MODS.formy.normal_(*philam)

    def wraplatlon(self, lat, lon):
        return wrap90(lat), wrap180(lon)

    latlon = wraplatlon  # default

    def latlon3(self, lon1, lat2, lon2, wrap):
        if wrap:
            lat2,  lon2 = self.latlon(lat2, lon2)
            lon21, lon2 = unroll180(lon1, lon2)
        else:
            lon21 = lon2 - lon1
        return lon21, lat2, lon2

    def _latlonop(self, wrap):
        if wrap and self._normal is not None:
            return self.latlon
        else:
            return _passargs

    def wraphilam(self, phi, lam):
        return wrapPI_2(phi), wrapPI(lam)

    philam = wraphilam  # default

    def philam3(self, lam1, phi2, lam2, wrap):
        if wrap:
            phi2,  lam2 = self.philam(phi2, lam2)
            lam21, lam2 = unrollPI(lam1, lam2)
        else:
            lam21 = lam2 - lam1
        return lam21, phi2, lam2

    def _philamop(self, wrap):
        if wrap and self._normal is not None:
            return self.philam
        else:
            return _passargs

    def point(self, ll, wrap=True):  # in .points._fractional, -.PointsIter.iterate, ...
        '''Return C{ll} or a copy, I{normalized} or I{wrap}'d.
        '''
        if wrap and self._normal is not None:
            lat, lon = ll.latlon
            if fabs(lon) > 180 or fabs(lat) > 90:
                _n = self.latlon
                ll = ll.copy(name=_n.__name__)
                ll.latlon = _n(lat, lon)
        return ll

_Wrap = _Wrap()  # PYCHOK singleton


# def _wrap(angle, wrap, modulo):
#     '''(INTERNAL) Angle wrapper M{((wrap-modulo)..+wrap]}.
#
#        @arg angle: Angle (C{degrees}, C{radians} or C{grades}).
#        @arg wrap: Range (C{degrees}, C{radians} or C{grades}).
#        @arg modulo: Upper limit (360 C{degrees}, PI2 C{radians} or 400 C{grades}).
#
#        @return: The B{C{angle}}, wrapped (C{degrees}, C{radians} or C{grades}).
#     '''
#     a = float(angle)
#     if not (wrap - modulo) <= a < wrap:
#         # math.fmod(-1.5, 3.14) == -1.5, but -1.5 % 3.14 == 1.64
#         # math.fmod(-1.5, 360)  == -1.5, but -1.5 % 360 == 358.5
#         a %= modulo
#         if a > wrap:
#             a -= modulo
#     return a


def wrap90(deg):
    '''Wrap degrees to M{[-90..+90]}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees90}).
    '''
    w = wrap180(deg)
    return (w - _180_0) if w > 90 else ((w + _180_0) if w < -90 else w)


def wrap180(deg):
    '''Wrap degrees to M{[-180..+180]}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees180}).
    '''
    d =  float(deg)
    w = _umod_360(d)
    if w > 180:
        w -= _360_0
    elif d < 0 and w == 180:
        w = -w
    return w


def wrap360(deg):  # see .streprs._umod_360
    '''Wrap degrees to M{[0..+360)}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees360}).
    '''
    return _umod_360(float(deg))


def wrapPI(rad):
    '''Wrap radians to M{[-PI..+PI]}.

       @arg rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI}).
    '''
    r =  float(rad)
    w = _umod_PI2(r)
    if w > PI:
        w -= PI2
    elif r < 0 and w == PI:
        w = -PI
    return w


def wrapPI2(rad):
    '''Wrap radians to M{[0..+2PI)}.

       @arg rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI2}).
    '''
    return _umod_PI2(float(rad))


def wrapPI_2(rad):
    '''Wrap radians to M{[-PI/2..+PI/2]}.

       @arg rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI_2}).
    '''
    w = wrapPI(rad)
    return (w - PI) if w > PI_2 else ((w + PI) if w < (-PI_2) else w)


# def wraplatlon(lat, lon):
#     '''Both C{wrap90(B{lat})} and C{wrap180(B{lon})}.
#     '''
#     return wrap90(lat), wrap180(lon)


def wrap_normal(*normal):
    '''Define the operation for the keyword argument C{B{wrap}=True},
       across L{pygeodesy}: I{wrap}, I{normalize} or I{no-op}.  For
       backward compatibility, the default is I{wrap}.

       @arg normal: If C{True}, I{normalize} lat- and longitude using
                    L{normal} or L{normal_}, if C{False}, I{wrap} the
                    lat- and longitude individually by L{wrap90} or
                    L{wrapPI_2} respectively L{wrap180}, L{wrapPI} or
                    if C{None}, leave lat- and longitude I{unchanged}.
                    Do not supply any value to get the current setting.

       @return: The previous L{wrap_normal} setting (C{bool} or C{None}).
    '''
    t = _Wrap.normal
    if normal:
        _Wrap.normal = normal[0]
    return t


# def wraphilam(phi, lam,):
#     '''Both C{wrapPI_2(B{phi})} and C{wrapPI(B{lam})}.
#     '''
#     return wrapPI_2(phi), wrapPI(lam)


def yard2m(yards):
    '''Convert I{UK} yards to meter.

       @arg yards: Value in yards (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{yards}}.
    '''
    return Float(yards=yards) * _M_YARD_UK

# **) MIT License
#
# Copyright (C) 2016-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
