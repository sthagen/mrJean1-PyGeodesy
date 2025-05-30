
# -*- coding: utf-8 -*-

u'''Ellipsoidal and spherical earth models.

Classes L{a_f2Tuple}, L{Ellipsoid} and L{Ellipsoid2}, an L{Ellipsoids} registry and
2 dozen functions to convert I{equatorial} radius, I{polar} radius, I{eccentricities},
I{flattenings} and I{inverse flattening}.

See module L{datums} for L{Datum} and L{Transform} information and other details.

Following is the list of predefined L{Ellipsoid}s, all instantiated lazily.

@var Ellipsoids.Airy1830: Ellipsoid(name='Airy1830', a=6377563.396, f=0.00334085, f_=299.3249646, b=6356256.90923729)
@var Ellipsoids.AiryModified: Ellipsoid(name='AiryModified', a=6377340.189, f=0.00334085, f_=299.3249646, b=6356034.44793853)
@var Ellipsoids.ATS1977: Ellipsoid(name='ATS1977', a=6378135, f=0.00335281, f_=298.257, b=6356750.30492159)
@var Ellipsoids.Australia1966: Ellipsoid(name='Australia1966', a=6378160, f=0.00335289, f_=298.25, b=6356774.71919531)
@var Ellipsoids.Bessel1841: Ellipsoid(name='Bessel1841', a=6377397.155, f=0.00334277, f_=299.1528128, b=6356078.962818)
@var Ellipsoids.BesselModified: Ellipsoid(name='BesselModified', a=6377492.018, f=0.00334277, f_=299.1528128, b=6356173.5087127)
@var Ellipsoids.CGCS2000: Ellipsoid(name='CGCS2000', a=6378137, f=0.00335281, f_=298.2572221, b=6356752.31414036)
@var Ellipsoids.Clarke1866: Ellipsoid(name='Clarke1866', a=6378206.4, f=0.00339008, f_=294.97869821, b=6356583.8)
@var Ellipsoids.Clarke1880: Ellipsoid(name='Clarke1880', a=6378249.145, f=0.00340756, f_=293.465, b=6356514.86954978)
@var Ellipsoids.Clarke1880IGN: Ellipsoid(name='Clarke1880IGN', a=6378249.2, f=0.00340755, f_=293.46602129, b=6356515)
@var Ellipsoids.Clarke1880Mod: Ellipsoid(name='Clarke1880Mod', a=6378249.145, f=0.00340755, f_=293.46630766, b=6356514.96639549)
@var Ellipsoids.CPM1799: Ellipsoid(name='CPM1799', a=6375738.7, f=0.00299052, f_=334.39, b=6356671.92557493)
@var Ellipsoids.Delambre1810: Ellipsoid(name='Delambre1810', a=6376428, f=0.00321027, f_=311.5, b=6355957.92616372)
@var Ellipsoids.Engelis1985: Ellipsoid(name='Engelis1985', a=6378136.05, f=0.00335282, f_=298.2566, b=6356751.32272154)
@var Ellipsoids.Everest1969: Ellipsoid(name='Everest1969', a=6377295.664, f=0.00332445, f_=300.8017, b=6356094.667915)
@var Ellipsoids.Everest1975: Ellipsoid(name='Everest1975', a=6377299.151, f=0.00332445, f_=300.8017255, b=6356098.14512013)
@var Ellipsoids.Fisher1968: Ellipsoid(name='Fisher1968', a=6378150, f=0.00335233, f_=298.3, b=6356768.33724438)
@var Ellipsoids.GEM10C: Ellipsoid(name='GEM10C', a=6378137, f=0.00335281, f_=298.2572236, b=6356752.31424783)
@var Ellipsoids.GPES: Ellipsoid(name='GPES', a=6378135, f=0, f_=0, b=6378135)
@var Ellipsoids.GRS67: Ellipsoid(name='GRS67', a=6378160, f=0.00335292, f_=298.24716743, b=6356774.51609071)
@var Ellipsoids.GRS80: Ellipsoid(name='GRS80', a=6378137, f=0.00335281, f_=298.2572221, b=6356752.31414035)
@var Ellipsoids.Helmert1906: Ellipsoid(name='Helmert1906', a=6378200, f=0.00335233, f_=298.3, b=6356818.16962789)
@var Ellipsoids.IAU76: Ellipsoid(name='IAU76', a=6378140, f=0.00335281, f_=298.257, b=6356755.28815753)
@var Ellipsoids.IERS1989: Ellipsoid(name='IERS1989', a=6378136, f=0.00335281, f_=298.257, b=6356751.30156878)
@var Ellipsoids.IERS1992TOPEX: Ellipsoid(name='IERS1992TOPEX', a=6378136.3, f=0.00335281, f_=298.25722356, b=6356751.61659215)
@var Ellipsoids.IERS2003: Ellipsoid(name='IERS2003', a=6378136.6, f=0.00335282, f_=298.25642, b=6356751.85797165)
@var Ellipsoids.Intl1924: Ellipsoid(name='Intl1924', a=6378388, f=0.003367, f_=297, b=6356911.94612795)
@var Ellipsoids.Intl1967: Ellipsoid(name='Intl1967', a=6378157.5, f=0.0033529, f_=298.24961539, b=6356772.2)
@var Ellipsoids.Krassovski1940: Ellipsoid(name='Krassovski1940', a=6378245, f=0.00335233, f_=298.3, b=6356863.01877305)
@var Ellipsoids.Krassowsky1940: Ellipsoid(name='Krassowsky1940', a=6378245, f=0.00335233, f_=298.3, b=6356863.01877305)
@var Ellipsoids.Maupertuis1738: Ellipsoid(name='Maupertuis1738', a=6397300, f=0.0052356, f_=191, b=6363806.28272251)
@var Ellipsoids.Mercury1960: Ellipsoid(name='Mercury1960', a=6378166, f=0.00335233, f_=298.3, b=6356784.28360711)
@var Ellipsoids.Mercury1968Mod: Ellipsoid(name='Mercury1968Mod', a=6378150, f=0.00335233, f_=298.3, b=6356768.33724438)
@var Ellipsoids.NWL1965: Ellipsoid(name='NWL1965', a=6378145, f=0.00335289, f_=298.25, b=6356759.76948868)
@var Ellipsoids.OSU86F: Ellipsoid(name='OSU86F', a=6378136.2, f=0.00335281, f_=298.2572236, b=6356751.51693008)
@var Ellipsoids.OSU91A: Ellipsoid(name='OSU91A', a=6378136.3, f=0.00335281, f_=298.2572236, b=6356751.6165948)
@var Ellipsoids.Plessis1817: Ellipsoid(name='Plessis1817', a=6376523, f=0.00324002, f_=308.64, b=6355862.93325557)
@var Ellipsoids.PZ90: Ellipsoid(name='PZ90', a=6378136, f=0.0033528, f_=298.2578393, b=6356751.36174571)
@var Ellipsoids.SGS85: Ellipsoid(name='SGS85', a=6378136, f=0.00335281, f_=298.257, b=6356751.30156878)
@var Ellipsoids.SoAmerican1969: Ellipsoid(name='SoAmerican1969', a=6378160, f=0.00335289, f_=298.25, b=6356774.71919531)
@var Ellipsoids.Sphere: Ellipsoid(name='Sphere', a=6371008.771415, f=0, f_=0, b=6371008.771415)
@var Ellipsoids.SphereAuthalic: Ellipsoid(name='SphereAuthalic', a=6371000, f=0, f_=0, b=6371000)
@var Ellipsoids.SpherePopular: Ellipsoid(name='SpherePopular', a=6378137, f=0, f_=0, b=6378137)
@var Ellipsoids.Struve1860: Ellipsoid(name='Struve1860', a=6378298.3, f=0.00339294, f_=294.73, b=6356657.14266956)
@var Ellipsoids.WGS60: Ellipsoid(name='WGS60', a=6378165, f=0.00335233, f_=298.3, b=6356783.28695944)
@var Ellipsoids.WGS66: Ellipsoid(name='WGS66', a=6378145, f=0.00335289, f_=298.25, b=6356759.76948868)
@var Ellipsoids.WGS72: Ellipsoid(name='WGS72', a=6378135, f=0.00335278, f_=298.26, b=6356750.52001609)
@var Ellipsoids.WGS84: Ellipsoid(name='WGS84', a=6378137, f=0.00335281, f_=298.25722356, b=6356752.31424518)
@var Ellipsoids.WGS84_NGS: Ellipsoid(name='WGS84_NGS', a=6378137, f=0.00335281, f_=298.2572221, b=6356752.31414035)
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

# from pygeodesy.albers import AlbersEqualAreaCylindrical  # _MODS
from pygeodesy.basics import copysign0, isbool, _isin, isint,  typename
from pygeodesy.constants import EPS, EPS0, EPS02, EPS1, INF, NINF, PI4, PI_2, PI_3, R_M, R_MA, R_FM, \
                               _EPSqrt, _EPStol as _TOL, _floatuple as _T, _isfinite, _over, \
                               _0_0s, _0_0, _0_5, _1_0, _1_EPS, _2_0, _4_0, _90_0, \
                               _0_25, _3_0  # PYCHOK used!
from pygeodesy.errors import _AssertionError, IntersectionError, _ValueError, _xattr, _xkwds_not
from pygeodesy.fmath import cbrt, cbrt2, fdot, Fhorner, fpowers, hypot, hypot_, \
                            hypot1, hypot2, sqrt3,  Fsum
# from pygeodesy.fsums import Fsum  # from .fmath
# from pygeodesy.internals import typename  # from .basics
from pygeodesy.interns import NN, _a_, _Airy1830_, _AiryModified_, _b_, _Bessel1841_, _beta_, \
                             _Clarke1866_, _Clarke1880IGN_, _DMAIN_, _DOT_, _f_, _GRS80_, \
                             _height_, _Intl1924_, _incompatible_, _invalid_, _Krassovski1940_, \
                             _Krassowsky1940_, _lat_, _meridional_, _negative_, _not_finite_, \
                             _prime_vertical_, _radius_, _Sphere_, _SPACE_, _vs_, _WGS72_, _WGS84_
# from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS  # from .named
from pygeodesy.named import _lazyNamedEnumItem as _lazy, _name__, _NamedEnum, \
                                _NamedEnumItem, _NamedTuple, _Pass,  _ALL_LAZY, _MODS
from pygeodesy.namedTuples import Distance2Tuple, Vector3Tuple, Vector4Tuple
from pygeodesy.props import deprecated_Property_RO, Property_RO, property_doc_, \
                            deprecated_property_RO, property_RO, property_ROver
from pygeodesy.streprs import Fmt, fstr, instr, strs, unstr
# from pygeodesy.triaxials import _hartzell3  # _MODS
from pygeodesy.units import Azimuth, Bearing, Distance, Float, Float_, Height, Lamd, Lat, \
                            Meter, Meter2, Meter3, Phi, Phid, Radius, Radius_, Scalar
from pygeodesy.utily import atan1, atan1d, atan2b, degrees90, m2radians, radians2m, sincos2d

from math import asinh, atan, atanh, cos, degrees, exp, fabs, radians, sin, sinh, sqrt, tan  # as _tan

__all__ = _ALL_LAZY.ellipsoids
__version__ = '25.05.12'

_f_0_0    = Float(f =_0_0)  # zero flattening
_f__0_0   = Float(f_=_0_0)  # zero inverse flattening
# see U{WGS84_f<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Constants.html>}
_f__WGS84 = Float(f_=_1_0 / (1000000000 / 298257223563))  # 298.257223562_999_97 vs 298.257223563


def _aux(lat, inverse, auxLat, clip=90):
    '''Return a named auxiliary latitude in C{degrees}.
    '''
    return Lat(lat, clip=clip, name=_lat_ if inverse else typename(auxLat))


def _s2_c2(phi):
    '''(INTERNAL) Return 2-tuple C{(sin(B{phi})**2, cos(B{phi})**2)}.
    '''
    if phi:
        s2 = sin(phi)**2
        if s2 > EPS:
            c2 = _1_0 - s2
            if c2 > EPS:
                if c2 < EPS1:
                    return s2, c2
            else:
                return _1_0, _0_0  # phi == PI_2
    return _0_0, _1_0  # phi == 0


class a_f2Tuple(_NamedTuple):
    '''2-Tuple C{(a, f)} specifying an ellipsoid by I{equatorial}
       radius C{a} in C{meter} and scalar I{flattening} C{f}.

       @see: Class L{Ellipsoid2}.
    '''
    _Names_ = (_a_,   _f_)  # name 'f' not 'f_'
    _Units_ = (_Pass, _Pass)

    def __new__(cls, a, f, **name):
        '''New L{a_f2Tuple} ellipsoid specification.

           @arg a: Equatorial radius (C{scalar} > 0).
           @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: An L{a_f2Tuple}C{(a, f)}.

           @raise UnitError: Invalid B{C{a}} or B{C{f}}.

           @note: C{abs(B{f}) < }L{EPS<pygeodesy.constants.EPS>} is
                  forced to C{B{f}=0}, I{spherical}.

           @note: Negative C{B{f}} produces a I{prolate} ellipsoid.
        '''
        a = Radius_(a=a)  # low=EPS, high=None
        f = Float_( f=f, low=None, high=EPS1)
        if fabs(f) < EPS:  # force spherical
            f = _f_0_0
        return _NamedTuple.__new__(cls, a, f, **name)

    @Property_RO
    def b(self):
        '''Get the I{polar} radius (C{meter}), M{a * (1 - f)}.
        '''
        return a_f2b(self.a, self.f)  # PYCHOK .a and .f

    def ellipsoid(self, **name):
        '''Return an L{Ellipsoid} for this 2-tuple C{(a, f)}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise NameError: A registered C{ellipsoid} with the
                             same B{C{name}} already exists.
        '''
        return Ellipsoid(self.a, f=self.f, name=self._name__(name))  # PYCHOK .a and .f

    @Property_RO
    def f_(self):
        '''Get the I{inverse} flattening (C{scalar}), M{1 / f} == M{a / (a - b)}.
        '''
        return f2f_(self.f)  # PYCHOK .f


class Circle4Tuple(_NamedTuple):
    '''4-Tuple C{(radius, height, lat, beta)} of the C{radius} and C{height},
       both conventionally in C{meter} of a parallel I{circle of latitude} at
       (geodetic) latitude C{lat} and the I{parametric (or reduced) auxiliary
       latitude} C{beta}, both in C{degrees90}.

       The C{height} is the (signed) distance along the z-axis between the
       parallel and the equator.  At near-polar C{lat}s, the C{radius} is C{0},
       the C{height} is the ellipsoid's (signed) polar radius and C{beta}
       equals C{lat}.
    '''
    _Names_ = (_radius_, _height_, _lat_, _beta_)
    _Units_ = ( Radius,   Height,   Lat,   Lat)


class Curvature2Tuple(_NamedTuple):
    '''2-Tuple C{(meridional, prime_vertical)} of radii of curvature, both in
       C{meter}, conventionally.
    '''
    _Names_ = (_meridional_, _prime_vertical_)
    _Units_ = ( Meter,        Meter)

    @property_RO
    def transverse(self):
        '''Get this I{prime_vertical}, aka I{transverse} radius of curvature.
        '''
        return self.prime_vertical


class Ellipsoid(_NamedEnumItem):
    '''Ellipsoid with I{equatorial} and I{polar} radii, I{flattening}, I{inverse
       flattening} and other, often used, I{cached} attributes, supporting
       I{oblate} and I{prolate} ellipsoidal and I{spherical} earth models.
    '''
    _a  = 0  # equatorial radius, semi-axis (C{meter})
    _b  = 0  # polar radius, semi-axis (C{meter}): a * (f - 1) / f
    _f  = 0  # (1st) flattening: (a - b) / a
    _f_ = 0  # inverse flattening: 1 / f = a / (a - b)

    _geodsolve  = NN  # means, use PYGEODESY_GEODSOLVE
    _KsOrder    =  8  # Krüger series order (4, 6 or 8)
    _rhumbsolve = NN  # means, use PYGEODESY_RHUMBSOLVE

    def __init__(self, a, b=None, f_=None, f=None, **name):
        '''New L{Ellipsoid} from the I{equatorial} radius I{and} either
           the I{polar} radius or I{inverse flattening} or I{flattening}.

           @arg a: Equatorial radius, semi-axis (C{meter}).
           @arg b: Optional polar radius, semi-axis (C{meter}).
           @arg f_: Inverse flattening: M{a / (a - b)} (C{float} >>> 1.0).
           @arg f: Flattening: M{(a - b) / a} (C{scalar}, near zero for
                   spherical).
           @kwarg name: Optional, unique C{B{name}=NN} (C{str}).

           @raise NameError: Ellipsoid with the same B{C{name}} already exists.

           @raise ValueError: Invalid B{C{a}}, B{C{b}}, B{C{f_}} or B{C{f}} or
                              B{C{f_}} and B{C{f}} are incompatible.

           @note: M{abs(f_) > 1 / EPS} or M{abs(1 / f_) < EPS} is forced
                  to M{1 / f_ = 0}, spherical.
        '''
        ff_ =  f, f_  # assertion below
        n   = _name__(**name) if name else NN
        try:
            a = Radius_(a=a)  # low=EPS
            if not _isfinite(a):
                raise ValueError(_SPACE_(_a_, _not_finite_))

            if b:  # not _isin(b, None, _0_0)
                b  = Radius_(b=b)  # low=EPS
                f  = a_b2f(a, b) if f is None else Float(f=f)
                f_ = f2f_(f) if f_ is None else Float(f_=f_)
            elif f is not None:
                f  = Float(f=f)
                b  = a_f2b(a, f)
                f_ = f2f_(f) if f_ is None else Float(f_=f_)
            elif f_:
                f_ = Float(f_=f_)
                b  = a_f_2b(a, f_)  # a * (f_ - 1) / f_
                f  = f_2f(f_)
            else:  # only a, spherical
                f_ = f = 0
                b  = a  # superfluous

            if not f < _1_0:  # sanity check, see .ecef.Ecef.__init__
                raise ValueError(_SPACE_(_f_, _invalid_))
            if not _isfinite(b):
                raise ValueError(_SPACE_(_b_, _not_finite_))

            if fabs(f) < EPS or a == b or not f_:  # spherical
                b  =  a
                f  = _f_0_0
                f_ = _f__0_0

        except (TypeError, ValueError) as x:
            d = _xkwds_not(None, b=b, f_=f_, f=f)
            t =  instr(self, a=a, name=n, **d)
            raise _ValueError(t, cause=x)

        self._a  = a
        self._b  = b
        self._f  = f
        self._f_ = f_

        self._register(Ellipsoids, n)

        if f and f_:  # see test/testEllipsoidal
            d = dict(eps=_TOL)
            if None in ff_:  # both f_ and f given
                d.update(Error=_ValueError, txt=_incompatible_)
            self._assert(_1_0 / f,  f_=f_, **d)
            self._assert(_1_0 / f_, f =f,  **d)
        self._assert(self.b2_a2, e21=self.e21, eps=EPS)

    def __eq__(self, other):
        '''Compare this and an other ellipsoid.

           @arg other: The other ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Ellipsoid) and
                                  self.a == other.a and
                                 (self.f == other.f or self.b == other.b))

    def __hash__(self):
        return self._hash  # memoized

    @Property_RO
    def a(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self._a

    equatoradius = a  # = Requatorial

    @Property_RO
    def a2(self):
        '''Get the I{equatorial} radius I{squared} (C{meter} I{squared}), M{a**2}.
        '''
        return Meter2(a2=self.a**2)

    @Property_RO
    def a2_(self):
        '''Get the inverse of the I{equatorial} radius I{squared} (C{meter} I{squared}), M{1 / a**2}.
        '''
        return Float(a2_=_1_0 / self.a2)

    @Property_RO
    def a_b(self):
        '''Get the ratio I{equatorial} over I{polar} radius (C{float}), M{a / b} == M{1 / (1 - f)}.
        '''
        return Float(a_b=self.a / self.b if self.f else _1_0)

    @Property_RO
    def a2_b(self):
        '''Get the I{polar} meridional (or polar) radius of curvature (C{meter}), M{a**2 / b}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
                 and U{Moritz, H. (1980), Geodetic Reference System 1980
                 <https://WikiPedia.org/wiki/Earth_radius#cite_note-Moritz-2>}.

           @note: Symbol C{c} is used by IUGG and IERS for the U{polar radius of curvature
                  <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}, see L{c2}
                  and L{R2} or L{Rauthalic}.
        '''
        return Radius(a2_b=self.a2 / self.b if self.f else self.a)  # = rocPolar

    @Property_RO
    def a2_b2(self):
        '''Get the ratio I{equatorial} over I{polar} radius I{squared} (C{float}),
           M{(a / b)**2} == M{1 / (1 - e**2)} == M{1 / (1 - e2)} == M{1 / e21}.
        '''
        return Float(a2_b2=self.a_b**2 if self.f else _1_0)

    @Property_RO
    def a_f(self):
        '''Get the I{equatorial} radius and I{flattening} (L{a_f2Tuple}), see method C{toEllipsoid2}.
        '''
        return a_f2Tuple(self.a, self.f, name=self.name)

    a_f2 = a_f  # synonym

    @Property_RO
    def A(self):
        '''Get the UTM I{meridional (or rectifying)} radius (C{meter}).

           @note: C{A * PI / 2} ≈= L{L<Ellipsoid.L>}, the I{quarter meridian}.

           @see: I{Meridian arc unit} U{Q<https://StudyLib.net/doc/7443565/>},
                 the mean, meridional length I{per radian}.
        '''
        A, n = self.a, self.n
        if n:
            d = (n + _1_0) * 1048576 / A
            if d:  # use 6 n**2 terms, half-way between the _KsOrder's 4, 6, 8
                # <https://GeographicLib.SourceForge.io/C++/doc/tmseries30.html>
                # <https://GeographicLib.SourceForge.io/C++/doc/transversemercator.html> and
                # <https://www.MyGeodesy.id.AU/documents/Karney-Krueger%20equations.pdf> (3)
                # A *= fhorner(n**2, 1048576, 262144, 16384, 4096, 1600, 784, 441) / 1048576) / (1 + n)
                A = Radius(A=Fhorner(n**2, 1048576, 262144, 16384, 4096, 1600, 784, 441).fover(d))
        return A
#       # Moritz, H. <https://Geodesy.Geology.Ohio-State.EDU/course/refpapers/00740128.pdf>
#       # q =  4 / self.rocPolar
#       # Q = (1 - 3 / 4 * e'2 + 45 / 64 * e'4 - 175 / 256 * e'6 + 11025 / 16384 * e'8) * rocPolar
#       #   = (4 + e'2 * (-3 + e'2 * (45 / 16 + e'2 * (-175 / 64 + e'2 * 11025 / 4096)))) / q
#       return Radius(Q=Fhorner(self.e22, 4, -3, 45 / 16, -175 / 64, 11025 / 4096).fover(q))

    @Property_RO
    def _albersCyl(self):
        '''(INTERNAL) Helper for C{auxAuthalic}.
        '''
        return _MODS.albers.AlbersEqualAreaCylindrical(datum=self, name=self.name)

    @Property_RO
    def AlphaKs(self):
        '''Get the I{Krüger} U{Alpha series coefficients<https://GeographicLib.SourceForge.io/C++/doc/tmseries30.html>} (C{KsOrder}C{-tuple}).
        '''
        return self._Kseries(  # XXX int/int quotients may require  from __future__ import division as _; del _  # noqa: E702 ;
            #   n    n**2   n**3      n**4         n**5            n**6                 n**7                     n**8
            _T(1/2, -2/3,   5/16,    41/180,    -127/288,       7891/37800,         72161/387072,        -18975107/50803200),
                 _T(13/48, -3/5,    557/1440,    281/630,   -1983433/1935360,       13769/28800,         148003883/174182400),     # PYCHOK unaligned
                        _T(61/240, -103/140,   15061/26880,   167603/181440,    -67102379/29030400,       79682431/79833600),      # PYCHOK unaligned
                               _T(49561/161280, -179/168,    6601661/7257600,       97445/49896,      -40176129013/7664025600),    # PYCHOK unaligned
                                            _T(34729/80640, -3418889/1995840,    14644087/9123840,      2605413599/622702080),     # PYCHOK unaligned
                                                        _T(212378941/319334400, -30705481/10378368,   175214326799/58118860800),   # PYCHOK unaligned
                                                                            _T(1522256789/1383782400, -16759934899/3113510400),    # PYCHOK unaligned
                                                                                                  _T(1424729850961/743921418240))  # PYCHOK unaligned

    @Property_RO
    def area(self):
        '''Get the ellipsoid's surface area (C{meter} I{squared}), M{4 * PI * c2}.

           @see: Properties L{areax}, L{c2} and L{R2} and functions
                 L{ellipsoidalExact.areaOf} and L{ellipsoidalKarney.areaOf}.
        '''
        return Meter2(area=self.c2 * PI4)

    @Property_RO
    def areax(self):
        '''Get the ellipsoid's surface area (C{meter} I{squared}), M{4 * PI * c2x}, more
           accurate for very I{oblate} ellipsoids.

           @see: Properties L{area}, L{c2x} and L{R2x}, class L{GeodesicExact} and
                 functions L{ellipsoidalExact.areaOf} and L{ellipsoidalKarney.areaOf}.
        '''
        return Meter2(areax=self.c2x * PI4)

    def _assert(self, val, eps=_TOL, f0=_0_0, Error=_AssertionError, txt=NN, **name_value):
        '''(INTERNAL) Assert a C{name=value} vs C{val}.
        '''
        for n, v in name_value.items():
            if fabs(v - val) > eps:  # PYCHOK no cover
                t = (v, _vs_, val)
                t = _SPACE_.join(strs(t, prec=12, fmt=Fmt.g))
                t =  Fmt.EQUAL(self._DOT_(n), t)
                raise Error(t, txt=txt or Fmt.exceeds_eps(eps))
            return Float(v if self.f else f0, name=n)
        raise Error(unstr(self._DOT_(typename(self._assert)), val,
                                eps=eps, f0=f0, **name_value))

    def auxAuthalic(self, lat, inverse=False):
        '''Compute the I{authalic} auxiliary latitude or the I{inverse} thereof.

           @arg lat: The geodetic (or I{authalic}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{authalic} and
                           return the geodetic latitude (C{bool}).

           @return: The I{authalic} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/AuthalicLatitude<https://GeographicLib.SourceForge.io/
                 C++/doc/classGeographicLib_1_1Ellipsoid.html>}, U{Authalic latitude
                 <https://WikiPedia.org/wiki/Latitude#Authalic_latitude>}, and
                 U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, p 16.
        '''
        if self.f:
            f   = self._albersCyl._tanf if inverse else \
                  self._albersCyl._txif  # PYCHOK attr
            lat = atan1d(f(tan(Phid(lat))))  # PYCHOK attr
        return _aux(lat, inverse, Ellipsoid.auxAuthalic)

    def auxConformal(self, lat, inverse=False):
        '''Compute the I{conformal} auxiliary latitude or the I{inverse} thereof.

           @arg lat: The geodetic (or I{conformal}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{conformal} and
                           return the geodetic latitude (C{bool}).

           @return: The I{conformal} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/ConformalLatitude<https://GeographicLib.SourceForge.io/
                 C++/doc/classGeographicLib_1_1Ellipsoid.html>}, U{Conformal latitude
                 <https://WikiPedia.org/wiki/Latitude#Conformal_latitude>}, and
                 U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 15-16.
        '''
        if self.f:
            f   = self.es_tauf if inverse else self.es_taupf  # PYCHOK attr
            lat = atan1d(f(tan(Phid(lat))))  # PYCHOK attr
        return _aux(lat, inverse, Ellipsoid.auxConformal)

    def auxGeocentric(self, lat, inverse=False):
        '''Compute the I{geocentric} auxiliary latitude or the I{inverse} thereof.

           @arg lat: The geodetic (or I{geocentric}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the geocentric and
                           return the I{geocentric} latitude (C{bool}).

           @return: The I{geocentric} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/GeocentricLatitude<https://GeographicLib.SourceForge.io/
                 C++/doc/classGeographicLib_1_1Ellipsoid.html>}, U{Geocentric latitude
                 <https://WikiPedia.org/wiki/Latitude#Geocentric_latitude>}, and
                 U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 17-18.
        '''
        if self.f:
            f   = self.a2_b2 if inverse else self.b2_a2
            lat = atan1d(tan(Phid(lat)) * f)
        return _aux(lat, inverse, Ellipsoid.auxGeocentric)

    def auxIsometric(self, lat, inverse=False):
        '''Compute the I{isometric} auxiliary latitude or the I{inverse} thereof.

           @arg lat: The geodetic (or I{isometric}) latitude (C{degrees}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{isometric} and
                           return the geodetic latitude (C{bool}).

           @return: The I{isometric} (or geodetic) latitude in C{degrees}.

           @note: The I{isometric} latitude for geodetic C{+/-90} is far
                  outside the C{[-90..+90]} range but the inverse
                  thereof is the original geodetic latitude.

           @see: U{Inverse-/IsometricLatitude<https://GeographicLib.SourceForge.io/
                 C++/doc/classGeographicLib_1_1Ellipsoid.html>}, U{Isometric latitude
                 <https://WikiPedia.org/wiki/Latitude#Isometric_latitude>}, and
                 U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 15-16.
        '''
        if self.f:
            r   = Phid(lat, clip=0)
            lat = degrees(atan1(self.es_tauf(sinh(r))) if inverse else
                          asinh(self.es_taupf(tan(r))))
        # clip=0, since auxIsometric(+/-90) is far outside [-90..+90]
        return _aux(lat, inverse, Ellipsoid.auxIsometric, clip=0)

    def auxParametric(self, lat, inverse=False):
        '''Compute the I{parametric} auxiliary latitude or the I{inverse} thereof.

           @arg lat: The geodetic (or I{parametric}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{parametric} and
                           return the geodetic latitude (C{bool}).

           @return: The I{parametric} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/ParametricLatitude<https://GeographicLib.SourceForge.io/
                 C++/doc/classGeographicLib_1_1Ellipsoid.html>}, U{Parametric latitude
                 <https://WikiPedia.org/wiki/Latitude#Parametric_(or_reduced)_latitude>},
                 and U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, p 18.
        '''
        if self.f:
            lat = self._beta(Lat(lat), inverse=inverse)
        return _aux(lat, inverse, Ellipsoid.auxParametric)

    auxReduced = auxParametric  # synonymous

    def auxRectifying(self, lat, inverse=False):
        '''Compute the I{rectifying} auxiliary latitude or the I{inverse} thereof.

           @arg lat: The geodetic (or I{rectifying}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{rectifying} and
                           return the geodetic latitude (C{bool}).

           @return: The I{rectifying} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/RectifyingLatitude<https://GeographicLib.SourceForge.io/
                 C++/doc/classGeographicLib_1_1Ellipsoid.html>}, U{Rectifying latitude
                 <https://WikiPedia.org/wiki/Latitude#Rectifying_latitude>}, and
                 U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 16-17.
        '''
        if self.f:
            lat = Lat(lat)
            if 0 < fabs(lat) < _90_0:
                if inverse:
                    e = self._elliptic_e22
                    d = degrees90(e.fEinv(e.cE * lat / _90_0))
                    lat =  self.auxParametric(d, inverse=True)
                else:
                    lat = _over(self.Llat(lat), self.L) * _90_0
        return _aux(lat, inverse, Ellipsoid.auxRectifying)

    @Property_RO
    def b(self):
        '''Get the I{polar} radius, semi-axis (C{meter}).
        '''
        return self._b

    polaradius = b  # = Rpolar

    @Property_RO
    def b_a(self):
        '''Get the ratio I{polar} over I{equatorial} radius (C{float}), M{b / a == f1 == 1 - f}.

           @see: Property L{f1}.
        '''
        return self._assert(self.b / self.a, b_a=self.f1, f0=_1_0)

    @Property_RO
    def b2(self):
        '''Get the I{polar} radius I{squared} (C{float}), M{b**2}.
        '''
        return Meter2(b2=self.b**2)

    @Property_RO
    def b2_a(self):
        '''Get the I{equatorial} meridional radius of curvature (C{meter}), M{b**2 / a}, see C{rocMeridional}C{(0)}.

           @see: U{Radii of Curvature<https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return Radius(b2_a=_over(self.b2, self.a) if self.f else self.b)

    @Property_RO
    def b2_a2(self):
        '''Get the ratio I{polar} over I{equatorial} radius I{squared} (C{float}), M{(b / a)**2}
           == M{(1 - f)**2} == M{1 - e**2} == C{e21}.
        '''
        return Float(b2_a2=self.b_a**2 if self.f else _1_0)

    def _beta(self, lat, inverse=False):
        '''(INTERNAL) Get the I{parametric (or reduced) auxiliary latitude} or inverse thereof.
        '''
        s, c = sincos2d(lat)  # like Karney's tand(lat)
        s *= self.a_b if inverse else self.b_a
        return atan1d(s, c)

    @Property_RO
    def BetaKs(self):
        '''Get the I{Krüger} U{Beta series coefficients<https://GeographicLib.SourceForge.io/C++/doc/tmseries30.html>} (C{KsOrder}C{-tuple}).
        '''
        return self._Kseries(  # XXX int/int quotients may require  from __future__ import division as _; del _  # noqa: E702 ;
            #   n    n**2  n**3     n**4        n**5            n**6                 n**7                   n**8
            _T(1/2, -2/3, 37/96,   -1/360,    -81/512,      96199/604800,     -5406467/38707200,      7944359/67737600),
                  _T(1/48, 1/15, -437/1440,    46/105,   -1118711/3870720,       51841/1209600,      24749483/348364800),      # PYCHOK unaligned
                       _T(17/480, -37/840,   -209/4480,      5569/90720,       9261899/58060800,     -6457463/17740800),       # PYCHOK unaligned
                              _T(4397/161280, -11/504,    -830251/7257600,      466511/2494800,     324154477/7664025600),     # PYCHOK unaligned
                                          _T(4583/161280, -108847/3991680,    -8005831/63866880,     22894433/124540416),      # PYCHOK unaligned
                                                      _T(20648693/638668800, -16363163/518918400, -2204645983/12915302400),    # PYCHOK unaligne
                                                                          _T(219941297/5535129600, -497323811/12454041600),    # PYCHOK unaligned
                                                                                              _T(191773887257/3719607091200))  # PYCHOK unaligned

    @deprecated_Property_RO
    def c(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{R2} or C{Rauthalic}.'''
        return self.R2

    @Property_RO
    def c2(self):
        '''Get the I{authalic} earth radius I{squared} (C{meter} I{squared}).

           @see: Properties L{c2x}, L{area}, L{R2}, L{Rauthalic}, I{Karney's} U{equation (60)
                 <https://Link.Springer.com/article/10.1007%2Fs00190-012-0578-z>} and C++ U{Ellipsoid.Area
                 <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Ellipsoid.html>},
                 U{Authalic radius<https://WikiPedia.org/wiki/Earth_radius#Authalic_radius>}, U{Surface area
                 <https://WikiPedia.org/wiki/Ellipsoid>} and U{surface area
                 <https://www.Numericana.com/answer/geometry.htm#oblate>}.
        '''
        return self._c2f(False)

    @Property_RO
    def c2x(self):
        '''Get the I{authalic} earth radius I{squared} (C{meter} I{squared}), more accurate for very I{oblate}
           ellipsoids.

           @see: Properties L{c2}, L{areax}, L{R2x}, L{Rauthalicx}, class L{GeodesicExact} and I{Karney}'s comments at C++
                 attribute U{GeodesicExact._c2<https://GeographicLib.SourceForge.io/C++/doc/GeodesicExact_8cpp_source.html>}.
        '''
        return self._c2f(True)

    def _c2f(self, c2x):
        '''(INTERNAL) Helper for C{.c2} and C{.c2x}.
        '''
        f, c2 = self.f, self.b2
        if f:
            e = self.e
            if e > EPS0:
                if f > 0:  # .isOblate
                    c2 *= (asinh(sqrt(self.e22abs)) if c2x else atanh(e)) / e
                elif f < 0:  # .isProlate
                    c2 *=  atan1(e) / e  # XXX asin?
            c2 = Meter2(c2=(self.a2 + c2) * _0_5)
        return c2

    def circle4(self, lat):
        '''Get the equatorial or a parallel I{circle of latitude}.

           @arg lat: Geodetic latitude (C{degrees90}, C{str}).

           @return: A L{Circle4Tuple}C{(radius, height, lat, beta)}.

           @raise RangeError: Latitude B{C{lat}} outside valid range and
                              L{rangerrors<pygeodesy.rangerrors>} is C{True}.

           @raise TypeError: Invalid B{C{lat}}.

           @raise ValueError: Invalid B{C{lat}}.

           @see: Definition of U{I{p} and I{z} under B{Parametric (or reduced) latitude}
                 <https://WikiPedia.org/wiki/Latitude>}, I{Karney's} C++ U{CircleRadius and CircleHeight
                 <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Ellipsoid.html>}
                 and method C{Rlat}.
        '''
        lat = Lat(lat)
        if lat:
            b = lat
            if fabs(lat) < _90_0:
                if self.f:
                    b = self._beta(lat)
                z, r = sincos2d(b)
                r *= self.a
                z *= self.b
            else:  # near-polar
                r, z = _0_0, copysign0(self.b, lat)
        else:  # equator
            r = self.a
            z = lat = b = _0_0
        return Circle4Tuple(r, z, lat, b)

    def degrees2m(self, deg, lat=0):
        '''Convert an angle to the distance along the equator or
           along a parallel of (geodetic) latitude.

           @arg deg: The angle (C{degrees}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Distance (C{meter}, same units as the equatorial
                    and polar radii) or C{0} for near-polar B{C{lat}}.

           @raise RangeError: Latitude B{C{lat}} outside valid range and
                              L{rangerrors<pygeodesy.rangerrors>} is C{True}.

           @raise ValueError: Invalid B{C{deg}} or B{C{lat}}.
        '''
        return self.radians2m(radians(deg), lat=lat)

    def distance2(self, lat0, lon0, lat1, lon1):
        '''I{Approximate} the distance and (initial) bearing between
           two points based on the U{local, flat earth approximation
           <https://www.EdWilliams.org/avform.htm#flat>} aka U{Hubeny
           <https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

           I{Suitable only for distances of several hundred Km or Miles
           and only between points not near-polar}.

           @arg lat0: From latitude (C{degrees}).
           @arg lon0: From longitude (C{degrees}).
           @arg lat1: To latitude (C{degrees}).
           @arg lon1: To longitude (C{degrees}).

           @return: A L{Distance2Tuple}C{(distance, initial)} with C{distance}
                    in same units as this ellipsoid's axes.

           @note: The meridional and prime_vertical radii of curvature are
                  taken and scaled I{at the initial latitude}, see C{roc2}.

           @see: Function L{pygeodesy.flatLocal}/L{pygeodesy.hubeny}.
        '''
        phi0 = Phid(lat0=lat0)
        m, n = self.roc2_(phi0, scaled=True)
        m *= Phid(lat1=lat1) - phi0
        n *= Lamd(lon1=lon1) - Lamd(lon0=lon0)
        return Distance2Tuple(hypot(m, n), atan2b(n, m))

    @Property_RO
    def e(self):
        '''Get the I{unsigned, (1st) eccentricity} (C{float}), M{sqrt(1 - (b / a)**2))}, see C{a_b2e}.

           @see: Property L{es}.
        '''
        return Float(e=sqrt(self.e2abs) if self.e2 else _0_0)

    @deprecated_Property_RO
    def e12(self):  # see property ._e12
        '''DEPRECATED, use property C{e21}.'''
        return self.e21

#   @Property_RO
#   def _e12(self):  # see property ._elliptic_e12
#       # (INTERNAL) until e12 above can be replaced with e21.
#       return self.e2 / (_1_0 - self.e2)  # see I{Karney}'s Ellipsoid._e12 = e2 / (1 - e2)

    @Property_RO
    def e2(self):
        '''Get the I{signed, (1st) eccentricity squared} (C{float}), M{f * (2 - f)
           == 1 - (b / a)**2}, see C{a_b2e2}.
        '''
        return self._assert(a_b2e2(self.a, self.b), e2=f2e2(self.f))

    @Property_RO
    def e2abs(self):
        '''Get the I{unsigned, (1st) eccentricity squared} (C{float}).
        '''
        return fabs(self.e2)

    @Property_RO
    def e21(self):
        '''Get 1 less I{1st eccentricity squared} (C{float}), M{1 - e**2}
           == M{1 - e2} == M{(1 - f)**2} == M{b**2 / a**2}, see C{b2_a2}.
        '''
        return self._assert((_1_0 - self.f)**2, e21=_1_0 - self.e2, f0=_1_0)

#   _e2m   = e21  # see I{Karney}'s Ellipsoid._e2m = 1 - _e2
    _1_e21 = a2_b2  # == M{1 / e21} == M{1 / (1 - e**2)}

    @Property_RO
    def e22(self):
        '''Get the I{signed, 2nd eccentricity squared} (C{float}), M{e2 / (1 - e2)
           == e2 / (1 - f)**2 == (a / b)**2 - 1}, see C{a_b2e22}.
        '''
        return self._assert(a_b2e22(self.a, self.b), e22=f2e22(self.f))

    @Property_RO
    def e22abs(self):
        '''Get the I{unsigned, 2nd eccentricity squared} (C{float}).
        '''
        return fabs(self.e22)

    @Property_RO
    def e32(self):
        '''Get the I{signed, 3rd eccentricity squared} (C{float}), M{e2 / (2 - e2)
           == (a**2 - b**2) / (a**2 + b**2)}, see C{a_b2e32}.
        '''
        return self._assert(a_b2e32(self.a, self.b), e32=f2e32(self.f))

    @Property_RO
    def e32abs(self):
        '''Get the I{unsigned, 3rd eccentricity squared} (C{float}).
        '''
        return fabs(self.e32)

    @Property_RO
    def e4(self):
        '''Get the I{unsignd, (1st) eccentricity} to 4th power (C{float}), M{e**4 == e2**2}.
        '''
        return Float(e4=self.e2**2 if self.e2 else _0_0)

    eccentricity     = e    # eccentricity
#   eccentricity2    = e2   # eccentricity squared
    eccentricity1st2 = e2   # first eccentricity squared, signed
    eccentricity2nd2 = e22  # second eccentricity squared, signed
    eccentricity3rd2 = e32  # third eccentricity squared, signed

    def ecef(self, Ecef=None):
        '''Return U{ECEF<https://WikiPedia.org/wiki/ECEF>} converter.

           @kwarg Ecef: ECEF class to use, default L{EcefKarney}.

           @return: An ECEF converter for this C{ellipsoid}.

           @raise TypeError: Invalid B{C{Ecef}}.

           @see: Module L{pygeodesy.ecef}.
        '''
        return _MODS.ecef._4Ecef(self, Ecef)

    @Property_RO
    def _elliptic_e12(self):  # see I{Karney}'s Ellipsoid._e12
        '''(INTERNAL) Elliptic helper for C{Rhumb}.
        '''
        e12 = _over(self.e2, self.e2 - _1_0)  # NOT DEPRECATED .e12!
        return _MODS.elliptic.Elliptic(e12)

    @Property_RO
    def _elliptic_e22(self):  # aka ._elliptic_ep2
        '''(INTERNAL) Elliptic helper for C{auxRectifying}, C{L}, C{Llat}.
        '''
        return _MODS.elliptic.Elliptic(-self.e22abs)  # complex

    equatoradius = a  # Requatorial

    def e2s(self, s):
        '''Compute norm M{sqrt(1 - e2 * s**2)}.

           @arg s: Sine value (C{scalar}).

           @return: Norm (C{float}).

           @raise ValueError: Invalid B{C{s}}.
        '''
        return sqrt(self.e2s2(s)) if self.e2 else _1_0

    def e2s2(self, s):
        '''Compute M{1 - e2 * s**2}.

           @arg s: Sine value (C{scalar}).

           @return: Result (C{float}).

           @raise ValueError: Invalid B{C{s}}.
        '''
        r = _1_0
        if self.e2:
            try:
                r -= self.e2 * Scalar(s=s)**2
                if r < 0:
                    raise ValueError(_negative_)
            except (TypeError, ValueError) as x:
                t = self._DOT_(typename(Ellipsoid.e2s2))
                raise _ValueError(t, s, cause=x)
        return r

    @Property_RO
    def es(self):
        '''Get the I{signed (1st) eccentricity} (C{float}).

           @see: Property L{e}.
        '''
        # note, self.e is always non-negative
        return Float(es=copysign0(self.e, self.f))  # see .ups

    def es_atanh(self, x):
        '''Compute M{es * atanh(es * x)} or M{-es * atan(es * x)}
           for I{oblate} respectively I{prolate} ellipsoids where
           I{es} is the I{signed} (1st) eccentricity.

           @raise ValueError: Invalid B{C{x}}.

           @see: Function U{Math::eatanhe<https://GeographicLib.SourceForge.io/
                 C++/doc/classGeographicLib_1_1Math.html>}.
        '''
        return self._es_atanh(Scalar(x=x)) if self.f else _0_0

    def _es_atanh(self, x):  # see .albers._atanhee, .AuxLat._atanhee
        '''(INTERNAL) Helper for .es_atanh, ._es_taupf2 and ._exp_es_atanh.
        '''
        es = self.es  # signOf(es) == signOf(f)
        return es * (atanh(es * x) if es > 0 else  # .isOblate
                    (-atan(es * x) if es < 0 else  # .isProlate
                     _0_0))  # .isSpherical

    @Property_RO
    def es_c(self):
        '''Get M{(1 - f) * exp(es_atanh(1))} (C{float}), M{b_a * exp(es_atanh(1))}.
        '''
        return Float(es_c=(self._exp_es_atanh_1 * self.b_a) if self.f else _1_0)

    def es_tauf(self, taup):
        '''Compute I{Karney}'s U{equations (19), (20) and (21)
           <https://ArXiv.org/abs/1002.1417>}.

           @see: I{Karney}'s C++ method U{Math::tauf<https://GeographicLib.
                 SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>} and
                 and I{Veness}' JavaScript method U{toLatLon<https://www.
                 Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>}.
        '''
        t = Scalar(taup=taup)
        if self.f:  # .isEllipsoidal
            a =  fabs(t)
            T = (self._exp_es_atanh_1 if a > 70 else self._1_e21) * t
            if fabs(T * _EPSqrt) < _2_0:  # handles +/- INF and NAN
                s = (a * _TOL) if a > _1_0 else _TOL
                for T, _, d in self._es_tauf3(t, T):  # max 2
                    if fabs(d) < s:
                        break
            t = Scalar(tauf=T)
        return t

    def _es_tauf3(self, taup, T, N=9):  # in .utm.Utm._toLLEB
        '''(INTERNAL) Yield a 3-tuple C{(τi, iteration, delta)} for at most
           B{C{N}} Newton iterations, converging rapidly except when C{delta}
           toggles on +/-1.12e-16 or +/-4.47e-16, see C{.utm.Utm._toLLEB}.
        '''
        e    = self._1_e21
        _F2_ = Fsum(T).fsum2f_  # τ0
        _tf2 = self._es_taupf2
        for i in range(1, N + 1):
            a, h = _tf2(T)
            # = (taup - a) / hypot1(a) / ((e + T**2) / h)
            d = _over((taup - a) * (T**2 + e), hypot1(a) * h)
            T, d = _F2_(d)  # τi, (τi - τi-1)
            yield T, i, d

    def es_taupf(self, tau):
        '''Compute I{Karney}'s U{equations (7), (8) and (9)
           <https://ArXiv.org/abs/1002.1417>}.

           @see: I{Karney}'s C++ method U{Math::taupf<https://GeographicLib.
                 SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>}.
        '''
        t = Scalar(tau=tau)
        if self.f:  # .isEllipsoidal
            t, _ = self._es_taupf2(t)
            t = Scalar(taupf=t)
        return t

    def _es_taupf2(self, tau):
        '''(INTERNAL) Return 2-tuple C{(es_taupf(tau), hypot1(tau))}.
        '''
        if _isfinite(tau):
            h = hypot1(tau)
            s = sinh(self._es_atanh(tau / h))
            a = hypot1(s) * tau - h * s
        else:
            a, h = tau, INF
        return a, h

    @Property_RO
    def _exp_es_atanh_1(self):
        '''(INTERNAL) Helper for .es_c and .es_tauf.
        '''
        return exp(self._es_atanh(_1_0)) if self.es else _1_0

    @Property_RO
    def f(self):
        '''Get the I{flattening} (C{scalar}), M{(a - b) / a}, C{0} for spherical, negative for prolate.
        '''
        return self._f

    @Property_RO
    def f_(self):
        '''Get the I{inverse flattening} (C{scalar}), M{1 / f} == M{a / (a - b)}, C{0} for spherical, see C{a_b2f_}.
        '''
        return self._f_

    @Property_RO
    def f1(self):
        '''Get the I{1 - flattening} (C{float}), M{f1 == 1 - f == b / a}.

           @see: Property L{b_a}.
        '''
        return Float(f1=_1_0 - self.f)

    @Property_RO
    def f2(self):
        '''Get the I{2nd flattening} (C{float}), M{(a - b) / b == f / (1 - f)}, C{0} for spherical, see C{a_b2f2}.
        '''
        return self._assert(self.a_b - _1_0, f2=f2f2(self.f))

    @deprecated_Property_RO
    def geodesic(self):
        '''DEPRECATED, use property C{geodesicw}.'''
        return self.geodesicw

    def geodesic_(self, exact=True):
        '''Get the an I{exact} C{Geodesic...} instance for this ellipsoid.

           @kwarg exact: If C{bool} return L{GeodesicExact}C{(exact=B{exact}, ...)},
                         otherwise a L{Geodesic}, L{GeodesicExact} or L{GeodesicSolve}
                         instance for I{this} ellipsoid.

           @return: The C{exact} geodesic (C{Geodesic...}).

           @raise TypeError: Invalid B{C{exact}}.

           @raise ValueError: Incompatible B{C{exact}} ellipsoid.
        '''
        if isbool(exact):  # for consistenccy with C{.rhumb_}
            g = _MODS.geodesicx.GeodesicExact(self, C4order=30 if exact else 24,
                                              name=self.name)
        else:
            g =  exact
            E = _xattr(g, ellipsoid=None)
            if not (E is self and isinstance(g, self._Geodesics)):
                raise _ValueError(exact=g, ellipsoid=E, txt_not_=self.name)
        return g

    @property_ROver
    def _Geodesics(self):
        '''(INTERNAL) Get all C{Geodesic...} classes, I{once}.
        '''
        t = (_MODS.geodesicx.GeodesicExact,
             _MODS.geodsolve.GeodesicSolve)
        try:
            t += (_MODS.geodesicw.Geodesic,
                  _MODS.geodesicw._wrapped.Geodesic)
        except ImportError:
            pass
        return t  # overwrite property_ROver

    @property_RO
    def geodesicw(self):
        '''Get this ellipsoid's I{wrapped} U{geodesicw.Geodesic
           <https://GeographicLib.SourceForge.io/Python/doc/code.html>}, provided
           I{Karney}'s U{geographiclib<https://PyPI.org/project/geographiclib>}
           package is installed.
        '''
        # if not self.isEllipsoidal:
        #     raise _IsnotError(_ellipsoidal_, ellipsoid=self)
        return _MODS.geodesicw.Geodesic(self)

    @property_RO
    def geodesicx(self):
        '''Get this ellipsoid's I{exact} L{GeodesicExact}.
        '''
        # if not self.isEllipsoidal:
        #     raise _IsnotError(_ellipsoidal_, ellipsoid=self)
        return _MODS.geodesicx.GeodesicExact(self, name=self.name)

    @property
    def geodsolve(self):
        '''Get this ellipsoid's L{GeodesicSolve}, the I{wrapper} around utility
           U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>},
           provided the path to the C{GeodSolve} executable is specified with env
           variable C{PYGEODESY_GEODSOLVE} or re-/set with this property..
        '''
        # if not self.isEllipsoidal:
        #     raise _IsnotError(_ellipsoidal_, ellipsoid=self)
        return _MODS.geodsolve.GeodesicSolve(self, path=self._geodsolve, name=self.name)

    @geodsolve.setter  # PYCHOK setter!
    def geodsolve(self, path):
        '''Re-/set the (fully qualified) path to the U{GeodSolve
           <https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>} executable,
           overriding env variable C{PYGEODESY_GEODSOLVE} (C{str}).
        '''
        self._geodsolve = path

    def hartzell4(self, pov, los=None):
        '''Compute the intersection of this ellipsoid's surface and a Line-Of-Sight
           from a Point-Of-View in space.

           @arg pov: Point-Of-View outside this ellipsoid (C{Cartesian}, L{Ecef9Tuple}
                     or L{Vector3d}).
           @kwarg los: Line-Of-Sight, I{direction} to this ellipsoid (L{Los}, L{Vector3d})
                       or C{True} for the I{normal, perpendicular, plumb} to the surface
                       of this ellipsoid or C{False} or C{None} to point to its center.

           @return: L{Vector4Tuple}C{(x, y, z, h)} with the cartesian coordinates C{x},
                    C{y} and C{z} of the projection on or the intersection with this
                    ellipsoid and the I{distance} C{h} from B{C{pov}} to C{(x, y, z)}
                    along B{C{los}}, all in C{meter}, conventionally.

           @raise IntersectionError: Null B{C{pov}} or B{C{los}} vector, or B{C{pov}}
                                     is inside this ellipsoid or B{C{los}} points
                                     outside this ellipsoid or in opposite direction.

           @raise TypeError: Invalid B{C{pov}} or B{C{los}}.

           @see: U{I{Satellite Line-of-Sight Intersection with Earth}<https://StephenHartzell.
                 Medium.com/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6>} and
                 methods L{Ellipsoid.height4} and L{Triaxial.hartzell4}.
        '''
        try:
            v, d, i = _MODS.triaxials._hartzell3(pov, los, self._triaxial)
        except Exception as x:
            raise IntersectionError(pov=pov, los=los, cause=x)
        return Vector4Tuple(v.x, v.y, v.z, d, iteration=i, name__=self.hartzell4)

    @Property_RO
    def _hash(self):
        return hash((self.a, self.f))

    def height4(self, xyz, normal=True):
        '''Compute the projection on and the height of a cartesian above or below
           this ellipsoid's surface.

           @arg xyz: The cartesian (C{Cartesian}, L{Ecef9Tuple}, L{Vector3d},
                     L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg normal: If C{True}, the projection is perpendicular to (the nearest
                          point on) this ellipsoid's surface, otherwise the C{radial}
                          line to this ellipsoid's center (C{bool}).

           @return: L{Vector4Tuple}C{(x, y, z, h)} with the cartesian coordinates C{x},
                    C{y} and C{z} of the projection on and the height C{h} above or
                    below this ellipsoid's surface, all in C{meter}, conventionally.

           @raise ValueError: Null B{C{xyz}}.

           @raise TypeError: Non-cartesian B{C{xyz}}.

           @see: U{Distance to<https://StackOverflow.com/questions/22959698/distance-from-given-point-to-given-ellipse>}
                 and U{intersection with<https://MathWorld.wolfram.com/Ellipse-LineIntersection.html>} an ellipse and
                 methods L{Ellipsoid.hartzell4} and L{Triaxial.height4}.
        '''
        v = _MODS.vector3d._otherV3d(xyz=xyz)
        r =  v.length

        a, b, i = self.a, self.b, None
        if r < EPS0:  # EPS
            v =  v.times(_0_0)
            h = -a

        elif self.isSpherical:
            v = v.times(a / r)
            h = r - a

        elif normal:  # perpendicular to ellipsoid
            x, y = hypot(v.x, v.y), fabs(v.z)
            if x < EPS0:  # PYCHOK no cover
                z = copysign0(b, v.z)
                v = Vector3Tuple(v.x, v.y, z)
                h = y - b  # polar
            elif y < EPS0:  # PYCHOK no cover
                t = a / r
                v = v.times_(t, t, 0)  # force z=0.0
                h = x - a  # equatorial
            else:  # normal in 1st quadrant
                x, y, i = _MODS.triaxials._plumbTo3(x, y, self)
                t, v = v, v.times_(x, x, y)
                h = t.minus(v).length

        else:  # radial to ellipsoid's center
            h =  hypot_(a * v.z, b * v.x, b * v.y)
            t = (a * b / h) if h > EPS0 else _0_0  # EPS
            v =  v.times(t)
            h =  r * (_1_0 - t)

        return Vector4Tuple(v.x, v.y, v.z, h, iteration=i, name__=self.height4)

    def _hubeny_2(self, phi2, phi1, lam21, scaled=True, squared=True):
        '''(INTERNAL) like function C{pygeodesy.flatLocal_}/C{pygeodesy.hubeny_},
           returning the I{angular} distance in C{radians squared} or C{radians}
        '''
        m, n = self.roc2_((phi2 + phi1) * _0_5, scaled=scaled)
        h, r = (hypot2, self.a2_) if squared else (hypot, _1_0 / self.a)
        return h(m * (phi2 - phi1), n * lam21) * r

    @Property_RO
    def isEllipsoidal(self):
        '''Is this model I{ellipsoidal} (C{bool})?
        '''
        return self.f != 0

    @Property_RO
    def isOblate(self):
        '''Is this ellipsoid I{oblate} (C{bool})?  I{Prolate} or
           spherical otherwise.
        '''
        return self.f > 0

    @Property_RO
    def isProlate(self):
        '''Is this ellipsoid I{prolate} (C{bool})?  I{Oblate} or
           spherical otherwise.
        '''
        return self.f < 0

    @Property_RO
    def isSpherical(self):
        '''Is this ellipsoid I{spherical} (C{bool})?
        '''
        return self.f == 0

    def _Kseries(self, *AB8Ks):
        '''(INTERNAL) Compute the 4-, 6- or 8-th order I{Krüger} Alpha
           or Beta series coefficients per I{Karney}'s U{equations (35)
           and (36)<https://ArXiv.org/pdf/1002.1417v3.pdf>}.

           @arg AB8Ks: 8-Tuple of 8-th order I{Krüger} Alpha or Beta series
                       coefficient tuples.

           @return: I{Krüger} series coefficients (L{KsOrder}C{-tuple}).

           @see: I{Karney}'s 30-th order U{TMseries30
                 <https://GeographicLib.SourceForge.io/C++/doc/tmseries30.html>}.
        '''
        k = self.KsOrder
        if self.n:
            ns = fpowers(self.n, k)
            ks = tuple(fdot(AB8Ks[i][:k-i], *ns[i:]) for i in range(k))
        else:
            ks = _0_0s(k)
        return ks

    @property_doc_(''' the I{Krüger} series' order (C{int}), see properties C{AlphaKs}, C{BetaKs}.''')
    def KsOrder(self):
        '''Get the I{Krüger} series' order (C{int} 4, 6 or 8).
        '''
        return self._KsOrder

    @KsOrder.setter  # PYCHOK setter!
    def KsOrder(self, order):
        '''Set the I{Krüger} series' order (C{int} 4, 6 or 8).

           @raise ValueError: Invalid B{C{order}}.
        '''
        if not (isint(order) and _isin(order, 4, 6, 8)):
            raise _ValueError(order=order)
        if self._KsOrder != order:
            Ellipsoid.AlphaKs._update(self)
            Ellipsoid.BetaKs._update(self)
            self._KsOrder = order

    @Property_RO
    def L(self):
        '''Get the I{quarter meridian} C{L}, aka the C{polar distance}
           along a meridian between the equator and a pole (C{meter}),
           M{b * Elliptic(-e2 / (1 - e2)).cE} or M{b * PI / 2}.
        '''
        r = self._elliptic_e22.cE if self.f else PI_2
        return Distance(L=self.b * r)

    def Llat(self, lat):
        '''Return the I{meridional length}, the distance along a meridian
           between the equator and a (geodetic) latitude, see C{L}.

           @arg lat: Geodetic latitude (C{degrees90}).

           @return: The meridional length at B{C{lat}}, negative on southern
                    hemisphere (C{meter}).
        '''
        r = self._elliptic_e22.fEd(self.auxParametric(lat)) if self.f else Phid(lat)
        return Distance(Llat=self.b * r)

    Lmeridian = Llat  # meridional distance

    @property_RO
    def _Lpd(self):
        '''Get the I{quarter meridian} per degree (C{meter}), M{self.L / 90}.
        '''
        return Meter(_Lpd=self.L / _90_0)

    @property_RO
    def _Lpr(self):
        '''Get the I{quarter meridian} per radian (C{meter}), M{self.L / PI_2}.
        '''
        return Meter(_Lpr=self.L / PI_2)

    @deprecated_Property_RO
    def majoradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{a} or C{Requatorial}.'''
        return self.a

    def m2degrees(self, distance, lat=0):
        '''Convert a distance to an angle along the equator or along
           a parallel of (geodetic) latitude.

           @arg distance: Distance (C{meter}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Angle (C{degrees}) or C{INF} for near-polar B{C{lat}}.

           @raise RangeError: Latitude B{C{lat}} outside valid range and
                              L{rangerrors<pygeodesy.rangerrors>} is C{True}.

           @raise ValueError: Invalid B{C{distance}} or B{C{lat}}.
       '''
        return degrees(self.m2radians(distance, lat=lat))

    def m2radians(self, distance, lat=0):
        '''Convert a distance to an angle along the equator or along
           a parallel of (geodetic) latitude.

           @arg distance: Distance (C{meter}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Angle (C{radians}) or C{INF} for near-polar B{C{lat}}.

           @raise RangeError: Latitude B{C{lat}} outside valid range and
                              L{rangerrors<pygeodesy.rangerrors>} is C{True}.

           @raise ValueError: Invalid B{C{distance}} or B{C{lat}}.
        '''
        r = self.circle4(lat).radius if lat else self.a
        return m2radians(distance, radius=r, lat=0)

    @deprecated_Property_RO
    def minoradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{b}, C{polaradius} or C{Rpolar}.'''
        return self.b

    @Property_RO
    def n(self):
        '''Get the I{3rd flattening} (C{float}), M{f / (2 - f) == (a - b) / (a + b)}, see C{a_b2n}.
        '''
        return self._assert(a_b2n(self.a, self.b), n=f2n(self.f))

    flattening    = f
    flattening1st = f
    flattening2nd = f2
    flattening3rd = n

    polaradius = b  # Rpolar

#   Q = A  # I{meridian arc unit} C{Q}, the mean, meridional length I{per radian}

    @deprecated_Property_RO
    def quarteradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{L} or method C{Llat}.'''
        return self.L

    @Property_RO
    def R1(self):
        '''Get the I{mean} earth radius per I{IUGG} (C{meter}), M{(2 * a + b) / 3 == a * (1 - f / 3)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}
                 and method C{Rgeometric}.
        '''
        r = Fsum(self.a, self.a, self.b).fover(_3_0) if self.f else self.a
        return Radius(R1=r)

    Rmean = R1

    @Property_RO
    def R2(self):
        '''Get the I{authalic} earth radius (C{meter}), M{sqrt(c2)}.

           @see: C{R2x}, C{c2}, C{area} and U{Earth radius
                 <https://WikiPedia.org/wiki/Earth_radius>}.
        '''
        return Radius(R2=sqrt(self.c2) if self.f else self.a)
#       # Moritz, H. <https://Geodesy.Geology.Ohio-State.EDU/course/refpapers/00740128.pdf>
#       # R2 = (1 - 2/3 * e'2 + 26/45 * e'4 - 100/189 * e'6 + 7034/14175 * e'8) * rocPolar
#       #    = (3 + e'2 * (-2 + e'2 * (26/15 + e'2 * (-100/63 + e'2 * 7034/4725)))) * rocPolar / 3
#       return Fhorner(self.e22, 3, -2, 26 / 15, -100 / 63, 7034 / 4725).fover(3 / self.rocPolar)

    Rauthalic = R2

    @Property_RO
    def R2x(self):
        '''Get the I{authalic} earth radius (C{meter}), M{sqrt(c2x)}.

           @see: C{R2}, C{c2x} and C{areax}.
        '''
        return Radius(R2x=sqrt(self.c2x) if self.f else self.a)

    Rauthalicx = R2x

    @Property_RO
    def R3(self):
        '''Get the I{volumetric} earth radius (C{meter}), M{(a * a * b)**(1/3)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>} and C{volume}.
        '''
        r = (cbrt(self.b_a) * self.a) if self.f else self.a
        return Radius(R3=r)

    Rvolumetric = R3

    def radians2m(self, rad, lat=0):
        '''Convert an angle to the distance along the equator or along
           a parallel of (geodetic) latitude.

           @arg rad: The angle (C{radians}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Distance (C{meter}, same units as the equatorial
                    and polar radii) or C{0} for near-polar B{C{lat}}.

           @raise RangeError: Latitude B{C{lat}} outside valid range and
                              L{rangerrors<pygeodesy.rangerrors>} is C{True}.

           @raise ValueError: Invalid B{C{rad}} or B{C{lat}}.
        '''
        r = self.circle4(lat).radius if lat else self.a
        return radians2m(rad, radius=r, lat=0)

    @Property_RO
    def Rbiaxial(self):
        '''Get the I{biaxial, quadratic} mean earth radius (C{meter}), M{sqrt((a**2 + b**2) / 2)}.

           @see: C{Rtriaxial}
        '''
        a, b = self.a, self.b
        if b < a:
            b  = sqrt(_0_5 + self.b2_a2 * _0_5) * a
        elif b > a:
            b *= sqrt(_0_5 + self.a2_b2 * _0_5)
        return Radius(Rbiaxial=b)

    Requatorial = a  # for consistent naming

    def Rgeocentric(self, lat):
        '''Compute the I{geocentric} earth radius of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Geocentric earth radius (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Geocentric Radius
                 <https://WikiPedia.org/wiki/Earth_radius#Geocentric_radius>}
        '''
        r, a = self.a, Phid(lat)
        if a and self.f:
            if fabs(a) < PI_2:
                s2, c2 = _s2_c2(a)
                b2_a2_s2 = self.b2_a2 * s2
                # R == sqrt((a2**2 * c2 + b2**2 * s2) / (a2 * c2 + b2 * s2))
                #   == sqrt(a2**2 * (c2 + (b2 / a2)**2 * s2) / (a2 * (c2 + b2 / a2 * s2)))
                #   == sqrt(a2 * (c2 + (b2 / a2)**2 * s2) / (c2 + (b2 / a2) * s2))
                #   == a * sqrt((c2 + b2_a2 * b2_a2 * s2) / (c2 + b2_a2 * s2))
                #   == a * sqrt((c2 + b2_a2 * b2_a2_s2) / (c2 + b2_a2_s2))
                r *= sqrt((c2 + b2_a2_s2 * self.b2_a2) / (c2 + b2_a2_s2))
            else:
                r  = self.b
        return Radius(Rgeocentric=r)

    @Property_RO
    def Rgeometric(self):
        '''Get the I{geometric} mean earth radius (C{meter}), M{sqrt(a * b)}.

           @see: C{R1}.
        '''
        g = sqrt(self.a * self.b) if self.f else self.a
        return Radius(Rgeometric=g)

    def rhumb_(self, exact=True):
        '''Get the an I{exact} C{Rhumb...} instance for this ellipsoid.

           @kwarg exact: If C{bool} or C{None} return L{Rhumb}C{(exact=B{exact}, ...)},
                         otherwise a L{Rhumb}, L{RhumbAux} or L{RhumbSolve} instance
                         for I{this} ellipsoid.

           @return: The C{exact} rhumb (C{Rhumb...}).

           @raise TypeError: Invalid B{C{exact}}.

           @raise ValueError: Incompatible B{C{exact}} ellipsoid.
        '''
        if isbool(exact):  # use Rhumb for backward compatibility
            r = _MODS.rhumb.ekx.Rhumb(self, exact=exact, name=self.name)
        else:
            r =  exact
            E = _xattr(r, ellipsoid=None)
            if not (E is self and isinstance(r, self._Rhumbs)):
                raise _ValueError(exact=r, ellipsosid=E, txt_not_=self.name)
        return r

    @property_RO
    def rhumbaux(self):
        '''Get this ellipsoid's I{Auxiliary} C{rhumb.RhumbAux}.
        '''
        # if not self.isEllipsoidal:
        #     raise _IsnotError(_ellipsoidal_, ellipsoid=self)
        return _MODS.rhumb.aux_.RhumbAux(self, name=self.name)

    @property_RO
    def rhumbekx(self):
        '''Get this ellipsoid's I{Elliptic, Krüger} C{rhumb.Rhumb}.
        '''
        # if not self.isEllipsoidal:
        #     raise _IsnotError(_ellipsoidal_, ellipsoid=self)
        return _MODS.rhumb.ekx.Rhumb(self, name=self.name)

    @property_ROver
    def _Rhumbs(self):
        '''(INTERNAL) Get all C{Rhumb...} classes, I{once}.
        '''
        r = _MODS.rhumb
        return (r.aux_.RhumbAux,  # overwrite property_ROver
                r.ekx.Rhumb, r.solve.RhumbSolve)

    @property
    def rhumbsolve(self):
        '''Get this ellipsoid's L{RhumbSolve}, the I{wrapper} around utility
           U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>},
           provided the path to the C{RhumbSolve} executable is specified with env
           variable C{PYGEODESY_RHUMBSOLVE} or re-/set with this property.
        '''
        # if not self.isEllipsoidal:
        #     raise _IsnotError(_ellipsoidal_, ellipsoid=self)
        return _MODS.rhumb.solve.RhumbSolve(self, path=self._rhumbsolve, name=self.name)

    @rhumbsolve.setter  # PYCHOK setter!
    def rhumbsolve(self, path):
        '''Re-/set the (fully qualified) path to the U{RhumbSolve
           <https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>} executable,
           overriding env variable C{PYGEODESY_RHUMBSOLVE} (C{str}).
        '''
        self._rhumbsolve = path

    @deprecated_property_RO
    def rhumbx(self):
        '''DEPRECATED on 2023.11.28, use property C{rhumbekx}. '''
        return self.rhumbekx

    def Rlat(self, lat):
        '''I{Approximate} the earth radius of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Approximate earth radius (C{meter}).

           @raise RangeError: Latitude B{C{lat}} outside valid range and
                              L{rangerrors<pygeodesy.rangerrors>} is C{True}.

           @raise TypeError: Invalid B{C{lat}}.

           @raise ValueError: Invalid B{C{lat}}.

           @note: C{Rlat(B{90})} equals C{Rpolar}.

           @see: Method C{circle4}.
        '''
        # r = a - (a - b) * |lat| / 90
        r = self.a
        if self.f and lat:  # .isEllipsoidal
            r -= (r - self.b) * fabs(Lat(lat)) / _90_0
            r  =  Radius(Rlat=r)
        return r

    Rpolar = b  # for consistent naming

    def roc1_(self, sa, ca=None):
        '''Compute the I{prime-vertical}, I{normal} radius of curvature
           of (geodetic) latitude, I{unscaled}.

           @arg sa: Sine of the latitude (C{float}, [-1.0..+1.0]).
           @kwarg ca: Optional cosine of the latitude (C{float}, [-1.0..+1.0])
                      to use an alternate formula.

           @return: The prime-vertical radius of curvature (C{float}).

           @note: The delta between both formulae with C{Ellipsoids.WGS84}
                  is less than 2 nanometer over the entire latitude range.

           @see: Method L{roc2_} and class L{EcefYou}.
        '''
        if not self.f:  # .isSpherical
            n = self.a
        elif ca is None:
            r = self.e2s2(sa)  # see .roc2_ and _EcefBase._forward
            n = sqrt(self.a2 / r) if r > EPS02 else _0_0
        elif ca:  # derived from EcefYou.forward
            h = hypot(ca, self.b_a * sa) if sa else fabs(ca)
            n = self.a / h
        elif sa:
            n = self.a2_b / fabs(sa)
        else:
            n = self.a
        return n

    def roc2(self, lat, scaled=False):
        '''Compute the I{meridional} and I{prime-vertical}, I{normal}
           radii of curvature of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).
           @kwarg scaled: Scale prime_vertical by C{cos(radians(B{lat}))} (C{bool}).

           @return: An L{Curvature2Tuple}C{(meridional, prime_vertical)} with
                    the radii of curvature.

           @raise ValueError: Invalid B{C{lat}}.

           @see: Methods L{roc2_} and L{roc1_}, U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and meridional and
                 prime vertical U{Radii of Curvature<https://WikiPedia.org/wiki/
                 Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phid(lat), scaled=scaled)

    def roc2_(self, phi, scaled=False):
        '''Compute the I{meridional} and I{prime-vertical}, I{normal} radii of
           curvature of (geodetic) latitude.

           @arg phi: Latitude (C{radians}).
           @kwarg scaled: Scale prime_vertical by C{cos(B{phi})} (C{bool}).

           @return: An L{Curvature2Tuple}C{(meridional, prime_vertical)} with the
                    radii of curvature.

           @raise ValueError: Invalid B{C{phi}}.

           @see: Methods L{roc2} and L{roc1_}, property L{rocEquatorial2}, U{Local,
                 flat earth approximation<https://www.EdWilliams.org/avform.htm#flat>}
                 and the meridional and prime vertical U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        a = fabs(Phi(phi))
        if self.f:
            r = self.e2s2(sin(a))
            if r > EPS02:
                n = self.a / sqrt(r)
                m = n *  self.e21 / r
            else:
                m = n = _0_0
        else:
            m = n = self.a
        if scaled and a:
            n *= cos(a) if a < PI_2 else _0_0
        return Curvature2Tuple(m, n)

    def rocAzimuth(self, lat, azimuth):
        '''Compute the I{directional} radius of curvature of (geodetic) latitude
           and C{azimuth} compass direction.

           @see: Method L{rocBearing<Ellipsoid.rocBearing>} for details, using C{azimuth} for C{bearing}.
        '''
        return Radius(rocAzimuth=self._rocDirectional(lat, Azimuth(azimuth)))

    def rocBearing(self, lat, bearing):
        '''Compute the I{directional} radius of curvature of (geodetic) latitude
           and C{bearing} compass direction.

           @arg lat: Latitude (C{degrees90}).
           @arg bearing: Direction (compass C{degrees360}).

           @return: Directional radius of curvature (C{meter}).

           @raise RangeError: Latitude B{C{lat}} outside valid range and
                              L{rangerrors<pygeodesy.rangerrors>} is C{True}.

           @raise ValueError: Invalid B{C{lat}} or B{C{bearing}}.

           @see: U{Radii of Curvature<https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
        '''
        return Radius(rocBearing=self._rocDirectional(lat, Bearing(bearing)))

    def _rocDirectional(self, lat, deg):
        '''(INTERNAL) Helper for C{rocAzimuth} and C{rocBearing}.
        '''
        if self.f:
            s2, c2 = _s2_c2(radians(deg))
            m, n = self.roc2_(Phid(lat))
            if n < m:  # == n / (c2 * n / m + s2)
                c2 *= n / m
            elif m < n:  # == m / (c2 + s2 * m / n)
                s2 *= m / n
                n   = m
            r = _over(n, c2 + s2)  # == 1 / (c2 / m + s2 / n)
        else:
            r =  self.b  # == self.a
        return r

    @Property_RO
    def rocEquatorial2(self):
        '''Get the I{meridional} and I{prime-vertical}, I{normal} radii of curvature
           at the equator as L{Curvature2Tuple}C{(meridional, prime_vertical)}.

           @see: Methods L{rocMeridional} and L{rocPrimeVertical}, properties L{b2_a},
                 L{a2_b}, C{rocPolar} and polar and equatorial U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        m = self.b2_a if self.f else self.a
        return Curvature2Tuple(m, self.a)

    def rocGauss(self, lat):
        '''Compute the I{Gaussian} radius of curvature of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Gaussian radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: Non-directional U{Radii of Curvature<https://WikiPedia.org/wiki/
                 Earth_radius#Radii_of_curvature>}
        '''
        # using ...
        #    m, n = self.roc2_(Phid(lat))
        #    return sqrt(m * n)
        # ... requires 1 or 2 sqrt
        g = self.b
        if self.f:
            s2, c2 = _s2_c2(Phid(lat))
            g = _over(g, c2 + self.b2_a2 * s2)
        return Radius(rocGauss=g)

    def rocMean(self, lat):
        '''Compute the I{mean} radius of curvature of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Mean radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: Non-directional U{Radii of Curvature<https://WikiPedia.org/wiki/
                 Earth_radius#Radii_of_curvature>}
        '''
        if self.f:
            m, n = self.roc2_(Phid(lat))
            m *= _over(n * _2_0, m + n)  # == 2 / (1 / m + 1 / n)
        else:
            m  =  self.a
        return Radius(rocMean=m)

    def rocMeridional(self, lat):
        '''Compute the I{meridional} radius of curvature of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Meridional radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: Methods L{roc2} and L{roc2_}, U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and U{Radii of
                 Curvature<https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        r = self.roc2_(Phid(lat)) if lat else self.rocEquatorial2
        return Radius(rocMeridional=r.meridional)

    rocPolar = a2_b  # synonymous

    def rocPrimeVertical(self, lat):
        '''Compute the I{prime-vertical}, I{normal} radius of curvature of
           (geodetic) latitude, aka the I{transverse} radius of curvature.

           @arg lat: Latitude (C{degrees90}).

           @return: Prime-vertical radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: Methods L{roc2}, L{roc2_} and L{roc1_}, U{Local, flat earth
                 approximation<https://www.EdWilliams.org/avform.htm#flat>} and
                 U{Radii of Curvature<https://WikiPedia.org/wiki/
                 Earth_radius#Radii_of_curvature>}.
        '''
        r = self.roc2_(Phid(lat)) if lat else self.rocEquatorial2
        return Radius(rocPrimeVertical=r.prime_vertical)

    rocTransverse = rocPrimeVertical  # synonymous

    @deprecated_Property_RO
    def Rquadratic(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Rbiaxial} or C{Rtriaxial}.'''
        return self.Rbiaxial

    @deprecated_Property_RO
    def Rr(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Rrectifying}.'''
        return self.Rrectifying

    @Property_RO
    def Rrectifying(self):
        '''Get the I{rectifying} earth radius (C{meter}), M{((a**(3/2) + b**(3/2)) / 2)**(2/3)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}.
        '''
        r = self.a
        if self.f:
            r *= cbrt2((sqrt3(self.b_a) + _1_0) * _0_5)
        return Radius(Rrectifying=r)

    @deprecated_Property_RO
    def Rs(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Rgeometric}.'''
        return self.Rgeometric

    @Property_RO
    def Rtriaxial(self):
        '''Get the I{triaxial, quadratic} mean earth radius (C{meter}), M{sqrt((3 * a**2 + b**2) / 4)}.

           @see: C{Rbiaxial}
        '''
        q, b = self.a, self.b
        if b < q:
            q *= sqrt((self.b2_a2 + _3_0)        * _0_25)
        elif b > q:
            q  = sqrt((self.a2_b2 * _3_0 + _1_0) * _0_25) * b
        return Radius(Rtriaxial=q)

    def toEllipsoid2(self, **name):
        '''Get a copy of this ellipsoid as an L{Ellipsoid2}.

           @kwarg name: Optional, unique C{B{name}=NN} (C{str}).

           @see: Property C{a_f}.
        '''
        return Ellipsoid2(self, None, **name)

    def toStr(self, prec=8, terse=4, **sep_name):  # PYCHOK expected
        '''Return this ellipsoid as a text string.

           @kwarg prec: Number of decimal digits, unstripped (C{int}).
           @kwarg terse: Limit the number of items (C{int}, 0...18),
                         use C{B{terse}=0} or C{=None} for all.
           @kwarg sep_name: Optional C{B{name}=NN} (C{str}) or C{None}
                      to exclude this ellipsoid's name and separator
                      C{B{sep}=", "} to join the items (C{str}).

           @return: This C{Ellipsoid}'s attributes (C{str}).
        '''
        E =  Ellipsoid
        t = (E.a, E.f, E.f_, E.b, E.f2, E.n, E.e,
                       E.e2, E.e21, E.e22, E.e32,
                       E.A, E.L, E.R1, E.R2, E.R3,
                       E.Rbiaxial, E.Rtriaxial)
        if terse:
            t = t[:terse]
        return self._instr(prec=prec, props=t, **sep_name)

    def toTriaxial(self, **name):
        '''Convert this ellipsoid to a L{Triaxial_}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{Triaxial_} or L{Triaxial} with the C{X} axis
                    pointing east and C{Z} pointing north.

           @see: Method L{Triaxial_.toEllipsoid}.
        '''
        T = self._triaxial
        return T.copy(**name) if name else T

    @property_RO
    def _triaxial(self):
        '''(INTERNAL) Get this ellipsoid's un-/ordered C{Triaxial/_}.
        '''
        a, b, m = self.a, self.b, _MODS.triaxials
        T = m.Triaxial if a > b else m.Triaxial_
        return T(a, a, b, name=self.name)

    @Property_RO
    def volume(self):
        '''Get the ellipsoid's I{volume} (C{meter**3}), M{4 / 3 * PI * R3**3}.

           @see: C{R3}.
        '''
        return Meter3(volume=self.a2 * self.b * PI_3 * _4_0)


class Ellipsoid2(Ellipsoid):
    '''An L{Ellipsoid} specified by I{equatorial} radius and I{flattening}.
    '''
    def __init__(self, a, f=None, **name):
        '''New L{Ellipsoid2}.

           @arg a: Equatorial radius, semi-axis (C{meter}) or a previous
                   L{Ellipsoid} instance.
           @arg f: Flattening: (C{float} < 1.0, negative for I{prolate}),
                   if B{C{a}} is in C{meter}.
           @kwarg name: Optional, unique C{B{name}=NN} (C{str}).

           @raise NameError: Ellipsoid with that B{C{name}} already exists.

           @raise ValueError: Invalid B{C{a}} or B{C{f}}.

           @note: C{abs(B{f}) < EPS} is forced to C{B{f}=0}, I{spherical}.
                  Negative C{B{f}} produces a I{prolate} ellipsoid.
        '''
        if f is None and isinstance(a, Ellipsoid):
            Ellipsoid.__init__(self, a.a,   f =a.f,
                                     b=a.b, f_=a.f_, **name)
        else:
            Ellipsoid.__init__(self, a, f=f, **name)


def _ispherical_a_b(a, b):
    '''(INTERNAL) C{True} for spherical or invalid C{a} or C{b}.
    '''
    return a < EPS0 or b < EPS0 or fabs(a - b) < EPS0


def _ispherical_f(f):
    '''(INTERNAL) C{True} for spherical or invalid C{f}.
    '''
    return f > EPS1 or fabs(f) < EPS


def _ispherical_f_(f_):
    '''(INTERNAL) C{True} for spherical or invalid C{f_}.
    '''
    f_ = fabs(f_)
    return f_ < EPS or f_ > _1_EPS


def a_b2e(a, b):
    '''Return C{e}, the I{1st eccentricity} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The I{unsigned}, (1st) eccentricity (C{float} or C{0}), M{sqrt(1 - (b / a)**2)}.

       @note: The result is always I{non-negative} and C{0} for I{near-spherical} ellipsoids.
    '''
    e2 = _a2b2e2(a, b, b2=False)
    return Float(e=sqrt(fabs(e2)) if e2 else _0_0)  # == sqrt(fabs((a - b) * (a + b))) / a


def a_b2e2(a, b):
    '''Return C{e2}, the I{1st eccentricity squared} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The I{signed}, (1st) eccentricity I{squared} (C{float} or C{0}), M{1 - (b / a)**2}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    return Float(e2=_a2b2e2(a, b, b2=False))


def a_b2e22(a, b):
    '''Return C{e22}, the I{2nd eccentricity squared} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The I{signed}, 2nd eccentricity I{squared} (C{float} or C{0}), M{(a / b)**2 - 1}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    return Float(e22=_a2b2e2(a, b, a2=False))


def a_b2e32(a, b):
    '''Return C{e32}, the I{3rd eccentricity squared} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The I{signed}, 3rd eccentricity I{squared} (C{float} or C{0}),
                M{(a**2 - b**2) / (a**2 + b**2)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    return Float(e32=_a2b2e2(a, b))


def _a2b2e2(a, b, a2=True, b2=True):
    '''(INTERNAL) Helper for C{a_b2e}, C{a_b2e2}, C{a_b2e22} and C{a_b2e32}.
    '''
    if _ispherical_a_b(a, b):
        e2   = _0_0
    else:  # a > 0, b > 0
        a, b = (_1_0, b / a) if a > b else (a / b, _1_0)
        a2b2 =  float(a - b) * (a + b)
        e2   = _over(a2b2, (a**2 if a2 else _0_0) +
                           (b**2 if b2 else _0_0)) if a2b2 else _0_0
    return e2


def a_b2f(a, b):
    '''Return C{f}, the I{flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The flattening (C{scalar} or C{0}), M{(a - b) / a}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    f = 0 if _ispherical_a_b(a, b) else _over(float(a - b), a)
    return _f_0_0 if _ispherical_f(f) else Float(f=f)


def a_b2f_(a, b):
    '''Return C{f_}, the I{inverse flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The inverse flattening (C{scalar} or C{0}), M{a / (a - b)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    f_ = 0 if _ispherical_a_b(a, b) else _over(a, float(a - b))
    return _f__0_0 if _ispherical_f_(f_) else Float(f_=f_)


def a_b2f2(a, b):
    '''Return C{f2}, the I{2nd flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The I{signed}, 2nd flattening (C{scalar} or C{0}), M{(a - b) / b}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    t = 0 if _ispherical_a_b(a, b) else float(a - b)
    return Float(f2=_0_0 if fabs(t) < EPS0 else _over(t, b))


def a_b2n(a, b):
    '''Return C{n}, the I{3rd flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The I{signed}, 3rd flattening (C{scalar} or C{0}), M{(a - b) / (a + b)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    t = 0 if _ispherical_a_b(a, b) else float(a - b)
    return Float(n=_0_0 if fabs(t) < EPS0 else _over(t, a + b))


def a_f2b(a, f):
    '''Return C{b}, the I{polar} radius for a given I{equatorial} radius and I{flattening}.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The polar radius (C{float}), M{a * (1 - f)}.
    '''
    b = a if _ispherical_f(f) else (a * (_1_0 - f))
    return Radius_(b=a if _ispherical_a_b(a, b) else b)


def a_f_2b(a, f_):
    '''Return C{b}, the I{polar} radius for a given I{equatorial} radius and I{inverse flattening}.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg f_: Inverse flattening (C{scalar} >>> 1).

       @return: The polar radius (C{float}), M{a * (f_ - 1) / f_}.
    '''
    b = a if _ispherical_f_(f_) else _over(a * (f_ - _1_0), f_)
    return Radius_(b=a if _ispherical_a_b(a, b) else b)


def b_f2a(b, f):
    '''Return C{a}, the I{equatorial} radius for a given I{polar} radius and I{flattening}.

       @arg b: Polar radius (C{scalar} > 0).
       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The equatorial radius (C{float}), M{b / (1 - f)}.
    '''
    t = _1_0 - f
    a =  b if fabs(t) < EPS0 else _over(b, t)
    return Radius_(a=b if _ispherical_a_b(a, b) else a)


def b_f_2a(b, f_):
    '''Return C{a}, the I{equatorial} radius for a given I{polar} radius and I{inverse flattening}.

       @arg b: Polar radius (C{scalar} > 0).
       @arg f_: Inverse flattening (C{scalar} >>> 1).

       @return: The equatorial radius (C{float}), M{b * f_ / (f_ - 1)}.
    '''
    t = f_ - _1_0
    a = b if _ispherical_f_(f_) or fabs(t)      < EPS0 \
                                or fabs(t - f_) < EPS0 else _over(b * f_, t)
    return Radius_(a=b if _ispherical_a_b(a, b) else a)


def e2f(e):
    '''Return C{f}, the I{flattening} for a given I{1st eccentricity}.

       @arg e: The (1st) eccentricity (0 <= C{float} < 1)

       @return: The flattening (C{scalar} or C{0}).

       @see: Function L{e22f}.
    '''
    return e22f(e**2)


def e22f(e2):
    '''Return C{f}, the I{flattening} for a given I{1st eccentricity squared}.

       @arg e2: The (1st) eccentricity I{squared}, I{signed} (L{NINF} < C{float} < 1)

       @return: The flattening (C{float} or C{0}), M{e2 / (sqrt(1 - e2) + 1)}.
    '''
    return Float(f=_over(e2, sqrt(_1_0 - e2) + _1_0)) if e2 and e2 < _1_0 else _f_0_0


def f2e2(f):
    '''Return C{e2}, the I{1st eccentricity squared} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The I{signed}, (1st) eccentricity I{squared} (C{float} < 1), M{f * (2 - f)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    return Float(e2=_0_0 if _ispherical_f(f) else (f * (_2_0 - f)))


def f2e22(f):
    '''Return C{e22}, the I{2nd eccentricity squared} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The I{signed}, 2nd eccentricity I{squared} (C{float} > -1 or C{INF}),
                M{f * (2 - f) / (1 - f)**2}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for near-spherical ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1Ellipsoid.html>}.
    '''
    # e2 / (1 - e2) == f * (2 - f) / (1 - f)**2
    t = (_1_0 - f)**2
    return Float(e22=INF if t < EPS0 else _over(f2e2(f), t))  # PYCHOK type


def f2e32(f):
    '''Return C{e32}, the I{3rd eccentricity squared} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The I{signed}, 3rd eccentricity I{squared} (C{float}),
                M{f * (2 - f) / (1 + (1 - f)**2)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1Ellipsoid.html>}.
    '''
    # e2 / (2 - e2) == f * (2 - f) / (1 + (1 - f)**2)
    e2 = f2e2(f)
    return Float(e32=_over(e2, _2_0 - e2))


def f_2f(f_):
    '''Return C{f}, the I{flattening} for a given I{inverse flattening}.

       @arg f_: Inverse flattening (C{scalar} >>> 1).

       @return: The flattening (C{scalar} or C{0}), M{1 / f_}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    f = 0 if _ispherical_f_(f_) else _over(_1_0, f_)
    return _f_0_0 if _ispherical_f(f) else Float(f=f)  # PYCHOK type


def f2f_(f):
    '''Return C{f_}, the I{inverse flattening} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The inverse flattening (C{scalar} or C{0}), M{1 / f}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    f_ = 0 if _ispherical_f(f) else _over(_1_0, f)
    return _f__0_0 if _ispherical_f_(f_) else Float(f_=f_)  # PYCHOK type


def f2f2(f):
    '''Return C{f2}, the I{2nd flattening} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The I{signed}, 2nd flattening (C{scalar} or C{INF}), M{f / (1 - f)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    t = _1_0 - f
    return Float(f2=_0_0 if _ispherical_f(f) else
                    (INF if  fabs(t) < EPS   else _over(f, t)))  # PYCHOK type


def f2n(f):
    '''Return C{n}, the I{3rd flattening} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The I{signed}, 3rd flattening (-1 <= C{float} < 1),
                M{f / (2 - f)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    return Float(n=_0_0 if _ispherical_f(f) else _over(f, float(_2_0 - f)))


def n2e2(n):
    '''Return C{e2}, the I{1st eccentricity squared} for a given I{3rd flattening}.

       @arg n: The 3rd flattening (-1 <= C{scalar} < 1).

       @return: The I{signed}, (1st) eccentricity I{squared} (C{float} or NINF),
                M{4 * n / (1 + n)**2}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.

       @see: U{Flattening<https://WikiPedia.org/wiki/Flattening>}.
    '''
    t = (n + _1_0)**2
    return Float(e2=_0_0 if fabs(n) < EPS0 else
                   (NINF if      t  < EPS0 else _over(_4_0 * n, t)))


def n2f(n):
    '''Return C{f}, the I{flattening} for a given I{3rd flattening}.

       @arg n: The 3rd flattening (-1 <= C{scalar} < 1).

       @return: The flattening (C{scalar} or NINF), M{2 * n / (1 + n)}.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    t = n + _1_0
    f = 0 if fabs(n) < EPS0 else (NINF if t < EPS0 else _over(_2_0 * n, t))
    return _f_0_0 if _ispherical_f(f) else Float(f=f)


def n2f_(n):
    '''Return C{f_}, the I{inverse flattening} for a given I{3rd flattening}.

       @arg n: The 3rd flattening (-1 <= C{scalar} < 1).

       @return: The inverse flattening (C{scalar} or C{0}), M{1 / f}.

       @see: L{n2f} and L{f2f_}.
    '''
    return f2f_(n2f(n))


class Ellipsoids(_NamedEnum):
    '''(INTERNAL) L{Ellipsoid} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    def _Lazy(self, a, b, f_, **kwds):
        '''(INTERNAL) Instantiate the L{Ellipsoid}.
        '''
        return Ellipsoid(a, b=b, f_=f_, **kwds)

Ellipsoids = Ellipsoids(Ellipsoid)  # PYCHOK singleton
'''Some pre-defined L{Ellipsoid}s, all I{lazily} instantiated.'''
# <https://www.GNU.org/software/gama/manual/html_node/Supported-ellipsoids.html>
# <https://GSSC.ESA.int/navipedia/index.php/Reference_Frames_in_GNSS>
# <https://kb.OSU.edu/dspace/handle/1811/77986>
# <https://www.IBM.com/docs/en/db2/11.5?topic=systems-supported-spheroids>
# <https://w3.Energistics.org/archive/Epicentre/Epicentre_v3.0/DataModel/LogicalDictionary/StandardValues/ellipsoid.html>
# <https://GitHub.com/locationtech/proj4j/blob/master/src/main/java/org/locationtech/proj4j/datum/Ellipsoid.java>
Ellipsoids._assert(  # <https://WikiPedia.org/wiki/Earth_ellipsoid>
    Airy1830       = _lazy(_Airy1830_,       *_T(6377563.396, _0_0,               299.3249646)),  # b=6356256.909
    AiryModified   = _lazy(_AiryModified_,   *_T(6377340.189, _0_0,               299.3249646)),  # b=6356034.448
#   APL4_9         = _lazy('APL4_9',         *_T(6378137.0,   _0_0,               298.24985392)),  # Appl. Phys. Lab. 1965
#   ANS            = _lazy('ANS',            *_T(6378160.0,   _0_0,               298.25)),  # Australian Nat. Spheroid
#   AN_SA96        = _lazy('AN_SA96',        *_T(6378160.0,   _0_0,               298.24985392)),  # Australian Nat. South America
    Australia1966  = _lazy('Australia1966',  *_T(6378160.0,   _0_0,               298.25)),  # b=6356774.7192
    ATS1977        = _lazy('ATS1977',        *_T(6378135.0,   _0_0,               298.257)),  # "Average Terrestrial System"
    Bessel1841     = _lazy(_Bessel1841_,     *_T(6377397.155,  6356078.962818,    299.152812797)),
    BesselModified = _lazy('BesselModified', *_T(6377492.018, _0_0,               299.1528128)),
#   BesselNamibia  = _lazy('BesselNamibia',  *_T(6377483.865, _0_0,               299.1528128)),
    CGCS2000       = _lazy('CGCS2000',       *_T(R_MA,        _0_0,               298.257222101)),  # BeiDou Coord System (BDC)
#   Clarke1858     = _lazy('Clarke1858',     *_T(6378293.639, _0_0,               294.260676369)),
    Clarke1866     = _lazy(_Clarke1866_,     *_T(6378206.4,    6356583.8,         294.978698214)),
    Clarke1880     = _lazy('Clarke1880',     *_T(6378249.145,  6356514.86954978,  293.465)),
    Clarke1880IGN  = _lazy(_Clarke1880IGN_,  *_T(6378249.2,    6356515.0,         293.466021294)),
    Clarke1880Mod  = _lazy('Clarke1880Mod',  *_T(6378249.145,  6356514.96639549,  293.466307656)),  # aka Clarke1880Arc
    CPM1799        = _lazy('CPM1799',        *_T(6375738.7,    6356671.92557493,  334.39)),  # Comm. des Poids et Mesures
    Delambre1810   = _lazy('Delambre1810',   *_T(6376428.0,    6355957.92616372,  311.5)),  # Belgium
    Engelis1985    = _lazy('Engelis1985',    *_T(6378136.05,   6356751.32272154,  298.2566)),
#   Everest1830    = _lazy('Everest1830',    *_T(6377276.345, _0_0,               300.801699997)),
#   Everest1948    = _lazy('Everest1948',    *_T(6377304.063, _0_0,               300.801699997)),
#   Everest1956    = _lazy('Everest1956',    *_T(6377301.243, _0_0,               300.801699997)),
    Everest1969    = _lazy('Everest1969',    *_T(6377295.664,  6356094.667915,    300.801699997)),
    Everest1975    = _lazy('Everest1975',    *_T(6377299.151,  6356098.14512013,  300.8017255)),
    Fisher1968     = _lazy('Fisher1968',     *_T(6378150.0,    6356768.33724438,  298.3)),
#   Fisher1968Mod  = _lazy('Fisher1968Mod',  *_T(6378155.0,   _0_0,               298.3)),
    GEM10C         = _lazy('GEM10C',         *_T(R_MA,         6356752.31424783,  298.2572236)),
    GPES           = _lazy('GPES',           *_T(6378135.0,    6356750.0,        _0_0)),  # "Gen. Purpose Earth Spheroid"
    GRS67          = _lazy('GRS67',          *_T(6378160.0,   _0_0,               298.247167427)),  # Lucerne b=6356774.516
#   GRS67Truncated = _lazy('GRS67Truncated', *_T(6378160.0,   _0_0,               298.25)),
    GRS80          = _lazy(_GRS80_,          *_T(R_MA,         6356752.314140347, 298.25722210088)),  # IUGG, ITRS, ETRS89
#   Hayford1924    = _lazy('Hayford1924',    *_T(6378388.0,    6356911.94612795,  None)),  # aka Intl1924 f_=297
    Helmert1906    = _lazy('Helmert1906',    *_T(6378200.0,    6356818.16962789,  298.3)),
#   Hough1960      = _lazy('Hough1960',      *_T(6378270.0,   _0_0,               297.0)),
    IAU76          = _lazy('IAU76',          *_T(6378140.0,   _0_0,               298.257)),  # Int'l Astronomical Union
    IERS1989       = _lazy('IERS1989',       *_T(6378136.0,   _0_0,               298.257)),  # b=6356751.302
    IERS1992TOPEX  = _lazy('IERS1992TOPEX',  *_T(6378136.3,    6356751.61659215,  298.257223563)),  # IERS/TOPEX/Poseidon/McCarthy
    IERS2003       = _lazy('IERS2003',       *_T(6378136.6,    6356751.85797165,  298.25642)),
    Intl1924       = _lazy(_Intl1924_,       *_T(6378388.0,   _0_0,               297.0)),  # aka Hayford b=6356911.9462795
    Intl1967       = _lazy('Intl1967',       *_T(6378157.5,    6356772.2,         298.24961539)),  # New Int'l
    Krassovski1940 = _lazy(_Krassovski1940_, *_T(6378245.0,    6356863.01877305,  298.3)),  # spelling
    Krassowsky1940 = _lazy(_Krassowsky1940_, *_T(6378245.0,    6356863.01877305,  298.3)),  # spelling
#   Kaula          = _lazy('Kaula',          *_T(6378163.0,   _0_0,               298.24)),  # Kaula 1961
#   Lerch          = _lazy('Lerch',          *_T(6378139.0,   _0_0,               298.257)),  # Lerch 1979
    Maupertuis1738 = _lazy('Maupertuis1738', *_T(6397300.0,    6363806.28272251,  191.0)),  # France
    Mercury1960    = _lazy('Mercury1960',    *_T(6378166.0,    6356784.28360711,  298.3)),
    Mercury1968Mod = _lazy('Mercury1968Mod', *_T(6378150.0,    6356768.33724438,  298.3)),
#   MERIT          = _lazy('MERIT',          *_T(6378137.0,   _0_0,               298.257)),  # MERIT 1983
#   NWL10D         = _lazy('NWL10D',         *_T(6378135.0,   _0_0,               298.26)),  # Naval Weapons Lab.
    NWL1965        = _lazy('NWL1965',        *_T(6378145.0,    6356759.76948868,  298.25)),  # Naval Weapons Lab.
#   NWL9D          = _lazy('NWL9D',          *_T(6378145.0,    6356759.76948868,  298.25)),  # NWL1965
    OSU86F         = _lazy('OSU86F',         *_T(6378136.2,    6356751.51693008,  298.2572236)),
    OSU91A         = _lazy('OSU91A',         *_T(6378136.3,    6356751.6165948,   298.2572236)),
#   Plessis1817    = _lazy('Plessis1817',    *_T(6397523.0,    6355863.0,         153.56512242)),  # XXX incorrect?
    Plessis1817    = _lazy('Plessis1817',    *_T(6376523.0,    6355862.93325557,  308.64)),  # XXX IGN France 1972
#   Prolate        = _lazy('Prolate',        *_T(6356752.3,    R_MA,             _0_0)),
    PZ90           = _lazy('PZ90',           *_T(6378136.0,   _0_0,               298.257839303)),  # GLOSNASS PZ-90 and PZ-90.11
#   SEAsia         = _lazy('SEAsia',         *_T(6378155.0,   _0_0,               298.3)),  # SouthEast Asia
    SGS85          = _lazy('SGS85',          *_T(6378136.0,    6356751.30156878,  298.257)),  # Soviet Geodetic System
    SoAmerican1969 = _lazy('SoAmerican1969', *_T(6378160.0,    6356774.71919531,  298.25)),  # South American
    Sphere         = _lazy(_Sphere_,         *_T(R_M,          R_M,              _0_0)),  # pseudo
    SphereAuthalic = _lazy('SphereAuthalic', *_T(R_FM,         R_FM,             _0_0)),  # pseudo
    SpherePopular  = _lazy('SpherePopular',  *_T(R_MA,         R_MA,             _0_0)),  # EPSG:3857 Spheroid
    Struve1860     = _lazy('Struve1860',     *_T(6378298.3,    6356657.14266956,  294.73)),
#   Walbeck        = _lazy('Walbeck',        *_T(6376896.0,   _0_0,               302.78)),
#   WarOffice      = _lazy('WarOffice',      *_T(6378300.0,   _0_0,               296.0)),
    WGS60          = _lazy('WGS60',          *_T(6378165.0,    6356783.28695944,  298.3)),
    WGS66          = _lazy('WGS66',          *_T(6378145.0,    6356759.76948868,  298.25)),
    WGS72          = _lazy(_WGS72_,          *_T(6378135.0,   _0_0,               298.26)),  # b=6356750.52
    WGS84          = _lazy(_WGS84_,          *_T(R_MA,        _0_0,           _f__WGS84)),  # GPS b=6356752.3142451793
#   U{NOAA/NOS/NGS/inverse<https://GitHub.com/noaa-ngs/inverse/blob/main/invers3d.f>}
    WGS84_NGS      = _lazy('WGS84_NGS',      *_T(R_MA,        _0_0,               298.257222100882711243162836600094))
)

_EWGS84 = Ellipsoids.WGS84  # (INTERNAL) shared

if __name__ == _DMAIN_:

    from pygeodesy.interns import _COMMA_, _NL_, _NLATvar_
    from pygeodesy import nameof, printf

    for E in (_EWGS84, Ellipsoids.GRS80,  # NAD83,
               Ellipsoids.Sphere, Ellipsoids.SpherePopular,
               Ellipsoid(_EWGS84.b, _EWGS84.a, name='_Prolate')):
        e = f2n(E.f) - E.n
        printf('# %s: %s', _DOT_('Ellipsoids', E.name), E.toStr(prec=10, terse=0), nl=1)
        printf('# e=%s, f_=%s, f=%s, n=%s (%s)', fstr(E.e,  prec=13, fmt=Fmt.e),
                                                 fstr(E.f_, prec=13, fmt=Fmt.e),
                                                 fstr(E.f,  prec=13, fmt=Fmt.e),
                                                 fstr(E.n,  prec=13, fmt=Fmt.e),
                                                 fstr(e,    prec=9,  fmt=Fmt.e))
        printf('# %s %s', Ellipsoid.AlphaKs.name, fstr(E.AlphaKs, prec=20))
        printf('# %s  %s', Ellipsoid.BetaKs.name, fstr(E.BetaKs,  prec=20))
        printf('# %s %s', nameof(Ellipsoid.KsOrder), E.KsOrder)  # property

    # __doc__ of this file, force all into registry
    t = [NN] + Ellipsoids.toRepr(all=True, asorted=True).split(_NL_)
    printf(_NLATvar_.join(i.strip(_COMMA_) for i in t))

# % python3.13 -m pygeodesy.ellipsoids

# Ellipsoids.WGS84: name='WGS84', a=6378137, f=0.0033528107, f_=298.257223563, b=6356752.3142451793, f2=0.0033640898, n=0.0016792204, e=0.0818191908, e2=0.00669438, e21=0.99330562, e22=0.0067394967, e32=0.0033584313, A=6367449.1458234144, L=10001965.7293127235, R1=6371008.7714150595, R2=6371007.1809184738, R3=6371000.7900091587, Rbiaxial=6367453.6345163295, Rtriaxial=6372797.5559594007
# e=8.1819190842622e-02, f_=2.98257223563e+02, f=3.3528106647475e-03, n=1.6792203863837e-03 (0.0e+00)
# AlphaKs 0.00083773182062446994, 0.00000076085277735725, 0.00000000119764550324, 0.00000000000242917068, 0.00000000000000571182, 0.0000000000000000148, 0.00000000000000000004, 0.0
# BetaKs  0.00083773216405794875, 0.0000000590587015222, 0.00000000016734826653, 0.00000000000021647981, 0.00000000000000037879, 0.00000000000000000072, 0.0, 0.0
# KsOrder 8

# Ellipsoids.GRS80: name='GRS80', a=6378137, f=0.0033528107, f_=298.2572221009, b=6356752.3141403468, f2=0.0033640898, n=0.0016792204, e=0.081819191, e2=0.00669438, e21=0.99330562, e22=0.0067394968, e32=0.0033584313, A=6367449.1457710434, L=10001965.7292304561, R1=6371008.7713801153, R2=6371007.1808835147, R3=6371000.7899741363, Rbiaxial=6367453.6344640013, Rtriaxial=6372797.5559332585
# e=8.1819191042833e-02, f_=2.9825722210088e+02, f=3.3528106811837e-03, n=1.6792203946295e-03 (0.0e+00)
# AlphaKs 0.00083773182472890429, 0.00000076085278481561, 0.00000000119764552086, 0.00000000000242917073, 0.00000000000000571182, 0.0000000000000000148, 0.00000000000000000004, 0.0
# BetaKs  0.0008377321681623882, 0.00000005905870210374, 0.000000000167348269, 0.00000000000021647982, 0.00000000000000037879, 0.00000000000000000072, 0.0, 0.0
# KsOrder 8

# Ellipsoids.Sphere: name='Sphere', a=6371008.7714149999, f=0, f_=0, b=6371008.7714149999, f2=0, n=0, e=0, e2=0, e21=1, e22=0, e32=0, A=6371008.7714149999, L=10007557.1761167478, R1=6371008.7714149999, R2=6371008.7714149999, R3=6371008.7714149999, Rbiaxial=6371008.7714149999, Rtriaxial=6371008.7714149999
# e=0.0e+00, f_=0.0e+00, f=0.0e+00, n=0.0e+00 (0.0e+00)
# AlphaKs 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
# BetaKs  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
# KsOrder 8

# Ellipsoids.SpherePopular: name='SpherePopular', a=6378137, f=0, f_=0, b=6378137, f2=0, n=0, e=0, e2=0, e21=1, e22=0, e32=0, A=6378137, L=10018754.171394622, R1=6378137, R2=6378137, R3=6378137, Rbiaxial=6378137, Rtriaxial=6378137
# e=0.0e+00, f_=0.0e+00, f=0.0e+00, n=0.0e+00 (0.0e+00)
# AlphaKs 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
# BetaKs  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
# KsOrder 8

# Ellipsoids._Prolate: name='_Prolate', a=6356752.3142451793, f=-0.0033640898, f_=-297.257223563, b=6378137, f2=-0.0033528107, n=-0.0016792204, e=0.0820944379, e2=-0.0067394967, e21=1.0067394967, e22=-0.00669438, e32=-0.0033584313, A=6367449.1458234144, L=10035500.5204500332, R1=6363880.5428301189, R2=6363878.9413582645, R3=6363872.5644020075, Rbiaxial=6367453.6345163295, Rtriaxial=6362105.2243882557
# e=8.2094437949696e-02, f_=-2.97257223563e+02, f=-3.3640898209765e-03, n=-1.6792203863837e-03 (0.0e+00)
# AlphaKs -0.00084149152514366627, 0.00000076653480614871, -0.00000000120934503389, 0.0000000000024576225, -0.00000000000000578863, 0.00000000000000001502, -0.00000000000000000004, 0.0
# BetaKs  -0.00084149187224351817, 0.00000005842735196773, -0.0000000001680487236, 0.00000000000021706261, -0.00000000000000038002, 0.00000000000000000073, -0.0, 0.0
# KsOrder 8

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
