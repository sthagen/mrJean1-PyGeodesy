
# -*- coding: utf-8 -*-

u'''I{Veness}' Universal Transverse Mercator (UTM) projection.

Classes L{Utm} and L{UTMError} and functions L{parseUTM5}, L{toUtm8} and
L{utmZoneBand5}.

Pure Python implementation of UTM / WGS-84 conversion functions using
an ellipsoidal earth model, transcoded from JavaScript originals by
I{(C) Chris Veness 2011-2024} published under the same MIT Licence**, see
U{UTM<https://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>} and
U{Module utm<https://www.Movable-Type.co.UK/scripts/geodesy/docs/module-utm.html>}.

The U{UTM<https://WikiPedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>}
system is a 2-dimensional Cartesian coordinate system providing another way
to identify locations on the surface of the earth.  UTM is a set of 60
transverse Mercator projections, normally based on the WGS-84 ellipsoid.
Within each zone, coordinates are represented as B{C{easting}}s and B{C{northing}}s,
measured in metres.

This module includes some of I{Charles Karney}'s U{'Transverse Mercator with an
accuracy of a few nanometers'<https://ArXiv.org/pdf/1002.1417v3.pdf>}, 2011
(building on Krüger's U{'Konforme Abbildung des Erdellipsoids in der Ebene'
<https://bib.GFZ-Potsdam.DE/pub/digi/krueger2.pdf>}, 1912) and C++ class
U{TransverseMercator<https://GeographicLib.SourceForge.io/C++/doc/
classGeographicLib_1_1TransverseMercator.html>}.

Some other references are U{Universal Transverse Mercator coordinate system
<https://WikiPedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>},
U{Transverse Mercator Projection<https://GeographicLib.SourceForge.io/tm.html>}
and Henrik Seidel U{'Die Mathematik der Gauß-Krueger-Abbildung'
<https://DE.WikiPedia.org/wiki/Gauß-Krüger-Koordinatensystem>}, 2006.
'''

from pygeodesy.basics import len2, map2, neg  # splice
from pygeodesy.constants import EPS, EPS0, _K0_UTM, _0_0, _0_0001
from pygeodesy.datums import _ellipsoidal_datum, _WGS84,  _under
from pygeodesy.dms import degDMS, parseDMS2
from pygeodesy.errors import MGRSError, RangeError, _ValueError, \
                            _xkwds_get, _xkwds_pop2
from pygeodesy.fmath import fdot_, fdot3, hypot, hypot1,  _operator
# from pygeodesy.internals import _under  # from .datums
from pygeodesy.interns import MISSING, NN, _by_, _COMMASPACE_, _N_, \
                             _NS_, _outside_, _range_, _S_, _scale0_, \
                             _SPACE_, _UTM_, _V_, _X_, _zone_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.namedTuples import EasNor2Tuple, UtmUps5Tuple, \
                                  UtmUps8Tuple, UtmUpsLatLon5Tuple
from pygeodesy.props import deprecated_method, property_doc_, \
                            Property_RO
from pygeodesy.streprs import Fmt, unstr
from pygeodesy.units import Band, Int, Lat, Lon, Meter, Zone
from pygeodesy.utily import atan1, atan2, degrees90, degrees180, sincos2
from pygeodesy.utmupsBase import _hemi, _LLEB, _parseUTMUPS5, _to4lldn, \
                                 _to3zBhp, _to3zll, _UPS_LATS, _UPS_ZONE, \
                                 _UTM_LAT_MAX, _UTM_ZONE_MAX, \
                                 _UTM_LAT_MIN, _UTM_ZONE_MIN, \
                                 _UTM_ZONE_OFF_MAX, UtmUpsBase

from math import asinh, atanh, cos, cosh, degrees, fabs, radians, \
                 sin, sinh, tan, tanh  # tan as _tan
# import operator as _operator  # from .fmath

__all__ = _ALL_LAZY.utm
__version__ = '24.11.26'

_Bands = 'CDEFGHJKLMNPQRSTUVWXX'  # UTM latitude bands C..X (no
# I|O) 8° each, covering 80°S to 84°N and X repeated for 80-84°N
_bandLat_      = 'bandLat'
_FalseEasting  =  Meter(  500e3)  # falsed offset (C{meter})
_FalseNorthing =  Meter(10000e3)  # falsed offset (C{meter})
_SvalbardXzone = {32: 9, 34: 21, 36: 33}  # [zone] longitude


class UTMError(_ValueError):
    '''Universal Transverse Mercator (UTM parse or other L{Utm} issue.
    '''
    pass


class _Kseries(object):
    '''(INTERNAL) Alpha or Beta Krüger series.

       Krüger series summations for B{C{eta}}, B{C{ksi}}, B{C{p}} and B{C{q}},
       caching the C{cos}, C{cosh}, C{sin} and C{sinh} values for
       the given B{C{eta}} and B{C{ksi}} angles (in C{radians}).
    '''
    def __init__(self, AB, x, y):
        '''(INTERNAL) New Alpha or Beta Krüger series

           @arg AB: Krüger Alpha or Beta series coefficients
                      (C{4-, 6- or 8-tuple}).
           @arg x: Eta angle (C{radians}).
           @arg y: Ksi angle (C{radians}).
        '''
        n,   j2 = len2(range(2, len(AB) * 2 + 1, 2))
        _m2, _x = map2, _operator.mul

        self._ab =  AB
        self._pq = _m2(_x, j2, AB)
#       assert len(self._ab) == len(self._pq) == n

        x2 = _m2(_x, j2, (x,) * n)
        self._chx = _m2(cosh, x2)
        self._shx = _m2(sinh, x2)
#       assert len(x2) == len(self._chx) == len(self._shx) == n

        y2 = _m2(_x, j2, (y,) * n)
        self._cy = _m2(cos, y2)
        self._sy = _m2(sin, y2)
        # self._sy, self._cy = splice(sincos2(*y2))  # PYCHOK false
#       assert len(y2) == len(self._cy) == len(self._sy) == n

    def xs(self, x0):
        '''(INTERNAL) Eta summation (C{float}).
        '''
        return fdot3(self._ab, self._cy, self._shx, start=x0)

    def ys(self, y0):
        '''(INTERNAL) Ksi summation (C{float}).
        '''
        return fdot3(self._ab, self._sy, self._chx, start=y0)

    def ps(self, p0):
        '''(INTERNAL) P summation (C{float}).
        '''
        return fdot3(self._pq, self._cy, self._chx, start=p0)

    def qs(self, q0):
        '''(INTERNAL) Q summation (C{float}).
        '''
        return fdot3(self._pq, self._sy, self._shx, start=q0)


class Utm(UtmUpsBase):
    '''Universal Transverse Mercator (UTM) coordinate.
    '''
#   _band   =  NN        # latitudinal band letter ('C'|..|'X', no 'I'|'O')
    _Bands  = _Bands     # latitudinal Band letters (C{tuple})
    _Error  =  UTMError  # or etm.ETMError
#   _scale  =  None      # grid scale factor (C{scalar}) or C{None}
    _scale0 = _K0_UTM    # central scale factor (C{scalar})
    _zone   =  0         # longitudinal zone (C{int} 1..60)

    def __init__(self, zone=31, hemisphere=_N_, easting=166022,  # PYCHOK expected
                                northing=0, band=NN, datum=_WGS84, falsed=True,
                                gamma=None, scale=None, **name_convergence):
        '''New L{Utm} UTM coordinate.

           @kwarg zone: Longitudinal UTM zone (C{int}, 1..60) or zone with/-out
                        I{latitudinal} Band letter (C{str}, '1C'|..|'60X').
           @kwarg hemisphere: Northern or southern hemisphere (C{str}, C{'N[orth]'}
                              or C{'S[outh]'}).
           @kwarg easting: Easting, see B{C{falsed}} (C{meter}).
           @kwarg northing: Northing, see B{C{falsed}} (C{meter}).
           @kwarg band: Optional, I{latitudinal} band (C{str}, 'C'|..|'X', no 'I'|'O').
           @kwarg datum: Optional, this coordinate's datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg falsed: If C{True}, both B{C{easting}} and B{C{northing}} are
                          falsed (C{bool}).
           @kwarg gamma: Optional meridian convergence, bearing off grid North,
                         clockwise from true North to save (C{degrees}) or C{None}.
           @kwarg scale: Optional grid scale factor to save (C{scalar}) or C{None}.
           @kwarg name_convergence: Optional C{B{name}=NN} (C{str}) and DEPRECATED
                       keyword argument C{B{convergence}=None}, use B{C{gamma}}.

           @raise TypeError: Invalid or near-spherical B{C{datum}}.

           @raise UTMError: Invalid B{C{zone}}, B{C{hemishere}}, B{C{easting}},
                            B{C{northing}}, B{C{band}}, B{C{convergence}} or
                            B{C{scale}}.
        '''
        if name_convergence:
            gamma, name = _xkwds_pop2(name_convergence, convergence=gamma)
            if name:
                self.name = name

        self._zone, B, _ = _to3zBlat(zone, band)

        h = str(hemisphere)[:1].upper()
        if h not in _NS_:
            raise self._Error(hemisphere=hemisphere)
        self._hemisphere = h

        e, n = easting, northing  # Easting(easting), ...
#       if not falsed:
#           e, n = _false2(e, n, h)
#       # check easting/northing (with 40km overlap
#       # between zones) - is this worthwhile?
#       @raise RangeError: If B{C{easting}} or B{C{northing}} outside
#                          the valid UTM range.
#       if 120e3 > e or e > 880e3:
#           raise RangeError(easting=easting)
#       if 0 > n or n > _FalseNorthing:
#           raise RangeError(northing=northing)

        UtmUpsBase.__init__(self, e, n, band=B, datum=datum, falsed=falsed,
                                                gamma=gamma, scale=scale)

    def __eq__(self, other):
        return isinstance(other, Utm) and other.zone       == self.zone \
                                      and other.hemisphere == self.hemisphere \
                                      and other.easting    == self.easting \
                                      and other.northing   == self.northing \
                                      and other.band       == self.band \
                                      and other.datum      == self.datum

    def __repr__(self):
        return self.toRepr(B=True)

    def __str__(self):
        return self.toStr()

    def _xcopy2(self, Xtm, **name):
        '''(INTERNAL) Make copy as an B{C{Xtm}} instance.

           @arg Xtm: Class to return the copy (C{Xtm=Etm}, C{Xtm=Utm} or
                     C{self.classof}).
        '''
        return Xtm(self.zone, self.hemisphere, self.easting, self.northing,
                   band=self.band, datum=self.datum, falsed=self.falsed,
                   gamma=self.gamma, scale=self.scale, name=self._name__(name))

    @property_doc_(''' the I{latitudinal} band.''')
    def band(self):
        '''Get the I{latitudinal} band (C{'C'|..|'X'}).
        '''
        if not self._band:
            self._toLLEB()
        return self._band

    @band.setter  # PYCHOK setter!
    def band(self, band):
        '''Set or reset the I{latitudinal} band letter (C{'C'|..|'X'})
           or C{None} or C{""} to reset.

           @raise TypeError: Invalid B{C{band}}.

           @raise ValueError: Invalid B{C{band}}.
        '''
        self._band1(band)

    @Property_RO
    def _etm(self):
        '''(INTERNAL) Cache for method L{toEtm}.
        '''
        return self._xcopy2(_MODS.etm.Etm)

    @Property_RO
    def falsed2(self):
        '''Get the easting and northing falsing (L{EasNor2Tuple}C{(easting, northing)}).
        '''
        e = n = 0
        if self.falsed:
            e = _FalseEasting  # relative to central meridian
            if self.hemisphere == _S_:  # relative to equator
                n = _FalseNorthing
        return EasNor2Tuple(e, n)

    def parse(self, strUTM, **name):
        '''Parse a string to a similar L{Utm} instance.

           @arg strUTM: The UTM coordinate (C{str}), see function L{parseUTM5}.
           @kwarg name: Optional instance name (C{str}), overriding this name.

           @return: The similar instance (L{Utm}).

           @raise UTMError: Invalid B{C{strUTM}}.

           @see: Functions L{pygeodesy.parseUPS5} and L{pygeodesy.parseUTMUPS5}.
        '''
        return parseUTM5(strUTM, datum=self.datum, Utm=self.classof,
                                 name=self._name__(name))

    @deprecated_method
    def parseUTM(self, strUTM):  # PYCHOK no cover
        '''DEPRECATED, use method L{Utm.parse}.'''
        return self.parse(strUTM)

    @Property_RO
    def pole(self):
        '''Get the top center of (stereographic) projection, C{""} always.
        '''
        return NN  # N/A for UTM

    def toEtm(self):
        '''Copy this UTM to an ETM coordinate.

           @return: The ETM coordinate (L{Etm}).
        '''
        return self._etm

    def toLatLon(self, LatLon=None, eps=EPS, unfalse=True, **LatLon_kwds):
        '''Convert this UTM coordinate to an (ellipsoidal) geodetic point.

           @kwarg LatLon: Optional, ellipsoidal class to return the geodetic
                          point (C{LatLon}) or C{None}.
           @kwarg eps: Optional convergence limit, L{EPS} or above (C{float}).
           @kwarg unfalse: Unfalse B{C{easting}} and B{C{northing}}
                           if falsed (C{bool}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: This UTM as (B{C{LatLon}}) or if B{C{LatLon}} is
                    C{None}, as L{LatLonDatum5Tuple}C{(lat, lon, datum,
                    gamma, scale)}.

           @raise TypeError: Invalid B{C{datum}} or B{C{LatLon}} is not ellipsoidal.

           @raise UTMError: Invalid meridional radius or H-value.

        '''
        if eps < EPS:
            eps = EPS  # less doesn't converge

        if self._latlon and self._latlon._toLLEB_args == (unfalse, eps):
            return self._latlon5(LatLon)
        else:
            self._toLLEB(unfalse=unfalse, eps=eps)
            return self._latlon5(LatLon, **LatLon_kwds)

    def _toLLEB(self, unfalse=True, eps=EPS):  # PYCHOK signature
        '''(INTERNAL) Compute (ellipsoidal) lat- and longitude.
        '''
        x, y = self.eastingnorthing2(falsed=not unfalse)

        E  = self.datum.ellipsoid
        # from Karney 2011 Eq 15-22, 36
        A0 = self.scale0 * E.A
        if A0 < EPS0:
            raise self._Error(meridional=A0)
        x =  x / A0  # /= chokes PyChecker
        y =  y / A0
        K = _Kseries(E.BetaKs, x, y)  # Krüger series
        x =  neg(K.xs(-x))  # η' eta
        y =  neg(K.ys(-y))  # ξ' ksi

        sy, cy = sincos2(y)
        shx = sinh(x)
        H   = hypot(shx, cy)
        if H < EPS0:
            raise self._Error(H=H)

        T =  sy / H  # τʹ == τ0
        p = _0_0  # previous d
        e = _0_0001 * eps
        for T, i, d in E._es_tauf3(T, T):  # 4-5 trips
            # d may toggle on +/-1.12e-16 or +/-4.47e-16,
            # see the references at C{Ellipsoid.es_tauf}
            if fabs(d) < eps or fabs(d + p) < e:
                break
            p = d
        else:
            t = unstr(self.toLatLon, eps=eps, unfalse=unfalse)
            raise self._Error(Fmt.no_convergence(d, eps), txt=t)

        a = atan1(T)  # phi, lat
        b = atan2(shx, cy)
        if unfalse and self.falsed:
            b += radians(_cmlon(self.zone))
        ll = _LLEB(degrees90(a), degrees180(b), datum=self.datum, name=self.name)

        # gamma and scale: Karney 2011 Eq 26, 27 and 28
        p = neg(K.ps(-1))
        q =     K.qs(0)
        g = degrees(atan1(tan(y) * tanh(x)) + atan2(q, p))
        k = hypot(p, q) * E.a / A0
        if k:
            k = E.e2s(sin(a)) * hypot1(T) * H / k  # INF?
        ll._iteration = i
        self._latlon5args(ll, g, k, _toBand, unfalse, eps)

    def toRepr(self, prec=0, fmt=Fmt.SQUARE, sep=_COMMASPACE_, B=False, cs=False, **unused):  # PYCHOK expected
        '''Return a string representation of this UTM coordinate.

           Note that UTM coordinates are rounded, not truncated (unlike
           MGRS grid references).

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg fmt: Enclosing backets format (C{str}).
           @kwarg sep: Optional separator between name:value pairs (C{str}).
           @kwarg B: Optionally, include latitudinal band (C{bool}).
           @kwarg cs: Optionally, include meridian convergence and grid
                      scale factor (C{bool} or non-zero C{int} to specify
                      the precison like B{C{prec}}).

           @return: This UTM as a string C{"[Z:09[band], H:N|S, E:meter,
                    N:meter]"} plus C{", C:degrees, S:float"} if C{B{cs}
                    is True} (C{str}).
        '''
        return self._toRepr(fmt, B, cs, prec, sep)

    def toStr(self, prec=0, sep=_SPACE_, B=False, cs=False):  # PYCHOK expected
        '''Return a string representation of this UTM coordinate.

           To distinguish from MGRS grid zone designators, a space is
           left between the zone and the hemisphere.

           Note that UTM coordinates are rounded, not truncated (unlike
           MGRS grid references).

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg sep: Optional separator to join (C{str}) or C{None}
                       to return an unjoined C{tuple} of C{str}s.
           @kwarg B: Optionally, include latitudinal band (C{bool}).
           @kwarg cs: Optionally, include meridian convergence and grid
                      scale factor (C{bool} or non-zero C{int} to specify
                      the precison like B{C{prec}}).

           @return: This UTM as a string with C{zone[band], hemisphere,
                    easting, northing, [convergence, scale]} in
                    C{"00 N|S meter meter"} plus C{" degrees float"} if
                    C{B{cs} is True} (C{str}).
        '''
        return self._toStr(self.hemisphere, B, cs, prec, sep)

    def toUps(self, pole=NN, eps=EPS, falsed=True, **unused):
        '''Convert this UTM coordinate to a UPS coordinate.

           @kwarg pole: Optional top/center of the UPS projection,
                        (C{str}, 'N[orth]'|'S[outh]').
           @kwarg eps: Optional convergence limit, L{EPS} or above
                       (C{float}), see method L{Utm.toLatLon}.
           @kwarg falsed: If C{True}, false both easting and northing
                          (C{bool}).

           @return: The UPS coordinate (L{Ups}).
        '''
        u = self._ups
        if u is None or u.pole != (pole or u.pole) or bool(falsed) != u.falsed:
            ll  =  self.toLatLon(LatLon=_LLEB, eps=eps, unfalse=True)
            ups = _MODS.ups
            self._ups = u = ups.toUps8(ll, Ups=ups.Ups, falsed=falsed,
                                           name=self.name, pole=pole)
        return u

    def toUtm(self, zone, eps=EPS, falsed=True, **unused):
        '''Convert this UTM coordinate to a different zone.

           @arg zone: New UTM zone (C{int}).
           @kwarg eps: Optional convergence limit, L{EPS} or above
                       (C{float}), see method L{Utm.toLatLon}.
           @kwarg falsed: If C{True}, fFalse both easting and northing
                          (C{bool}).

           @return: The UTM coordinate (L{Utm}).
        '''
        if zone == self.zone and falsed == self.falsed:
            return self.copy()
        elif zone:
            u = self._utm
            if u is None or u.zone != zone or bool(falsed) != u.falsed:
                ll = self.toLatLon(LatLon=_LLEB, eps=eps, unfalse=True)
                self._utm = u = toUtm8(ll, Utm=self.classof, falsed=falsed,
                                           name=self.name, zone=zone)
            return u
        raise self._Error(zone=zone)

    @Property_RO
    def zone(self):
        '''Get the (longitudinal) zone (C{int}, 1..60).
        '''
        return self._zone


def _cmlon(zone):
    '''(INTERNAL) Central meridian longitude (C{degrees180}).
    '''
    return (zone * 6) - 183


def _false2(e, n, h):
    '''(INTERNAL) False easting and northing.
    '''
    # Karney, "Test data for the transverse Mercator projection (2009)"
    # <https://GeographicLib.SourceForge.io/C++/doc/transversemercator.html>
    # and <https://Zenodo.org/record/32470#.W4LEJS2ZON8>
    e += _FalseEasting  # make e relative to central meridian
    if h == _S_:
        n += _FalseNorthing  # make n relative to equator
    return e, n


def _parseUTM5(strUTM, datum, Xtm, falsed, Error=UTMError, **name):  # imported by .etm
    '''(INTERNAL) Parse a string representing a UTM coordinate,
       consisting of C{"zone[band] hemisphere easting northing"},
       see L{pygeodesy.parseETM5} and L{parseUTM5}.
    '''
    z, h, e, n, B = _parseUTMUPS5(strUTM, None, Error=Error)
    if _UTM_ZONE_MIN > z or z > _UTM_ZONE_MAX or (B and B not in _Bands):
        raise Error(strUTM=strUTM, zone=z, band=B)

    return UtmUps5Tuple(z, h, e, n, B, Error=Error, **name) if Xtm is None else \
                    Xtm(z, h, e, n, band=B, datum=datum, falsed=bool(falsed), **name)


def parseUTM5(strUTM, datum=_WGS84, Utm=Utm, falsed=True, **name):
    '''Parse a string representing a UTM coordinate, consisting
       of C{"zone[band] hemisphere easting northing"}.

       @arg strUTM: A UTM coordinate (C{str}).
       @kwarg datum: Optional datum to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg Utm: Optional class to return the UTM coordinate (L{Utm})
                   or C{None}.
       @kwarg falsed: Use C{B{falsed}=False} if both easting and northing
                      are I{not falsed} (C{bool}).
       @kwarg name: Optional C{B{name}=NN} (C{str}).

       @return: The UTM coordinate (B{C{Utm}}) or if C{B{Utm} is None}, a
                L{UtmUps5Tuple}C{(zone, hemipole, easting, northing, band)}.
                The C{hemipole} is the C{'N'|'S'} hemisphere.

       @raise UTMError: Invalid B{C{strUTM}}.

       @raise TypeError: Invalid B{C{datum}}.
    '''
    return _parseUTM5(strUTM, datum, Utm, falsed, **name)


def toUtm8(latlon, lon=None, datum=None, Utm=Utm, falsed=True,
                             strict=True, zone=None, **name_cmoff):
    '''Convert a lat-/longitude point to a UTM coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal)
                    geodetic C{LatLon} point.
       @kwarg lon: Longitude (C{degrees}), required if B{C{latlon}} is
                   C{degrees}, ignored otherwise.
       @kwarg datum: Optional datum for this UTM coordinate, overriding
                     B{C{latlon}}'s datum (L{Datum}, L{Ellipsoid}, L{Ellipsoid2}
                     or L{a_f2Tuple}).
       @kwarg Utm: Optional class to return the UTM coordinate (L{Utm}) or C{None}.
       @kwarg falsed: If C{True}, false both easting and northing (C{bool}).
       @kwarg strict: If C{True}, restrict B{C{lat}} to UTM ranges (C{bool}).
       @kwarg zone: Optional UTM zone to enforce (C{int} or C{str}).
       @kwarg name_cmoff: Optional C{B{name}=NN} (C{str}) and DEPRECATED keyword
                   argument C{B{cmoff}=True} to offset the longitude from the zone's
                   central meridian (C{bool}), use C{B{falsed}} instead.

       @return: The UTM coordinate (B{C{Utm}}) or if C{B{Utm} is None} or C{B{falsed}
                is False}, a L{UtmUps8Tuple}C{(zone, hemipole, easting, northing, band,
                datum, gamma, scale)} where C{hemipole} is the C{'N'|'S'} hemisphere.

       @raise RangeError: If B{C{lat}} outside the valid UTM bands or if B{C{lat}}
                          or B{C{lon}} outside the valid range and L{rangerrors
                          <pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{datum}} or B{C{latlon}} not ellipsoidal.

       @raise UTMError: Invalid B{C{zone}}.

       @raise ValueError: If B{C{lon}} is missing or B{C{latlon}} is invalid.

       @note: Implements Karney’s method, using 8-th order Krüger series, giving
              results accurate to 5 nm (or better) for distances up to 3,900 Km
              from the central meridian.
    '''
    z, B, lat, lon, d, f, n = _to7zBlldfn(latlon, lon, datum,
                                          falsed, zone, strict,
                                          UTMError, **name_cmoff)
    d = _ellipsoidal_datum(d, name=n)
    E = d.ellipsoid

    a, b = radians(lat), radians(lon)
    # easting, northing: Karney 2011 Eq 7-14, 29, 35
    sb, cb = sincos2(b)

    T   = tan(a)
    T12 = hypot1(T)
    S   = sinh(E.e * atanh(E.e * T / T12))

    T_ = fdot_(T, hypot1(S), -S, T12)
    H  = hypot(T_, cb)

    y = atan2(T_, cb)  # ξ' ksi
    x = asinh(sb / H)  # η' eta

    A0 = E.A * getattr(Utm, _under(_scale0_), _K0_UTM)  # Utm is class or None

    K = _Kseries(E.AlphaKs, x, y)  # Krüger series
    y =  K.ys(y) * A0  # ξ
    x =  K.xs(x) * A0  # η

    # convergence: Karney 2011 Eq 23, 24
    p_ = K.ps(1)
    q_ = K.qs(0)
    g  = degrees(atan2(T_ * tan(b), hypot1(T_)) + atan2(q_, p_))
    # scale: Karney 2011 Eq 25
    k  = E.e2s(sin(a)) * T12 / H * (A0 / E.a * hypot(p_, q_))

    return _toXtm8(Utm, z, lat, x, y,
                        B, d, g, k, f, n, latlon, EPS)


def _toBand(lat, *unused, **strict_Error):  # see ups._toBand
    '''(INTERNAL) Get the I{latitudinal} Band (row) letter.
    '''
    if _UTM_LAT_MIN <= lat < _UTM_LAT_MAX:  # [-80, 84) like Veness
        return _Bands[int(lat - _UTM_LAT_MIN) >> 3]
    elif _xkwds_get(strict_Error, strict=True):
        r = _range_(_UTM_LAT_MIN, _UTM_LAT_MAX, ropen=True)
        t = _SPACE_(_outside_, _UTM_, _range_, r)
        E = _xkwds_get(strict_Error, Error=RangeError)
        raise E(lat=degDMS(lat), txt=t)
    else:
        return NN  # None


def _toXtm8(Xtm, z, lat, x, y, B, d, g, k, f,  # PYCHOK 13+ args
                 n, latlon, eps_other, Error=UTMError):
    '''(INTERNAL) Helper for methods L{toEtm8} and L{toUtm8}.
    '''
    h = _hemi(lat)
    if f:
        x, y = _false2(x, y, h)
    if Xtm is None:  # DEPRECATED
        r = UtmUps8Tuple(z, h, x, y, B, d, g, k, Error=Error, name=n)
    else:
        r = Xtm(z, h, x, y, band=B, datum=d, falsed=f,
                                    gamma=g, scale=k, name=n)
        if isinstance(latlon, _LLEB) and d is latlon.datum:  # see ups.toUtm8
            r._latlon5args(latlon, g, k, _toBand, f, eps_other)  # XXX weakref(latlon)?
        elif not r._band:
            r._band = _toBand(lat)
    return r


def _to3zBlat(zone, band, Error=UTMError):  # in .mgrs
    '''(INTERNAL) Check and return zone, Band and band latitude.

       @arg zone: Zone number or string.
       @arg band: Band letter.
       @arg Error: Exception to raise (L{UTMError}).

       @return: 3-Tuple (zone, Band, latitude).
    '''
    z, B, _ = _to3zBhp(zone, band, Error=Error)
    if not (_UTM_ZONE_MIN <= z <= _UTM_ZONE_MAX or
           (_UPS_ZONE == z and Error is MGRSError)):
        raise Error(zone=zone)

    b = None
    if B:
        if z == _UPS_ZONE:  # polar
            try:
                b = Lat(_UPS_LATS[B], name=_bandLat_)
            except KeyError:
                raise Error(band=band or B, zone=z)
        else:  # UTM
            b = _Bands.find(B)
            if b < 0:
                raise Error(band=band or B, zone=z)
            b = Int((b << 3) - 80, name=_bandLat_)
        B = Band(B)
    elif Error is not UTMError:
        raise Error(band=band, txt=MISSING)

    return Zone(z), B, b


def _to4zBll(lat, lon, cmoff=True, strict=True, Error=RangeError):
    '''(INTERNAL) Return zone, Band and lat- and (central) longitude in degrees.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg cmoff: If C{True}, offset B{C{lon}} from the zone's central meridian.
       @kwarg strict: Restrict B{C{lat}} to the UTM ranges (C{bool}).
       @kwarg Error: Error for out of UTM range B{C{lat}}s.

       @return: 4-Tuple (zone, Band, lat, lon).
    '''
    z, lat, lon = _to3zll(lat, lon)  # in .utmupsBase

    x = lon - _cmlon(z)  # z before Norway/Svalbard
    if fabs(x) > _UTM_ZONE_OFF_MAX:
        t = _SPACE_(_outside_, _UTM_, _zone_, str(z), _by_, degDMS(x, prec=6))
        raise Error(lon=degDMS(lon), txt=t)

    B = _toBand(lat, strict=strict, Error=Error)
    if B == _X_:  # and 0 <= lon < 42: z = int(lon + 183) // 6 + 1
        x = _SvalbardXzone.get(z, None)
        if x:  # Svalbard/Spitsbergen archipelago
            z += 1 if lon >= x else -1
    elif B == _V_ and z == 31 and lon >= 3:
        z += 1  # SouthWestern Norway

    if cmoff:  # lon off central meridian
        lon -= _cmlon(z)  # z I{after} Norway/Svalbard
    return Zone(z), (Band(B) if B else None), Lat(lat), Lon(lon)


def _to7zBlldfn(latlon, lon, datum, falsed, zone, strict, Error, **name_cmoff):
    '''(INTERNAL) Determine 7-tuple (zone, band, lat, lon, datum,
        falsed, name) for methods L{toEtm8} and L{toUtm8}.
    '''
    f, name = _xkwds_pop2(name_cmoff, cmoff=falsed)  # DEPRECATED
    lat, lon, d, n = _to4lldn(latlon, lon, datum, name)
    z, B, lat, lon = _to4zBll(lat, lon, cmoff=f, strict=strict)
    if zone:  # re-zone for ETM/UTM
        r, _, _ = _to3zBhp(zone, B)
        if r != z:
            if not _UTM_ZONE_MIN <= r <= _UTM_ZONE_MAX:
                raise Error(zone=zone)
            if f:  # re-offset from central meridian
                lon += _cmlon(z) - _cmlon(r)
            z = r
    return z, B, lat, lon, d, f, n


def utmZoneBand5(lat, lon, cmoff=False, **name):
    '''Return the UTM zone number, Band letter, hemisphere and
       (clipped) lat- and longitude for a given location.

       @arg lat: Latitude in degrees (C{scalar} or C{str}).
       @arg lon: Longitude in degrees (C{scalar} or C{str}).
       @kwarg cmoff: If C{True}, offset longitude from the zone's central
                     meridian (C{bool}).
       @kwarg name: Optional C{B{name}=NN} (C{str}).

       @return: A L{UtmUpsLatLon5Tuple}C{(zone, band, hemipole, lat, lon)}
                where C{hemipole} is the C{'N'|'S'} UTM hemisphere.

       @raise RangeError: If B{C{lat}} outside the valid UTM bands or if
                          B{C{lat}} or B{C{lon}} outside the valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise ValueError: Invalid B{C{lat}} or B{C{lon}}.
    '''
    lat, lon = parseDMS2(lat, lon)
    z, B, lat, lon = _to4zBll(lat, lon, cmoff=cmoff)
    return UtmUpsLatLon5Tuple(z, B, _hemi(lat), lat, lon, **name)

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
