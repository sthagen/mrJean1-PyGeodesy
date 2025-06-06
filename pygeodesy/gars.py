
# -*- coding: utf-8 -*-

u'''I{Global Area Reference System} (GARS) en-/decoding.

Class L{Garef} and several functions to encode, decode and inspect
GARS references.

Transcoded from I{Charles Karney}'s C++ class U{GARS
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1GARS.html>}.

@see: U{Global Area Reference System
      <https://WikiPedia.org/wiki/Global_Area_Reference_System>} and U{NGA
      (GARS)<https://Earth-Info.NGA.mil/GandG/coordsys/grids/gars.html>}.
'''

from pygeodesy.basics import isstr, _splituple,  typename
from pygeodesy.constants import _off90, _1_over, _0_5
from pygeodesy.errors import _ValueError, _xkwds, _xStrError
# from pygeodesy.internals import typename  # from .basics
from pygeodesy.interns import NN, _0to9_, _AtoZnoIO_, _INV_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _name__,  Fmt, Property_RO
from pygeodesy.namedTuples import LatLon2Tuple, LatLonPrec3Tuple
# from pygeodesy.props import Property_RO  # from .named
# from pygeodesy.streprs import Fmt  # from .named
from pygeodesy.units import Int_, Lat, Lon, Precision_, Scalar_, Str

from math import floor

__all__ = _ALL_LAZY.gars
__version__ = '25.05.07'

_Digits  = _0to9_
_LatLen  =    2
_LatOrig =  -90
_Letters = _AtoZnoIO_
_LonLen  =    3
_LonOrig = -180
_MaxPrec =    2

_MinLen = _LonLen + _LatLen
_MaxLen = _MinLen + _MaxPrec

_M1 = _M2 = 2
_M3 =  3
_M4 = _M1 * _M2 * _M3

_LatOrig_M4   = _LatOrig * _M4
_LatOrig_M1   = _LatOrig * _M1
_LonOrig_M4   = _LonOrig * _M4
_LonOrig_M1_1 = _LonOrig * _M1 - 1

_Resolutions  = _1_over(_M1), _1_over(_M1 * _M2), _1_over(_M4)


def _2divmod2(ll, _Orig_M4):
    x = int(floor(ll * _M4)) - _Orig_M4
    i = (x * _M1) // _M4
    x -= i * _M4 // _M1
    return i, x


def _2fll(lat, lon, *unused):
    '''(INTERNAL) Convert lat, lon.
    '''
    # lat, lon = parseDMS2(lat, lon)
    return (Lat(lat, Error=GARSError),
            Lon(lon, Error=GARSError))


# def _2Garef(garef):
#     '''(INTERNAL) Check or create a L{Garef} instance.
#     '''
#     if not isinstance(garef, Garef):
#         try:
#             garef = Garef(garef)
#         except (TypeError, ValueError):
#             raise _xStrError(Garef, Str, garef=garef)
#     return garef


def _2garstr2(garef):
    '''(INTERNAL) Check a garef string.
    '''
    try:
        n, garstr = len(garef), garef.upper()
        if n < _MinLen or n > _MaxLen \
                       or garstr.startswith(_INV_) \
                       or not garstr.isalnum():
            raise ValueError()
        return garstr, _2Precision(n - _MinLen)

    except (AttributeError, TypeError, ValueError) as x:
        raise GARSError(typename(Garef), garef, cause=x)


def _2Precision(precision):
    '''(INTERNAL) Return a L{Precision_} instance.
    '''
    return Precision_(precision, Error=GARSError, low=0, high=_MaxPrec)


class GARSError(_ValueError):
    '''Global Area Reference System (GARS) encode, decode or other L{Garef} issue.
    '''
    pass


class Garef(Str):
    '''Garef class, a named C{str}.
    '''
    # no str.__init__ in Python 3
    def __new__(cls, lat_gll, lon=None, precision=1, **name):
        '''New L{Garef} from an other L{Garef} instance or garef C{str}
           or from a lat- and longitude.

           @arg lat_gll: Latitude (C{degrees90}), a garef (L{Garef},
                         C{str}) or a location (C{LatLon}, C{LatLon*Tuple}).
           @kwarg lon: Logitude (C{degrees180)}, required if B{C{lat_gll}}
                       is C{degrees90}, ignored otherwise.
           @kwarg precision: The desired garef resolution and length (C{int}
                             0..2), see L{encode<pygeodesy.gars.encode>}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: New L{Garef}.

           @raise GARSError: INValid B{C{lat_gll}}.

           @raise RangeError: Invalid B{C{lat_gll}} or B{C{lon}}.

           @raise TypeError: Invalid B{C{lat_gll}} or B{C{lon}}.
        '''
        if lon is None:
            if isinstance(lat_gll, Garef):
                g, ll, p = str(lat_gll), lat_gll.latlon, lat_gll.precision
            elif isstr(lat_gll):
                ll = _splituple(lat_gll)
                if len(ll) > 1:
                    g, ll, p = _encode3(ll[0], ll[1], precision)
                else:
                    g, ll =  lat_gll.upper(), None
                    _, p  = _2garstr2(g)  # validate
            else:  # assume LatLon
                try:
                    g, ll, p = _encode3(lat_gll.lat, lat_gll.lon, precision)
                except AttributeError:
                    raise _xStrError(Garef, gll=lat_gll, Error=GARSError)
        else:
            g, ll, p = _encode3(lat_gll, lon, precision)

        self = Str.__new__(cls, g, name=_name__(name, _or_nameof=lat_gll))
        self._latlon    = ll
        self._precision = p
        return self

    @Property_RO
    def decoded3(self):
        '''Get this garef's attributes (L{LatLonPrec3Tuple}).
        '''
        lat, lon = self.latlon
        return LatLonPrec3Tuple(lat, lon, self.precision, name=self.name)

    @Property_RO
    def _decoded3(self):
        '''(INTERNAL) Initial L{LatLonPrec3Tuple}.
        '''
        return decode3(self)

    @Property_RO
    def latlon(self):
        '''Get this garef's (center) lat- and longitude (L{LatLon2Tuple}).
        '''
        lat, lon = self._latlon or self._decoded3[:2]
        return LatLon2Tuple(lat, lon, name=self.name)

    @Property_RO
    def precision(self):
        '''Get this garef's precision (C{int}).
        '''
        p = self._precision
        return self._decoded3.precision if p is None else p

    def toLatLon(self, LatLon=None, **LatLon_kwds):
        '''Return (the center of) this garef cell as an instance
           of the supplied C{LatLon} class.

           @kwarg LatLon: Class to use (C{LatLon} or C{None}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                               keyword arguments.

           @return: This garef location as B{C{LatLon}} or if
                    C{B{LatLon} is None} as L{LatLonPrec3Tuple}.
        '''
        return self.decoded3 if LatLon is None else LatLon(
              *self.latlon, **_xkwds(LatLon_kwds, name=self.name))


def decode3(garef, center=True, **name):
    '''Decode a C{garef} to lat-, longitude and precision.

       @arg garef: To be decoded (L{Garef} or C{str}).
       @kwarg center: If C{True}, use the garef's center, otherwise
                      the south-west, lower-left corner (C{bool}).

       @return: A L{LatLonPrec3Tuple}C{(lat, lon, precision)}.

       @raise GARSError: Invalid B{C{garef}}, INValid, non-alphanumeric
                         or bad length B{C{garef}}.
    '''
    def _Error(i):
        return GARSError(garef=Fmt.SQUARE(repr(garef), i))

    def _ll(chars, g, i, j, lo, hi):
        ll, b = 0, len(chars)
        for i in range(i, j):
            d = chars.find(g[i])
            if d < 0:
                raise _Error(i)
            ll = ll * b + d
        if ll < lo or ll > hi:
            raise _Error(j)
        return ll

    def _ll2(lon, lat, g, i, m):
        d = _Digits.find(g[i])
        if d < 1 or d > m * m:
            raise _Error(i)
        d, r = divmod(d - 1, m)
        lon = lon * m + r
        lat = lat * m + (m - 1 - d)
        return lon, lat

    g, precision = _2garstr2(garef)

    lon = _ll(_Digits,  g,       0, _LonLen, 1, 720) + _LonOrig_M1_1
    lat = _ll(_Letters, g, _LonLen, _MinLen, 0, 359) + _LatOrig_M1
    if precision > 0:
        lon, lat = _ll2(lon, lat, g, _MinLen, _M2)
        if precision > 1:
            lon, lat = _ll2(lon, lat, g, _MinLen + 1, _M3)

    if center:  # ll = (ll * 2 + 1) / 2
        lon += _0_5
        lat += _0_5

    n = _name__(name, _or_nameof=garef)
    r = _Resolutions[precision]  # == 1.0 / unit
    return LatLonPrec3Tuple(Lat(lat * r, Error=GARSError),
                            Lon(lon * r, Error=GARSError),
                            precision, name=n)


def encode(lat, lon, precision=1):
    '''Encode a lat-/longitude as a C{garef} of the given precision.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg precision: Optional, the desired C{garef} resolution
                         and length (C{int} 0..2).

       @return: The C{garef} (C{str}).

       @raise RangeError: Invalid B{C{lat}} or B{C{lon}}.

       @raise GARSError: Invalid B{C{precision}}.

       @note: The C{garef} length is M{precision + 5} and the C{garef}
              resolution is B{30′} for B{C{precision}} 0, B{15′} for 1
              and B{5′} for 2, respectively.
    '''
    g, _, _ = _encode3(lat, lon, precision)
    return g


def _encode3(lat, lon, precision):  # MCCABE 14
    '''Return 3-tuple C{(garef, (lat, lon), p)}.
    '''
    def _digit(x, y, m):
        return _Digits[m * (m - y - 1) + x + 1],

    def _str(chars, x, n):
        s, b = [], len(chars)
        for i in range(n):
            x, i = divmod(x, b)
            s.append(chars[i])
        return tuple(reversed(s))

    p = _2Precision(precision)

    lat, lon = _2fll(lat, lon)

    ix, x = _2divmod2(       lon,  _LonOrig_M4)
    iy, y = _2divmod2(_off90(lat), _LatOrig_M4)

    g = _str(_Digits,  ix + 1, _LonLen) + \
        _str(_Letters, iy,     _LatLen)
    if p > 0:
        ix, x = divmod(x, _M3)
        iy, y = divmod(y, _M3)
        g += _digit(ix, iy, _M2)
        if p > 1:
            g += _digit(x, y, _M3)

    return NN.join(g), (lat, lon), p


def precision(res):
    '''Determine the L{Garef} precision to meet a required (geographic)
       resolution.

       @arg res: The required resolution (C{degrees}).

       @return: The L{Garef} precision (C{int} 0..2).

       @raise ValueError: Invalid B{C{res}}.

       @see: Function L{gars.encode} for more C{precision} details.
    '''
    r = Scalar_(res=res)
    for p in range(_MaxPrec):
        if resolution(p) <= r:
            return p
    return _MaxPrec


def resolution(prec):
    '''Determine the (geographic) resolution of a given L{Garef} precision.

       @arg prec: The given precision (C{int}).

       @return: The (geographic) resolution (C{degrees}).

       @raise GARSError: Invalid B{C{prec}}.

       @see: Function L{gars.encode} for more C{precision} details.
    '''
    p = Int_(prec=prec, Error=GARSError, low=-1, high=_MaxPrec + 1)
    return _Resolutions[max(0, min(p, _MaxPrec))]


__all__ += _ALL_DOCS(decode3,  # functions
                     encode, precision, resolution)

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
