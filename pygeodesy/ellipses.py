
# -*- coding: utf-8 -*-

u'''Class C{Ellipse} for 2-D ellipse attributes, like perimeter, area, etc.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.constants import EPS, EPS_2, INT0, PI, PI_2, PI2, \
                               _0_0, _1_0, _4_0, _over, _1_over
from pygeodesy.constants import _0_5, _3_0, _10_0, MANT_DIG as _DIG53  # PYCHOK used!
# from pygeodesy.ellipsoids import Ellipsoid  # _MODS
from pygeodesy.errors import _ConvergenceError, _ValueError
from pygeodesy.fmath import fhorner, hypot
from pygeodesy.fsums import _fsum  # PYCHOK used!
from pygeodesy.internals import typename,  _DOT_
# from pygeodesy.interns import _DOT_  # from .internals
# from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS  # from .utily
from pygeodesy.named import _Named,  unstr
from pygeodesy.props import Property_RO, property_RO, property_ROnce
# from pygeodesy.streprs import unstr  # from .named
# from pygeodesy.triaxials import Triaxial_  # _MODS
from pygeodesy.utily import atan2, sincos2,  _ALL_LAZY, _MODS

from math import fabs, radians, sqrt
# import operator as _operator  # from .fmath

__all__ = _ALL_LAZY.ellipses
__version__ = '26.02.12'

_TOL53    =  sqrt(EPS_2)     # sqrt(pow(_0_5, _DIG53))
_TOL53_53 = _TOL53 / _DIG53  # "flat" b/a tolerance, 1.9e-10
# assert _DIG53 == 53


class Ellipse(_Named):
    '''Class to compute attributes of a 2-D ellipse, like perimeter, area and arcs.
    '''
    _flat  =  False
    _maxit = _DIG53

    def __init__(self, a, b, **name):
        '''New L{Ellipse} with semi-axes B{C{a}} and B{C{b}}.

           The ellipse is oblate if C{a > b}, prolate if C{a < b}
           circular if C{a == b} and "flat" if C{min(a, b) near 0}.

           @arg a: X semi-axis length (C{meter}, conventionally).
           @arg b: Y semi-axis length (C{meter}, conventionally).

           @raise ValueError: Invalid B{C{a}} or B{C{b}}.
        '''
        if name:
            self.name = name
        self._a = a
        self._b = b
        r = a < b
        if r:  # prolate
            a, b = b, a
        if b < 0:  # PYCHOK no cover
            raise self._Error(None)
        if a > b:
            if _isFlat(a, b):
                self._flat = True
                P = a * _4_0
            else:  # pro-/oblate
                P = None
        else:  # circle
            P = a * PI2
        self._Pab4 = r, P, a, b

    @Property_RO
    def a(self):
        '''Return semi-axis C{B{a}} of this ellipse (C{meter}, conventionally).
        '''
        return self._a

    def arc(self, deg2, deg1=0):
        '''Compute the length of U{elliptic arc<https://www.JohnDCook.com/blog/2022/11/02/elliptic-arc-length/>}
           C{(B{deg2} - B{deg1})}, both counter-clockwise from semi-axis B{C{a}} to B{C{b}} of the ellipse.

           @arg deg2: End angle of the elliptic arc (C{degrees}).
           @kwarg deg1: Start angle of the elliptic arc (C{degrees}).

           @return: Arc length, signed (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).
        '''
        return self.arc_(radians(deg2), (radians(deg1) if deg1 else _0_0))

    def arc_(self, rad2, rad1=0):
        '''Compute the length of U{elliptic arc<https://www.JohnDCook.com/blog/2022/11/02/elliptic-arc-length/>}
           C{(B{rad2} - B{rad1})}, both counter-clockwise from semi-axis B{C{a}} to B{C{b}} of the ellipse.

           @arg rad2: End angle of the elliptic arc (C{radians}).
           @kwarg rad1: Start angle of the elliptic arc (C{radians}).

           @return: Arc length, signed (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).
        '''
        r, L, a, _ = self._Pab4
        if L is None:
            _e = self._ellipe or self._ellipE
            k  = self._k
            r  = PI_2 if r else _0_0
            L  = self._arc(_e, k, r + rad2)
            r += rad1
            if r:
                L -= self._arc(_e, k, r)
            L *= a
        else:
            L *= (rad2 - rad1) / PI2
        return L

    def _arc(self, _e, k, r):
        '''(INTERNAL) Helper for method C{.arc_}.
        '''
        t, r = divmod(r, PI2)
        L = _e(k, r)  # phi=r
        if t:  # + t * perimeter
            t *= _e(k) * _4_0
            L +=  t
        return L

    @Property_RO
    def area(self):
        '''Return the area of this ellipse (C{meter**2}, conventionally).
        '''
        return self.a * self.b * PI

    @Property_RO
    def b(self):
        '''Return semi-axis C{B{b}} of this ellipse (C{meter}, conventionally).
        '''
        return self._b

    @Property_RO
    def _Ek(self):
        '''(INTERNAL) Get the C{Elliptic(k)} instance.
        '''
        return _MODS.elliptic._Ek(self._k)

    def _ellipE(self, k, phi=None):  # PYCHOK k
        '''(INTERNAL) Get the in-/complete integral of the 2nd kind.
        '''
        # assert k == self._Ek.k2
        return self._Ek.cE if phi is None else self._Ek.fE(phi)

    @property_ROnce
    def _ellipe(self):
        '''(INTERNAL) Wrap functions C{scipy.special.ellipe} and C{-.ellipeinc}, I{once}.
        '''
        try:
            from scipy.special import ellipe, ellipeinc

            def _ellipe(k, phi=None):
                r = ellipe(k) if phi is None else ellipeinc(phi, k)
                return float(r)

        except (AttributeError, ImportError):
            _ellipe = None
        return _ellipe  # overwrite property_ROnce

    def _Error(self, where, **cause):
        '''(INTERNAL) Build an error.
        '''
        t = self.named3
        u = unstr(t, a=self.a, b=self.b)
        if where:
            t =  typename(where, where)
            u = _DOT_(u, t)
        return _ValueError(u, **cause)

    @Property_RO
    def foci(self):
        '''Get the foci lengths (C{meter}, conventionally), C{positive} if
           the ellipse is oblate with foci on semi-axis C{a}, C{negative}
           if prolate with foci on semi-axis C{b} or C{0} if circular with
           coincident foci in the ellipse' center.
        '''
        c = self._k
        if c:
            r, _, a, _ = self._Pab4
            c = sqrt(c) * a
            if r:  # prolate
                c = -c
        return c

    @property_ROnce
    def _GKs(self):
        '''(INTERNAL) Compute the coefficients for property C{.perimeterGK}, I{once}.
        '''
        # U{numerators<https://OEIS.org/A056981>}, U{denominators<https://OEIS.org/A056982>}
        return (1, 1 / 4, 1 / 64, 1 / 256, 25 / 16384, 49 / 65536,
                441 / 1048576, 1089 / 4194304)  # overwrite property_ROnce

    def _HGKs(self, h, maxit):
        '''(INTERNAL) Yield the terms for property C{.perimeterHGK}.
        '''
        s = t = _1_0
        yield s
        for u in range(-1, maxit * 2, 2):
            t *= u / (u + 3) * h
            t2 = t**2
            yield t2
            p  = s
            s += t2
            if s == p:  # 44 trips
                break
        else:  # PYCHOK no cover
            raise _ConvergenceError(maxit, s, p)

    @property_RO
    def isCircular(self):
        '''Is this ellipse circular? (C{bool})
        '''
        return self.a == self.b

    @property_RO
    def isFlat(self):
        '''Is this ellipse "flat", too pro-/oblate? (C{bool})
        '''
        return self._flat

    @property_RO
    def isOblate(self):
        '''Is this ellipse oblate? (C{bool})
        '''
        return self.a > self.b

    @property_RO
    def isProlate(self):
        '''Is this ellipse prolate? (C{bool})
        '''
        return self.a < self.b

    @Property_RO
    def _k(self):
        '''(INTERNAL) Get C{0 <= k <= 1}.
        '''
        # C{k} is aka Elliptic C{k2} and SciPy's C{m}
        _, _, a, b = self._Pab4
        return ((_1_0 - (b / a)**2) if 0 < b else _1_0) if b < a else _0_0

    @Property_RO
    def perimeterAGM(self):
        '''Compute the perimeter of this ellipse using the U{Arithmetic-Geometric Mean
           <https://PaulBourke.net/geometry/ellipsecirc>} formula (C{meter}, same units
           as semi-axes B{C{a}} and B{C{b}}).
        '''
        _, P, a, b = self._Pab4
        if P is None:
            t  = _TOL53
            m  = -1
            c  = a + b
            ds = [c**2]
            _d = ds.append
            for _ in range(self._maxit):  # 4..5 trips
                b  = sqrt(a * b)
                a  = c * _0_5
                c  = a + b
                d  = a - b
                m *= 2
                _d(d**2 * m)
                if d <= (b * t):
                    break
            else:  # PYCHOK no cover
                raise _ConvergenceError(self._maxit, _over(d, b), t)
            P = _over(_fsum(ds) * PI, c)  # nonfinites=True
        return P

    @Property_RO
    def perimeter4Arc3(self):
        '''Compute the perimeter (and arcs) of this ellipse using the U{4-Arc
           <https://PaulBourke.net/geometry/ellipsecirc>} approximation as a
           3-Tuple C{(P, Ra, Rb)} with perimeter C{P}, arc radii C{Ra} and
           C{Rb} at the respective semi-axes (all in C{meter}, same units as
           semi-axes B{C{a}} and B{C{b}}).
        '''
        r, P, a, b = self._Pab4
        if P is None:
            h = hypot(a, b)
            t = atan2(b, a)
            s, c = sincos2(t)
            L  = (h - (a - b)) * _0_5
            Ra = _over(L, c)
            Rb = _over(h - L, s)
            P  = (t * Rb + (PI_2 - t) * Ra) * _4_0
        elif a > b:  # flat
            Ra, Rb = _0_0, _1_over(b)  # INF
        else:  # circle
            Ra, Rb =  a, b
        return (P, Rb, Ra) if r else (P, Ra, Rb)

#   @Property_RO
#   def perimeterCR(self):
#       '''Compute the perimeter of this ellipse using U{Rackauckas'
#          <https://www.ChrisRackauckas.com/assets/Papers/ChrisRackauckas-The_Circumference_of_an_Ellipse.pdf>}
#          approximation, also U{here<https://ExtremeLearning.com.AU/a-formula-for-the-perimeter-of-an-ellipse>}
#          (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).
#       '''
#       _, P, a, b = self._Pab4
#       if P is None:
#           P  =   a + b
#           h  = ((a - b) / P)**2
#           P *= (fhorner(h, 135168, -85760,  -5568, 3867) /
#                 fhorner(h, 135168, -119552, 22208, 345)) * PI
#       return P

    @Property_RO
    def perimeterGK(self):
        '''Compute the perimeter of this ellipse using the U{Gauss-Kummer
           <https://www.JohnDCook.com/blog/2023/05/28/approximate-ellipse-perimeter>} series,
           C{B{b / a} > 0.75} (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).
        '''
        P, h = self._Ph2
        if h:
            P *= fhorner(h**2, *self._GKs) * PI
        return P

    @Property_RO
    def perimeter2k(self):
        '''Compute the perimeter of this ellipse using the complete integral of the 2nd
           kind, C{Elliptic.cE} (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).
        '''
        return self._perimeter2k(self._ellipE)

    @Property_RO
    def perimeter2k_(self):
        '''Compute the perimeter of this ellipse using U{SciPy's ellipe
           <https://www.JohnDCook.com/perimeter_ellipse.html>} function
           if available, otherwise use property C{perimeter2k} (C{meter},
           same units as semi-axes B{C{a}} and B{C{b}}).
        '''
        return self._perimeter2k(self._ellipe or self._ellipE)

    def _perimeter2k(self, _ellip):
        '''(INTERNAL) Helper for methods C{.PE2k} and C{.Pe2k}.
        '''
        _, P, a, _ = self._Pab4
        if P is None:  # see .ellipsoids.Ellipsoid.L
            k =  self._k
            P = _ellip(k) * a * _4_0
        return P

    @Property_RO
    def _Ph2(self):
        _, P, a, b = self._Pab4
        if P is None:
            P =  a + b
            h = (a - b) / P
        else:
            h =  None
        return P, h

    @Property_RO
    def perimeterHGK(self):
        '''Compute the perimeter of this ellipse using the U{Hypergeometric Gauss-Kummer
           <https://web.Tecnico.ULisboa.PT/~mcasquilho/compute/com/,ellips/PerimeterOfEllipse.pdf>}
           series (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).
        '''
        P, h = self._Ph2
        if h:
            hs =  self._HGKs(h, self._maxit)
            P *= _fsum(hs) * PI  # nonfinites=True
        return P

#   @Property_RO
#   def perimeterLS(self):
#       '''Compute the perimeter of this ellipse using the U{Linderholm-Segal
#          <https://www.JohnDCook.com/blog/2021/03/24/perimeter-of-an-ellipse>}
#          formula, aka C{3/2 norm} (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).
#       '''
#       _, P, a, b = self._Pab4
#       if P is None:
#           n = pow(a, _1_5) + pow(b, _1_5)
#           P = pow(n * _0_5, _2_3rd) * PI2
#       return P

    @Property_RO
    def perimeter2R(self):
        '''Compute the perimeter of this ellipse using U{Ramanujan's 2nd
           <https://PaulBourke.net/geometry/ellipsecirc>} approximation,
           C{B{b / a} > 0.9} (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).
        '''
        P, h = self._Ph2
        if h:
            h *= _3_0 * h
            h /=  sqrt(_4_0 - h) + _10_0  # /= chokes PyChecker?
            P *= (h + _1_0) * PI
        return P

    @Property_RO
    def R2(self):
        '''Compute the I{authalic} radius of this ellipse, C{sqrt(B{a} * B{b})}
           (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).
        '''
        a, b = self.a, self.b
        return sqrt(a * b) if a != b else float(a)

    def Roc_(self, rad_x, *y):
        '''Compute the U{radius of curvature<https://WikiPedia.org/wiki/Radius_of_curvature>}
           at angle C{B{rad} radians} or at point C{(B{x}, B{y})} I{on} this ellipse.

           @return: Curvature (C{meter}, same units as semi-axes B{C{a}} and B{C{b}}).

           @raise ValueError: C{(B{x}, B{y})} not located on the ellipse.
        '''
        try:
            a, b = self.a, self.b
            if b != a:
                r = float(rad_x)
                if y:
                    x, y = r, float(y[0])
                    if self._sideOf(x, y, EPS):
                        raise _ValueError(x=x, y=y)
                    r = atan2(y, x)
                s, c = sincos2(r)
                r = _over(hypot(s * a, c * b)**3, a * b)
            else:  # circle
                r = float(a)
        except Exception as x:
            raise self._Error(self.Roc_, cause=x)
        return r

    def _sideOf(self, x, y, eps):
        '''(INTERNAL) Helper for methods C{.Roc} and C{.sideOf}.
        '''
        a, b = self.a, self.b
        s = hypot(x * b, y * a) - (a * b)
        return INT0 if fabs(s) < eps else s

    def sideOf(self, x, y, eps=EPS):
        '''Return a C{positive}, C{negative} or C{0} scalar if point C{(B{x}, B{y})}
           is C{outside}, C{inside} respectively C{on} this ellipse.
        '''
        try:
            return self._sideOf(x, y, eps)
        except Exception as X:
            raise self._Error(self.sideOf, x=x, y=y, cause=X)

    def toEllipsoid(self):
        '''Return an L{Ellipsoid<pygeodesy.Ellipsoid>} from this ellipse's
           C{a} and C{b} semi-axes.
        '''
        return _MODS.ellipsoids.Ellipsoid(self.a, b=self.b, name=self.name)

    def toTriaxial_(self, c=EPS):
        '''Return a L{Triaxial_<pygeodesy.Triaxial_>} from this ellipse's
           C{a} and C{b} semi-axes with B{C{c}} as minor semi-axis.
        '''
        return _MODS.triaxials.Triaxial_(self.a, b=self.b, c=c, name=self.name)  # 'NN'


def _isFlat(a, b):  # in .triaxials.bases
    '''(INTERNAL) Is C{b <<< a}?
    '''
    return b < (a * _TOL53_53)

# **) MIT License
#
# Copyright (C) 2026-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
