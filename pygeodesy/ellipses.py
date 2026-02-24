
# -*- coding: utf-8 -*-

u'''Class C{Ellipse} for 2-D ellipse attributes, like perimeter, area, etc.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

# from pygeodesy.basics import islistuple  # _MODS
from pygeodesy.constants import EPS, EPS_2, INT0, NEG0, PI, PI_2, PI3_2, PI2, \
                               _0_0, _1_0, _4_0, _isfinite, _over, _1_over  # _N_1_0
from pygeodesy.constants import _0_5, _3_0, _10_0, MANT_DIG as _DIG53  # PYCHOK used!
# from pygeodesy.ellipsoids import Ellipsoid  # _MODS
from pygeodesy.errors import _ConvergenceError, _ValueError
from pygeodesy.fmath import fhorner, hypot
from pygeodesy.fsums import _fsum  # PYCHOK used!
from pygeodesy.internals import typename,  _DOT_
# from pygeodesy.interns import _DOT_  # from .internals
# from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS  # from .utily
from pygeodesy.named import _NamedBase,  unstr
from pygeodesy.props import Property_RO, property_RO, property_ROnce
# from pygeodesy.streprs import unstr  # from .named
# from pygeodesy.triaxials import Triaxial_ , TriaxialError # _MODS
from pygeodesy.units import Degrees, Meter, Meter2, Radians, Radius, Scalar
from pygeodesy.utily import atan2, sincos2, sincos2d,  _ALL_LAZY, _MODS
# from pygeodesy.vector3d import Vector3d  # _MODS

from math import degrees, fabs, radians, sqrt
# import operator as _operator  # from .fmath

__all__ = _ALL_LAZY.ellipses
__version__ = '26.02.23'

_TOL53    =  sqrt(EPS_2)     # sqrt(pow(_0_5, _DIG53))
_TOL53_53 = _TOL53 / _DIG53  # "flat" b/a tolerance, 1.9e-10
# assert _DIG53 == 53


class Ellipse(_NamedBase):
    '''Class to compute various attributes of a 2-D ellipse.
    '''
#   _ab3   = (a, b, a * b)  # unordered
    _flat  =  False
    _maxit = _DIG53
#   _Pab4  = (r, P, a, b)  # a >= b, ordered

    def __init__(self, a, b, **name):
        '''New L{Ellipse} with semi-axes B{C{a}} and B{C{b}}.

           The ellipse is C{oblate} if C{B{a} > B{b}}, C{prolate} if
           C{B{a} < B{b}}, C{circular} if C{B{a} == B{b}} and C{"flat"}
           if C{min(B{a}, B{b}) <<< max(B{a}, B{b})}.

           @arg a: X semi-axis length (C{meter}, conventionally).
           @arg b: Y semi-axis length (C{meter}, conventionally).

           @raise ValueError: Invalid B{C{a}} or B{C{b}}.
        '''
        if name:
            self.name = name
        self._ab3 = a, b, (a * b)  # unordered

        r = a < b
        if r:  # prolate
            a, b = b, a
        if b < 0 or not _isfinite(a):  # PYCHOK no cover
            raise self._Error(None)
        if a > b:
            if _isFlat(a, b):
                self._flat = True
                P = a * _4_0
            else:  # pro-/oblate
                P = None
        else:  # circular
            P = a * PI2
        self._Pab4 = r, P, a, b  # ordered

    @Property_RO
    def a(self):
        '''Get semi-axis C{B{a}} of this ellipse (C{meter}, conventionally).
        '''
        a, _, _ = self._ab3
        return Meter(a=a)

    def arc(self, deg2, deg1=0):
        '''Compute the length of U{elliptic arc<https://www.JohnDCook.com/blog/2022/11/02/elliptic-arc-length/>}
           C{(B{deg2} - B{deg1})}, both counter-clockwise from semi-axis B{C{a}} to B{C{b}} of the ellipse.

           @arg deg2: End angle of the elliptic arc (C{degrees}).
           @kwarg deg1: Start angle of the elliptic arc (C{degrees}).

           @return: Arc length, signed (C{meter}, conventionally).
        '''
        return self.arc_(radians(deg2), (radians(deg1) if deg1 else _0_0))

    def arc_(self, rad2, rad1=0):
        '''Compute the length of U{elliptic arc<https://www.JohnDCook.com/blog/2022/11/02/elliptic-arc-length/>}
           C{(B{rad2} - B{rad1})}, both counter-clockwise from semi-axis B{C{a}} to B{C{b}} of the ellipse.

           @arg rad2: End angle of the elliptic arc (C{radians}).
           @kwarg rad1: Start angle of the elliptic arc (C{radians}).

           @return: Arc length, signed (C{meter}, conventionally).
        '''
        r, L, a, _ = self._Pab4
        if L is None:
            _e = self._ellipe or self._ellipE
            k  = self.e2
            r  = PI_2 if r else _0_0
            L  = self._arc(_e, k, r + rad2)
            r += rad1
            if r:
                L -= self._arc(_e, k, r)
            L *= a
        else:
            L *= (rad2 - rad1) / PI2
        return Meter(arc=L)

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
        '''Get the area of this ellipse (C{meter**2}, conventionally).
        '''
        _, _, ab = self._ab3
        return Meter2(area=ab * PI)

    @Property_RO
    def b(self):
        '''Get semi-axis C{B{b}} of this ellipse (C{meter}, conventionally).
        '''
        _, b, _ = self._ab3
        return Meter(b=b)

    @Property_RO
    def c(self):
        '''Get the C{linear excentricity B{c}}, I{unsigned} (C{meter}, conventionally).
        '''
        return Meter(c=fabs(self.foci))

    @Property_RO
    def e(self):
        '''Get the excentricity (C{scalar, 0 <= B{e} <= 1}).
        '''
        e = self.e2
        return Scalar(e=sqrt(e) if 0 < e < 1 else e)

    @Property_RO
    def e2(self):
        '''Get the excentricity I{squared} (C{scalar, 0 <= B{e2} <= 1}).
        '''
        # C{e2} is aka C{k}, Elliptic C{k2} and SciPy's C{m}
        _, _, a, b = self._Pab4
        return Scalar(e2=((_1_0 - (b / a)**2) if 0 < b else _1_0) if b < a else _0_0)

    @Property_RO
    def _Ek(self):
        '''(INTERNAL) Get the C{Elliptic(k)} instance.
        '''
        return _MODS.elliptic._Ek(self.e2)

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

    def _Error(self, where, **cause):  # PYCHOK no cover
        '''(INTERNAL) Build an L{EllipseError}.
        '''
        t = self.named3
        u = unstr(t, a=self.a, b=self.b)
        if where:
            t =  typename(where, where)
            u = _DOT_(u, t)
        return EllipseError(u, **cause)

    @Property_RO
    def foci(self):
        '''Get the U{linear excentricity<https://WikiPedia.org/wiki/Ellipse#Standard_equation>},
           I{signed} (C{meter}, conventionally), C{positive} if this ellipse is oblate, C{negative}
           if prolate or C{0} if circular.  See also property L{Ellipse.c}.
        '''
        c = float(self.e)
        if c:
            r, _, a, _ = self._Pab4
            c  *= a
            if r:  # prolate
                c = -c
        return Meter(foci=c)  # signed

    @property_ROnce
    def _GKs(self):
        '''(INTERNAL) Compute the coefficients for property C{.perimeterGK}, I{once}.
        '''
        # U{numerators<https://OEIS.org/A056981>}, U{denominators<https://OEIS.org/A056982>}
        return (1, 1 / 4, 1 / 64, 1 / 256, 25 / 16384, 49 / 65536,
                441 / 1048576, 1089 / 4194304)  # overwrite property_ROnce

    def hartzell4(self, x, y, los=False):
        '''Compute the intersection of this ellipse with a Line-Of-Sight from Point-Of-View
           C{(B{x}, B{y})} I{outside} this ellipse.

           @kwarg los: Line-Of-Sight, I{direction} to the ellipse (L{Los}, L{Vector3d},
                       L{Vector2Tuple} or 2-tuple C{(dx, dy)}) or C{True} for the I{normal,
                       perpendicular, plumb} to this ellipse or C{False} or C{None} to
                       point to its center.

           @return: L{Vector4Tuple}C{(x, y, z, h)} with coordinates C{x}, C{y} and C{z=0}
                    of the intersection and C{h} the distance to "Point-Of-View" C{(B{x},
                    B{y})} I{along the} B{C{los}}, all in C{meter}, conventionally.

           @raise EllipseError: Invalid B{C{x}}, B{C{y}} or B{C{los}} or B{C{los}} points
                                outside or away from this ellipse.

           @see: Function L{hartzell4<triaxials.triaxial5.hartzell4>} for further details.
        '''
        V3d = _MODS.vector3d.Vector3d
        if los not in (True, False, None):
            try:
                los = V3d(los.x, los.y, 0)
            except (AttributeError, TypeError):
                if _MODS.basics.islistuple(los, minum=2):
                    los = V3d(*map(float, los[:2]))
        return self._triaxialX(self.hartzell4, V3d(x, y, 0), los=los)

    def height4(self, x, y, **normal_eps):
        '''Compute the projection on and distance to this ellipse from a point C{(B{x}, B{y})}
           in- or outside this ellipse.

           @kwarg normal_eps: With default C{B{normal}=True} the projection is I{perpendicular,
                         plumb} to this ellipse, otherwise C{radially} to its center (C{bool}).
                         Tolerance C{B{eps}=EPS} for root finding and validation (C{scalar}),
                         use a negative value to skip validation.

           @return: L{Vector4Tuple}C{(x, y, z, h)} with coordinates C{x}, C{y} and C{z=0}
                    of the projection on or the intersection with the ellipse and C{h} the
                    I{signed, normal distance} to the ellipse in C{meter}, conventionally.
                    Positive C{h} indicates, C{x} and/or C{y} are outside the ellipse,
                    negative C{h} means inside.

           @raise EllipseError: Invalid B{C{x}}, B{C{y}} or B{C{eps}}, no convergence in
                  root finding or validation failed.

           @see: Methods L{Ellipse.normal3d}, L{Ellipse.normal4} and function L{height4
                 <triaxials.triaxial5.height4>}.
        '''
        return self._triaxialX(self.height4, x, y, 0, **normal_eps)

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
        '''Is this ellipse oblate (foci on semi-axis C{a})? (C{bool})
        '''
        return self.a > self.b

    @property_RO
    def isProlate(self):
        '''Is this ellipse prolate (foci on semi-axis C{b})? (C{bool})
        '''
        return self.a < self.b

    @Property_RO
    def lati(self):
        '''Get the U{semi-latus rectum<https://WikiPedia.org/wiki/Ellipse#Standard_equation>},
           I{signed} (C{meter}, conventionally), C{positive} if this ellipse is oblate or
           circular, C{0} if "flat" and oblate, C{negative} if prolate or C{NEG0} if "flat"
           and prolate.  See also property L{Ellipse.p}.
        '''
        r, _, a, p = self._Pab4
        if 0 < p < a:
            p *= p / a
        if r:
            p = -p if p else NEG0
        return Meter(lati=p)  # signed

    def normal3d(self, deg_x, y=None, **length):
        '''Get a 3-D vector I{perpendicular to} this ellipse from point C{(B{x}, B{y})}
           C{on} this ellipse or at C{B{deg} degrees} along this ellipse.

           @kwarg length: Optional, signed C{B{length}=1} in out-/inward direction
                          (C{scalar}).

           @return: A C{Vector3d(x_, y_, z_=0)} normalized to B{C{length}}, pointing
                    out- or inward for postive respectively negative B{C{length}}.

           @raise EllipseError: Invalid B{C{x}} and/or B{C{y}}.

           @see: Methods L{Ellipse.height4}, L{Ellipse.normal4}, L{Ellipse.sideOf} and
                 C{Triaxial_.normal3d}.
        '''
        return self._triaxialX(self.normal3d, *self._xy03(deg_x, y), **length)

    def normal4(self, deg_x, y=None, **height_normal):
        '''Compute a point at B{C{height}} above or below this ellipse point C{(B{x},
           B{y})} C{on} this ellipse or at C{B{deg} degrees} along this ellipse.

           @kwarg height_normal: The desired distance C{B{height}=0} in- or outside this
                         ellipse (C{meter}, conventionally) and C{B{normal}=True},  If
                         C{B{normal}=True}, the B{C{height}} is I{perpendicular, plumb}
                         to this ellipse, otherwise C{radially} to its center (C{bool}).

           @return: L{Vector4Tuple}C{(x, y, z, h)} with coordinates C{x}, C{y} and C{z=0}
                    and C{h} the I{signed, normal distance} to the ellipse in C{meter},
                    conventionally.  Positive C{h} indicates, C{x} and/or C{y} are outside
                    the ellipse, negative C{h} means inside.

           @raise EllipseError: Invalid B{C{x}} and/or B{C{y}}.

           @see: Methods L{Ellipse.height4}, L{Ellipse.normal3d}, L{Ellipse.sideOf} and
                 C{Triaxial_.normal4}.
        '''
        return self._triaxialX(self.normal4, *self._xy03(deg_x, y), **height_normal)

    @Property_RO
    def p(self):
        '''Get the C{semi-latus rectum B{p} (aka B{𝓁}, script-small-l)}, I{unsigned}
           (C{meter}, conventionally).
        '''
        return Meter(p=fabs(self.lati))

    @Property_RO
    def perimeterAGM(self):
        '''Compute the perimeter of this ellipse using the U{Arithmetic-Geometric Mean
           <https://PaulBourke.net/geometry/ellipsecirc>} formula (C{meter}, conventionally).
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
        return Meter(perimeterAGM=P)

    @Property_RO
    def perimeter4Arc3(self):
        '''Compute the perimeter (and arcs) of this ellipse using the U{4-Arc
           <https://PaulBourke.net/geometry/ellipsecirc>} approximation as a
           3-Tuple C{(P, Ra, Rb)} with perimeter C{P}, arc radii C{Ra} and C{Rb}
           at the respective semi-axes (all in C{meter}, conventionally).
        '''
        r, P, a, b = self._Pab4
        if P is None:
            h = hypot(a, b)
            t = atan2(b, a)
            s, c = sincos2(t)
            L = (h - (a - b)) * _0_5
            a = _over(L, c)
            b = _over(h - L, s)
            P = (t * b + (PI_2 - t) * a) * _4_0
        elif a > b:  # flat
            a, b = _0_0, _1_over(b)  # INF
#       else:  # circular
#           pass
        if r:
            a, b = b, a
        return Meter(perimeter4Arc=P), Radius(Ra=a), Radius(Rb=b)

#   @Property_RO
#   def perimeterCR(self):
#       '''Compute the perimeter of this ellipse using U{Rackauckas'
#          <https://www.ChrisRackauckas.com/assets/Papers/ChrisRackauckas-The_Circumference_of_an_Ellipse.pdf>}
#          approximation, also U{here<https://ExtremeLearning.com.AU/a-formula-for-the-perimeter-of-an-ellipse>}
#          (C{meter}, conventionally).
#       '''
#       _, P, a, b = self._Pab4
#       if P is None:
#           P  =   a + b
#           h  = ((a - b) / P)**2
#           P *= (fhorner(h, 135168, -85760,  -5568, 3867) /
#                 fhorner(h, 135168, -119552, 22208, 345)) * PI
#       return Meter(perimeterCR=P)

    @Property_RO
    def perimeterGK(self):
        '''Compute the perimeter of this ellipse using the U{Gauss-Kummer
           <https://www.JohnDCook.com/blog/2023/05/28/approximate-ellipse-perimeter>}
           series, C{B{b / a} > 0.75} (C{meter}, conventionally).
        '''
        P, h = self._Ph2
        if h:
            P *= fhorner(h**2, *self._GKs) * PI
        return Meter(perimeterGK=P)

    @Property_RO
    def perimeter2k(self):
        '''Compute the perimeter of this ellipse using the complete integral
           of the 2nd kind, C{Elliptic.cE} (C{meter}, conventionally).
        '''
        return self._perimeter2k(self._ellipE)

    @Property_RO
    def perimeter2k_(self):
        '''Compute the perimeter of this ellipse using U{SciPy's ellipe
           <https://www.JohnDCook.com/perimeter_ellipse.html>} function
           if available, otherwise use property C{perimeter2k} (C{meter},
           conventionally).
        '''
        return self._perimeter2k(self._ellipe or self._ellipE)

    def _perimeter2k(self, _ellip):
        '''(INTERNAL) Helper for methods C{.PE2k} and C{.Pe2k}.
        '''
        _, P, a, _ = self._Pab4
        if P is None:  # see .ellipsoids.Ellipsoid.L
            k =  self.e2
            P = _ellip(k) * a * _4_0
        return Meter(perimeter2k=P)

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
           series (C{meter}, conventionally).
        '''
        P, h = self._Ph2
        if h:
            hs =  self._HGKs(h, self._maxit)
            P *= _fsum(hs) * PI  # nonfinites=True
        return Meter(perimeterHGK=P)

#   @Property_RO
#   def perimeterLS(self):
#       '''Compute the perimeter of this ellipse using the U{Linderholm-Segal
#          <https://www.JohnDCook.com/blog/2021/03/24/perimeter-of-an-ellipse>}
#          formula, aka C{3/2 norm} (C{meter}, conventionally).
#       '''
#       _, P, a, b = self._Pab4
#       if P is None:
#           n = pow(a, _1_5) + pow(b, _1_5)
#           P = pow(n * _0_5, _2_3rd) * PI2
#       return Meter(perimeterLS=P)

    @Property_RO
    def perimeter2R(self):
        '''Compute the perimeter of this ellipse using U{Ramanujan's 2nd
           <https://PaulBourke.net/geometry/ellipsecirc>} approximation,
           C{B{b / a} > 0.9} (C{meter}, conventionally).
        '''
        P, h = self._Ph2
        if h:
            h *= _3_0 * h
            h /=  sqrt(_4_0 - h) + _10_0  # /= chokes PyChecker?
            P *= (h + _1_0) * PI
        return Meter(perimeter2R=P)

    @Property_RO
    def R2(self):
        '''Get the I{authalic} radius of this ellipse, C{sqrt(B{a} * B{b})}
           (C{meter}, conventionally).
        '''
        a, b, ab = self._ab3
        return Radius(R2=sqrt(ab) if a != b else float(a))

    def Roc(self, deg_x, y=None, **eps):
        '''Compute the U{radius of curvature<https://WikiPedia.org/wiki/Radius_of_curvature>} at
           point C{(B{x}, B{y})} I{on} this ellipse or at C{B{deg} degrees} along this ellipse.

           @see: Method L{Roc_<Ellipse.Roc_>} for ruther details.
        '''
        x = radians(deg_x) if y is None else deg_x
        return self.Roc_(x, y, **eps)

    def Roc_(self, rad_x, y=None, **eps):
        '''Compute the U{radius of curvature<https://WikiPedia.org/wiki/Radius_of_curvature>} at
           point C{(B{x}, B{y})} I{on} this ellipse or at C{B{rad} radians} along this ellipse.

           @kwarg eps: See method C{sideOf}, use C{B{eps}=0} to permit any points.

           @return: Curvature (C{meter}, conventionally).

           @raise ValueError: Point C{(B{x}, B{y})} not near this ellipse, unless C{B{eps}=0}.
        '''
        try:
            a, b, ab = self._ab3
            if b != a:
                s, c, _ = self._scr3(rad_x, y, **eps)
                r = _over(hypot(s * a, c * b)**3, ab)
            else:  # circular
                r = float(a)
        except Exception as x:
            raise self._Error(self.Roc_, cause=x)
        return Radius(Roc=r)

    def _scr3(self, rad_x, y, eps=EPS):
        '''(INTERNAL) Helper for methods C{.Roc_} and C{.slope_}.
        '''
        r = float(rad_x)
        if y is not None:
            x, y = r, float(y)
            if eps and eps > 0:
                s = self._sideOf(x, y, eps)
                if s:
                    raise _ValueError(x=x, y=y, eps=eps, sideOf=s)
            r = atan2(y, x)
        s, c = sincos2(r)
        return s, c, r

    def _sideOf(self, x, y, eps):
        '''(INTERNAL) Helper for methods C{._scr3} and C{.sideOf}.
        '''
        a, b, ab = self._ab3
        s = ab or max(a, b)
        if s:
            s = (hypot(x * b, y * a) - s) / s
#           s = max(_N_1_0, min(_1_0, s))
        else:  # dot
            s = _1_0 if x or y else _0_0
        return INT0 if fabs(s) < eps else s

    def sideOf(self, x, y, eps=EPS):
        '''Return a C{positive}, C{negative} or C{0} fraction if point C{(B{x}, B{y})}
           is C{outside}, C{inside} respectively C{on} this ellipse.
        '''
        try:
            return Scalar(sideOf=self._sideOf(x, y, eps))
        except Exception as X:
            raise self._Error(self.sideOf, x=x, y=y, cause=X)

    def slope(self, deg_x, y=None, **eps):
        '''Compute the slope of the tangent at point C{(B{x}, B{y})} I{on} this ellipse
           or at C{B{deg} degrees} along this ellipse.

           @see: Method L{slope_<Ellipse.slope_>} for ruther details.
        '''
        x = radians(deg_x) if y is None else deg_x
        return Degrees(slope=degrees(self.slope_(x, y, **eps)))

    def slope_(self, rad_x, y=None, **eps):
        '''Compute the slope of the tangent at point C{(B{x}, B{y})} I{on} this ellipse
           or at C{B{rad} radians} along this ellipse.

           @kwarg eps: See method C{sideOf}, use C{B{eps}=0} to permit any points.

           @return: Slope (C{radians}), negative for C{0 <= B{rad} < PI/2}.

           @raise ValueError: C{(B{x}, B{y})} not near this ellipse, unless C{B{eps}=0}.
        '''
        try:
            s, c, _ = self._scr3(rad_x, y, **eps)
            a, b, _ = self._ab3
            r = atan2(-b**2 * c, a**2 * s)
            if r < 0:
                r += PI2
            if r >= PI3_2:
                r -= PI2
        except Exception as x:
            raise self._Error(self.slope_, cause=x)
        return Radians(slope=r or _0_0)  # no -0.0

    def toEllipsoid(self):
        '''Return an L{Ellipsoid<pygeodesy.Ellipsoid>} from this ellipse'
           C{a} and C{b} semi-axes.
        '''
        return _MODS.ellipsoids.Ellipsoid(self.a, b=self.b, name=self.name)

    def toStr(self, prec=8, terse=2, **sep_name):  # PYCHOK signature
        '''Return this ellipse as a text string.

           @kwarg prec: Number of decimal digits, unstripped (C{int}).
           @kwarg terse: Limit the number of items (C{int}, 0...9),
                         use C{B{terse}=0} or C{=None} for all.
           @kwarg sep_name: Optional C{B{name}=NN} (C{str}) or C{None}
                      to exclude this ellipse' name and separator
                      C{B{sep}=", "} to join the items (C{str}).

           @return: This C{Ellipse}' attributes (C{str}).
        '''
        E = Ellipse
        t = E.a, E.b
        if (terse or 0) != 2:
            t += E.c, E.e, E.e2, E.p, E.area, E.perimeter2k, E.R2
            if terse:
                t = t[:terse]
        return self._instr(prec=prec, props=t, **sep_name)

    def toTriaxial_(self, c=EPS):
        '''Return a L{Triaxial_<pygeodesy.Triaxial_>} from this ellipse'
           C{a} and C{b} semi-axes with B{C{c}} as minor semi-axis.
        '''
        return _MODS.triaxials.Triaxial_(self.a, b=self.b, c=c, name=self.name)  # 'NN'

    def _triaxialX(self, method, *args, **kwds):
        '''(INTERNAL) Map triaxial exceptions to L{EllipseError}s.
        '''
        try:
            _m = getattr(self.toTriaxial_(), method.__name__)
            return _m(*args, **kwds)
        except Exception as x:
            raise self._Error(method, cause=x)

    def _xy03(self, deg_x, y):
        if y is None:
            y, x = sincos2d(deg_x)
            y *= self.b
            x *= self.a
        else:
            x  = float(deg_x)
            y  = float(y)
        return x, y, 0


class EllipseError(_ValueError):
    '''Raised for any L{Ellipse} or C{ellipses} issue.
    '''
    pass  # ...


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
