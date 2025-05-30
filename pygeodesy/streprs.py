
# -*- coding: utf-8 -*-

u'''Floating point and other formatting utilities.
'''

from pygeodesy.basics import isint, islistuple, isscalar, isstr, itemsorted, \
                            _zip, _0_0,  typename
# from pygeodesy.constants import _0_0
from pygeodesy.errors import _or, _IsnotError, _TypeError, _ValueError, \
                             _xkwds_get, _xkwds_item2
# from pygeodesy.internals import typename  # from .basics
from pygeodesy.interns import NN, _0_, _0to9_, MISSING, _BAR_, _COMMASPACE_, \
                             _DOT_, _E_, _ELLIPSIS_, _EQUAL_, _H_, _LR_PAIRS, \
                             _N_, _name_, _not_scalar_, _PERCENT_, _SPACE_, \
                             _STAR_, _UNDER_
from pygeodesy.interns import _convergence_, _distant_, _e_, _exceeds_, \
                              _EQUALSPACED_, _f_, _F_, _g_, _limit_, \
                              _no_, _tolerance_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS

from math import fabs, log10 as _log10

__all__ = _ALL_LAZY.streprs
__version__ = '25.04.14'

_at_        = 'at'         # PYCHOK used!
_EN_PREC    =  6           # max MGRS/OSGR precision, 1 micrometer
_EN_WIDE    =  5           # number of MGRS/OSGR units, log10(_100km)
_eps_       = 'eps'        # PYCHOK used!
_OKd_       = '._-'        # acceptable name characters
_PAREN_g    = '(%g)'       # PYCHOK used!
_RESIDUAL_  = 'RESIDUAL'   # PYCHOK used!
_threshold_ = 'threshold'  # PYCHOK used!


class _Fmt(str):
    '''(INTERNAL) Callable formatting.
    '''
    name = NN

    def __call__(self, *name_value_, **name_value):
        '''Format a C{name=value} pair or C{name, value} pair or
           just a single C{value}.
        '''
        for n, v in name_value.items():
            break
        else:
            n, v = name_value_[:2] if len(name_value_) > 1 else \
                  (NN, (name_value_ or MISSING))
        t = str.__mod__(self, v)
        return NN(n, t) if n else t

#   def __mod__(self, arg, **unused):
#       '''Regular C{%} operator.
#       '''
#       return str.__mod__(self, arg)


class Fstr(str):
    '''(INTERNAL) C{float} format.
    '''
    name = NN

    def __call__(self, flt, prec=None, ints=False):
        '''Format the B{C{flt}} like function L{fstr}.
        '''
        # see also function C{fstr} if isscalar case below
        t = str.__mod__(_pct(self), flt) if prec is None else next(
           _streprs(prec, (flt,), self, ints, True, None))
        return t

    def __mod__(self, arg, **unused):
        '''Regular C{%} operator.

           @arg arg: A C{scalar} value to be formatted (either
                     the C{scalar}, or a 1-tuple C{(scalar,)},
                     or 2-tuple C{(prec, scalar)}.

           @raise TypeError: Non-scalar B{C{arg}} value.

           @raise ValueError: Invalid B{C{arg}}.
        '''
        def _error(arg):
            n = _DOT_(typename(Fstr), self.name or self)
            return _SPACE_(n, _PERCENT_, repr(arg))

        prec = 6  # default std %f and %F
        if islistuple(arg):
            n = len(arg)
            if n == 1:
                arg = arg[0]
            elif n == 2:
                prec, arg = arg
            else:
                raise _ValueError(_error(arg))

        if not isscalar(arg):
            raise _TypeError(_error(arg))
        return self(arg, prec=prec)  # Fstr.__call__(self, arg, prec=prec)


class _Sub(str):
    '''(INTERNAL) Class list formatter.
    '''
    # see .ellipsoidalNvector.LatLon.deltaTo
    def __call__(self, *Classes):
        t = _or(*map(typename, Classes))
        return  str.__mod__(self, t or MISSING)


class Fmt(object):
    '''Formatting options.
    '''
    ANGLE         = _Fmt('<%s>')
    COLON         = _Fmt(':%s')
#   COLONSPACE    = _Fmt(': %s')  # == _COLONSPACE_(n, v)
#   COMMASPACE    = _Fmt(', %s')  # == _COMMASPACE_(n, v)
    convergence   = _Fmt(_convergence_(_PAREN_g))
    CURLY         = _Fmt('{%s}')  # BRACES
    distant       = _Fmt(_distant_('(%.3g)'))
    DOT           = _Fmt('.%s')   # == NN(_DOT_, n)
    e             =  Fstr(_e_)
    E             =  Fstr(_E_)
    EQUAL         = _Fmt(_EQUAL_(NN, '%s'))
    EQUALg        = _Fmt(_EQUAL_(NN, '%g'))
    EQUALSPACED   = _Fmt(_EQUALSPACED_(NN, '%s'))
    exceeds_eps   = _Fmt(_exceeds_(_eps_, _PAREN_g))
    exceeds_limit = _Fmt(_exceeds_(_limit_, _PAREN_g))
    exceeds_R     = _Fmt(_exceeds_(_RESIDUAL_, _PAREN_g))
    f             =  Fstr(_f_)
    F             =  Fstr(_F_)
    g             =  Fstr(_g_)
    G             =  Fstr('G')
    h             =  Fstr('%+.*f')  # height, .streprs.hstr
    limit         = _Fmt(' %s limit')  # .units
    LOPEN         = _Fmt('(%s]')  # left-open range (L, R]
    PAREN         = _Fmt('(%s)')
    PAREN_g       = _Fmt(_PAREN_g)
    PARENSPACED   = _Fmt(' (%s)')
    QUOTE2        = _Fmt('"%s"')
    ROPEN         = _Fmt('[%s)')  # right-open range [L, R)
#   SPACE         = _Fmt(' %s')   # == _SPACE_(n, v)
    SQUARE        = _Fmt('[%s]')  # BRACKETS
    TAG           =  ANGLE
    TAGEND        = _Fmt('</%s>')
    tolerance     = _Fmt(_tolerance_(_PAREN_g))
    zone          = _Fmt('%02d')  # .epsg, .mgrs, .utmupsBase

    def __init__(self):
        for n, a in type(self).__dict__.items():
            if isinstance(a, (Fstr, _Fmt)):
                setattr(a, _name_, n)
        self.__name__ = typename(type(self))

    def __call__(self, obj, prec=9):
        '''Return C{str(B{obj})} or C{repr(B{obj})}.
        '''
        return str(obj) if isint(obj) else next(
              _streprs(prec, (obj,), Fmt.g, False, False, repr))

    def INDEX(self, name=NN, i=None, **name_i):
        '''Return C{"B{name}" if B{i} is None else "B{name}[B{i}]"}.
        '''
        if name_i:
            name, i = _xkwds_item2(name_i)
        return name if i is None else self.SQUARE(name, i)

    def no_convergence(self, _d, *tol, **thresh):
        '''Return C{"no convergence (B{_d})"}, C{"no convergence
           (B{_d}), tolerance (B{tol})"} or C{"no convergence
           (B{_d}), threshold (B{tol})"}.
        '''
        t = Fmt.convergence(fabs(_d))
        if tol:
            t = _COMMASPACE_(t, Fmt.tolerance(tol[0]))
            if thresh and _xkwds_get(thresh, thresh=False):
                t = t.replace(_tolerance_, _threshold_)
        return _no_(t)

    def repr_at(self, inst, text=NN):
        '''Return a C{repr} string C{"<B{text} at B{hex_id}>"}.
        '''
        return self.ANGLE(_SPACE_((text or inst), _at_, hex(id(inst))))

Fmt = Fmt()  # PYCHOK singleton

_DOTSTAR_ = Fmt.DOT(_STAR_)  # in _streprs above
# formats %G and %g drop all trailing zeros and the
# decimal point, making the float appear as an int
_Gg       = (Fmt.G, Fmt.g)
_FfEeGg   = (Fmt.F, Fmt.f, Fmt.E, Fmt.e) + _Gg  # float formats
_Fspec_   = NN('[%[<flags>][<width>]', _DOTSTAR_, ']', _BAR_.join(_FfEeGg))  # in testStreprs

del _convergence_, _distant_, _e_, _eps_, _exceeds_, _EQUALSPACED_,\
    _f_, _F_, _g_, _limit_, _PAREN_g, _RESIDUAL_


def anstr(name, OKd=_OKd_, sub=_UNDER_):
    '''Make a valid name of alphanumeric and OKd characters.

       @arg name: The original name (C{str}).
       @kwarg OKd: Other acceptable characters (C{str}).
       @kwarg sub: Substitute for invalid charactes (C{str}).

       @return: The modified name (C{str}).

       @note: Leading and trailing whitespace characters are removed,
              intermediate whitespace characters are coalesced and
              substituted.
    '''
    s = n = str(name).strip()
    for c in n:
        if not (c.isalnum() or c in OKd or c in sub):
            s = s.replace(c, _SPACE_)
    return sub.join(s.strip().split())


def attrs(inst, *names, **Nones_True__pairs_kwds):  # prec=6, fmt=Fmt.F, ints=False, Nones=True, sep=_EQUAL_
    '''Get instance attributes as I{name=value} strings, with C{float}s
       formatted by function L{fstr}.

       @arg inst: The instance (any C{type}).
       @arg names: The attribute names, all other positional (C{str}).
       @kwarg Nones_True__pairs_kwds: Keyword argument for function L{pairs}, except
              C{B{Nones}=True} to in-/exclude missing or C{None}-valued attributes.

       @return: A C{tuple(B{sep}.join(t) for t in zip(B{names}, reprs(values, ...)))}
                of C{str}s.
    '''
    def _items(inst, names, Nones):
        for n in names:
            v = getattr(inst, n, None)
            if Nones or v is not None:
                yield n, v

    def _Nones_kwds(Nones=True, **kwds):
        return Nones, kwds

    Nones, kwds = _Nones_kwds(**Nones_True__pairs_kwds)
    return pairs(_items(inst, names, Nones), **kwds)


def enstr2(easting, northing, prec, *extras, **wide_dot):
    '''Return an MGRS/OSGR easting, northing string representations.

       @arg easting: Easting from false easting (C{meter}).
       @arg northing: Northing from from false northing (C{meter}).
       @arg prec: Precision, the number of I{decimal} digits (C{int}) or if
                  negative, the number of I{units to drop}, like MGRS U{PRECISION
                  <https://GeographicLib.SourceForge.io/C++/doc/GeoConvert.1.html#PRECISION>}.
       @arg extras: Optional leading items (C{str}s).
       @kwarg wide_dot: Optional keword argument C{B{wide}=%d} for the number of I{unit digits}
                        (C{int}) and C{B{dot}=False} (C{bool}) to insert a decimal point.

       @return: B{C{extras}} + 2-tuple C{(str(B{easting}), str(B{northing}))} or
                + 2-tuple C{("", "")} for C{B{prec} <= -B{wide}}.

       @raise ValueError: Invalid B{C{easting}}, B{C{northing}} or B{C{prec}}.

       @note: The B{C{easting}} and B{C{northing}} values are I{truncated, not rounded}.
    '''
    t = extras
    try:  # like .dms.compassPoint
        p = min(int(prec), _EN_PREC)
        w = p + _xkwds_get(wide_dot, wide=_EN_WIDE)
        if w > 0:
            f  =  10**p  # truncate
            d  = (-p) if p > 0 and _xkwds_get(wide_dot, dot=False) else 0
            t += (_0wdot(w, int(easting  * f), d),
                  _0wdot(w, int(northing * f), d))
        else:  # prec <= -_EN_WIDE
            t += (NN, NN)
    except (TypeError, ValueError) as x:
        raise _ValueError(easting=easting, northing=northing, prec=prec, cause=x)
    return t

if enstr2.__doc__:  # PYCHOK expected
    enstr2.__doc__  %= (_EN_WIDE,)


def _enstr2m3(estr, nstr, wide=_EN_WIDE):  # in .mgrs, .osgr
    '''(INTERNAL) Convert east- and northing C{str}s to meter and resolution.
    '''
    def _s2m2(s, m):  # e or n str to float meter
        if _DOT_ in s:
            m  = 1  # meter
        else:
            s += _0_ * wide
            s  = _DOT_(s[:wide], s[wide:wide+_EN_PREC])
        return float(s), m

    e, m = _s2m2(estr, 0)
    n, m = _s2m2(nstr, m)
    if not m:
        p = max(len(estr), len(nstr))  # 2 = Km, 5 = m, 7 = cm
        m = 10**max(-_EN_PREC, wide - p)  # resolution, meter
    return e, n, m


def fstr(floats, prec=6, fmt=Fmt.F, ints=False, sep=_COMMASPACE_, strepr=None, force=True):
    '''Convert one or more floats to string, optionally stripped of trailing zero decimals.

       @arg floats: Single or a list, sequence, tuple, etc. (C{scalar}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.  In
                    addition, trailing decimal zeros are stripped for U{alternate,
                    form '#'<https://docs.Python.org/3/library/stdtypes.html
                    #printf-style-string-formatting>}.
       @kwarg fmt: Optional C{float} format (C{letter}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).
       @kwarg sep: Separator joining the B{C{floats}} (C{str}).
       @kwarg strepr: Optional callable to format non-C{floats} (typically
                      C{repr}, C{str}) or C{None} to raise a TypeError and used
                      only if C{B{force} is not True}.
       @kwarg force: If C{True}, format all B{C{floats}} using B{C{fmt}},
                     otherwise use B{C{strepr}} for non-C{floats} (C{bool}).

       @return: The C{sep.join(strs(floats, ...)} joined (C{str}) or single
                C{strs((floats,), ...)} (C{str}) if B{C{floats}} is C{scalar}.
    '''
    if isscalar(floats):  # see Fstr.__call__ above
        return next(_streprs(prec, (floats,), fmt, ints, force, strepr))
    else:
        return sep.join(_streprs(prec, floats, fmt, ints, force, strepr))


def _fstrENH2(inst, prec, m, fmt=Fmt.F):  # in .css, .lcc, .utmupsBase
    # (INTERNAL) For C{Css.} and C{Lcc.} C{toRepr} and C{toStr} and C{UtmUpsBase._toStr}.
    t = inst.easting, inst.northing
    t = tuple(_streprs(prec, t, fmt, False, True, None))
    T = _E_, _N_
    if m is not None and fabs(inst.height):  # fabs(self.height) > EPS
        t +=  hstr(inst.height, prec=-2, m=m),
        T += _H_,
    return t, T


def _fstrLL0(inst, prec, toRepr):  # in .azimuthal, .css
    # (INTERNAL) For C{_AlbersBase.}, C{_AzimuthalBase.} and C{CassiniSoldner.}
    t = tuple(_streprs(prec, inst.latlon0, Fmt.F, False, True, None))
    if toRepr:
        n = inst.name
        if n:
            t += Fmt.EQUAL(_name_, repr(n)),
        t = Fmt.PAREN(inst.classname, _COMMASPACE_.join(t))
    return t


def fstrzs(efstr, ap1z=False):
    '''Strip trailing zero decimals from a C{float} string.

       @arg efstr: Float with or without exponent (C{str}).
       @kwarg ap1z: Append the decimal point and one zero decimal
                    if the B{C{efstr}} is all digits (C{bool}).

       @return: Float (C{str}).
    '''
    s = efstr.find(_DOT_)
    if s >= 0:
        e = efstr.rfind(Fmt.e)
        if e < 0:
            e = efstr.rfind(Fmt.E)
            if e < 0:
                e = len(efstr)
        s += 2  # keep 1st _DOT_ + _0_
        if s < e and efstr[e-1] == _0_:
            efstr = NN(efstr[:s], efstr[s:e].rstrip(_0_), efstr[e:])

    elif ap1z:
        # %.G and %.g formats may drop the decimal
        # point and all trailing zeros, ...
        if efstr.isdigit():
            efstr += _DOT_ + _0_  # ... append or ...
        else:  # ... insert one dot and zero
            e = efstr.rfind(Fmt.e)
            if e < 0:
                e = efstr.rfind(Fmt.E)
            if e > 0:
                efstr = NN(efstr[:e], _DOT_, _0_, efstr[e:])

    return efstr


def hstr(height, prec=2, fmt=Fmt.h, ints=False, m=NN):
    '''Return a string for the height value.

       @arg height: Height value (C{float}).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional C{float} format (C{letter}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).
       @kwarg m: Optional unit of the height (C{str}).
    '''
    h = next(_streprs(prec, (height,), fmt, ints, True, None))
    return NN(h, str(m)) if m else h


def instr(inst, *args, **kwds):
    '''Return the string representation of an instantiation.

       @arg inst: The instance (any C{type}).
       @arg args: Optional positional arguments.
       @kwarg kwds: Optional keyword arguments.

       @return: Representation (C{str}).
    '''
    return unstr(_MODS.named.classname(inst), *args, **kwds)


def lrstrip(txt, lrpairs=_LR_PAIRS):
    '''Left- I{and} right-strip parentheses, brackets, etc. from a string.

       @arg txt: String to be stripped (C{str}).
       @kwarg lrpairs: Parentheses, etc. to remove (C{dict} of one or several
                       C{(Left, Right)} pairs).

       @return: Stripped B{C{txt}} (C{str}).
    '''
    _e, _s, _n = str.endswith, str.startswith, len
    while _n(txt) > 2:
        for L, R in lrpairs.items():
            if _e(txt, R) and _s(txt, L):
                txt = txt[_n(L):-_n(R)]
                break  # restart
        else:
            return txt


def pairs(items, prec=6, fmt=Fmt.F, ints=False, sep=_EQUAL_):
    '''Convert items to I{name=value} strings, with C{float}s handled like L{fstr}.

       @arg items: Name-value pairs (C{dict} or 2-{tuple}s of any C{type}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional C{float} format (C{letter}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).
       @kwarg sep: Separator joining I{names} and I{values} (C{str}).

       @return: A C{tuple(B{sep}.join(t) for t in B{items}))} of C{str}s.
    '''
    try:
        if isinstance(items, dict):
            items = itemsorted(items)
        elif not islistuple(items):
            items = tuple(items)
        # can't unzip empty items tuple, list, etc.
        n, v = _zip(*items) if items else ((), ())  # strict=True
    except (TypeError, ValueError):
        raise _IsnotError(dict, '2-tuples', items=items)
    v = _streprs(prec, v, fmt, ints, False, repr)
    return tuple(sep.join(t) for t in _zip(map(str, n), v))  # strict=True


def _pct(fmt):
    '''(INTERNAL) Prefix C{%} if needed.
    '''
    return fmt if _PERCENT_ in fmt else NN(_PERCENT_, fmt)


def reprs(objs, prec=6, fmt=Fmt.F, ints=False):
    '''Convert objects to C{repr} strings, with C{float}s handled like L{fstr}.

       @arg objs: List, sequence, tuple, etc. (any C{type}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional C{float} format (C{letter}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).

       @return: A C{tuple(map(fstr|repr, objs))} of C{str}s.
    '''
    return tuple(_streprs(prec, objs, fmt, ints, False, repr)) if objs else ()


def _resolution10(resolution, Error=ValueError):  # in .mgrs, .osgr
    '''(INTERNAL) Validate C{resolution} in C{meter}.
    '''
    try:
        r = int(_log10(resolution))
        if _EN_WIDE < r or r < -_EN_PREC:
            raise ValueError
    except (ValueError, TypeError):
        raise Error(resolution=resolution)
    return _MODS.units.Meter(resolution=10**r)


def _streprs(prec, objs, fmt, ints, force, strepr):
    '''(INTERNAL) Helper for C{fstr}, C{pairs}, C{reprs} and C{strs}
    '''
    # <https://docs.Python.org/3/library/stdtypes.html#printf-style-string-formatting>
    if fmt in _FfEeGg:
        fGg = fmt in _Gg
        fmt = NN(_PERCENT_, _DOT_, abs(prec), fmt)

    elif fmt.startswith(_PERCENT_):
        fGg = False
        try:  # to make sure fmt is valid
            f = fmt.replace(_DOTSTAR_, Fmt.DOT(abs(prec)))
            _ = f % (_0_0,)
        except (TypeError, ValueError):
            raise _ValueError(fmt=fmt, txt_not_=repr(_DOTSTAR_))
        fmt = f

    else:
        raise _ValueError(fmt=fmt, txt_not_=repr(_Fspec_))

    for i, o in enumerate(objs):
        if force or isinstance(o, float):
            t = fmt % (float(o),)
            if ints and t.rstrip(_0to9_ if isint(o, both=True) else
                                 _0_).endswith(_DOT_):
                t = t.split(_DOT_)[0]
            elif prec > 1:
                t = fstrzs(t, ap1z=fGg)
        elif strepr:
            t = strepr(o)
        else:
            t = Fmt.PARENSPACED(Fmt.SQUARE(objs=i), o)
            raise TypeError(_SPACE_(t, _not_scalar_))
        yield t


def strs(objs, prec=6, fmt=Fmt.F, ints=False):
    '''Convert objects to C{str} strings, with C{float}s handled like L{fstr}.

       @arg objs: List, sequence, tuple, etc. (any C{type}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional C{float} format (C{letter}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).

       @return: A C{tuple(map(fstr|str, objs))} of C{str}s.
    '''
    return tuple(_streprs(prec, objs, fmt, ints, False, str)) if objs else ()


def unstr(where, *args, **kwds_):
    '''Return the string representation of an invokation.

       @arg where: Class, function, method (C{type}) or name (C{str}).
       @arg args: Optional positional arguments.
       @kwarg kwds_: Optional keyword arguments, except C{B{_Cdot}=None},
                     C{B{_ELLIPSIS}=False} and C{B{_fmt}=Fmt.g}.

       @return: Representation (C{str}).
    '''
    def _C_e_g_kwds3(_Cdot=None, _ELLIPSIS=0, _fmt=Fmt.g, **kwds):
        return _Cdot, _ELLIPSIS, _fmt, kwds

    C, e, g, kwds = _C_e_g_kwds3(**kwds_)
    if e and len(args) > (e + 1):
        t  =  reprs(args[:e], fmt=g)
        t += _ELLIPSIS_,
        t +=  reprs(args[-1:], fmt=g)
    else:
        t  =  reprs(args, fmt=g) if args else ()
    if kwds:
        t += pairs(itemsorted(kwds), fmt=g)
    n = where if isstr(where) else typename(where)  # _NN_
    if C and hasattr(C, n):
        try:  # bound method of class C?
            where = where.__func__
        except AttributeError:
            pass  # method of C?
        if getattr(C, n, None) is where:
            n = _DOT_(typename(C), n)
    return Fmt.PAREN(n, _COMMASPACE_.join(t))


def _0wd(*w_i):  # in .osgr, .wgrs
    '''(INTERNAL) Int formatter'.
    '''
    return '%0*d' % w_i


def _0wdot(w, f, dot=0):
    '''(INTERNAL) Int and Float formatter'.
    '''
    s = _0wd(w, int(f))
    if dot:
        s = _DOT_(s[:dot], s[dot:])
    return s


def _0wpF(*w_p_f):  # in .dms, .osgr
    '''(INTERNAL) Float deg, min, sec formatter'.
    '''
    return '%0*.*f' % w_p_f  # XXX was F


def _xzipairs(names, values, sep=_COMMASPACE_, fmt=NN, pair_fmt=Fmt.COLON):
    '''(INTERNAL) Zip C{names} and C{values} into a C{str}, joined and bracketed.
    '''
    try:
        t = sep.join(pair_fmt(*t) for t in _zip(names, values))  # strict=True
    except Exception as x:
        raise _ValueError(names=names, values=values, cause=x)
    return (fmt % (t,)) if fmt else t  # enc


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
