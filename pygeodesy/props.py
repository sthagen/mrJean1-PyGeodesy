
# -*- coding: utf-8 -*-

u'''Mutable, immutable and caching/memoizing properties and
deprecation decorators.

To enable C{DeprecationWarning}s from C{PyGeodesy}, set env var
C{PYGEODESY_WARNINGS} to a non-empty string I{AND} run C{python}
with command line option C{-X dev} or with one of the C{-W}
choices, see callable L{DeprecationWarnings} below.
'''

from pygeodesy.basics import isclass as _isclass
from pygeodesy.errors import _AssertionError, _AttributeError, \
                             _xcallable, _xkwds_get
from pygeodesy.internals import _tailof, typename
from pygeodesy.interns import MISSING, NN, _an_, _COMMASPACE_, \
                             _DEPRECATED_, _DOT_, _EQUALSPACED_, \
                             _immutable_, _invalid_, _module_, \
                             _N_A_, _NL_, _not_, _SPACE_, _UNDER_
# from pygeodesy.named import callname  # _MODS, avoid circular
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, \
                             _FOR_DOCS, _WARNINGS_X_DEV
# from pygeodesy.streprs import Fmt  # _MODS

from functools import wraps as _wraps

__all__ = _ALL_LAZY.props
__version__ = '25.05.26'

_class_       = 'class'
_DNL_         = _NL_ * 2  # PYCHOK used!
_dont_use_    = _DEPRECATED_ + ", don't use."
_function_    = 'function'
_has_been_    = 'has been'  # PYCHOK used!
_method_      = 'method'
_not_an_inst_ = _not_(_an_, 'instance')


def _allPropertiesOf(Clas_or_inst, *Bases, **excls):
    '''(INTERNAL) Yield all C{R/property/_RO}s at C{Clas_or_inst}
       as specified in the C{Bases} arguments, except C{excls}.
    '''
    if _isclass(Clas_or_inst):
        S = Clas_or_inst,  # just this Clas
    else:  # class and super-classes of inst
        try:
            S = type(Clas_or_inst).__mro__[:-1]  # not object
        except AttributeError:
            raise
            S = ()  # not an inst
    B =  Bases or _PropertyBase
    P = _property_RO___
    for C in S:
        for n, p in C.__dict__.items():
            if isinstance(p, B) and p.name == n and not (
               isinstance(p, P) or n in excls):
                yield p


def _allPropertiesOf_n(n, Clas_or_inst, *Bases, **excls):
    '''(INTERNAL) Assert the number of C{R/property/_RO}s at C{Clas_or_inst}.
    '''
    t = tuple(p.name for p in _allPropertiesOf(Clas_or_inst, *Bases, **excls))
    if len(t) != n:
        raise _AssertionError(_COMMASPACE_.join(t), Clas_or_inst,
                          txt=_COMMASPACE_(len(t), _not_(n)))
    return t


def _hasProperty(inst, name, *Classes):  # in .named._NamedBase._update
    '''(INTERNAL) Check whether C{inst} has a C{P/property/_RO} by this C{name}.
    '''
    p = getattr(inst.__class__, name, None)  # walks __class__.__mro__
    return bool(p and isinstance(p, Classes or _PropertyBase)
                  and p.name == name)


# def _isclass(obj):
#     '''(INTERNAL) Get and overwrite C{_isclass}.
#     '''
#     _MODS.getmodule(__name__)._isclass = f = _MODS.basics.isclass
#     return f(obj)


def _update_all(inst, *attrs, **Base_needed):
    '''(INTERNAL) Zap all I{cached} L{property_RO}s, L{Property}s,
       L{Property_RO}s and the named C{attrs} of an instance.

       @return: The number of updates (C{int}), if any.
    '''
    if _isclass(inst):
        raise _AssertionError(inst, txt=_not_an_inst_)
    try:
        d = inst.__dict__
    except AttributeError:
        return 0
    u = len(d)
    if u >  _xkwds_get(Base_needed, needed=0):
        B = _xkwds_get(Base_needed, Base=_PropertyBase)
        for p in _allPropertiesOf(inst, B):
            p._update(inst)  # d.pop(p.name, None)

        if attrs:
            _update_attrs(inst, *attrs)  # remove attributes from inst.__dict__
        u -= len(d)
    return u  # updates


# def _update_all_from(inst, other, **Base):
#     '''(INTERNAL) Update all I{cached} L{Property}s and
#        L{Property_RO}s of instance C{inst} from C{other}.
#
#        @return: The number of updates (C{int}), if any.
#     '''
#     if _isclass(inst):
#         raise _AssertionError(inst, txt=_not_an_inst_)
#     try:
#         d = inst.__dict__
#         f = other.__dict__
#     except AttributeError:
#         return 0
#     u = len(f)
#     if u:
#         u = len(d)
#         B = _xkwds_get(Base, Base=_PropertyBase)
#         for p in _allPropertiesOf(inst, B):
#             p._update_from(inst, other)
#         u -= len(d)
#     return u  # number of updates


def _update_attrs(inst, *attrs):
    '''(INTERNAL) Zap all named C{attrs} of an instance.

       @return: The number of updates (C{int}), if any.
    '''
    try:
        d = inst.__dict__
    except AttributeError:
        return 0
    u = len(d)
    if u:  # zap attrs from inst.__dict__
        _p = d.pop
        for a in attrs:
            _ = _p(a, MISSING)
#           if _ is MISSING and not hasattr(inst, a):
#               n = _MODS.named.classname(inst, prefixed=True)
#               a = _DOT_(n, _SPACE_(a, _invalid_))
#               raise _AssertionError(a, txt=repr(inst))
#           _ = _p(a, None)  # redo: hasattr side effect
        u -= len(d)
        # assert u >= 0
    return u  # number of named C{attrs} zapped


class _PropertyBase(property):
    '''(INTERNAL) Base class for C{P/property/_RO}.
    '''
    def __init__(self, method, fget, fset, doc=NN):

        _xcallable(getter=method, fget=fget)

        self.method = method
        self.name   = typename(method)
        d = doc or method.__doc__
        if _FOR_DOCS and d:
            self.__doc__ = d   # PYCHOK no cover

        property.__init__(self, fget, fset, self._fdel, d or _N_A_)

    def _Error(self, kind, nameter, farg):
        '''(INTERNAL) Return an C{AttributeError} instance.
        '''
        if farg:
            n = _DOT_(self.name, typename(nameter))
            n = _SPACE_(n, typename(farg))
        else:
            n = nameter
        e = _SPACE_(kind, _MODS.named.classname(self))
        return _AttributeError(e, txt=n)

    def _fdel(self, inst):
        '''Zap the I{cached/memoized} C{property} value.
        '''
        self._update(inst, None)   # PYCHOK no cover

    def _fget(self, inst):
        '''Get and I{cache/memoize} the C{property} value.
        '''
        try:  # to get the value cached in instance' __dict__
            return inst.__dict__[self.name]
        except KeyError:
            # cache the value in the instance' __dict__
            inst.__dict__[self.name] = val = self.method(inst)
            return val

    def _fset_error(self, inst, val):
        '''Throws an C{AttributeError}, always.
        '''
        n = _MODS.named.classname(inst)
        n = _DOT_(n, self.name)
        n = _EQUALSPACED_(n, repr(val))
        raise self._Error(_immutable_, n, None)

    def _update(self, inst, *unused):
        '''(INTERNAL) Zap the I{cached/memoized} C{inst.__dict__[name]} item.
        '''
        inst.__dict__.pop(self.name, None)  # name, NOT _name

    def _update_from(self, inst, other):
        '''(INTERNAL) Copy a I{cached/memoized} C{inst.__dict__[name]} item
           from C{other.__dict__[name]} if present, otherwise zap it.
        '''
        n = self.name  # name, NOT _name
        v = other.__dict__.get(n, MISSING)
        if v is MISSING:
            inst.__dict__.pop(n, None)
        else:
            inst.__dict__[n] = v

    def deleter(self, fdel):
        '''Throws an C{AttributeError}, always.
        '''
        raise self._Error(_invalid_, self.deleter, fdel)

    def getter(self, fget):
        '''Throws an C{AttributeError}, always.
        '''
        raise self._Error(_invalid_, self.getter, fget)

    def setter(self, fset):
        '''Throws an C{AttributeError}, always.
        '''
        raise self._Error(_immutable_, self.setter, fset)


class Property_RO(_PropertyBase):
    # No __doc__ on purpose
    def __init__(self, method, doc=NN):  # PYCHOK expected
        '''New I{immutable}, I{caching}, I{memoizing} C{property} I{Factory}
           to be used as C{decorator}.

           @arg method: The callable being decorated as this C{property}'s C{getter},
                        to be invoked only once.
           @kwarg doc: Optional property documentation (C{str}).

           @note: Like standard Python C{property} without a C{setter}, but with
                  a more descriptive error message when set.

           @see: Python 3's U{functools.cached_property<https://docs.Python.org/3/
                 library/functools.html#functools.cached_property>} and U{-.cache
                 <https://Docs.Python.org/3/library/functools.html#functools.cache>}
                 to I{cache} or I{memoize} the property value.

           @see: Luciano Ramalho, "Fluent Python", O'Reilly, 2016 p. 636
                 Example 19-24 or 2022 p. 870 Example 22-28 and U{class
                 Property<https://docs.Python.org/3/howto/descriptor.html>}.
        '''
        _fget = method if _FOR_DOCS else self._fget  # XXX force method.__doc__ to epydoc
        _PropertyBase.__init__(self, method, _fget, self._fset_error)

    def __get__(self, inst, *unused):  # PYCHOK 2 vs 3 args
        if inst is None:
            return self
        try:  # to get the cached value immediately
            return inst.__dict__[self.name]
        except (AttributeError, KeyError):
            return self._fget(inst)


class Property(Property_RO):
    # No __doc__ on purpose
    __init__ = Property_RO.__init__
    '''New I{mutable}, I{caching}, I{memoizing} C{property} I{Factory}
       to be used as C{decorator}.

       @see: L{Property_RO} for more details.

       @note: Unless and until the C{setter} is defined, this L{Property} behaves
              like an I{immutable}, I{caching}, I{memoizing} L{Property_RO}.
    '''

    def setter(self, method):
        '''Make this C{Property} I{mutable}.

           @arg method: The callable being decorated as this C{Property}'s C{setter}.

           @note: Setting a new property value always clears the previously I{cached}
                  or I{memoized} value I{after} invoking the B{C{method}}.
        '''
        def _fset(inst, val):
            '''Set and I{cache}, I{memoize} the C{property} value.
            '''
            _ = method(inst, val)
            self._update(inst)  # un-cache this item

        return self._setters(method, _fset)

    def setter_(self, method):
        '''Make this C{Property} I{mutable}.

           @arg method: The callable being decorated as this C{Property}'s C{setter}
                        and returning the new property value to be I{cached} or
                        I{memoized}.
        '''
        def _fset(inst, val):
            '''Set and I{cache}, I{memoize} the C{property} value.
            '''
            val = method(inst, val)
            inst.__dict__[self.name] = val

        return self._setters(method, _fset)

    def _setters(self, method, _fset):
        _xcallable(setter=method, fset=_fset)
        if _FOR_DOCS:  # XXX force method.__doc__ into epydoc
            _PropertyBase.__init__(self, self.method, self.method, method)
        else:  # class Property <https://docs.Python.org/3/howto/descriptor.html>
            _PropertyBase.__init__(self, self.method, self._fget, _fset)
        return self


class property_RO(_PropertyBase):
    # No __doc__ on purpose
    _uname = NN

    def __init__(self, method, doc=NN):  # PYCHOK expected
        '''New I{immutable}, standard C{property} to be used as C{decorator}.

           @arg method: The callable being decorated as C{property}'s C{getter}.
           @kwarg doc: Optional property documentation (C{str}).

           @note: Like standard Python C{property} without a setter, but with
                  a more descriptive error message when set.

           @see: L{Property_RO}.
        '''
        _PropertyBase.__init__(self, method, method, self._fset_error, doc=doc)
        self._uname = NN(_UNDER_, self.name)  # actual attr UNDER<name>

    def _update(self, inst, *Clas):  # PYCHOK signature
        '''(INTERNAL) Zap the I{cached} C{B{inst}.__dict__[_name]} item.
        '''
        uname = self._uname
        if uname in inst.__dict__:
            if Clas:  # overrides inst.__class__
                d = Clas[0].__dict__.get(uname, MISSING)
            else:
                d = getattr(type(inst), uname, MISSING)
#               if d is MISSING:  # XXX superfluous
#                   for c in type(inst).__mro__[:-1]:
#                       if uname in c.__dict__:
#                           d = c.__dict__[uname]
#                           break
            if d is None:  # remove inst value
                inst.__dict__.pop(uname)


class _property_RO___(_PropertyBase):
    # No __doc__ on purpose

    def __init__(self, method, doc=NN):  # PYCHOK expected
        '''New C{property_ROnce} or C{property_ROver}, holding a singleton value as
           class attribute for all instances of that class and overwriting C{self},
           the C{property_ROver} instance in the 1st invokation.

           @see: L{property_RO} for further details.
        '''
        _PropertyBase.__init__(self, method, self._fget, self._fset_error, doc=doc)

    def _fdel(self, *unused):  # PYCHOK no cover
        '''Silently ignored, always.
        '''
        pass

    def _update(self, *unused):  # PYCHOK signature
        '''(INTERNAL) No-op, ignore updates.
        '''
        pass


class property_ROnce(_property_RO___):
    # No __doc__ on purpose

    def _fget(self, inst):
        '''Get the C{property} value, only I{once} and memoize/cache it.
        '''
        try:
            v = self._val
        except AttributeError:
            v = self._val = self.method(inst)
        return v


class property_ROver(_property_RO___):
    # No __doc__ on purpose

    def _fget(self, inst):
        '''Get the C{property} value I{once} and overwrite C{self},
           this C{property} instance in the (super-)class of C{self}
           where this property is define as a L{property_ROver}.
        '''
        v = self.method(inst)
        n = self.name
        C = type(inst)
        for c in C.__mro__:  # [:-1]
            if getattr(c, n, None) is self:
                setattr(c, n, v)  # overwrite property_ROver
                break
        else:
            n = _DOT_(typename(C), n)
            raise _AssertionError(_EQUALSPACED_(n, v))
        return v


class _NamedProperty(property):  # in .named
    '''Class C{property} with a C{.name} attribute.
    '''
    @Property_RO
    def name(self):
        '''Get the name of this C{property} (C{str}).
        '''
        return typename(self.fget)


def property_doc_(doc):
    '''Decorator for a standard C{property} with basic documentation.

       @arg doc: The property documentation (C{str}).

       @example:

        >>>class Clas(object):
        >>>
        >>>     @property_doc_("documentation text.")
        >>>     def name(self):
        >>>         ...
        >>>
        >>>     @name.setter
        >>>     def name(self, value):
        >>>         ...
    '''
    # See Luciano Ramalho, "Fluent Python", O'Reilly, 2016 p. 212+
    # Example 7-23 or 2022 p. 331+ Example 9-22 and <https://
    # Python-3-Patterns-Idioms-Test.ReadTheDocs.io/en/latest/PythonDecorators.html>

    def _documented_property(method):
        '''(INTERNAL) Return the documented C{property}.
        '''
        t = 'get and set' if doc.startswith(_SPACE_) else NN
        return _NamedProperty(method, None, None, NN('Property to ', t, doc))

    return _documented_property


def _deprecated(call, kind, qual_d):
    '''(INTERNAL) Decorator for DEPRECATED functions, methods, etc.

       @see: Brett Slatkin, "Effective Python", 2019 page 105, 2nd
             ed, Addison-Wesley.
    '''
    doc = _DEPRECATEDof(call)

    @_wraps(call)  # PYCHOK self?
    def _deprecated_call(*args, **kwds):
        if qual_d:  # function
            q =  qual_d
        elif args:  # method
            q = _qualified(args[0], typename(call))
        else:  # PYCHOK no cover
            q =  typename(call)
        _throwarning(kind, q, doc)
        return call(*args, **kwds)

    return _deprecated_call


def deprecated_class(cls_or_class):
    '''Use inside __new__ or __init__ of a DEPRECATED class.

       @arg cls_or_class: The class (C{cls} or C{Class}).

       @note: NOT a decorator!
    '''
    if _WARNINGS_X_DEV:
        q = _DOT_(cls_or_class.__module__, typename(cls_or_class))
        _throwarning(_class_, q, cls_or_class.__doc__)  # _DDOC_


def deprecated_function(call):
    '''Decorator for a DEPRECATED function.

       @arg call: The deprecated function (C{callable}).

       @return: The B{C{call}} DEPRECATED.
    '''
    return _deprecated(call, _function_, _DOT_(
                       call.__module__, typename(call))) if \
           _WARNINGS_X_DEV else call


def deprecated_method(call):
    '''Decorator for a DEPRECATED method.

       @arg call: The deprecated method (C{callable}).

       @return: The B{C{call}} DEPRECATED.
    '''
    return _deprecated(call, _method_, NN) if _WARNINGS_X_DEV else call


def _deprecated_module(name):  # PYCHOK no cover
    '''(INTERNAL) Callable within a DEPRECATED module.
    '''
    if _WARNINGS_X_DEV:
        _throwarning(_module_, name, _dont_use_)


def _DEPRECATEDof(obj):
    '''(INTERNAL) Get uniform DEPRECATED __doc__ string.
    '''
    try:
        d = obj.__doc__.strip()
        i = d.find(_DEPRECATED_)
    except AttributeError:
        i = -1
    return _DOT_(_DEPRECATED_, NN) if i < 0 else d[i:]


if _WARNINGS_X_DEV:
    class deprecated_property(_PropertyBase):
        '''Decorator for a DEPRECATED C{property} or C{Property}.
        '''
        def __init__(self, method):
            '''Decorator for a DEPRECATED C{property} or C{Property} getter.
            '''
            doc = _DEPRECATEDof(method)

            def _fget(inst):  # PYCHOK no cover
                '''Get the C{property} or C{Property} value.
                '''
                q = _qualified(inst, self.name)
                _throwarning(typename(property), q, doc)
                return self.method(inst)  # == method

            _PropertyBase.__init__(self, method, _fget, None, doc=doc)

        def setter(self, method):
            '''Decorator for a DEPRECATED C{property} or C{Property} setter.

               @arg method: The callable being decorated as this C{Property}'s C{setter}.

               @note: Setting a new property value always clears the previously I{cached}
                      or I{memoized} value I{after} invoking the B{C{method}}.
            '''
            if not callable(method):
                _PropertyBase.setter(self, method)  # PYCHOK no cover

            if _FOR_DOCS:  # XXX force method.__doc__ into epydoc
                _PropertyBase.__init__(self, self.method, self.method, method)
            else:

                def _fset(inst, val):
                    '''Set the C{property} or C{Property} value.
                    '''
                    q = _qualified(inst, self.name)
                    _throwarning(typename(property), q, _DEPRECATEDof(method))
                    method(inst, val)
                    # self._update(inst)  # un-cache this item

                # class Property <https://docs.Python.org/3/howto/descriptor.html>
                _PropertyBase.__init__(self, self.method, self._fget, _fset)
            return self

else:  # PYCHOK no cover
    class deprecated_property(property):  # PYCHOK expected
        '''Decorator for a DEPRECATED C{property} or C{Property}.
        '''
        pass

deprecated_Property = deprecated_property


def deprecated_Property_RO(method):
    '''Decorator for a DEPRECATED L{Property_RO}.

       @arg method: The C{Property_RO.fget} method (C{callable}).

       @return: The B{C{method}} DEPRECATED.
    '''
    return _deprecated_RO(method, Property_RO)


def deprecated_property_RO(method):
    '''Decorator for a DEPRECATED L{property_RO}.

       @arg method: The C{property_RO.fget} method (C{callable}).

       @return: The B{C{method}} DEPRECATED.
    '''
    return _deprecated_RO(method, property_RO)


def _deprecated_RO(method, _RO):
    '''(INTERNAL) Create a DEPRECATED C{property_RO} or C{Property_RO}.
    '''
    doc = _DEPRECATEDof(method)

    if _WARNINGS_X_DEV:

        class _Deprecated_RO(_PropertyBase):
            __doc__ = doc

            def __init__(self, method):
                _PropertyBase.__init__(self, method, self._fget, self._fset_error, doc=doc)

            def _fget(self, inst):  # PYCHOK no cover
                q = _qualified(inst, self.name)
                _throwarning(typename(_RO), q, doc)
                return self.method(inst)

        return _Deprecated_RO(method)
    else:  # PYCHOK no cover
        return _RO(method, doc=doc)


def _qualified(inst, name):
    '''(INTERNAL) Fully qualify a name.
    '''
    # _DOT_(inst.classname, name), not _DOT_(inst.named4, name)
    t =  type(inst)
    q = _DOT_(t.__module__, typename(t), name)  # _DMODULE_
    return q


class DeprecationWarnings(object):
    '''(INTERNAL) Handle C{DeprecationWaring}s.
    '''
    _Warnings = 0

    def __call__(self):  # for backward compatibility
        '''Have any C{DeprecationWarning}s been reported or raised?

           @return: The number of C{DeprecationWarning}s (C{int}) so
                    far or C{None} if not enabled.

           @note: To get C{DeprecationWarning}s if any, run C{python}
                  with env var C{PYGEODESY_WARNINGS} set to a non-empty
                  string I{AND} use C{python[3]} command line option
                  C{-X dev}, C{-W always} or C{-W error}, etc.
        '''
        return self.Warnings

    @property_ROver
    def _Fmt(self):
        '''Get C{streprs.Fmt}, I{once}.
        '''
        return _MODS.streprs.Fmt

    @property_ROver
    def _stacklevel3(self):
        '''Get C{dict(stacklevel=3)}, I{once}.
        '''
        return dict(stacklevel=3)

    def throw(self, kind, name, doc, **stacklevel):  # stacklevel=3
        '''Report or raise a C{DeprecationWarning}.

           @arg kind: Warning kind (C{str}), C{"method"}, C{"funtion"}, ...
           @arg name: Qualified name (C{str}) of B{C{kind}}.
           @arg doc: The __doc__ (C{str}) of B{C{kind}}, C{"DEPRECATED ...}.
        '''
        link = _tailof(name) or name
        if link is not name:  # make "link<name>"
            link = self._Fmt.ANGLE(link, name)
        link =  self._Fmt.CURLY(L=link)  # "L{link}"
        text =  doc.split(_DNL_, 1)[0].strip()
        text = _SPACE_(kind, link, _has_been_, *text.split())
        kwds =  stacklevel if stacklevel else self._stacklevel3
        # XXX invoke warn or raise DeprecationWarning(text)
        self._warn(text, category=DeprecationWarning, **kwds)
        self._Warnings += 1

    @property_ROver
    def _warn(self):
        '''Get Python's C{warnings.warn} function, I{once}.
        '''
        from warnings import warn as w
        return w

    @property_RO
    def Warnings(self):
        '''Get the number of C{DeprecationWarning}s (C{int}) so
           far or C{None} if not enabled.
        '''
        return self._Warnings if _WARNINGS_X_DEV else None

DeprecationWarnings = DeprecationWarnings()  # PYCHOK singleton
_throwarning        = DeprecationWarnings.throw

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
