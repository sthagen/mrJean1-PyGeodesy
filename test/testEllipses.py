
# -*- coding: utf-8 -*-

# Test module L{ellipses}.

__all__ = ('Tests',)
__version__ = '26.02.14'

from bases import startswith, TestsBase

from pygeodesy import Ellipse, Ellipsoids, EPS, PI_4, elliperim


class Tests(TestsBase):

    def testEllipse(self):

        def _a(name, n=14):
            return Ellipse.__name__ + '.' + (name + ' ' * 12)[:n]

        def _n(name):
            return _a('perimeter' + name)

        a = 6378172.0
        for b, x in ((6378102.0, '40075016.6858801'),
                     (a * 0.9,   '38097844.6222377'),
                     (a / 2,     '30897294.5'),
                     (a / 4,     '273573'),
                     (a / 8,     '26106'),
                     (EPS,       '25512')):
            m = 1. - (b / a)**2
            self.test('a, b, b/a, m', (a, b, (b / a), m), '(6378172.0, ', known=startswith, nl=1)
            E = Ellipse(a, b)
            for p, n in ((E.perimeter2k,  '2k'),
                         (E.perimeter2k_, '2k_'),
                         (E.perimeterAGM, 'AGM'),
                         (E.perimeterHGK, 'HGK'),
                         (E.perimeterGK,  'GK'),
                         (E.perimeter2R,  '2R')):
                if p is not None:
                    self.test(_n(n), p, x, known=startswith, prec=9)

            p = E.perimeter4Arc3
            t = str(p).lstrip('(').rstrip(')')
            self.test(_n('4Arc3'), t, x[:1], known=startswith)

            n = elliperim.__name__ + ' DEPRECATED  '  # for backward compatibility
            x = x.split('.')[0]  # str(int(float(x)))
            self.test(n, elliperim(a, b), x, known=startswith, nl=1)
            self.test(n, elliperim(a, a), 40075236.597, prec=3)
            self.test(n, elliperim(a, 0), a * 4, prec=1)
            self.test(n, elliperim(0, b), b * 4, prec=1)
            self.test(n, elliperim(0, 0), '0.0')

            m =  Ellipse(b, a).arc
            n = _a('arc')
            x =  E.arc(45)
            self.test(n, x, x, prec=6)
            self.test(n, m(90, 45), x, prec=6)
            x = E.arc(270)
            self.test(n, x, x, prec=6)
            self.test(n, m(360, 90), x, prec=6)

            x = E.foci
            self.test(_a('foci'), x, x)
            self.test(_a('R2'), E.R2, E.R2)
            x = E.Roc_(PI_4)
            self.test(_a('Roc_'), x, x)
            x = E.sideOf(0, 1)
            self.test(_a('sideOf'), x, x)

            x = E.toEllipsoid().toStr()
            self.test(_a('toEllipsoid', n=-12), x, x)
            x = E.toTriaxial_().toStr()
            self.test(_a('toTriaxial_', n=-12), x, x, nt=1)

            for n, W in Ellipsoids.items(all=True, asorted=True):
                self.test(n + '.polarimeter/2k', W.polarimeter, W.toEllipse().perimeter2k, prec=6)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testEllipse()
    t.results()
    t.exit()
