import six
import json
import logging
import numpy as np
#try:
#    import matplotlib.pyplot as plt
#    HAS_MATPLOTLIB = True
#except ImportError:
#    HAS_MATPLOTLIB = False


# initialize logger
logger = logging.getLogger(__name__)


def to_radians(angle, force_degrees=False):
    '''Convert arbitrary angular input to radians

    Supports slopes in radians, degrees or tangent notation ``y:x``,
    for example ``1:3``. If the input is a float larger than 2 pi it
    is assumed to be an angle in degrees. This behavior can be forced
    by an optional argument.

    Parameters
    ----------
    angle : float or str
      Float or string notation of a slope definition.
    force_degrees : bool
      Flag indicating whether to assume input is in degrees.

    Returns
    float
      Angle definition in radians.

    '''

    if isinstance(angle, six.string_types+(six.text_type, str)):
        if ':' in angle:
            dy, dx = angle.split(':')
            return np.arctan2(float(dy), float(dx))
        else:
            angle = float(angle)
            
    if np.abs(angle) > 2*np.pi or force_degrees:
        return np.radians(angle)

    return angle


class Profile:
    '''Class to define arbitrary cross-shore profiles consisting of linear
    segments. The class operates more or less as a list containing
    :class:`ProfileSegment` objects. The geometry of the profile is
    computed from the segments upon request. The geometry and position
    of individual segments are not necessarily defined, but can be
    derived from the interdependency between segments. The geometry
    and position of the profile as a whole should always be
    determined.

    Examples
    --------
    >>> p = Profile()
    >>> p.append(slope='1:3', z1=0., z2=10.)
    >>> p.append(slope=0, dx=10.)
    >>> p.append(slope='-1:3', z2=0.)
    >>> p.plot()

    '''

    
    def __init__(self):

        self.segments = []


    def __repr__(self):
        s = ''
        s += 'Profile:\n'
        for segment in self.segments:
            s += '  %s\n' % '\n  '.join(repr(segment).strip().split('\n'))
        return s

    
    def __iter__(self):
        for segment in self.segments:
            yield segment


    def __len__(self):
        return self.length

    
    def prepend(self, *args, **kwargs):
        self.insert(0, *args, **kwargs)

        
    def append(self, *args, **kwargs):
        self.segments.append(self._segment(*args, **kwargs))


    def insert(self, i, *args, **kwargs):
        self.segments.insert(i, self._segment(*args, **kwargs))


    def pop(self, i=-1):
        return self.segments.pop(i)


    def split(self, x=None, z=None): # FIXME
        xp, zp = self.get_geometry()
        
        if x is not None:
            ix = (x >= xp[:-1]) & (x < xp[1:])
            if any(ix):
                ix = np.where(ix)[0]
                segment = self.segments.pop(ix)
                segments = segment.split(
                    x=(x - xp[ix]) / (xp[ix+1] - xp[ix]))
                self.segments.insert(ix, segments[1])
                self.segments.insert(ix, segments[0])
            else:
                raise ValueError('Coordinate outside profile scope.')

        
    def get_geometry(self):

        n = len(self.segments)
        x = [None] * (n+1)
        z = [None] * (n+1)
        dx = [None] * n
        dz = [None] * n

        # collect known positions
        for i, segment in enumerate(self.segments):
            if segment.x1 is not None:
                if x[i] is None:
                    x[i] = segment.x1
                elif x[i] != segment.x1:
                    raise ValueError('Geometry overdetermined in x-direction: '
                                     '%0.2f and %0.2f.' % (x[i], segment.x1))
            if segment.x2 is not None:
                if x[i+1] is None:
                    x[i+1] = segment.x2
                elif x[i+1] != segment.x2:
                    raise ValueError('Geometry overdetermined in x-direction: '
                                     '%0.2f and %0.2f.' % (x[i+1], segment.x2))
            if segment.z1 is not None:
                if z[i] is None:
                    z[i] = segment.z1
                elif z[i] != segment.z1:
                    raise ValueError('Geometry overdetermined in z-direction: '
                                     '%0.2f and %0.2f.' % (z[i], segment.z1))
            if segment.z2 is not None:
                if z[i+1] is None:
                    z[i+1] = segment.z2
                elif z[i+1] != segment.z2:
                    raise ValueError('Geometry overdetermined in x-direction: '
                                     '%0.2f and %0.2f.' % (z[i+1], segment.z2))
            
            dx[i], dz[i] = segment.get_geometry()

        # assume origin if underdetermined
        if all([xi is None for xi in x]):
            x[0] = 0.
        if all([zi is None for zi in z]):
            z[0] = 0.

        # forward sweep
        for i in range(n):
            if x[i+1] is None and x[i] is not None and dx[i] is not None:
                x[i+1] = x[i] + dx[i]
            if z[i+1] is None and z[i] is not None and dz[i] is not None:
                z[i+1] = z[i] + dz[i]
            if i < n-1 and (dx[i+1] is None or dz[i+1] is None):
                geom = self.segments[i+1].get_geometry(x1=x[i+1], z1=z[i+1])
                if dx[i+1] is None:
                    dx[i+1] = geom[0]
                if dz[i+1] is None:
                    dz[i+1] = geom[1]

        # backward sweep
        for i in range(n-1, -1, -1):
            if x[i-1] is None and x[i] is not None and dx[i-1] is not None:
                x[i-1] = x[i] - dx[i-1]
            if z[i-1] is None and z[i] is not None and dz[i-1] is not None:
                z[i-1] = z[i] - dz[i-1]
        
        if any([xi is None for xi in x]) or any([zi is None for zi in z]):
            logger.debug('Current profile coordinates: %s' % list(zip(x,z)))
            logger.debug('Current section geometries: %s' % list(zip(dx,dz)))
            raise ValueError('Geometry underdetermined.')
                
        return x, z


    def to_json(self, **kwargs):
        return json.dumps([dict(s) for s in self.segments], **kwargs).replace(' ','')

    
#    def plot(self, ax=None, **kwargs):
#        if not HAS_MATPLOTLIB:
#            raise ImportError('Matplotlib not available.')
#        if ax is None:
#            fig, ax = plt.subplots()
#        else:
#            fig = ax.get_figure()
#        x, z = self.get_geometry()
#        ax.plot(x, z, **kwargs)
#        return fig, ax

    
    @staticmethod
    def _segment(*args, **kwargs):
        if len(args) > 0 and isinstance(args[0], ProfileSegment):
            return args[0]
        else:
            return ProfileSegment(*args, **kwargs)


    @property
    def length(self):
        x = self.get_geometry()[0]
        if len(x) > 2:
            return x[-1] - x[0]
        else:
            return np.nan
        

    @property
    def height(self):
        z = self.get_geometry()[1]
        if len(z) > 2:
            return z[-1] - z[0]
        else:
            return np.nan


class ProfileSegment:


    def __init__(self, x1=None, z1=None, x2=None, z2=None, dx=None,
                 dz=None, slope=None):
        
        self.x1 = x1
        self.z1 = z1
        self.x2 = x2
        self.z2 = z2
        self.dx = dx
        self.dz = dz
        self.slope = slope


    def __repr__(self):
        s = ''
        s += 'ProfileSegment:\n'
        for k, v in dict(self).items():
            if isinstance(v, six.string_types+(six.text_type, str)):
                s += '  %s = %s\n' % (k, v)
            else:
                s += '  %s = %0.2f\n' % (k, v)
        return s


    def __iter__(self):
        if self.x1 is not None:
            yield 'x1', self.x1
        if self.z1 is not None:
            yield 'z1', self.z1
        if self.x2 is not None:
            yield 'x2', self.x2
        if self.z2 is not None:
            yield 'z2', self.z2
        if self.dx is not None:
            yield 'dx', self.dx
        if self.dz is not None:
            yield 'dz', self.dz
        if self.slope is not None:
            yield 'slope', self.slope


    def split(self, fx=None, fz=None): # FIXME
        segments = [ProfileSegment(**dict(self)),
                    ProfileSegment(**dict(self))]

        if x is not None:
            if self.dx is not None:
                segments[0].dx *= fx
                segments[1].dx *= 1. - fx
            
        return segments

        
    def get_geometry(self, x1=None, z1=None):

        dx = set()
        if self.dx is not None:
            dx.add(self.dx)
        if self.x1 is not None and self.x2 is not None:
            dx.add(self.x2 - self.x1)
        if self.dz is not None and self.slope is not None:
            dx.add(self.dz / np.tan(to_radians(self.slope)))
        if self.z1 is not None and self.z2 is not None and self.slope is not None:
            dx.add((self.z2 - self.z1) / np.tan(to_radians(self.slope)))

        dz = set()
        if self.dz is not None:
            dz.add(self.dz)
        if self.z1 is not None and self.z2 is not None:
            dz.add(self.z2 - self.z1)
        if self.dx is not None and self.slope is not None:
            dz.add(self.dx * np.tan(to_radians(self.slope)))
        if self.x1 is not None and self.x2 is not None and self.slope is not None:
            dz.add((self.x2 - self.x1) * np.tan(to_radians(self.slope)))

        if len(dx) == 0:
            logger.debug('Geometry underdetermined in x-direction, '
                         'using relative position...')
        if len(dz) == 0:
            logger.debug('Geometry underdetermined in z-direction, '
                         'using relative position...')

        if len(dx) == 0 and z1 is not None and self.z2 is not None and self.slope is not None:
            dx.add((self.z2 - z1) / np.tan(to_radians(self.slope)))
        if len(dz) == 0 and z1 is not None and self.z2 is not None:
            dz.add(self.z2 - z1)
        if len(dx) == 0 and x1 is not None and self.x2 is not None:
            dx.add(self.x2 - x1)
        if len(dz) == 0 and x1 is not None and self.x2 is not None and self.slope is not None:
            dz.add((self.x2 - x1) * np.tan(to_radians(self.slope)))

        if len(dx) == 0:
            logger.debug('Geometry still underdetermined in x-direction.')
            return None, None
        if len(dz) == 0:
            logger.debug('Geometry still underdetermined in z-direction.')
            return None, None

        if len(dx) > 1:
            raise ValueError('Geometry overdetermined in x-direction: %s' % str(dx))
        if len(dz) > 1:
            raise ValueError('Geometry overdetermined in z-direction: %s' % str(dz))

        return dx.pop(), dz.pop()


    def to_json(self, **kwargs):
        return json.dumps(dict(self), **kwargs).replace(' ','')
