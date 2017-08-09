import six
import json
import logging
import numpy as np
import matplotlib.pyplot as plt


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

    if isinstance(angle, six.string_types+(six.text_type,)):
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

    Notes
    -----
    Geometry solver is not finalized. The first segment should be
    determined in terms of geometry and position otherwise an error is
    raised. The algorithm should find a fully determined segment and
    work both ways from there.

    '''

    
    def __init__(self):

        self.x0 = None
        self.z0 = None
        self.segments = []


    def __repr__(self):
        s = ''
        s += 'Profile:\n'
        for segment in self.segments:
            s += '  %s\n' % '\n  '.join(repr(segment).strip().split('\n'))
        return s

    
    def __iter__(self):
        for segment in segments:
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
        x = [0]
        z = [0]
        for i, segment in enumerate(self.segments):
            logger.debug('Positioning section #%d...' % (i+1))
            dx, dz = segment.get_geometry(x1=x[i], z1=z[i])
            if dx is not None and dz is not None:
                x.append(x[i] + dx)
                z.append(z[i] + dz)
                logger.debug('Current profile coordinates: %s.' % str(list(zip(x, z))))
            
                if segment.x1 is not None:
                    self.set_reference_position(x=segment.x1-x[i])
                if segment.x2 is not None:
                    self.set_reference_position(x=segment.x2-x[i+1])
                if segment.z1 is not None:
                    self.set_reference_position(z=segment.z1-z[i])
                if segment.z2 is not None:
                    self.set_reference_position(z=segment.z2-z[i+1])
                logger.debug('Current origin: %s.' % str((self.x0, self.z0)))
            else:
                x.append(None)
                z.append(None)
                raise ValueError('Geometry underdetermined.')

        # reverse sweep
        # for i in range(len(self.segments))[:-1:-1]:
        #     segment = self.segments[i]
        #     if x[i+1] is None or z[i+1] is None:
        #         logger.debug('Positioning section #%d...' % (i+1))
        #         dx, dz = segment.get_geometry(x2=x[i+2], z2=z[i+2])
        #         if dx is not None and dz is not None:
        #             x[i+1] = x[i+2] - dx
        #             z[i+1] = z[i+2] - dz
        #             logger.debug('Current profile coordinates: %s.' % str(list(zip(x, z))))
        #         else:
        #             raise ValueError('Geometry underdetermined.')

        x = np.asarray(x)
        z = np.asarray(z)
        
        if self.x0 is not None:
            x += self.x0
        if self.z0 is not None:
            z += self.z0
                
        return x, z


    def set_reference_position(self, x=None, z=None):
        if x is not None:
            if self.x0 is None:
                self.x0 = x
            elif self.x0 != x:
                raise ValueError('Position overdetermined in x-direction: '
                                 'x=%0.2f and x=%0.2f.' % (self.x0, x))

        if z is not None:
            if self.z0 is None:
                self.z0 = z
            elif self.z0 != z:
                raise ValueError('Position overdetermined in z-direction: '
                                 'z=%0.2f and z=%0.2f.' % (self.z0, z))


    def to_json(self, **kwargs):
        return json.dumps([dict(s) for s in self.segments], **kwargs).replace(' ','')

    
    def plot(self, ax=None, **kwargs):
        
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        x, z = self.get_geometry()
        ax.plot(x, z, **kwargs)
        return fig, ax

    
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
            if isinstance(v, six.string_types+(six.text_type,)):
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
