from __future__ import absolute_import

import numpy as np
from profilebuilder.profilebuilder import *


class DikeProfile(Profile):


    def __init__(self, length=100., crest_height=10., crest_width=10.,
                 slope='1:3'):

        super(DikeProfile, self).__init__()

        self.append(slope=slope, z1=0., z2=crest_height)
        self.append(slope=0., dx=crest_width)
        self.append(slope=-to_radians(slope), z2=0)

        dL = np.maximum(0., length - self.length)

        self.prepend(slope=0., dx=dL/2.)
        self.append(slope=0., dx=dL/2.)


class DikeWithBermProfile(DikeProfile):


    def __init__(self, length=100., berm_height=5., berm_width=5.,
                 slope='1:3', **kwargs):

        super(DikeWithBermProfile, self).__init__(length=length,
                                                  slope=slope, **kwargs)

        self.insert(1, slope=slope, z1=0., z2=berm_height)
        self.insert(2, slope=0., dx=berm_width)
        self.segments[3].z1 = None

        dL = np.maximum(0., self.length - length)

        self.segments[0].dx -= dL/2.
        self.segments[-1].dx -= dL/2.
