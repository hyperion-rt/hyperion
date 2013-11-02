from __future__ import print_function, division

from ..util.functions import FreezableClass


class Density(FreezableClass):
    """
    Base class for density structures
    """
    pass


class Disk(Density):
    """
    Base class for disk density structures
    """
    pass


class Envelope(Density):
    """
    Base class for envelope density structures
    """
    pass
