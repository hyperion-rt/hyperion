from distutils.version import LooseVersion

import yt

if LooseVersion(yt.__version__) >= LooseVersion('3'):
    from .yt3_wrappers import *
else:
    from .yt2_wrappers import *