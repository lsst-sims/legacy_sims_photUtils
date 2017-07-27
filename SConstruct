# -*- python -*-
import os
from lsst.sconsUtils import scripts
scripts.BasicSConstruct("sims_photUtils")
base_env = Environment(ENV=os.environ)
for key in os.environ:
    base_env = os.environ[key]
