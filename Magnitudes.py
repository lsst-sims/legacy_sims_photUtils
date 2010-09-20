from Sed import Sed
from Bandpass import *

# Routines specific for instance catalog operations.
def loadSeds(sedList, dataDir = "./"):
    """Generate dictionary of SEDs required for generating magnitudes

    Given a dataDir and a list of seds return a dictionary with sedName and sed as key, value
    """
    sedDict={}
    for sedName in sedList:
        if sedName in sedDict:
            continue
        else:
            sed = Sed.Sed()
            sed.readSED_flambda(dataDir+"data/seds/"+ sedName)
            if sed.needResample():
                sed.resampleSED()             
            sedDict[sedName] = sed

    return sedDict

# Routines for InstanceCatalog
def loadBandpasses(bandpassList, dataDir="./"):
    """ Generate dictionary of bandpasses for the LSST nominal throughputs

    Given a list of of filter throughputs return a dictionary of filteNames and bandpass key, values
    """

    bandpassDict = {}
    for filter in bandpassList:
        bandpass = Bandpass()
        bandpass.readThroughput(dataDir  + filter + ".dat")
        bandpassDict[filter] = bandpass
    return bandpassDict

