
__all__ = ["PhysicalParameters"]

class PhysicalParameters(object):
    """
    A class to store physical constants and other immutable parameters
    used by the sims_photUtils code
    """

    #the quantities below are in nanometers
    minwavelen = 300.0
    maxwavelen = 1150.0
    wavelenstep = 0.1

    lightspeed = 299792458.0      # speed of light, = 2.9979e8 m/s
    planck = 6.626068e-27        # planck's constant, = 6.626068e-27 ergs*seconds
    nm2m = 1.00e-9               # nanometers to meters conversion = 1e-9 m/nm
    ergsetc2jansky = 1.00e23     # erg/cm2/s/Hz to Jansky units (fnu)

