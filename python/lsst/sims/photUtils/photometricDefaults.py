import numpy

__all__ = ["PhotometricDefaults"]

class PhotometricDefaults(object):
    """
    This class exists to store default values of parameters characterizing
    noise due to the telescope in one place.
    """

    # The following *wavelen* parameters are default values for gridding wavelen/sb/flambda.

    minwavelen = 300.0
    maxwavelen = 1150.0
    wavelenstep = 0.1

    lightspeed = 299792458.0      # speed of light, = 2.9979e8 m/s
    planck = 6.626068e-27        # planck's constant, = 6.626068e-27 ergs*seconds
    nm2m = 1.00e-9               # nanometers to meters conversion = 1e-9 m/nm
    ergsetc2jansky = 1.00e23     # erg/cm2/s/Hz to Jansky units (fnu)

    exptime = 15.0                    # Default exposure time. (option for method calls).
    nexp = 2                          # Default number of exposures. (option for methods).
    effarea = numpy.pi*(6.5*100/2.0)**2   # Default effective area of primary mirror. (option for methods).
    gain = 2.3                        # Default gain. (option for method call).
    rdnoise = 5                       # Default value - readnoise electrons or adu per pixel (per exposure)
    darkcurrent = 0.2                 # Default value - dark current electrons or adu per pixel per second
    othernoise = 4.69                 # Default value - other noise electrons or adu per pixel per exposure
    platescale = 0.2                  # Default value - "/pixel
    seeing = {'u': 0.77, 'g':0.73, 'r':0.70, 'i':0.67, 'z':0.65, 'y':0.63}  # Default seeing values (in ")
