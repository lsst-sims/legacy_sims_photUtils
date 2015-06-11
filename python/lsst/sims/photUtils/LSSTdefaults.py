import numpy

__all__ = ["LSSTdefaults"]

class LSSTdefaults(object):
    """
    This class exists to store default values of seeing, m5, and gamma taken from the over
    view paper (arXiv 0805.2366, Table 2, 29 August 2014 version)
    """

    def __init__(self):
        self._seeing = {'u': 0.77, 'g':0.73, 'r':0.70, 'i':0.67, 'z':0.65, 'y':0.63}
        self._m5 = {'u':23.68, 'g':24.89, 'r':24.43, 'i':24.00, 'z':24.45, 'y':22.60}
        self._gamma = {'u':0.037, 'g':0.038, 'r':0.039, 'i':0.039, 'z':0.040, 'y':0.040}


    def m5(self, tag):
        """
        From arXiv 0805.2366 29 August 2014 version  (Table 2):

        Typical 5-sigma depth for point sources at zenith, assuming
        exposure time of 2 x 15 seconds and observing conditions as listed

        @param [in] the name of a filter i.e. 'u', 'g', 'r', 'i', 'z', or 'y'

        @param [out] the corresponding m5 value
        """
        return self._m5[tag]


    def seeing(self, tag):
        """
        From arXiv 0805.2366 29 August 2014 version (Table 2):

        The expected delivered median zenith seeing in arcsec.  For larger
        airmass, X, seeing is proportional to X^0.6.

        @param [in] the name of a filter i.e. 'u', 'g', 'r', 'i', 'z', or 'y'

        @param [out] the corresponding seeing
        """

        return self._seeing[tag]


    def gamma(self, tag):
        """
        See Table 2 and Equaiton 5 of arXiv 0805.2366 29 August 2014 version.

        @param [in] the name of a filter i.e. 'u', 'g', 'r', 'i', 'z', or 'y'

        @param [out] the corresponding value of gamma as defined in the
        reference above
        """

        return self._gamma[tag]
