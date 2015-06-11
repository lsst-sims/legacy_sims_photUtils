import numpy

__all__ = ["PhotometricParameters"]

class PhotometricParameters(object):

    def __init__(self, exptime=15.0,
                 nexp=2,
                 effarea=numpy.pi*(6.5*100./2.0)**2,
                 gain=2.3,
                 readnoise=5.,
                 darkcurrent=0.2,
                 othernoise=4.69,
                 platescale=0.2):

        """
        @param [in] exptime exposure time in seconds (default 15)

        @param [in] nexp number of expusres (default 2)

        @param [in] effarea effective area in cm^2 (default for 6.5 meter diameter)

        @param [in] gain electrons per ADU (default 2.3)

        @param [in] readnoise electrons per pixel per exposure (default 5.0)

        @param [in] darkcurrent electons per pixel per second (default 0.2)

        @param [in] othernoise electrons per pixel per exposure (default 4.69)

        @param [in] platescale arcseconds per pixel (default 0.2)
        """

        self._exptime = exptime
        self._nexp = nexp
        self._effarea = effarea
        self._gain = gain
        self._platescale = platescale

        #The quantities below are measured in electrons.
        #This is taken from the specifications document LSE-30 on Docushare
        #Section 3.4.2.3 states that the total noise per pixel shall be 12.7 electrons
        #which the defaults sum to (remember to multply darkcurrent by the number
        #of seconds in an exposure=15).
        self._readnoise = readnoise
        self._darkcurrent = darkcurrent
        self._othernoise = othernoise


    @property
    def exptime(self):
        """
        exposure time in seconds
        """
        return self._exptime

    @exptime.setter
    def exptime(self, value):
        raise RuntimeError("You should not be setting exptime on the fly; " +
                           "Just instantiate a new case of PhotometricParameters")


    @property
    def nexp(self):
        """
        number of exposures
        """
        return self._nexp

    @nexp.setter
    def nexp(self, value):
        raise RuntimeError("You should not be setting nexp on the fly; " +
                           "Just instantiate a new case of PhotometricParameters")


    @property
    def effarea(self):
        """
        effective area in cm^2
        """
        return self._effarea

    @effarea.setter
    def effarea(self, value):
        raise RuntimeError("You should not be setting effarea on the fly; " +
                           "Just instantiate a new case of PhotometricParameters")


    @property
    def gain(self):
        """
        electrons per ADU
        """
        return self._gain

    @gain.setter
    def gain(self, value):
        raise RuntimeError("You should not be setting gain on the fly; " +
                           "Just instantiate a new case of PhotometricParameters")


    @property
    def platescale(self):
        """
        arcseconds per pixel
        """
        return self._platescale

    @platescale.setter
    def platescale(self, value):
        raise RuntimeError("You should not be setting platescale on the fly; " +
                           "Just instantiate a new case of PhotometricParameters")


    @property
    def readnoise(self):
        """
        electrons per pixel per exposure
        """
        return self._readnoise

    @readnoise.setter
    def readnoise(self, value):
        raise RuntimeError("You should not be setting readnoise on the fly; " +
                           "Just instantiate a new case of PhotometricParameters")


    @property
    def darkcurrent(self):
        """
        electrons per pixel per second
        """
        return self._darkcurrent

    @darkcurrent.setter
    def darkcurrent(self, value):
        raise RuntimeError("You should not be setting darkcurrent on the fly; " +
                           "Just instantiate a new case of PhotometricParameters")


    @property
    def othernoise(self):
        """
        electrons per pixel per exposure
        """
        return self._othernoise

    @othernoise.setter
    def othernoise(self,value):
        raise RuntimeError("You should not be setting othernoise on the fly; " +
                           "Just instantiate a new case of PhotometricParameters")

