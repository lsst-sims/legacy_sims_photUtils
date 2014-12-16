import numpy as np

__all__ = ["Spectrum"]

class Spectrum(object):
    """
    This class holds the spectral data of each model.
    """

    def __init__(self, fileName):

        """
        @param [in] fileName is the name of the SED file you want to load.                                                                                             """

        self.wave, self.flux = self.readSpectrum(fileName)
        self.name = None
        self.type = None
        self.age = None
        self.metallicity = None

    def readSpectrum(self, fileName):
        """ 
        Loads the spectrum and saves it wavelengths and flux

        @param [in] fileName is the name of the SED file you want to load.

        @param [out] wave is the wavelength array from the file.

        @param [out] flux is the flux (flambda) array from the file.
        """
        wave, flux = np.genfromtxt(fileName, unpack=True)
        return wave, flux
