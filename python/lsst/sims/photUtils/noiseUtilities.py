import numpy
from .Sed import Sed
from .Bandpass import Bandpass
from lsst.sims.photUtils import PhotometricDefaults

__all__ = ["calcM5", "setM5"]

def setM5(m5target, skysed, totalBandpass, hardware,
          expTime=PhotometricDefaults.exptime,
          nexp=PhotometricDefaults.nexp,
          readnoise=PhotometricDefaults.rdnoise,
          darkcurrent=PhotometricDefaults.darkcurrent,
          othernoise=PhotometricDefaults.othernoise,
          seeing=PhotometricDefaults.seeing['r'],
          platescale=PhotometricDefaults.platescale,
          gain=PhotometricDefaults.gain,
          effarea=PhotometricDefaults.effarea):
    """
    Take an SED representing the sky and normalize it so that
    m5 (the magnitude at which an object is detected in this
    bandpass at 5-sigma) is set to some specified value.

    @param [in] the desired value of m5

    @param [in] skysed is an instantiation of the Sed class representing
    sky emission

    @param [in] totalBandpass is an instantiation of the Bandpass class
    representing the total throughput of the telescope (instrumentation
    plus atmosphere)

    @param [in] hardware is an instantiation of the Bandpass class representing
    the throughput due solely to instrumentation.

    @param [in] expTime is the duration of a single exposure in seconds

    @param [in] nexp is the number of exposures being combined

    @param [in] readnoise

    @param [in] darkcurrent

    @param [in] othernoise

    @param [in] seeing in arcseconds

    @param [in] platescale in arcseconds per pixel

    @param [in] gain in electrons per ADU

    @param [in] effarea is the effective area of the primary mirror in square centimeters

    @param [out] returns an instantiation of the Sed class that is the skysed renormalized
    so that m5 has the desired value
    """

    #This is based on the LSST SNR document (v1.2, May 2010)
    #www.astro.washington.edu/users/ivezic/Astr511/LSST_SNRdoc.pdf

    #instantiate a flat SED
    flatSed = Sed()
    flatSed.setFlatSED()

    skySedOut = Sed(wavelen=numpy.copy(skysed.wavelen),
                    flambda=numpy.copy(skysed.flambda))

    #normalize the SED so that it has a magnitude equal to the desired m5
    fNorm = flatSed.calcFluxNorm(m5target, totalBandpass)
    flatSed.multiplyFluxNorm(fNorm)
    counts = flatSed.calcADU(totalBandpass, expTime=expTime*nexp, effarea=effarea, gain=gain)

    #calculate the effective number of pixels for a double-Gaussian PSF
    neff = flatSed.calcNeff(seeing, platescale)

    #calculate the square of the noise due to the instrument
    noise_instr_sq = flatSed.calcInstrNoiseSq(readnoise, darkcurrent, expTime, nexp, othernoise)

    #now solve equation 41 of the SNR document for the neff * sigma_total^2 term
    #given snr=5 and counts as calculated above
    nSigmaSq = (counts*counts)/25.0 - counts/gain

    skyNoiseTarget = nSigmaSq/neff - noise_instr_sq
    skyCountsTarget = skyNoiseTarget*gain
    skyCounts = skySedOut.calcADU(hardware, expTime=expTime*nexp, effarea=effarea, gain=gain) \
                    * platescale * platescale
    skySedOut.multiplyFluxNorm(skyCountsTarget/skyCounts)

    return skySedOut

def calcM5(skysed, totalBandpass, hardware, expTime=PhotometricDefaults.exptime,
           nexp=PhotometricDefaults.nexp, readnoise=PhotometricDefaults.rdnoise,
           darkcurrent=PhotometricDefaults.darkcurrent,
           othernoise=PhotometricDefaults.othernoise,
           seeing=PhotometricDefaults.seeing['r'], platescale=PhotometricDefaults.platescale,
           gain=PhotometricDefaults.gain, effarea=PhotometricDefaults.effarea):
    """
    Calculate the AB magnitude of a 5-sigma above sky background source.

    Pass into this function the bandpass, hardware only of bandpass, and sky sed objects.
    The exposure time, nexp, readnoise, darkcurrent, gain,
    seeing and platescale are also necessary.
    """
    #This comes from equation 45 of the SNR document (v1.2, May 2010)
    #www.astro.washington.edu/users/ivezic/Astr511/LSST_SNRdoc.pdf

    #create a flat fnu source
    flatsource = Sed()
    flatsource.setFlatSED()
    snr = 5.0
    v_n, noise_instr_sq, \
    noise_sky_sq, noise_skymeasurement_sq, \
    skycounts, neff = flatsource.calcNonSourceNoiseSq(skysed, hardware, readnoise,
                                                      darkcurrent, othernoise, seeing,
                                                      effarea, expTime, nexp, platescale,
                                                      gain)

    counts_5sigma = (snr**2)/2.0/gain + numpy.sqrt((snr**4)/4.0/gain + (snr**2)*v_n)

    #renormalize flatsource so that it has the required counts to be a 5-sigma detection
    #given the specified background
    counts_flat = flatsource.calcADU(totalBandpass, expTime=expTime*nexp, effarea=effarea, gain=gain)
    flatsource.multiplyFluxNorm(counts_5sigma/counts_flat)

    # Calculate the AB magnitude of this source.
    mag_5sigma = flatsource.calcMag(totalBandpass)
    return mag_5sigma
