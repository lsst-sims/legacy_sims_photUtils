import numpy
from .Sed import Sed
from .Bandpass import Bandpass
from lsst.sims.photUtils import PhotometricDefaults

__all__ = ["calcM5", "setM5", "calcGamma", "calcSNR_gamma"]

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

    The 5-sigma limiting magnitude (m5) for an observation is
    determined by a combination of the telescope and camera parameters
    (such as diameter of the mirrors and the readnoise) together with the
    sky background. This method (setM5) scales a provided sky background
    Sed so that an observation would have a target m5 value, for the
    provided hardware parameters. Using the resulting Sed in the
    'calcM5' method will return this target value for m5.

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

    The 5-sigma limiting magnitude (m5) for an observation is determined by
    a combination of the telescope and camera parameters (such as diameter
    of the mirrors and the readnoise) together with the sky background. This
    method (calcM5) calculates the expected m5 value for an observation given
    a sky background Sed and hardware parameters.

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

    @param [out] returns the value of m5 for the given bandpass and sky SED
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

def calcGamma(bandpass, m5,
              expTime=PhotometricDefaults.exptime,
              nexp=PhotometricDefaults.nexp,
              gain=PhotometricDefaults.gain,
              effarea=PhotometricDefaults.effarea):

    """
    Calculate the gamma parameter used for determining photometric
    signal to noise in equation 5 of the LSST overview paper
    (arXiv:0805.2366)

    @param [in] bandpass is an instantiation of the Bandpass class
    representing the bandpass for which you desire to calculate the
    gamma parameter

    @param [in] m5 is the magnitude at which a 5-sigma detection occurs
    in this Bandpass

    @param [in] expTime is the duration of a single exposure in seconds

    @param [in] nexp is the number of exposures being combined

    @param [in] gain is the number of electrons per ADU

    @param [in] effarea is the effective area of the primary mirror
    in square centimeters

    @param [out] gamma
    """
    #This is based on the LSST SNR document (v1.2, May 2010)
    #www.astro.washington.edu/users/ivezic/Astr511/LSST_SNRdoc.pdf
    #as well as equations 4-6 of the overview paper (arXiv:0805.2366)

    #instantiate a flat SED
    flatSed = Sed()
    flatSed.setFlatSED()

    #normalize the SED so that it has a magnitude equal to the desired m5
    fNorm = flatSed.calcFluxNorm(m5, bandpass)
    flatSed.multiplyFluxNorm(fNorm)
    counts = flatSed.calcADU(bandpass, expTime=expTime*nexp, effarea=effarea, gain=gain)

    #The expression for gamma below comes from:
    #
    #1) Take the approximation N^2 = N0^2 + alpha S from footnote 88 in the overview paper
    #where N is the noise in flux of a source, N0 is the noise in flux due to sky brightness
    #and instrumentation, S is the number of counts registered from the source and alpha
    #is some constant
    #
    #2) Divide by S^2 and demand that N/S = 0.2 for a source detected at m5. Solve
    #the resulting equation for alpha in terms of N0 and S5 (the number of counts from
    #a source at m5)
    #
    #3) Substitute this expression for alpha back into the equation for (N/S)^2
    #for a general source.  Re-factor the equation so that it looks like equation
    #5 of the overview paper (note that x = S5/S).  This should give you gamma = (N0/S5)^2
    #
    #4) Solve equation 41 of the SNR document for the neff * sigma_total^2 term
    #given snr=5 and counts as calculated above.  Note that neff * sigma_total^2
    #is N0^2 in the equation above
    #
    #This should give you

    gamma = 0.04 - 1.0/(counts*gain)

    return gamma

def calcSNR_gamma(fluxes, bandpasses, m5, gamma=None, sig2sys=None,
                 expTime=PhotometricDefaults.exptime,
                 nexp=PhotometricDefaults.nexp,
                 gain=PhotometricDefaults.gain,
                 effarea=PhotometricDefaults.effarea):
    """
    Calculate signal to noise in flux using the model from equation (5) of arXiv:0805.2366

    @param [in] fluxes is a numpy array of fluxes.  Each row is a different bandpass.
    Each column is a different object, i.e. fluxes[i][j] is the flux of the jth object
    in the ith bandpass.

    @param [in] bandpasses is a list of Bandpass objects corresponding to the
    bandpasses in which fluxes have been calculated

    @param [in] m5 is a numpy.array of 5-sigma limiting magnitudes, one for each bandpass.

    @param [in] gamma (optional) is the gamma parameter from equation(5) of
    arXiv:0805.2366.  If not provided, this method will calculate it.

    @param [in] sig2sys is the square of the systematic signal to noise ratio.

    @param [in] expTime (optional) is the duration of a single exposure in seconds

    @param [in] nexp (optional) is the number of exposures being combined

    @param [in] gain (optional) is the number of electrons per ADU

    @param [in] effarea (optional) is the effective area of the primary mirror
    in square centimeters

    @param [out] snr is a numpy array of the signal to noise ratio corresponding to
    fluxes.

    @param [out] gamma is a numpy array of the calculated gamma parameters for the
    bandpasses used here (in case the user wants to call this method again.
    """

    if fluxes.shape[0] != len(bandpasses):
        raise RuntimeError("Passed %d magnitudes to " % fluxes.shape[0] + \
                            " calcSNR_gamma; " + \
                            "passed %d bandpasses" % len(bandpasses))

    if gamma is not None and len(gamma) != len(bandpasses):
        raise RuntimeError("Passed %d bandpasses to " % len(bandpasses) + \
                           " calcSNR_gamma; " + \
                           "passed %d gamma parameters" % len(gamma))

    if len(m5) != len(bandpasses):
        raise RuntimeError("Passed %d bandpasses to " % len(bandpasses) + \
                           " calcSNR_gamma; " + \
                           "passed %d m5 values" % len(m5))

    if gamma is None:
        gg = []
        for b, m in zip(bandpasses, m5):
            gg.append(calcGamma(b, m, expTime=expTime, nexp=nexp, gain=gain, effarea=effarea))

        gamma = numpy.array(gg)

    m5Fluxes = numpy.array(numpy.power(10.0, -0.4*m5))

    noise = []
    for (gg, mf, ff) in zip(gamma, m5Fluxes, fluxes):
        fluxRatio = mf/ff

        if sig2sys is not None:
            sigmaSq = (0.04-gg)*fluxRatio+gg*fluxRatio*fluxRatio + sig2sys
        else:
            sigmaSq = (0.04-gg)*fluxRatio+gg*fluxRatio*fluxRatio

        noise.append(numpy.sqrt(sigmaSq))

    return 1.0/numpy.array(noise), gamma
