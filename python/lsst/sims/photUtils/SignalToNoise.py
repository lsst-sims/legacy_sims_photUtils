import numpy
from .Sed import Sed
from .Bandpass import Bandpass
from lsst.sims.photUtils import LSSTdefaults

__all__ = ["calcNeff", "calcInstrNoiseSq", "calcTotalNonSourceNoiseSq", "calcSNR_sed",
          "calcM5", "calcSkyCountsPerPixelForM5", "calcSkyCountsForM5", "calcGamma", "calcSNR_m5",
          "calcAstrometricError", "magErrorFromSNR", "calcMagError_m5", "calcMagError_sed"]


def calcNeff(seeing, platescale):
    """
    Calculate the effective number of pixels in a double Gaussian PSF

    @param [in] seeing in arcseconds

    @param [in] platescale in arcseconds per pixel

    @param [out] the effective number of pixels contained by a double
    Gaussian PSF

    see equation 31 of the SNR document
    https://docushare.lsstcorp.org/docushare/dsweb/ImageStoreViewer/LSE-40
    """
    return 2.436*(seeing/platescale)**2


def calcInstrNoiseSq(photParams):
    """
    Combine all of the noise due to intrumentation into one value

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [out] The noise due to all of these sources added in quadrature
    in ADU counts
    """

    return (photParams.nexp*photParams.readnoise**2 + \
           photParams.darkcurrent*photParams.exptime*photParams.nexp + \
           photParams.nexp*photParams.othernoise**2)/(photParams.gain*photParams.gain)



def calcTotalNonSourceNoiseSq(skySed, hardwarebandpass, photParams, seeing):
    """
    Calculate the noise due to things that are not the source being observed
    (i.e. intrumentation and sky background)

    @param [in] skySed -- an instantiation of the Sed class representing the sky
    (normalized so that skySed.calcMag() gives the sky brightness in magnitudes
    per square arcsecond)

    @param [in] hardwarebandpass -- an instantiation of the Bandpass class representing
    just the instrumentation throughputs

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] seeing in arcseconds

    @param [out] total non-source noise squared (in ADU counts)
    (this is simga^2_tot * neff in equation 41 of the SNR document
    https://docushare.lsstcorp.org/docushare/dsweb/ImageStoreViewer/LSE-40 )

    @param [out] noise squared due just to the instrument (in ADU counts)

    @param [out] noise squared due to the sky (in ADU counts)

    @param [out] noise squared due to sky measurement (in ADU counts, presently set to zero)

    @param [out] the effective number of pixels in a double Gaussian PSF
    """

    # This method outputs all of the parameters calculated along the way
    # so that the verbose version of calcSNR_sed still works

    # Calculate the effective number of pixels for double-Gaussian PSF
    neff = calcNeff(seeing, photParams.platescale)

    # Calculate the counts from the sky.
    # We multiply by two factors of the platescale because we expect the
    # skySed to be normalized such that calcADU gives counts per
    # square arc second, and we need to convert to counts per pixel.

    skycounts = skySed.calcADU(hardwarebandpass, photParams=photParams) \
                * photParams.platescale * photParams.platescale

    # Calculate the square of the noise due to instrumental effects.
    # Include the readout noise as many times as there are exposures

    noise_instr_sq = calcInstrNoiseSq(photParams=photParams)

    # Calculate the square of the noise due to sky background poisson noise
    noise_sky_sq = skycounts/photParams.gain

    # Discount error in sky measurement for now
    noise_skymeasurement_sq = 0

    total_noise_sq = neff*(noise_sky_sq + noise_instr_sq + noise_skymeasurement_sq)
    return total_noise_sq, noise_instr_sq, noise_sky_sq, noise_skymeasurement_sq, skycounts, neff


def calcSkyCountsPerPixelForM5(m5target, totalBandpass, photParams, seeing=None):

    """
    Calculate the number of sky counts per pixel expected for a given
    value of the 5-sigma limiting magnitude (m5)

    The 5-sigma limiting magnitude (m5) for an observation is
    determined by a combination of the telescope and camera parameters
    (such as diameter of the mirrors and the readnoise) together with the
    sky background.

    @param [in] the desired value of m5

    @param [in] totalBandpass is an instantiation of the Bandpass class
    representing the total throughput of the telescope (instrumentation
    plus atmosphere)

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] seeing in arcseconds

    @param [out] returns the expected number of sky counts per pixel
    """

    if seeing is None:
        seeing = LSSTdefaults().seeing('r')

    # instantiate a flat SED
    flatSed = Sed()
    flatSed.setFlatSED()

    # normalize the SED so that it has a magnitude equal to the desired m5
    fNorm = flatSed.calcFluxNorm(m5target, totalBandpass)
    flatSed.multiplyFluxNorm(fNorm)
    sourceCounts = flatSed.calcADU(totalBandpass, photParams=photParams)

    # calculate the effective number of pixels for a double-Gaussian PSF
    neff = calcNeff(seeing, photParams.platescale)

    # calculate the square of the noise due to the instrument
    noise_instr_sq = calcInstrNoiseSq(photParams=photParams)

    # now solve equation 41 of the SNR document for the neff * sigma_total^2 term
    # given snr=5 and counts as calculated above
    # SNR document can be found at
    # https://docushare.lsstcorp.org/docushare/dsweb/ImageStoreViewer/LSE-40

    nSigmaSq = (sourceCounts*sourceCounts)/25.0 - sourceCounts/photParams.gain

    skyNoiseTarget = nSigmaSq/neff - noise_instr_sq
    skyCountsTarget = skyNoiseTarget*photParams.gain

    # TODO:
    # This method should throw an error if skyCountsTarget is negative
    # unfortunately, that currently happens for default values of
    # m5 as taken from arXiv:0805.2366, table 2.  Adding the error
    # should probably wait for a later issue in which we hash out what
    # the units are for all of the parameters stored in PhotometricDefaults.

    return skyCountsTarget


def calcSkyCountsForM5(m5target, totalBandpass, photParams, seeing=None):

    """
    Calculate the number of sky counts (total; not per pixel) expected for a given
    value of the 5-sigma limiting magnitude (m5)

    The 5-sigma limiting magnitude (m5) for an observation is
    determined by a combination of the telescope and camera parameters
    (such as diameter of the mirrors and the readnoise) together with the
    sky background.

    @param [in] the desired value of m5

    @param [in] totalBandpass is an instantiation of the Bandpass class
    representing the total throughput of the telescope (instrumentation
    plus atmosphere)

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] seeing in arcseconds

    @param [out] returns the expected number of sky counts
    """

    if seeing is None:
        seeing = LSSTdefaults().seeing('r')

    neff = calcNeff(seeing=seeing, platescale=photParams.platescale)

    return neff*calcSkyCountsPerPixelForM5(m5target, totalBandpass, photParams, seeing=seeing)



def calcM5(skysed, totalBandpass, hardware, photParams, seeing=None):
    """
    Calculate the AB magnitude of a 5-sigma above sky background source.

    The 5-sigma limiting magnitude (m5) for an observation is determined by
    a combination of the telescope and camera parameters (such as diameter
    of the mirrors and the readnoise) together with the sky background. This
    method (calcM5) calculates the expected m5 value for an observation given
    a sky background Sed and hardware parameters.

    @param [in] skysed is an instantiation of the Sed class representing
    sky emission, normalized so that skysed.calcMag gives the sky brightness
    in magnitudes per square arcsecond.

    @param [in] totalBandpass is an instantiation of the Bandpass class
    representing the total throughput of the telescope (instrumentation
    plus atmosphere)

    @param [in] hardware is an instantiation of the Bandpass class representing
    the throughput due solely to instrumentation.

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] seeing in arcseconds

    @param [out] returns the value of m5 for the given bandpass and sky SED
    """
    # This comes from equation 45 of the SNR document (v1.2, May 2010)
    # https://docushare.lsstcorp.org/docushare/dsweb/ImageStoreViewer/LSE-40

    if seeing is None:
        seeing = LSSTdefaults().seeing('r')

    # create a flat fnu source
    flatsource = Sed()
    flatsource.setFlatSED()
    snr = 5.0
    v_n, noise_instr_sq, \
    noise_sky_sq, noise_skymeasurement_sq, \
    skycounts, neff = calcTotalNonSourceNoiseSq(skysed, hardware, photParams, seeing)

    counts_5sigma = (snr**2)/2.0/photParams.gain + \
                     numpy.sqrt((snr**4)/4.0/photParams.gain + (snr**2)*v_n)

    # renormalize flatsource so that it has the required counts to be a 5-sigma detection
    # given the specified background
    counts_flat = flatsource.calcADU(totalBandpass, photParams=photParams)
    flatsource.multiplyFluxNorm(counts_5sigma/counts_flat)

    # Calculate the AB magnitude of this source.
    mag_5sigma = flatsource.calcMag(totalBandpass)
    return mag_5sigma


def magErrorFromSNR(snr):
    """
    convert flux signal to noise ratio to an error in magnitude

    @param [in] snr is the signal to noise ratio in flux

    @param [out] the resulting error in magnitude
    """

    #see www.ucolick.org/~bolte/AY257/s_n.pdf section 3.1
    return 2.5*numpy.log10(1.0+1.0/snr)


def calcGamma(bandpass, m5, photParams):

    """
    Calculate the gamma parameter used for determining photometric
    signal to noise in equation 5 of the LSST overview paper
    (arXiv:0805.2366)

    @param [in] bandpass is an instantiation of the Bandpass class
    representing the bandpass for which you desire to calculate the
    gamma parameter

    @param [in] m5 is the magnitude at which a 5-sigma detection occurs
    in this Bandpass

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [out] gamma
    """
    # This is based on the LSST SNR document (v1.2, May 2010)
    # https://docushare.lsstcorp.org/docushare/dsweb/ImageStoreViewer/LSE-40
    # as well as equations 4-6 of the overview paper (arXiv:0805.2366)

    # instantiate a flat SED
    flatSed = Sed()
    flatSed.setFlatSED()

    # normalize the SED so that it has a magnitude equal to the desired m5
    fNorm = flatSed.calcFluxNorm(m5, bandpass)
    flatSed.multiplyFluxNorm(fNorm)
    counts = flatSed.calcADU(bandpass, photParams=photParams)

    # The expression for gamma below comes from:
    #
    # 1) Take the approximation N^2 = N0^2 + alpha S from footnote 88 in the overview paper
    # where N is the noise in flux of a source, N0 is the noise in flux due to sky brightness
    # and instrumentation, S is the number of counts registered from the source and alpha
    # is some constant
    #
    # 2) Divide by S^2 and demand that N/S = 0.2 for a source detected at m5. Solve
    # the resulting equation for alpha in terms of N0 and S5 (the number of counts from
    # a source at m5)
    #
    # 3) Substitute this expression for alpha back into the equation for (N/S)^2
    # for a general source.  Re-factor the equation so that it looks like equation
    # 5 of the overview paper (note that x = S5/S).  This should give you gamma = (N0/S5)^2
    #
    # 4) Solve equation 41 of the SNR document for the neff * sigma_total^2 term
    # given snr=5 and counts as calculated above.  Note that neff * sigma_total^2
    # is N0^2 in the equation above
    #
    # This should give you

    gamma = 0.04 - 1.0/(counts*photParams.gain)

    return gamma

def calcSNR_m5(magnitudes, bandpasses, m5, photParams, gamma=None):
    """
    Calculate signal to noise in flux using the model from equation (5) of arXiv:0805.2366

    @param [in] magnitudes is a numpy array.  Each row is a different bandpass.
    Each column is a different object, i.e. magnitudes[i][j] is the magnitude of the jth object
    in the ith bandpass.

    @param [in] bandpasses is a list of Bandpass objects corresponding to the
    bandpasses in which magnitudes have been calculated

    @param [in] m5 is a numpy.array of 5-sigma limiting magnitudes, one for each bandpass.

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] gamma (optional) is the gamma parameter from equation(5) of
    arXiv:0805.2366.  If not provided, this method will calculate it.

    @param [out] snr is a numpy array of the signal to noise ratio corresponding to
    the input magnitudes.

    @param [out] gamma is a numpy array of the calculated gamma parameters for the
    bandpasses used here (in case the user wants to call this method again.
    """

    if magnitudes.shape[0] != len(bandpasses):
        raise RuntimeError("Passed %d magnitudes to " % magnitudes.shape[0] + \
                            " calcSNR_m5; " + \
                            "passed %d bandpasses" % len(bandpasses))

    if gamma is not None and len(gamma) != len(bandpasses):
        raise RuntimeError("Passed %d bandpasses to " % len(bandpasses) + \
                           " calcSNR_m5; " + \
                           "passed %d gamma parameters" % len(gamma))

    if len(m5) != len(bandpasses):
        raise RuntimeError("Passed %d bandpasses to " % len(bandpasses) + \
                           " calcSNR_m5; " + \
                           "passed %d m5 values" % len(m5))

    if gamma is None:
        gg = []
        for b, m in zip(bandpasses, m5):
            gg.append(calcGamma(b, m, photParams=photParams))

        gamma = numpy.array(gg)

    dummySed = Sed()
    m5Fluxes = dummySed.fluxFromMag(m5)
    sourceFluxes = dummySed.fluxFromMag(magnitudes)

    noise = []
    for (gg, mf, ff) in zip(gamma, m5Fluxes, sourceFluxes):
        fluxRatio = mf/ff

        noiseSq = (0.04-gg)*fluxRatio+gg*fluxRatio*fluxRatio

        noise.append(numpy.sqrt(noiseSq))

    return 1.0/numpy.array(noise), gamma


def calcMagError_m5(magnitudes, bandpasses, m5, photParams, gamma=None):
    """
    Calculate magnitude error using the model from equation (5) of arXiv:0805.2366

    @param [in] magnitudes is a numpy array.  Each row is a different bandpass.
    Each column is a different object, i.e. magnitudes[i][j] is the magnitude of the jth object
    in the ith bandpass.

    @param [in] bandpasses is a list of Bandpass objects corresponding to the
    bandpasses in which magnitudes have been calculated

    @param [in] m5 is a numpy.array of 5-sigma limiting magnitudes, one for each bandpass.

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] gamma (optional) is the gamma parameter from equation(5) of
    arXiv:0805.2366.  If not provided, this method will calculate it.

    @param [out] is a numpy array of errors in magnitude
    """

    snr, gamma = calcSNR_m5(magnitudes, bandpasses, m5, photParams, gamma=gamma)

    if photParams.sigmaSys is not None:
        return numpy.sqrt(numpy.power(magErrorFromSNR(snr),2) + numpy.power(photParams.sigmaSys,2))
    else:
        return magErrorFromSNR(snr)


def calcSNR_sed(spectrum, totalbandpass, skysed, hardwarebandpass,
                    photParams, seeing, verbose=False):
    """
    Calculate the signal to noise ratio for a source, given the bandpass(es) and sky SED.

    For a given source, sky sed, total bandpass and hardware bandpass, as well as
    seeing / expTime, calculates the SNR with optimal PSF extraction
    assuming a double-gaussian PSF.

    @param [in] spectrum is an instantiation of the Sed class containing the SED of
    the object whose signal to noise ratio is being calculated

    @param [in] totalbandpass is an instantiation of the Bandpass class
    representing the total throughput (system + atmosphere)

    @param [in] skysed is an instantiation of the Sed class representing
    the sky emission per square arcsecond.

    @param [in] hardwarebandpass is an instantiation of the Bandpass class
    representing just the throughput of the system hardware.

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] seeing in arcseconds

    @param [in] verbose is a boolean

    @param [out] signal to noise ratio
    """

    # Calculate the counts from the source.
    sourcecounts = spectrum.calcADU(totalbandpass, photParams=photParams)

    # Calculate the (square of the) noise due to signal poisson noise.
    noise_source_sq = sourcecounts/photParams.gain

    non_source_noise_sq, \
    noise_instr_sq, \
    noise_sky_sq, \
    noise_skymeasurement_sq, \
    skycounts, neff = calcTotalNonSourceNoiseSq(skysed, hardwarebandpass, photParams, seeing)

    # Calculate total noise
    noise = numpy.sqrt(noise_source_sq + non_source_noise_sq)
    # Calculate the signal to noise ratio.
    snr = sourcecounts / noise
    if verbose:
        print "For Nexp %.1f of time %.1f: " % (photParams.nexp, photParams.expTime)
        print "Counts from source: %.2f  Counts from sky: %.2f" %(sourcecounts, skycounts)
        print "Seeing: %.2f('')  Neff pixels: %.3f(pix)" %(seeing, neff)
        print "Noise from sky: %.2f Noise from instrument: %.2f" \
            %(numpy.sqrt(noise_sky_sq), numpy.sqrt(noise_instr_sq))
        print "Noise from source: %.2f" %(numpy.sqrt(noise_source_sq))
        print " Total Signal: %.2f   Total Noise: %.2f    SNR: %.2f" %(sourcecounts, noise, snr)
        # Return the signal to noise value.
    return snr


def calcMagError_sed(spectrum, totalbandpass, skysed, hardwarebandpass,
                    photParams, seeing, verbose=False):
    """
    Calculate the magnitudeError for a source, given the bandpass(es) and sky SED.

    For a given source, sky sed, total bandpass and hardware bandpass, as well as
    seeing / expTime, calculates the SNR with optimal PSF extraction
    assuming a double-gaussian PSF.

    @param [in] spectrum is an instantiation of the Sed class containing the SED of
    the object whose signal to noise ratio is being calculated

    @param [in] totalbandpass is an instantiation of the Bandpass class
    representing the total throughput (system + atmosphere)

    @param [in] skysed is an instantiation of the Sed class representing
    the sky emission per square arcsecond.

    @param [in] hardwarebandpass is an instantiation of the Bandpass class
    representing just the throughput of the system hardware.

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] seeing in arcseconds

    @param [in] verbose is a boolean

    @param [out] magnitude error
    """

    snr = calcSNR_sed(spectrum, totalbandpass, skysed, hardwarebandpass,
                      photParams, seeing, verbose=verbose)

    if photParams.sigmaSys is not None:
        return numpy.sqrt(numpy.power(magErrorFromSNR(snr),2) + numpy.power(photParams.sigmaSys,2))
    else:
        return magErrorFromSNR(snr)


def calcAstrometricError(mag, m5, nvisit=1):
    """
    Calculate the astrometric error, for object catalog purposes.

    Returns astrometric error for a given SNR, in mas.
    """
    # The astrometric error can be applied to parallax or proper motion (for nvisit>1).
    # If applying to proper motion, should also divide by the # of years of the survey.
    # This is also referenced in the astroph/0805.2366 paper.
    # D. Monet suggests sqrt(Nvisit/2) for first 3 years, sqrt(N) for longer, in reduction of error
    # because of the astrometric measurement method, the systematic and random error are both reduced.
    # Zeljko says 'be conservative', so removing this reduction for now.
    rgamma = 0.039
    xval = numpy.power(10, 0.4*(mag-m5))
    # The average seeing is 0.7" (or 700 mas).
    error_rand = 700.0 * numpy.sqrt((0.04-rgamma)*xval + rgamma*xval*xval)
    error_rand = error_rand / numpy.sqrt(nvisit)
    # The systematic error floor in astrometry:
    error_sys = 10.0
    # These next few lines are the code removed due to Zeljko's 'be conservative' requirement.
    #if (nvisit<30):
    #    error_sys = error_sys/numpy.sqrt(nvisit/2.0)
    #if (nvisit>30):
    #    error_sys = error_sys/numpy.sqrt(nvisit)
    astrom_error = numpy.sqrt(error_sys * error_sys + error_rand*error_rand)
    return astrom_error
