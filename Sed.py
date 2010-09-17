""" 
sed - 

$Id$

4/13/2010  ljones@astro.washington.edu

Class data: 
wavelen (nm)
flambda (ergs/cm^s/s/nm)
fnu (Jansky)
zp  (units of fnu = -8.9 (if Janskys) or 48.6 (ergs/cm^2/s/hz)

It is important to note the units are NANOMETERS, not ANGSTROMS. It is possible to rig this so you can
use angstroms instead of nm, but you should know what you're doing and understand the grid limits. 

Methods: 
 Because of how these methods will be applied for catalog generation, (taking one base SED and then 
  applying various dust extinctions and redshifts), many of the methods will either work on,
  and update self, OR they can be given a set of lambda/flambda arrays and then will return 
  new versions of these arrays. 
 In general, the philosophy of Sed.py is to not define the boundaries of the wavelength (and thus flambda)
  until it is necessary for magnitude calculation, or other gridded-required conversions (such as
  calculating fnu). At those times, the min/max/step wavelengths are defined and the entire SED is
  resampled onto that grid. The default values are 300-1200 nanometers however this is user-definable.
  Please note that if you define your own grid, you will have to define it EVERY time you call a method
  that takes these arguments. 

 When considering whether to use the internal wavelen/flambda (self) values, versus input values:
  For consistency, anytime self.wavelen/flambda is used, it will be updated if the values are changed
  (except in the special case of calculating magnitudes), and if self.wavelen/flambda is updated, 
  self.fnu will be recalculated or set to None.
  If arrays are passed into a method, they will not be altered and the arrays which are returned will be
  allocated new memory. 

Method include: 
  setSED / setFlatSED / readSED_flambda / readSED_fnu -- to input information into Sed wavelen/flambda.
  getSED_flambda / getSED_fnu -- to return wavelen / flambda or fnu to the user.
  clearSED -- set everything to 0.
  synchronizeSED -- to grid wavelen/flambda/fnu onto the desired grid and calculate fnu.
  checkUseSelf / needResample -- not expected to be useful to the user, rather intended for internal use.
  resampleSED -- primarily internal use, but may be useful to user. Resamples SED onto specified grid.
  flambdaTofnu / fnuToflambda -- conversion methods, will resample/grid flambda and fnu in the process.
  redshiftSED -- as it says. 
  setupCCMab / addCCMDust -- separated into two components, so that a_x/b_x can be reused between SEDS
if the wavelength range and grid is the same for each SED (calculate a_x/b_x with setupCCMab). 
  multiplySED -- multiply two SEDS together.
  calcADU / calcMag / calcFlux -- with a Bandpass, calculate the ADU/magnitude/flux of a SED.
  calcFluxNorm / multiplyFluxNorm -- handle fluxnorm parameters (from UW LSST database) properly. These 
methods are intended to give a user an easy way to scale an SED to match an expected magnitude. 
  renormalizeSED  -- intended for rescaling SEDS to a common flambda or fnu level. 
  writeSED -- keep a file record of your SED. 
  calcSNR_psf / calcSNR_mag -- two methods to calculate the SNR of a SED. (_psf is more accurate, but 
requires knowing the sky count backgrounds. _mag assumes you know the m5 already). 
  calcMagError / calcAstrometricError -- currently only very rough values.
  manyMagCalc -- given a list of bandpasses, this will return an array of magnitudes (in the same 
order as the bandpasses) of this SED in each of those bandpasses. 

"""

import warnings as warning
import numpy as n

# The following *wavelen* parameters are default values for gridding wavelen/sb/flambda.
MINWAVELEN = 300
MAXWAVELEN = 1200
WAVELENSTEP = 0.1

LIGHTSPEED = 299792458       # speed of light, = 2.9979e8 m/s
PLANCK = 6.626068e-27        # planck's constant, = 6.626068e-27 ergs*seconds
NM2M = 1.00e-9               # nanometers to meters conversion = 1e-9 m/nm
ERGSETC2JANSKY = 1.00e23     # erg/cm2/s/Hz to Jansky units (fnu)

EXPTIME = 15                      # Default exposure time. (option for method calls).
NEXP = 2                          # Default number of exposures. (option for methods).
EFFAREA = n.pi*(6.5*100/2.0)**2   # Default effective area of primary mirror. (option for methods).
GAIN = 2.3                        # Default gain. (option for method call).
RDNOISE = 5                       # Default value - readnoise electrons per pixel (per exposure)
DARKCURRENT = 0.2                 # Default value - dark current electrons per pixel per second
OTHERNOISE = 4.69                 # Default value - other noise electrons or adu per pixel per exposure
PLATESCALE = 0.2                  # Default value - "/pixel
SEEING = {'u': 0.77, 'g':0.73, 'r':0.70, 'i':0.67, 'z':0.65, 'y':0.63}  # Default seeing values (in ")


class Sed: 
    """Class for holding and utilizing spectral energy distributions (SEDs)"""
    def __init__(self, wavelen=None, flambda=None, fnu=None):
        """Initialize sed object by giving filename or lambda/flambda array.

        Note that this does *not* regrid flambda and leaves fnu undefined."""
        self.fnu = None
        self.wavelen = None
        self.flambda = None
        #self.zp = -8.9  # default units, Jansky.
        self.zp = -2.5*n.log10(3631)
        # If init was given data to initialize class, use it.
        if (wavelen!= None) & ((flambda!=None) | (fnu!=None)):
            self.setSED(wavelen, flambda=flambda, fnu=fnu)
        return

    ### Methods for getters and setters.

    def setSED(self, wavelen, flambda=None, fnu=None):
        """Populate wavelen/flambda fields in sed by giving lambda/flambda or lambda/fnu array.

        If flambda present, this overrides fnu. Method sets fnu=None unless only fnu is given.
        Sets wavelen/flambda or wavelen/flambda/fnu over wavelength array given. """
        # Check wavelen array for type matches.
        if isinstance(wavelen, n.ndarray)==False:
            raise ValueError("Wavelength must be a numpy array")
        # Wavelen type ok - make new copy of data for self.
        self.wavelen = n.copy(wavelen)
        self.flambda=None
        self.fnu=None
        # Check if given flambda or fnu.
        if flambda!=None:
            # Check flambda data type and length.
            if (isinstance(flambda, n.ndarray)==False) | (len(flambda)!=len(self.wavelen)):
                raise ValueError("Flambda must be a numpy array of same length as Wavelen.")
            # Flambda ok, make a new copy of data for self.
            self.flambda = n.copy(flambda)
        else:
            # Were passed fnu instead : check fnu data type and length.
            if fnu==None:
                raise ValueError("Both fnu and flambda are 'None', cannot set the SED.")
            elif (isinstance(fnu, n.ndarray)==False) | (len(fnu)!=len(self.wavelen)):
                raise ValueError("(No Flambda) - Fnu must be numpy array of same length as Wavelen.")
            # Convert fnu to flambda.
            self.wavelen, self.flambda = self.fnuToflambda(wavelen, fnu)
        return

    def setFlatSED(self, wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP):
        """Populate the wavelength/flambda/fnu fields in sed according to a flat fnu source."""
        self.wavelen = n.arange(wavelen_min, wavelen_max+wavelen_step, wavelen_step, dtype='float')
        self.fnu = n.ones(len(self.wavelen), dtype='float') * 3631 #jansky
        self.fnuToflambda()
        return

    def readSED_flambda(self, filename):
        """Read a file containing [lambda Flambda] (lambda in nm) (Flambda erg/cm^2/s/nm).
        
        Does not resample wavelen/flambda onto grid; leave fnu=None. """
        # Try to open data file.
        try:
            f=open(filename, 'r')
        except IOError:
            raise IOError("The file %s for this sed does not exist" %(filename))
        # Read source SED from file - lambda, flambda should be first two columns in the file.
        # lambda should be in nm and flambda should be in ergs/cm2/s/nm
        sourcewavelen = []
        sourceflambda = []
        for line in f:
            if line.startswith("#"):
                continue
            values = line.split()
            sourcewavelen.append(float(values[0]))
            sourceflambda.append(float(values[1]))
        f.close()
        self.wavelen = n.array(sourcewavelen)
        self.flambda = n.array(sourceflambda)
        self.fnu = None 
        return

    def readSED_fnu(self, filename):
        """Read a file containing [lambda Fnu] (lambda in nm) (Fnu in Jansky).

        Does not resample wavelen/fnu/flambda onto a grid; leaves fnu set."""
        # Try to open the data file.
        try:
            f=open(filename, 'r')
        except IOError:
            raise IOError("The file %f for this sed does not exist" %(filename))
        # Read source SED from file - lambda, fnu should be first two columns in the file.
        # lambda should be in nm and fnu should be in Jansky.
        sourcewavelen = []
        sourcefnu = []
        for line in f:
            if line.startswith("#"):
                continue
            values = line.split()
            sourcewavelen.append(float(values[0]))
            sourcefnu.append(float(values[1]))
        f.close()
        # Convert to numpy arrays.
        sourcewavelen = n.array(sourcewavelen)
        sourcefnu = n.array(sourcefnu)
        # Convert fnu to flambda (regrids wavelen/flambda in the process).
        self.wavelen, self.flambda = self.fnuToflambda(sourcewavelen, sourcefnu)
        return

    def getSED_flambda(self):
        """Return copy of wavelen/flambda."""
        # Get new memory copies of the arrays.
        wavelen = n.copy(self.wavelen)
        flambda = n.copy(self.flambda)
        return wavelen, flambda

    def getSED_fnu(self):
        """Return copy of wavelen/fnu, without altering self."""
        wavelen = n.copy(self.wavelen)
        # Check if fnu currently set.
        if fnu!=None:
            # Get new memory copy of fnu.
            fnu = n.copy(self.fnu)
        else:
            # Fnu was not set .. grab copy fnu without changing self. 
            wavelen, fnu = self.flambdaTofnu(self.wavelen, self.flambda)
            # Now wavelen/fnu (new mem) are gridded evenly, but self.wavelen/flambda/fnu remain unchanged.
        return wavelen, fnu

    ## Methods that update or change self.

    def clearSED(self):
        """Reset all data in sed to None."""
        self.wavelen = None
        self.fnu = None
        self.flambda = None
        self.zp = -8.9
        return
    
    def synchronizeSED(self, wavelen_min=None, wavelen_max=None, wavelen_step=None):
        """Set all wavelen/flambda/fnu values, potentially on min/max/step grid.

        Uses flambda to recalculate fnu. If wavelen min/max/step are given, resamples
        wavelength/flambda/fnu onto an even grid with these values. """
        # Grid wavelength/flambda/fnu if desired.
        if ((wavelen_min!=None) & (wavelen_max!=None) & (wavelen_step!=None)):
            self.resampleSED(wavelen_min=wavelen_min, wavelen_max=wavelen_max,
                             wavelen_step=wavelen_step)
        # Reset or set fnu.
        self.flambdaTofnu()
        return
        
    ## Utilities common to several later methods.

    def checkUseSelf(self, wavelen, flux):
        """Simple utility to check if should be using self's data or passed arrays.
        
        Also does data integrity check on wavelen/flux if not self."""
        update_self = False
        if (wavelen==None) | (flux==None): 
            # Then one of the arrays was not passed - check if this is true for both arrays.
            if (wavelen!=None) | (flux!=None):
                # Then one of the arrays was passed - raise exception.
                raise ValueError("Must either pass *both* wavelen/flux pair, or use defaults.")
            update_self = True
        else:
            # Both of the arrays were passed in - check their validity. 
            if (isinstance(wavelen, n.ndarray)==False) | (isinstance(flux, n.ndarray)==False):
                raise ValueError("Must pass wavelen/flux as numpy arrays.")
            if len(wavelen)!=len(flux): 
                raise ValueError("Must pass equal length wavelen/flux arrays.")
        return update_self

    def needResample(self, wavelen_match=None, wavelen=None, 
                     wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP):
        """Check if wavelen or self.wavelen matches wavelen or wavelen_min/max/step grid."""
        # Check if should use self or passed wavelen.
        if wavelen==None:
            wavelen = self.wavelen
        # Check if wavelength arrays are equal, if wavelen_match passed. 
        if wavelen_match != None:
            if n.shape(wavelen_match) != n.shape(wavelen):
                need_regrid=True
            else:
                # check the elements to see if any vary
                need_regrid = n.any(abs(wavelen_match-wavelen)>1e-10)
        else:
            need_regrid = True
            # Check if wavelen_min/max/step are set - if ==None, then return (no regridding).
            # It's possible (writeSED) to call this routine, even with no final grid in mind.
            if ((wavelen_min == None) & (wavelen_max == None) & (wavelen_step==None)):
                need_regrid = False
                return need_regrid
            # Okay, now look at comparison of wavelen to the grid.
            wavelen_max_in = wavelen[len(wavelen)-1]
            wavelen_min_in = wavelen[0]
            wavelen_step_in = wavelen[1]-wavelen[0]
            # First check match to minimum/maximum and first step in array.
            if ((wavelen_min_in == wavelen_min) & (wavelen_max_in == wavelen_max)
                & (wavelen_step_in == wavelen_step)):
                # Then do a check to see if number of elements consistent with even step size.
                if len(wavelen) == (wavelen_max_in + wavelen_step_in - wavelen_min_in)/wavelen_step_in:
                    # Then now decent chance data is gridded properly.
                    need_regrid = False
        # At this point, need_grid=True unless it's proven to be False, so return value.
        return need_regrid
                    
            
    
    def resampleSED(self, wavelen=None, flux=None, wavelen_match=None,
                    wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP):
        """Resample flux onto grid defined by min/max/step OR another wavelength array.

        Give method wavelen/flux OR default to self.wavelen/self.flambda.
        Method either returns wavelen/flambda (if given those arrays) or updates wavelen/flambda in self. 
         If updating self, resets fnu to None. """
        # Is method acting on self.wavelen/flambda or passed in wavelen/flux arrays? Sort it out.
        update_self=self.checkUseSelf(wavelen, flux)
        if update_self:
            wavelen=self.wavelen
            flux=self.flambda
            self.fnu = None
        # Now, on with the resampling. 
        # The user should check if need this routine by trying needResample first.
        # In here, resampling will be done regardless if necessary or not.
        #  (this simplifies memory management, as if you call this funtion, you will be 
        #   getting new copies of data). 
        # Set up gridded wavelength or copy of wavelen array to match.
        if wavelen_match == None:
            wavelen_grid = n.arange(wavelen_min, wavelen_max+wavelen_step,
                                    wavelen_step, dtype='float')
        else:
            wavelen_grid = n.copy(wavelen_match)
        flux_grid = n.empty(len(wavelen), dtype='float')
        # Do the interpolation of wavelen/flux onto grid. (type/len failures will die here).
        flux_grid = n.interp(wavelen_grid, wavelen, flux, left=0.0, right=0.0)
        # Update self values if necessary.
        if update_self:
            self.wavelen = wavelen_grid
            self.flambda = flux_grid
        return wavelen_grid, flux_grid

    def flambdaTofnu(self, wavelen=None, flambda=None):
        """Convert flambda into fnu. 

        This routine assumes that flambda is in ergs/cm^s/s/nm and produces fnu in Jansky.
        Can act on self or user can provide wavelen/flambda and get back wavelen/fnu """
        # Change Flamda to Fnu by multiplying Flambda * lambda^2 = Fv
        # Fv dv = Fl dl .. Fv = Fl dl / dv = Fl dl / (dl*c/l/l) = Fl*l*l/c
        # Check - Is the method acting on self.wavelen/flambda/fnu or passed wavelen/flambda arrays? 
        update_self = self.checkUseSelf(wavelen, flambda)
        if update_self:
            wavelen = self.wavelen
            flambda=self.flambda
            self.fnu = None
        # Now on with the calculation. 
        # Calculate fnu.
        fnu = flambda * wavelen * wavelen * NM2M / LIGHTSPEED
        fnu = fnu * ERGSETC2JANSKY
        # If are using/updating self, then *all* wavelen/flambda/fnu will be gridded. 
        # This is so wavelen/fnu AND wavelen/flambda can be kept in sync.
        if update_self:
            self.wavelen = wavelen
            self.flambda = flambda
            self.fnu = fnu
        # Return wavelen, fnu.
        return wavelen, fnu

    def fnuToflambda(self, wavelen=None, fnu=None):
        """Convert fnu into flambda.
        
        Assumes fnu in units of Jansky and flambda in ergs/cm^s/s/nm.
        Can act on self or user can give wavelen/fnu and get wavelen/flambda returned"""
        # Fv dv = Fl dl .. Fv = Fl dl / dv = Fl dl / (dl*c/l/l) = Fl*l*l/c
        # Is method acting on self or passed arrays?
        update_self = self.checkUseSelf(wavelen, fnu)
        if update_self:
            wavelen = self.wavelen
            fnu = self.fnu
        # On with the calculation.
        # Calculate flambda.
        flambda = fnu / wavelen / wavelen * LIGHTSPEED / NM2M
        flambda = flambda / ERGSETC2JANSKY
        # If updating self, then *all of wavelen/fnu/flambda will be gridded and updated,
        # this is so wavelen/fnu AND wavelen/flambda can be kept in sync.
        if update_self: 
            self.wavelen = wavelen
            self.flambda = flambda
            self.fnu = fnu
        # Return wavelen/flambda.
        return wavelen, flambda

    ## methods to alter the sed

    def redshiftSED(self, redshift, dimming=False, wavelen=None, flambda=None):
        """Redshift an SED, optionally adding cosmological dimming. 
        
        Pass wavelen/flambda or redshift/update self.wavelen/flambda (unsets fnu)"""
        # Updating self or passed arrays?
        update_self = self.checkUseSelf(wavelen, flambda)
        if update_self:
            wavelen=self.wavelen
            flambda = self.flambda
            self.fnu = None
        else:
            # Make a copy of input data, because will change its values.
            wavelen = n.copy(wavelen)
            flambda = n.copy(flambda)
        # Okay, move onto redshifting the wavelen/flambda pair.
        wavelen = wavelen * (1+redshift)
        # Flambda now just has different wavelength for each value.
        # Add cosmological dimming if required.
        if dimming:
            flambda = flambda / (1+redshift)
        # Update self, if required - but just flambda (still no grid required).
        if update_self:
            self.wavelen = wavelen
            self.flambda = flambda
        return wavelen, flambda

    def setupCCMab(self, wavelen=None):
        """Calculate a(x) and b(x) for CCM dust model. (x=1/wavelen).
        
        Returns a(x) and b(x) can be common to many seds, wavelen is the same. """
        # This extinction law taken from Cardelli, Clayton and Mathis ApJ 1989.
        # The general form is A_l / A(V) = a(x) + b(x)/R_V  (where x=1/lambda in microns),
        # then different values for a(x) and b(x) depending on wavelength regime.
        # Also, the extinction is parametrized as R_v = A_v / E(B-V).
        # Magnitudes of extinction (A_l) translates to flux by a_l = -2.5log(f_red / f_nonred).
        if wavelen == None:
            wavelen = n.copy(self.wavelen)
        a_x = n.zeros(len(wavelen), dtype='float')
        b_x = n.zeros(len(wavelen), dtype='float')
        # Convert wavelength to x (in inverse microns).
        x = n.empty(len(wavelen), dtype=float)
        nm_to_micron = 1/1000.0
        x = 1.0 / (wavelen * nm_to_micron)  
        # Dust in infrared 0.3 /mu < x < 1.1 /mu (inverse microns).
        condition = (x>=0.3) & (x<=1.1)
        if len(a_x[condition]) > 0 :
            y = x[condition]
            a_x[condition] = 0.574 * y**1.61
            b_x[condition] = -0.527 * y**1.61
        # Dust in optical/NIR 1.1 /mu < x < 3.3 /mu region.
        condition = (x >=1.1) & (x<=3.3)
        if len(a_x[condition])>0:
            y = x[condition] - 1.82
            a_x[condition] = 1 + 0.104*y - 0.609*y**2 + 0.701*y**3 + 1.137*y**4 
            a_x[condition] = a_x[condition] - 1.718*y**5 - 0.827*y**6 + 1.647*y**7 - 0.505*y**8
            b_x[condition] = 1.952*y + 2.908*y**2 - 3.989*y**3 - 7.985*y**4
            b_x[condition] = b_x[condition] + 11.102*y**5 + 5.491*y**6 - 10.805*y**7 + 3.347*y**8
        # Dust in ultraviolet and UV (if needed for high-z) 3.3 /mu< x< 8 /mu.
        condition = (x>=3.3) & (x<5.9)
        if len(a_x[condition])>0:
            y = x[condition]
            a_x[condition] = 1.752 - 0.316*y - 0.104/((y-4.67)**2 + 0.341)
            b_x[condition] = -3.090 + 1.825*y + 1.206/((y-4.62)**2 + 0.263)
        condition = (x>5.9) & (x<8)
        if len(a_x[condition])>0:
            y = x[condition]
            Fa_x = n.empty(len(a_x[condition]), dtype=float)
            Fb_x = n.empty(len(a_x[condition]), dtype=float)
            Fa_x = -0.04473*(y-5.9)**2 - 0.009779*(y-5.9)**3
            Fb_x = 0.2130*(y-5.9)**2 + 0.1207*(y-5.9)**3
            a_x[condition] = 1.752 - 0.316*y - 0.104/((y-4.67)**2 + 0.341) + Fa_x
            b_x[condition] = -3.090 + 1.825*y + 1.206/((y-4.62)**2 + 0.263) + Fb_x
        # Dust in far UV (if needed for high-z) 8 /mu < x < 10 /mu region.
        condition = (x >= 8) & (x<= 11.)
        if len(a_x[condition])>0: 
            y = x[condition]-8.0
            a_x[condition] = -1.073 - 0.628*(y) + 0.137*(y)**2 - 0.070*(y)**3
            b_x[condition] = 13.670 + 4.257*(y) - 0.420*(y)**2 + 0.374*(y)**3
        return a_x, b_x

    def addCCMDust(self, a_x, b_x, A_v=None, ebv=None, R_v=3.1, wavelen=None, flambda=None):
        """Add CCM dust model extinction to the SED, modifying flambda and fnu.
        
        Specify any two of A_V, E(B-V) or R_V (=3.1 default) """
        # The extinction law taken from Cardelli, Clayton and Mathis ApJ 1989.
        # The general form is A_l / A(V) = a(x) + b(x)/R_V  (where x=1/lambda in microns).
        # Then, different values for a(x) and b(x) depending on wavelength regime.
        # Also, the extinction is parametrized as R_v = A_v / E(B-V).
        # The magnitudes of extinction (A_l) translates to flux by a_l = -2.5log(f_red / f_nonred).
        #
        # Figure out if updating self or passed arrays.
        update_self = self.checkUseSelf(wavelen, flambda)
        if update_self:
            wavelen = self.wavelen
            flambda = self.flambda
            self.fnu = None
        else:
            wavelen = n.copy(wavelen)
            flambda = n.copy(flambda)
        # Input parameters for reddening can include any of 3 parameters; only 2 are independent.
        # Figure out what parameters were given, and see if self-consistent.
        if R_v == 3.1:
            if A_v == None:
                A_v = R_v * ebv
            elif (A_v != None) & (ebv != None): 
                # Specified A_v and ebv, so R_v should be nondefault.
                R_v = A_v / ebv
        if (R_v != 3.1):
            if (A_v != None) & (ebv != None): 
                calcRv = A_v / ebv
                if calcRv != R_v: 
                    raise ValueError("CCM parametrization expects R_v = A_v / E(B-V);",
                                     "Please check input values, because values are inconsistent.")
            elif A_v == None:
                A_v = R_v * ebv
        # R_v and A_v values are specified or calculated.
        A_lambda = n.empty(len(wavelen), dtype=float)
        dust = n.empty(len(wavelen), dtype=float)
        A_lambda = (a_x + b_x / R_v) * A_v
        # dmag_red(dust) = -2.5 log10 (f_red / f_nored) : (f_red / f_nored) = 10**-0.4*dmag_red
        dust = n.power(10.0, -0.4*A_lambda)
        flambda = flambda * dust
        # Update self if required.
        if update_self:
            self.flambda = flambda 
        return wavelen, flambda

    def multiplySED(self, other_sed,
                    wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP):
        """Multiply two SEDs together - flambda * flambda - and return a new sed object.
        
        Does not alter self or other_sed"""
        # Set up wavelen/flambda of first object, on grid.
        wavelen_1 = self.wavelen
        flambda_1 = self.flambda
        if self.needResample(wavelen=wavelen_1, wavelen_min=wavelen_min,
                             wavelen_max=wavelen_max, wavelen_step=wavelen_step):
            wavelen_1, flambda_1 = self.resampleSED(wavelen_1, flambda_1,
                                                    wavelen_min=wavelen_min,
                                                    wavelen_max=wavelen_max,
                                                    wavelen_step=wavelen_step)
        # Set up wavelen/flambda of second object, on grid.  
        wavelen_2 = other_sed.wavelen
        flambda_2 = other_sed.flambda
        if self.needResample(wavelen=wavelen_2, wavelen_min=wavelen_min,
                             wavelen_max=wavelen_max, wavelen_step=wavelen_step):
            wavelen_2, flambda_2 = self.resampleSED(wavelen_2, flambda_2,
                                                    wavelen_min=wavelen_min, 
                                                    wavelen_max=wavelen_max,
                                                    wavelen_step=wavelen_step)
        # Multiply the two flambda together.
        wavelen = wavelen_1
        flambda = flambda_1 * flambda_2
        # Instantiate new sed object.
        new_sed = Sed(wavelen, flambda)
        return new_sed

    ## routines related to magnitudes

    def calcADU(self, bandpass, wavelen=None, fnu=None,
                expTime=EXPTIME, effarea=EFFAREA, gain=GAIN):
        """Calculate the number of adu from camera, using sb and fnu.

        Given wavelen/fnu arrays or use self. Passed wavelen/fnu arrays will be unchanged. 
         If method uses self, checks if fnu set; if not, does *not* permanently set fnu.
        Calculating the AB mag requires the wavelen/fnu pair to be on the same grid as bandpass; 
         (temporary values of these are used). """
        update_self = self.checkUseSelf(wavelen, fnu)
        if update_self:
            # Calculate fnu if required.
            if self.fnu == None:
                # If fnu not present, create temporary copy - on bandpass grid.
                wavelen, fnu = self.flambdaTofnu(self.wavelen, self.flambda)
            else:
                wavelen = self.wavelen
                fnu = self.fnu
        # Check bandpass/fnu (even if not self) are on same grid.
        if self.needResample(wavelen=wavelen, wavelen_match=bandpass.wavelen):
            # Here, not on the same grid so resample to match wavelen/fnu to bandpass.
            # Note that resampleSED allocates new memory for wavelen/fnu return values.
            wavelen, fnu = self.resampleSED(wavelen, fnu, wavelen_match=bandpass.wavelen)
        # Calculate the number of photons.
        dlambda = wavelen[1] - wavelen[0]
        # Nphoton in units of 10^-23 ergs/cm^s/nm. 
        nphoton = (fnu / wavelen * bandpass.sb).sum()
        adu = nphoton * (expTime * effarea/gain) * (1/ERGSETC2JANSKY) * (1/PLANCK) * dlambda 
        return adu


    def calcMag(self, bandpass, wavelen=None, fnu=None):
        """Calculate the AB magnitude of an object, using phi the normalized system response.

        Can pass wavelen/fnu arrays or use self. Passed wavelen/fnu arrays will be unchanged.
        If method uses self, checks if fnu set; if not, does *not* permanently set fnu. 
        Calculating the AB mag requires the wavelen/fnu pair to be on the same grid as bandpass; 
         (temporary values of these are used)"""
        # Note - the behavior in this first section might be considered a little odd. 
        # However, I felt calculating a magnitude should not (unexpectedly) regrid your 
        # wavelen/flambda information if you were using self., as this is not obvious from the "outside".
        # To preserve 'user logic', the wavelen/flambda of self are left untouched. Unfortunately 
        # this means, this method can be used inefficiently if calculating many magnitudes with
        # the same sed and same bandpass region - in that case, use self.flambdaTofnu() first 
        # on all of the seds in question, then hop to calculating magnitudes much more efficiently! 
        # (also, if bandpass and mag are using the same grid, this is more efficient). 
        update_self = self.checkUseSelf(wavelen, fnu)
        if update_self:
            # Calculate fnu if required.
            if self.fnu == None:
                # If fnu is not present, create temporary copy on bandpass grid.
                wavelen, fnu = self.flambdaTofnu(self.wavelen, self.flambda)
            else:
                wavelen = self.wavelen
                fnu = self.fnu
        # Continue with magnitude calculation.
        # Check if bandpass and wavelen/fnu are on the same grid. 
        if self.needResample(wavelen=wavelen, wavelen_match=bandpass.wavelen):                             
            # Here - not on the same grid, so resample to match wavelen/fnu to bandpass.
            # Note that resampleSED allocates new memory for wavelen/fnu return values.
            wavelen, fnu = self.resampleSED(wavelen, fnu, wavelen_match=bandpass.wavelen)
        # Calculate bandpass phi value if required.
        if bandpass.phi == None:
            bandpass.sbTophi()
        # Calculate flux in bandpass and then AB magnitude.
        dlambda = wavelen[1] - wavelen[0]
        mag = -2.5 * n.log10((fnu*bandpass.phi).sum()*dlambda) - self.zp
        return mag

    def calcFlux(self, bandpass, wavelen=None, fnu=None):
        """Calculate the F_b (integrated flux of an object, above the atmosphere), using phi.
        
        Passed wavelen/fnu arrays will be unchanged, but if uses self will check if fnu is set.
           If fnu not set, does *not* permanently set fnu. 
        Calculating the AB mag requires the wavelen/fnu pair to be on the same grid as bandpass; 
           (temporary values of these are used)"""
        update_self = self.checkUseSelf(wavelen, fnu)
        if update_self:
            # Calculate fnu if required.
            if self.fnu == None:
                # If fnu not present, create temporary copy - on bandpass grid.
                wavelen, fnu = self.flambdaTofnu(self.wavelen, self.flambda) 
            else:
                wavelen = self.wavelen
                fnu = self.fnu
        # Go on with magnitude calculation.
        # Check bandpass and wavelen/fnu are on the same grid.
        if self.needResample(wavelen=wavelen, wavelen_match=bandpass.wavelen):
            # Here - not on the same grid so resample to match wavelen/fnu to bandpass.
            # Note that resampleSED allocates new memory for wavelen/fnu return values.
            wavelen, fnu = self.resampleSED(wavelen, fnu, wavelen_match=bandpass.wavelen)
        # Calculate bandpass phi value if required.
        if bandpass.phi == None:
            bandpass.sbTophi()
        # Calculate flux in bandpass and return this value.
        dlambda = wavelen[1] - wavelen[0]        
        flux = (fnu*bandpass.phi).sum() * dlambda
        return flux

    def calcFluxNorm(self, magmatch, bandpass, wavelen=None, fnu=None,
                     wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP):
        """Calculate the fluxNorm (SED normalization value for a given mag) for a sed.
        
        Equivalent to adjusting a particular f_nu to Jansky's appropriate for the desired mag.
        Can pass wavelen/fnu or apply to self. Requires a gridded wavelen/fnu on min/max/step grid.
        Note that calcFluxNorm does not regrid self.wavelen/flambda/fnu permanently, so need
        to use the same grid here as you use in multiplyFluxNorm. """
        update_self = self.checkUseSelf(wavelen, fnu)
        if update_self:
            wavelen = self.wavelen
            fnu = self.fnu
            # Check possibility that fnu is not calculated yet.
            if fnu==None:
                # This only temporarily calculates fnu, so we don't permanently resample self here.
                wavelen, fnu = self.flambdaTofnu(wavelen=self.wavelen, flambda=self.flambda) 
        # Fluxnorm gets applied to f_nu (fluxnorm * SED(f_nu) * PHI = mag - 8.9 (AB zeropoint).
        # FluxNorm * SED => correct magnitudes for this object.
        # check if wavelen/fnu are on same grid as bandpass
        if self.needResample(wavelen=wavelen, wavelen_match=bandpass.wavelen):
            # Here - not on the same grid so resample to match wavelen/fnu to bandpass, 
            #  but don't store in self. 
            # Note that resampleSED allocates new memory for wavelen/fnu return values.
            wavelen, fnu = self.resampleSED(wavelen, fnu, wavelen_match=bandpass.wavelen)
        # Calculate fluxnorm. 
        dmag = magmatch - self.calcMag(bandpass, wavelen, fnu)
        fluxnorm = n.power(10, (-0.4*dmag))  
        return fluxnorm 
   
    def multiplyFluxNorm(self, fluxNorm, wavelen=None, fnu=None,
                         wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP):
        """Multiply wavelen/fnu (or self.wavelen/fnu) by fluxnorm.
        
        Returns wavelen/fnu arrays (or updates self). 
        Note that multiplyFluxNorm *does* regrid self.wavelen/flambda/fnu permanently""" 
        # Note that fluxNorm is intended to be applied to f_nu,
        # so that fluxnorm*fnu*phi = mag (expected magnitude).
        update_self = self.checkUseSelf(wavelen, fnu)
        if update_self:
            # Make sure fnu is defined.
            if self.fnu==None:
                self.flambdaTofnu()
            wavelen = self.wavelen
            fnu = self.fnu
            # Make sure on desired grid. 
            if self.needResample(wavelen, wavelen_min=wavelen_min,
                                 wavelen_max=wavelen_max, wavelen_step=wavelen_step):
                wavelen, fnu = self.resampleSED(wavelen, fnu, wavelen_min=wavelen_min,
                                                wavelen_max=wavelen_max,
                                                wavelen_step=wavelen_step)
        else:
            # Check wavelen/fnu are on desired grid:
            if self.needResample(wavelen, wavelen_min=wavelen_min,
                                 wavelen_max=wavelen_max, wavelen_step=wavelen_step):
                # Automatically get new copy of data here.
                wavelen, fnu = self.resampleSED(wavelen, fnu, wavelen_min=wavelen_min,
                                                wavelen_max=wavelen_max,
                                                wavelen_step=wavelen_step)
            else:
                # Require new copy of the data.
                wavelen = n.copy(wavelen)
                fnu = n.copy(fnu)
        # Apply fluxnorm.
        fnu = fnu * fluxNorm
        # Update self.
        if update_self:
            self.wavelen = wavelen
            self.fnu = fnu
            self.fnuToflambda()
        return wavelen, fnu

    def renormalizeSED(self, lambdanorm=500, normvalue=1, gap=0, normflux='flambda',
                       wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP):
        """Renormalize sed in flambda to have normflux=normvalue @ lambdanorm or averaged over gap.
        
        Can normalized in flambda or fnu values.
        Return a new sed object; only  uses self.wavelen/flambda as input wavelen/flambda"""
        # Normalizes the fnu/flambda SED at one wavelength or average value over small range (gap).
        # This is useful for generating SED catalogs, mostly, to make them match schema.
        # Do not use this for calculating specific magnitudes -- use calcfluxNorm and multiplyFluxNorm.
        wavelen = n.copy(self.wavelen)
        flambda = n.copy(self.flambda)
        # Start normalizing wavelen/flambda.
        if normflux=='flambda':
            if self.needResample(wavelen, wavelen_min, wavelen_max, wavelen_step):
                wavelen, flambda = self.resampleSED(wavelen, flambda, 
                                                    wavelen_min=wavelen_min,
                                                    wavelen_max=wavelen_max,
                                                    wavelen_step=wavelen_step)
            # "standard" schema have flambda = 1 at 500 nm
            if gap==0:
                lambdapt = n.arange(lambdanorm-wavelen_step/2.0, lambdanorm+wavelen_step/2.0, wavelen_step, dtype=float)
                flambda_atpt = n.zeros(len(lambdapt), dtype='float')
                flambda_atpt = n.interp(lambdapt, wavelen, flambda, left=None, right=None)
                gapval = flambda_atpt[0]
            else:
                lambdapt = n.arange(lambdanorm-gap, lambdanorm+gap, wavelen_step, dtype=float)
                flambda_atpt = n.zeros(len(lambdapt), dtype='float')
                flambda_atpt = n.interp(lambdapt, wavelen, flambda, left=None, right=None)
                gapval = flambda_atpt.sum()/len(lambdapt)
            # Now renormalize fnu and flambda in the case of normalizing flambda.
            konst = normvalue/gapval
            flambda = flambda * konst
        if normflux=='fnu':  
            # Get fnu and make sure on grid, in one step. 
            wavelen, fnu = self.flambdaTofnu(wavelen, flambda)
            if gap==0:
                lambdapt = n.arange(lambdanorm, lambdanorm+wavelen_step, wavelen_step, dtype=float)
                fnu_atpt = n.zeros(len(lambdapt), dtype='float')
                fnu_atpt = n.interp(lambdapt, wavelen, fnu, left=None, right=None)
                gapval = fnu_atpt[0]
            else:
                lambdapt = n.arange(lambdanorm-gap, lambdanorm+gap, wavelen_step, dtype=float)
                fnu_atpt = n.zeros(len(lambdapt), dtype='float')
                fnu_atpt = n.interp(lambdapt, wavelen, fnu, left=None, right=None)
                gapval = fnu_atpt.sum()/len(lambdapt)
            # Now renormalize fnu and flambda in the case of normalizing fnu.
            konst = normvalue/gapval
            fnu = fnu * konst
            flambda = self.fnutoflambda(wavelen,fnu)
        new_sed = Sed(wavelen=wavelen, flambda=flambda)
        return new_sed


    def writeSED(self, filename, print_header=None, print_fnu=False, 
                 wavelen_min=None, wavelen_max=None, wavelen_step=None):
        """Write SED (wavelen, flambda, optional fnu) out to file.
        
        Option of adding a header line (such as version info) to output file.
        Does not alter self, regardless of grid or presence/absence of fnu"""
        # This can be useful for debugging or recording an SED.
        f = open(filename, 'w')
        wavelen = self.wavelen
        flambda = self.flambda
        # See if need to regrid data (if regrid, new memory copy).
        if self.needResample(wavelen=wavelen, wavelen_min=wavelen_min, 
                             wavelen_max=wavelen_max, wavelen_step=wavelen_step):
            wavelen, flambda = self.resampleSED(wavelen, flambda, wavelen_min=wavelen_min,
                                                wavelen_max=wavelen_max,
                                                wavelen_step=wavelen_step)
        # Then just use this gridded wavelen/flambda to calculate fnu.
        # Print header.
        if print_header != None:
            print >>f, "#", print_header
        # Print standard header info.
        if print_fnu:
            wavelen, fnu = self.flambdaTofnu(wavelen, flambda)
            print >>f, "# Wavelength(nm)  Flambda(ergs/cm^s/s/nm)   Fnu(Jansky)"
        else:
            print >>f, "# Wavelength(nm)  Flambda(ergs/cm^s/s/nm)"
        for i in range(0, len(wavelen), 1):
            if print_fnu:
                fnu = self.flambdaTofnu(wavelen=wavelen, flambda=flambda)
                print >> f, wavelen[i], flambda[i], fnu[i]
            else:
                #print >> f, self.wavelen[i], self.flambda[i]
                print >> f, "%.2f %.7g" %(wavelen[i], flambda[i])
        # Done writing, close file.
        f.close()       
        return

    def calcSNR_psf(self, totalbandpass, skysed, hardwarebandpass, 
                    readnoise=RDNOISE, darkcurrent=DARKCURRENT, 
                    othernoise=OTHERNOISE, seeing=SEEING['r'], 
                    effarea=EFFAREA, expTime=EXPTIME, nexp=1, 
                    platescale=PLATESCALE, gain=1, verbose=False):
        """Calculate the signal to noise ratio for a source, given the bandpass(es) and sky SED.
        
        For a given source, sky sed, total bandpass and hardware bandpass, as well as 
        seeing / expTime, calculates the SNR with optimal PSF extraction 
        assuming a double-gaussian PSF. Assumes that all values (readnoise/othernoise
        /darkcurrent) are given in appropriate units for gain. (gain=1 -> values are in electrons,
        while if gain!=1 then values given should be in adu, including readnoise/darkcurrent)."""
        # Calculate the counts from the source.
        sourcecounts = self.calcADU(totalbandpass, expTime=expTime*nexp, effarea=effarea, gain=gain)
        # Calculate the counts from the sky.
        skycounts = (skysed.calcADU(hardwarebandpass, expTime=expTime*nexp, effarea=effarea, gain=gain)
                     * platescale * platescale)
        # Calculate the effective number of pixels for double-Gaussian PSF. 
        neff = 2.436*(seeing/platescale)**2
        # Calculate the (square of the) noise due to instrumental effects.
        # Include the readout noise twice because 30 seconds is two exposures.
        noise_instr_sq = nexp*readnoise**2 + darkcurrent*expTime*nexp + nexp*othernoise**2
        # Calculate the (square of the) noise due to sky background poisson noise.
        noise_sky_sq = skycounts/gain
        # Discount error in sky measurement for now.
        noise_skymeasurement_sq = 0
        # Calculate the (square of the) noise due to signal poisson noise.
        noise_source_sq = sourcecounts/gain
        # Calculate total noise
        noise = n.sqrt(noise_source_sq + neff*(noise_sky_sq+noise_instr_sq+noise_skymeasurement_sq))
        # Calculate the signal to noise ratio.
        snr = sourcecounts / noise
        if verbose:
            print "For Nexp %.1f of time %.1f: " % (nexp, expTime)
            print "Counts from source: %.2f  Counts from sky: %.2f" %(sourcecounts, skycounts)
            print "Seeing: %.2f('')  Neff pixels: %.3f(pix)" %(seeing, neff)
            print "Noise from sky: %.2f Noise from instrument: %.2f" \
                %(n.sqrt(noise_sky_sq), n.sqrt(noise_instr_sq))
            print "Noise from source: %.2f" %(n.sqrt(noise_source_sq))
            print " Total Signal: %.2f   Total Noise: %.2f    SNR: %.2f" %(sourcecounts, noise, snr)
            # Return the signal to noise value.
        return snr

    ## Below here, methods that are appropriate for sed, but don't require SED.

    def calcSNR_mag(self, mag, m5):
        """Calculate the signal to noise of an object, given only the 5-sigma limiting mag"""  
        flux_ratio = n.power(10, 0.4*(m5-mag))
        snr_obj = 5 * (flux_ratio)
        return snr_obj

    def calcMagError(self, mag, m5, nvisit=1):
        """Calculate the photometric error, for object catalog purposes.
        
        Returns photometric error for a given SNR"""
        # Calculate mag errors - see eqns 4-6 astroph/0805.2366.
        # Can either do single image (nvisit=1) or coadded (nvisit>1).
        # m5 should be the 5-sigma limiting magnitude of one image.
        error_sys = 0.005     # systematic error - 0.005
        if (nvisit>50):
            # Allowing for reduction in systematic error after many visits
            error_sys = 0.0025 
        # Rgamma value comes from overview paper.
        rgamma = 0.039        
        xval = n.power(10, 0.4*(mag-m5))
        error_rand = n.sqrt((0.04-rgamma)*xval + rgamma*xval*xval)/n.sqrt(nvisit)
        mag_error = n.sqrt(error_sys*error_sys + error_rand*error_rand)     
        return mag_error

    def calcAstrometricError(self, mag, m5, nvisit=1):
        """Calculate the astrometric error, for object catalog purposes.
        
        Returns astrometric error for a given SNR, in mas."""
        # The astrometric error can be applied to parallax or proper motion (for nvisit>1).
        # If applying to proper motion, should also divide by the # of years of the survey.
        # This is also referenced in the astroph/0805.2366 paper.
        # D. Monet suggests sqrt(Nvisit/2) for first 3 years, sqrt(N) for longer, in reduction of error
        # because of the astrometric measurement method, the systematic and random error are both reduced.
        # Zeljko says 'be conservative', so removing this reduction for now.
        rgamma = 0.039
        xval = n.power(10, 0.4*(mag-m5))
        # The average seeing is 0.7" (or 700 mas).
        error_rand = 700.0 * n.sqrt((0.04-rgamma)*xval + rgamma*xval*xval)
        error_rand = error_rand / n.sqrt(nvisit)
        # The systematic error floor in astrometry: 
        error_sys = 10.0
        # These next few lines are the code removed due to Zeljko's 'be conservative' requirement.
        #if (nvisit<30):
        #    error_sys = error_sys/n.sqrt(nvisit/2.0)
        #if (nvisit>30):
        #    error_sys = error_sys/n.sqrt(nvisit)
        astrom_error = n.sqrt(error_sys * error_sys + error_rand*error_rand)
        return astrom_error


## Bonus, many-magnitude calculation for many SEDs with a single bandpass
    
    def manyMagCalc(self, bandpasslist):
        """Calculate many magnitudes for many bandpasses using a single sed.

        This is LESS STABLE than calculating the magnitude independently, as it
        takes several things for granted. For example, it assumes that each SED is
        resampled onto the same wavelength array and that fnu has been calculated.
        Also that each bandpass in the bandpasslist has been sampled onto
        the same wavelength array and already has phi calculated."""
        dlambda = bandpasslist[0].wavelen[1] - bandpasslist[0].wavelen[0]
        # Calculate phis and resample onto same wavelength grid
        phi = n.empty((len(bandpasslist), len(bandpasslist[0].phi)), dtype='float')
        mags = n.empty(len(bandpasslist), dtype='float')
        i = 0
        for bandpass in bandpasslist:
            phi[i] = bandpass.phi
            i = i+1
        mags = -2.5*n.log10(n.sum(phi*self.fnu, axis=1)*dlambda) - self.zp
        return mags
