import astropy.cosmology as cosmology
import astropy.units as units

__all__ = ["CosmologyWrapper"]

class CosmologyWrapper(object):

    cosmologyInitialized = False

    def set_units(self):
        """
        This method specifies the units in which various outputs from the wrapper are expected
        (this is because the latest version of astropy.cosmology outputs quantities such as
        the Hubble parameter and luminosity distance with units attached; the version of
        astropy.cosmology that comes within anaconda does not do this as of 30 October 2014)
        """
        if not self.cosmologyInitialized:
            raise RuntimeError("Cannot call set_units; cosmology is not initialized")

        H = self.activeCosmology.H(0.0)
        if 'unit' in dir(H):
            self.hUnits = units.Unit("km / (Mpc s)")
        else:
            self.hunits = None

        dd = self.activeCosmology.comoving_distance(0.0)
        if 'unit' in dir(dd):
            self.distanceUnits = units.Mpc
        else:
            self.distanceUnits = None

        mm = self.activeCosmology.distmod(1.0)
        if 'unit' in dir(mm):
            self.modulusUnits = units.mag
        else:
            self.modulusUnits = None

    def set_current(self, universe):
        """
        Take the cosmology indicated by 'universe' and set it as the current/default
        cosmology (depending on the API of the version of astropy being run)

        universe is also assigned to self.activeCosmology, which is the cosmology that
        this wrapper's methods use for calculations.
        """

        if 'set_current' in dir(cosmology):
            cosmology.set_current(universe)
        elif 'default_cosmology' in dir(cosmology):
            cosmology.default_cosmology.set(universe)
        else:
            raise RuntimeError("CosmologyWrapper.set_current does not know how to handle this version of astropy")

        self.cosmologyInitialized = True
        self.activeCosmology = universe
        self.set_units()

    def get_current(self):
        """
        Return the cosmology currently stored as the current cosmology

        This is for users who want direct access to all of astropy.cosmology's methods,
        not just those wrapped by this class.

        documentation for astropy.cosmology can be found at the URL below (be sure to check which version of
        astropy you are running; as of 30 October 2014, the anaconda distributed with the stack
        comes with version 0.2.5)

        https://astropy.readthedocs.org/en/v0.2.5/cosmology/index.html
        """
        if not self.cosmologyInitialized:
            raise RuntimeError("Should not call CosmologyWrapper.get_current(); you have not set_current()")

        if 'get_current' in dir(cosmology):
            return cosmology.get_current()
        elif 'default_cosmology' in dir(cosmology):
            return cosmology.default_cosmology.get()
        else:
            raise RuntimeError("CosmologyWrapper.get_current does not know how to handle this version of astropy")

    def initializeCosmology(self, H0=72.0, Om0=0.25, Ode0=None, w0=None, wa=None):
        """
        Initialize the cosmology wrapper with the parameters specified
        (e.g. does not account for massive neutrinos)

        param [in] H0 is the Hubble parameter at the present epoch in km/s/Mpc

        param [in] Om0 is the current matter density paramter (fraction of critical density)

        param [in] Ode0 is the current dark energy density parameter

        param [in] w0 is the current dark energy equation of state w0 paramter

        param[in] wa is the current dark energy equation of state wa paramter

        The total dark energy equation of state as a function of z is
        w = w0 + wa z/(1+z)

        Currently, this wrapper class expects you to specify either a LambdaCDM (flat or non-flat) cosmology
        or a w0, wa (flat or non-flat) cosmology.
        """
        self.activeCosmology = None
        
        if w0 is not None and wa is None:
            wa = 0.0

        if w0 is None and wa is None and (Ode0 is None or Om0+Ode0==1.0):
            universe = cosmology.FlatLambdaCDM(H0=H0, Om0=Om0)
        elif w0 is None and wa is None:
            universe = cosmology.LambdaCDM(H0=H0, Om0=Om0, Ode0=Ode0)
        elif Ode0 is None or Om0+Ode0==1.0:
            universe = cosmology.Flatw0waCDM(H0=H0, Om0=Om0, w0=w0, wa=wa)
        else:
            universe = cosmology.w0waCDM(H0=H0, Om0=Om0, Ode0=Ode0,
                                         w0=w0, wa=wa)

        self.set_current(universe)

    def loadDefaultCosmology(self):
        self.initializeCosmology(H0=72.0, Om0=0.23)

    def H(self, redshift=0.0):
        """return the Hubble paramter in km/s/Mpc at the specified redshift"""
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call H; cosmology has not been initialized")

        H = self.activeCosmology.H(redshift)

        if 'value' in dir(H):
            if H.unit == self.hUnits:
                return H.value
            else:
                return H.to(self.hUnits).value
        else:
            return H


    def OmegaMatter(self, redshift=0.0):
        """return the matter density paramter (fraction of critical density) at the specified redshift"""
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaMatter; cosmology has not been initialized")

        return self.activeCosmology.Om(redshift)

    def OmegaDarkEnergy(self, redshift=0.0):
        """return the dark energy density paramter (fraction of critical density) at the specified redshift"""
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaDarkEnergy; cosmology has not been initialized")

        return self.activeCosmology.Ode(redshift)

    def OmegaPhotons(self, redshift=0.0):
        """return the photon density paramter (fraction of critical density) at the specified redshift"""
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaPhotons; cosmology has not been initialized")

        return self.activeCosmology.Ogamma(redshift)

    def OmegaNeutrinos(self, redshift=0.0):
        """
        return the neutrino density paramter (fraction of critical density) at the specified redshift

        assumes neutrinos are massless
        """
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaNeutrinos; cosmology has not been initialized")

        return self.activeCosmology.Onu(redshift)

    def OmegaCurvature(self, redshift=0.0):
        """
        return the effective curvature density paramter (fraction of critical density) at the
        specified redshift.

        Positive means the universe is open.

        Negative means teh universe is closed.

        Zero means the universe is flat.
        """
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaCurvature; cosmology has not been initialized")

        return self.activeCosmology.Ok(redshift)

    def w(self, redshift=0.0):
        """return the dark energy equation of state at the specified redshift"""
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call w; cosmology has not been initialized")

        return self.activeCosmology.w(redshift)

    def comovingDistance(self, redshift=0.0):
        """
        return the comoving distance to the specified redshift in Mpc

        note, this comoving distance is X in the FRW metric

        ds^2 = -c^2 dt^2 + a^2 dX^2 + a^2 sin^2(X) dOmega^2

        i.e. the curvature of the universe is folded into the sin()/sinh() function.
        This distande just integrates dX = c dt/a
        """
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call comovingDistance; cosmology has not been initialized")

        dd = self.activeCosmology.comoving_distance(redshift)

        if 'value' in dir(dd):
            if dd.unit == self.distanceUnits:
                return dd.value
            else:
                return dd.to(self.distanceUnits).value
        else:
            return dd

    def luminosityDistance(self, redshift=0.0):
        """
        the luminosity distance to the specified redshift in Mpc

        accounts for spatial curvature
        """
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call luminosityDistance; cosmology has not been initialized")

        dd = self.activeCosmology.luminosity_distance(redshift)

        if 'value' in dir(dd):
            if dd.unit == self.distanceUnits:
                return dd.value
            else:
                return dd.to(self.distanceUnits).value
        else:
            return dd

    def angularDiameterDistance(self, redshift=0.0):
        """angular diameter distance to the specified redshift in Mpc"""
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call angularDiameterDistance; cosmology has not been initialized")

        dd = self.activeCosmology.angular_diameter_distance(redshift)

        if 'value' in dir(dd):
            if dd.unit == self.distanceUnits:
                return dd.value
            else:
                return dd.to(self.distanceUnits).value
        else:
            return dd

    def distanceModulus(self, redshift=0.0):
        """distance modulus to the specified redshift"""
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call distanceModulus; cosmology has not been initialized")

        mm = self.activeCosmology.distmod(redshift)
        if 'unit' in dir(mm):
            if mm.unit == self.modulusUnits:
                return mm.value
            else:
                return mm.to(self.modulusUnits).value
        else:
            return mm
