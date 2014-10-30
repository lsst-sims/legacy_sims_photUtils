import astropy.cosmology as cosmology
import astropy.units as units

__all__ = ["CosmologyWrapper"]

class CosmologyWrapper(object):

    cosmologyInitialized = False

    def set_units(self):
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
        """
        if not self.cosmologyInitialized:
            raise RuntimeError("Should not call CosmologyWrapper.get_current(); you have not set_current()")

        if 'get_current' in dir(cosmology):
            return cosmology.get_current()
        elif 'default_cosmology' in dir(cosmology):
            return cosmology.default_cosmology.get()
        else:
            raise RuntimeError("CosmologyWrapper.get_current does not know how to handle this version of astropy")

    def Initialize(self, H0=72.0, Om0=0.25, Ode0=0.75, w0=-1.0, wa=0.0):
        self.activeCosmology = None

        self.H0 = H0
        self.Om0 = Om0
        self.Ode0 = Ode0
        self.w0 = w0
        self.wa = wa

        if self.w0==-1.0 and self.wa==0.0 and self.Om0+self.Ode0==1.0:
            universe = cosmology.FlatLambdaCDM(H0=self.H0, Om0=self.Om0)
        elif self.w0==-1.0 and self.wa==0.0:
            universe = cosmology.LambdaCDM(H0=self.H0, Om0=self.Om0, Ode0=self.Ode0)
        elif self.Om0+self.Ode0==1.0:
            universe = cosmology.Flatw0waCDM(H0=self.H0, Om0=self.Om0, w0=self.w0, wa=self.wa)
        else:
            universe = cosmology.w0waCDM(H0=self.H0, Om0=self.Om0, Ode0=self.Ode0,
                                         w0=self.w0, wa=self.wa)

        self.set_current(universe)

    def H(self, redshift=0.0):

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

        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaMatter; cosmology has not been initialized")

        return self.activeCosmology.Om(redshift)

    def OmegaDarkEnergy(self, redshift=0.0):
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaDarkEnergy; cosmology has not been initialized")

        return self.activeCosmology.Ode(redshift)

    def OmegaPhotons(self, redshift=0.0):
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaPhotons; cosmology has not been initialized")

        return self.activeCosmology.Ogamma(redshift)

    def OmegaNeutrinos(self, redshift=0.0):
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaNeutrinos; cosmology has not been initialized")

        return self.activeCosmology.Onu(redshift)

    def OmegaCurvature(self, redshift=0.0):
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call OmegaCurvature; cosmology has not been initialized")

        return self.activeCosmology.Ok(redshift)

    def w(self, redshift=0.0):
        if not self.cosmologyInitialized:
            raise RuntimeError("cannot call w; cosmology has not been initialized")

        return self.activeCosmology.w(redshift)

    def comovingDistance(self, redshift=0.0):
        """in Mpc"""
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
        """in Mpc"""
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
        """in Mpc"""
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
