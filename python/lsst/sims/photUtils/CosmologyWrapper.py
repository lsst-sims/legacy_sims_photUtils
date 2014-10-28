import astropy.cosmology as cosmology

__all__ = ["CosmologyWrapper"]

class CosmologyWrapper(object):

    cosmologyInitialized = False

    def _determine_API(self):
        if 'set_current' in dir(cosmology):
            self.cosmologyVersion='0.2.5'
        elif 'default_cosmology' in dir(cosmology):
            self.cosmologyVersion='1.0'
        else:
            raise RuntimeError("CosmologyWrapper does not know how to handle this version of astropy")

    def set_current(self, universe):
        """
        Take the cosmology indicated by 'universe' and set it as the current/default
        cosmology (depending on the API of the version of astropy being run)
        """
        
        if not hasattr(self, 'cosmologyVersion'):
            self._determineAPI()

        if self.cosmologyVersion == '0.2.5':
            cosmology.set_current(universe)
        elif self.cosmologyVersion == '1.0':
            cosmology.default_cosmology.set(universe)

        self.cosmologyInitialized = True
        self.activeCosmology = universe

    def get_current(self):
        """
        Return the cosmology currently stored as the current cosmology
        """
        if not self.cosmologyInitialized:
            raise RuntimeError("Should not call CosmologyWrapper.get_current(); you have not set_current()")

        if self.cosmologyVersion == '0.2.5':
            return cosmology.get_current()
        elif self.cosmologyVersion == '1.0':
            return cosmology.default_cosmology.get()

    def Initialize(self, H0=72.0, Om0=0.25, Ode0=0.75, w0=-1.0, wa=0.0):
        self.activeCosmology = None
        self._determine_API()

        print 'version ',self.cosmologyVersion

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
            return H.value
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
