import copy
import numpy
from collections import OrderedDict
from lsst.sims.photUtils import Bandpass, Sed

__all__ = ["CatSimBandpassDict"]

class CatSimBandpassDict(object):
    """
    This class will wrap an OrderedDict of Bandpass instantiations.

    Upon instantiation, this class's constructor will resample
    the input Bandpasses to be on the same wavelength grid (defined
    by the first input Bandpass).  The constructor will then calculate
    the 2-D phiArray for quick calculation of magnitudes in all
    Bandpasses simultaneously (see the member method calcMagListFromSed).
    """

    def __init__(self, bandpassList, bandpassNameList):
        """
        @param [in] bandpassList is a list of Bandpass instantiations

        @param [in] bandpassNameList is a list of tags to be associated
        with those Bandpasses
        """
        self._bandpassDict = OrderedDict()
        for bandpassName, bandpass in zip(bandpassNameList, bandpassList):

            if bandpassName in self._bandpassDict:
                raise RuntimeError("The bandpass %s occurs twice in your input " % bandpassName \
                                   + "to CatSimBandpassDict")

            self._bandpassDict[bandpassName] = copy.deepcopy(bandpass)

        dummySed = Sed()
        self._phiArray, self._wavelenStep = dummySed.setupPhiArray(self._bandpassDict.values())
        self._wavelen_match = self._bandpassDict.values()[0].wavelen
        self._nBandpasses = len(self._bandpassDict)


    def __getitem__(self, bandpass):
        return self._bandpassDict[bandpass]


    def __len__(self):
        return len(self._bandpassDict)


    def __iter__(self):
        for val in self._bandpassDict:
            yield val


    def values(self):
        return self._bandpassDict.values()


    def keys(self):
        return self._bandpassDict.keys()


    def calcMagListFromSed(self, sedobj, indices=None):
        """
        Return a list of magnitudes for a single Sed object.

        @param [in] sedobj is an Sed object

        @param [in] indices is an optional list of indices indicating which bandpasses to actually
        calculate magnitudes for.  Other magnitudes will be listed as 'None' (i.e. this method will
        return as many magnitudes as were loaded with the loadBandpassesFromFiles methods; it will
        just return nonsense for magnitudes you did not actually ask for)

        @param [out] magList is a list of magnitudes in the bandpasses stored in self.bandpassDict
        """

        if sedobj.wavelen is not None:

            # If the Sed's wavelength grid agrees with self._wavelen_match to one part in
            # 10^6, just use the Sed as-is.  Otherwise, copy it and resample it onto
            # self._wavelen_match
            if len(sedobj.wavelen)!=len(self._wavelen_match) or \
            not numpy.allclose(sedobj.wavelen, self._wavelen_match, atol=0.0, rtol=1.0e-6):
                dummySed = Sed(wavelen=sedobj.wavelen, flambda=sedobj.flambda)
                dummySed.resampleSED(wavelen_match=self._bandpassDict.values()[0].wavelen)
            else:
                dummySed = sedobj


            #for some reason, moving this call to flambdaTofnu()
            #to a point earlier in the
            #process results in some SEDs having 'None' for fnu.
            #
            #I looked more carefully at the documentation in Sed.py
            #Any time you update flambda in any way, fnu gets set to 'None'
            #This is to prevent the two arrays from getting out synch
            #(e.g. renormalizing flambda but forgettint to renormalize fnu)
            #
            dummySed.flambdaTofnu()

            if indices is not None:
                magList = [numpy.NaN]*self._nBandpasses

                magArray = dummySed.manyMagCalc(self._phiArray, self._wavelenStep, observedBandPassInd=indices)
                for i,ix in enumerate(indices):
                    magList[ix] = magArray[i]
            else:
                magList = dummySed.manyMagCalc(self._phiArray, self._wavelenStep)

            return numpy.array(magList)

        else:
            return numpy.array([numpy.NaN]*self._nBandpasses)



    def calcMagListFromSedList(self, sedList, indices=None):
        """
        Return a 2-D array of magnitudes from a CatSimSedList.
        Each row will correspond to a different Sed, each column
        will correspond to a different bandpass.

        @param [in] sedList is a CatSimSedList containing the Seds
        whose magnitudes are desired.

        @param [in] indices is an optional list of indices indicating which bandpasses to actually
        calculate magnitudes for.  Other magnitudes will be listed as 'None' (i.e. this method will
        return as many magnitudes as were loaded with the loadBandpassesFromFiles methods; it will
        just return nonsense for magnitudes you did not actually ask for)

        @param [out] output_list is a 2-D numpy array containing the magnitudes
        of each Sed (the rows) in each bandpass contained in this CatSimBandpassDict
        (the columns)
        """

        one_at_a_time = False
        if sedList.wavelenMatch is None:
            one_at_a_time = True
        elif len(sedList.wavelenMatch) != len(self._wavelen_match):
            one_at_a_time = True
        elif not numpy.allclose(sedList.wavelenMatch, self._wavelen_match, atol=0.0, rtol=1.0e-6):
            one_at_a_time = True

        output_list = []
        if one_at_a_time:
            for sed_obj in sedList:
                sub_list = self.calcMagListFromSed(sed_obj, indices=indices)
                output_list.append(sub_list)
        else:
            if indices is not None:
                for sed_obj in sedList:
                    sub_list = numpy.array([numpy.Nan]*self._nBandpasses)
                    if sed_obj.wavelen is not None:
                        sed_obj.flambdaTofnu()
                        mag_list = sed_obj.manyMagCalc(self._phiArray, self._wavelenStep, observedBandpassInd=indices)
                        for i,ix in enumerate(indices):
                            sub_list[ix] = mag_list[i]
                    output_list.append(sub_list)
            else:
                for sed_obj in sedList:
                    if sed_obj.wavelen is None:
                        sub_list = numpy.array([numpy.Nan]*self._nBandpasses)
                    else:
                        sed_obj.flambdaTofnu()
                        sub_list = sed_obj.manyMagCalc(self._phiArray, self._wavelenStep)
                    output_list.append(sub_list)

        return numpy.array(output_list)


    @property
    def phiArray(self):
        """
        A 2-D numpy array storing the values of phi (see eqn 2.3 of the science
        book) for all of the bandpasses in this dict.
        """
        return self._phiArray

    @phiArray.setter
    def phiArray(self, value):
        raise RuntimeError("You should not be setting phiArray on the fly " \
                           + "in a CatSimBandpassDict")


    @property
    def wavelenStep(self):
        """
        The step size of the wavelength grid for all of the bandpasses
        stored in this dict.
        """
        return self._wavelenStep

    @wavelenStep.setter
    def wavelenStep(self, value):
        raise RuntimeError("You should not be setting wavelenStep on the fly " \
                          + "in a CatSimBandpassDict")


    @property
    def nBandpasses(self):
        """
        The number of bandpasses stored in this dict.
        """
        return self._nBandpasses

    @nBandpasses.setter
    def nBandpasses(self, value):
        raise RuntimeError("You should not be setting nBandpasses on the fly " \
                           + "in a CatSimBandpassDict")


    @property
    def wavelenMatch(self):
        """
        The wavelength grid (in nm) on which all of the bandpass
        throughputs have been sampled.
        """
        return self._wavelen_match

    @wavelenMatch.setter
    def wavelenMatch(self, value):
        raise RuntimeError("You should not be setting wavelenMatch on the fly " \
                           + "in a CatSimBandpassDict")
