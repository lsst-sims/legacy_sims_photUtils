import numpy
from lsst.sims.photUtils import CosmologyWrapper
from lsst.sims.catalogs.measures.instance import compound

__all__ = ["ExampleCosmologyMixin"]

class ExampleCosmologyMixin(CosmologyWrapper):
    
    def get_distanceModulus(self):
        zz = self.column_by_name(redshift);
        return numpy.array([self.distanceModulus(redshift=zz)])
    
    @compound(
    'uAbs', 'gAbs', 'rAbs', 'iAbs', 'zAbs', 'yAbs',
    'uBulgeAbs', 'gBulgeAbs', 'rBulgeAbs', 'iBulgeAbs', 'zBulgeAbs', 'yBulgeAbs',
    'uDiskAbs', 'gDiskAbs', 'rDiskAbs', 'iDiskAbs', 'zDiskAbs', 'yDiskAbs',
    'uAgnAbs', 'gAgnAbs', 'rAgnAbs', 'iAgnAbs', 'zAgnAbs', 'yAgnAbs'
    )
    def absolute_magnitude_getter(self):
        uu = self.column_by_name('uRecalc')
        gg = self.column_by_name('gRecalc')
        rr = self.column_by_name('rRecalc')
        ii = self.column_by_name('iRecalc')
        zz = self.column_by_name('zRecalc')
        yy = self.column_by_name('yRecalc')
        
        ubulge = self.column_by_name('uBulge')
        gbulge = self.column_by_name('gBulge')
        rbulge = self.column_by_name('rBulge')
        ibulge = self.column_by_name('iBulge')
        zbulge = self.column_by_name('zBulge')
        ybulge = self.column_by_name('yBulge')
        
        udisk = self.column_by_name('uDisk')
        gdisk = self.column_by_name('gDisk')
        rdisk = self.column_by_name('rDisk')
        idisk = self.column_by_name('iDisk')
        zdisk = self.column_by_name('zDisk')
        ydisk = self.column_by_name('yDisk')
        
        uagn = self.column_by_name('uAgn')
        gagn = self.column_by_name('gAgn')
        ragn = self.column_by_name('rAgn')
        iagn = self.column_by_name('iAgn')
        zagn = self.column_by_name('zAgn')
        yagn = self.column_by_name('yAgn')
        
        redshift = self.column_by_name('redshift')
        modulus = self.column_by_name('distanceModulus')
        
        #now we need to undo cosmological dimming, since applyAvAndRedshift in
        #the PhotometryGalaxies mixin applies cosmological dimming, and that will
        #be redundant with the (1+z) factor in the luminosity distance
        #(which is used to calculate the cosmological distance modulus)
        
        brightening = -2.5*numpy.log10(1.0+redshift)
       
        brightening += modulus
        
        uu += brightening
        gg += brightening
        rr += brightening
        ii += brightening
        zz += brightening
        yy += brightening
        
        ubulge += brightening
        gbulge += brightening
        rbulge += brightening
        ibulge += brightening
        zbulge += brightening
        ybulge += brightening
        
        udisk += brighteting
        gdisk += brighteting
        rdisk += brightening
        idisk += brighteting
        zdisk += brightening
        ydisk += brightening
        
        uagn += brightening
        gagn += brightening
        ragn += brightening
        iagn += brightening
        zagn += brightening
        yagn += brightening
        
        return numpy.array([
                          uu, gg, rr, ii, zz, yy,
                          ubulge, gbulge, rbulge, ibulge, zbulge, ybulge,
                          udisk, gdisk, rdisk, idisk, zdisk, ydisk,
                          uagn, gagn, ragn, iagn, zagn, yagn
                          ])
