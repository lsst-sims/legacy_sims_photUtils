import numpy
from lsst.sims.photUtils import CosmologyObject
from lsst.sims.catalogs.measures.instance import compound

__all__ = ["ExampleCosmologyMixin"]

class ExampleCosmologyMixin():
    """
    This is an example mixin which applies cosmological distance modulus
    to galaxy magnitudes
    """

    cosmology = None

    def get_distanceModulus(self):
        if self.cosmology is None:
            self.cosmology = CosmologyObject()

        zz = self.column_by_name('redshift');
        return self.cosmology.distanceModulus(redshift=zz)

    @compound(
    'uAbs', 'gAbs', 'rAbs', 'iAbs', 'zAbs', 'yAbs',
    'uBulgeAbs', 'gBulgeAbs', 'rBulgeAbs', 'iBulgeAbs', 'zBulgeAbs', 'yBulgeAbs',
    'uDiskAbs', 'gDiskAbs', 'rDiskAbs', 'iDiskAbs', 'zDiskAbs', 'yDiskAbs',
    'uAgnAbs', 'gAgnAbs', 'rAgnAbs', 'iAgnAbs', 'zAgnAbs', 'yAgnAbs'
    )
    def get_cosmological_absolute_magnitudes(self):
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

        uu += modulus
        gg += modulus
        rr += modulus
        ii += modulus
        zz += modulus
        yy += modulus

        ubulge += modulus
        gbulge += modulus
        rbulge += modulus
        ibulge += modulus
        zbulge += modulus
        ybulge += modulus

        udisk += modulus
        gdisk += modulus
        rdisk += modulus
        idisk += modulus
        zdisk += modulus
        ydisk += modulus

        uagn += modulus
        gagn += modulus
        ragn += modulus
        iagn += modulus
        zagn += modulus
        yagn += modulus

        return numpy.array([
                          uu, gg, rr, ii, zz, yy,
                          ubulge, gbulge, rbulge, ibulge, zbulge, ybulge,
                          udisk, gdisk, rdisk, idisk, zdisk, ydisk,
                          uagn, gagn, ragn, iagn, zagn, yagn
                          ])
