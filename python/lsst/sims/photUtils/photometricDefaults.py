import numpy

__all__ = ["PhotometricDefaults"]

class PhotometricDefaults(object):
    """
    This class exists to store default values of parameters characterizing
    noise due to the telescope in one place.
    """

    # The following *wavelen* parameters are default values for gridding wavelen/sb/flambda.


    seeing = {'u': 0.77, 'g':0.73, 'r':0.70, 'i':0.67, 'z':0.65, 'y':0.63}  # Default seeing values (in ")

   #taken from table 2 of arxiv:0805.2366 (note that m5 is for 2 15 second exposures)
    m5 = {'u':23.68, 'g':24.89, 'r':24.43, 'i':24.00, 'z':24.45, 'y':22.60}
    gamma = {'u':0.037, 'g':0.038, 'r':0.039, 'i':0.039, 'z':0.040, 'y':0.040}
