__all__ = ['tbavg','extract_dynspec','__plot__','custom_colormap']

from numpy import load

def load_dict(savefile):
    # load a dictionary that has been saved using numpy.save
    return load(savefile).reshape(1)[0]
