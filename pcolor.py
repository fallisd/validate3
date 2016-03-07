import brewer2mpl
from discrete_cmap import discrete_cmap
from colormaps import viridis

def default_pcolor_args(data, anom=False):
    """Returns a dict with default pcolor params as key:value pairs

    Parameters
    ----------
    data : numpy array
    anom : boolean
           True if positive/negative display is wanted

    Returns
    -------
    dictionary
    """
    # Set 3-std range for colorbar to exclude outliers.
    if anom:
        # For anomalies, center range around 0
        anom_max = abs(data).mean() + abs(data).std() * 3.0
        vmin = -1 * anom_max
        vmax = anom_max
        # Anomaly cmap
        cmap = anom_cmap()

    else:
        # otherwise, center around the mean
        vmin = data.mean() - data.std() * 3.0
        vmax = data.mean() + data.std() * 3.0
        # Use true min/max if they are closer to the mean than the 3-std spread.
        if vmax > data.max():
            vmax = data.max()
        if vmin < data.min():
            vmin = data.min()
        # New mpl, colorblind friendly, continuously varying, default cmap
        cmap = viridis

    d = {'vmin': vmin,
         'vmax': vmax,
         'cmap': cmap,
         'rasterized': True
         }

    return d
    
def anom_cmap():
    """return a discrete blue-red cmap from colorbrewer"""
    ncols = 11
    cmap_anom = brewer2mpl.get_map('RdBu', 'diverging', ncols,
                                   reverse=True).mpl_colormap
    cmap_anom = discrete_cmap(ncols, cmap_anom)
    return cmap_anom
