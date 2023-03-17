import numpy as np
import scipy as scp

def smooth2level(y, x, level, span, method=''):
    '''Description:
    Smooths level of derivative of datapoints y,x with span number/percent of datapoints using method. If level=0, only datapoints are smoothed, if level > 0, all derivatives are smootzed as well.
    Function first mooths, then computes the derivative using gradient. This is done until level is reached. In the end, starting from the highest order derivative, smoothed values are integrated back to reach the previous derivatives and the datapoints themselves. The integration constant is chosen to be the first value of the previous derivative /datapoint
    
    NOTICE: smoothing to arbitrary level is not a good idea. 
    You can try with Er.dat using level=1 and level=2 or with Te.dat using level=2 and level=3. Datapoints themselves may get off the original ones if you use too high level.
    
    input:
    y       ... y datapoints
    x       ... x datapoints
    level   ... level of derivative to smooth (0 = function itself)
    span    ... number of datapoints for calculating smoothed values (from smooth function). or percent of datapoints for calculating smoothed values
    method  ... method for smoothing (optional, e.g. nearest)
    
    output: numpy array with dimension of level + 1 containing smoothed datapoints and derivatives'''

    out = np.zeros_like(y)
    for i in range(0, level):
        out = np.array([out, np.zeros_like(y)])
    out[:,0] = y

    # compute all derivatives up to level and smooth
    for j in range(0,level+1):
        if j > 0:
            out[:,j] = np.gradient(out[:,j-1], x)
        # Savitzky Golay filter
        if method == '':
            out[:,j] = scp.signal.savgol_filter(out[:,j], span, 3)
        else:
            out[:,j] = scp.signal.savgol_filter(out[:,j], span, method)
    
    # integrate back up to quantity
    for j in range(1,level+1):
        ind = level + 1 - j
        old = out[:,ind]
        out[:,ind] = np.trapz(out[:,ind+1], x=x)
        const = np.mean(out[:,ind] - old)
        out[:,ind] = out[:,ind] - const

    return out