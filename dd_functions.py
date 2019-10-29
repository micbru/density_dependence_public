'''
This will contain functions needed to make it easy to work with data and theoretical predictions.
I will add to it as I need more functions.
I can then call this into a jupyter notebook to make things easy and neat.
'''
import numpy as np
import pandas as pd
import scipy.stats as st

def dd_prob(n0,alpha):
    '''This gives the density dependent probability function as a vector of length n0+1 
        for a given n0 and alpha, where n0 is the species abundance.'''
    n = np.arange(n0+1)
    # First get the unnormed probability, using binom for large numbers
    bi_test = st.binom.pmf(n, n0, p=0.5)
    # If this is 0, we will get a warning, so use log to generate
    if np.any(bi_test==0.):
        logbinom = st.binom.logpmf(n, n0, p=0.5)
        binom = np.exp((alpha-1)*logbinom)
    # Otherwise it's fine to just take the exponent
    else:
        binom = bi_test**(alpha-1)
    unnormed = ((n/n0)**alpha+((n0-n)/n0)**alpha)*binom
    return unnormed/unnormed.sum()

# Create bisection function just using indexes
def bisect(df,xmax,ymax,level=1,xkey='gx',ykey='gy',skey='sp'):
    '''
    df = dataframe, where spatial and species information is stored. 
        Must be stored as a pandas dataframe with location data for each species, as well as species
        information so we can keep the species seperate when bisecting.
    level = int, this is how many bisections we want to draw. By default = 1 (one bisection)
    The maximum x and y coordinates of the dataset
        xmax, ymax
    The keys to access the spatial data and species data in the dataframe. By default work with BCI dataset
        xkey = 'gx' 
        ykey = 'gy'
        skey = 'sp'
    Note that this will bisect first in y, then alternate.
    '''
    # Setup a collection of arrays
    df_coll = []
    # Really we just need to get the 4 arguments bounding the box:
    #     arg_botom, arg_top, arg_left, arg_right.
    for l in np.arange(2**np.ceil(level/2)):
        # Get top and bottom args of box
        arg_bottom = df[ykey]>=ymax/(2**np.ceil(level/2))*l # Need >=
        arg_top = df[ykey]<ymax/(2**np.ceil(level/2))*(l+1) 
        # Careful about choosing ymax so we get the boundary. 
        # Choose the maximum data point + some small amount.
        if level == 1: # If we only bisect once we DON't use arg_left and arg_right
            df_coll.append(df[skey].loc[arg_bottom&arg_top].value_counts())
        else:
            for ll in np.arange(2**np.floor(level/2)):
                arg_left = df[xkey]>=xmax/(2**np.floor(level/2))*ll # Need >=
                arg_right = df[xkey]<xmax/(2**np.floor(level/2))*(ll+1)
                # Careful about choosing xmax so we get the boundary.
                # Choose maximum data point + some small amount
                # Now that we've bisected, add this to the list!
                df_coll.append(df[skey].loc[arg_left&arg_right&arg_bottom&arg_top].value_counts())
    # Concatenate to a pandas frame and return
    df_ret = pd.concat(df_coll,axis=1,sort=True)
    df_ret.columns=[i for i in range(2**level)]
    return df_ret

def create_f(df,thresh=0):
    '''
    Take a data frame of bisections (as created by the bisect function, or directly from bisection data)
        and create fraction plots.
    The fractions will be taken in the x direction by default, but that can be changed later if I want.
        Would have to sum over different indices...
    There is also "thresh" which is the threshold for number of individuals in a species
    I used to remove is there were too many NaN, but that is just removing high and low fractions, which is bad!
    '''
    # If the threshold is nonzero, drop the species where there are not enough individuals
    if thresh>0:
        df_data = df[df.T.sum()>thresh].fillna(0.)
    else:
        df_data = df.fillna(0.)

    # Get number of points
    nsp = len(df_data.columns)//2
    # Set up a numpy array to return
    fr_data = np.zeros((nsp,len(df_data)))
    # Set up n0_data
    n0_data = np.zeros((nsp,len(df_data)),dtype=int)
    # Get index data to put in dataframe as column later
    ids = np.tile(df_data.index,nsp) # Use tile instead of repeat to get proper repetition
    # Create the fractions by summing adjacent cells
    for i in range(nsp):
        n0_data[i] = df_data[2*i]+df_data[2*i+1]
        fr_data[i] = df_data[2*i]/n0_data[i]
    df_ret = pd.DataFrame({'sp': ids,'frac': fr_data.flatten(),'n0': n0_data.flatten()},columns=['sp','frac','n0'])
    df_ret.dropna(axis=0,how='any',inplace=True) # Drop if n0=0 (there will be a NaN)
    return df_ret

def loglikelihood(alpha,n,n0):
    '''
    Pass in a set of n's and n0's and see the likelihood for a given alpha. 
    Make sure n is integer.
    Note that depending on n0, we may have to have VERY tight bounds to minimize this.
    It's still possible to get 0's in dd_prob, which won't allow us to take the log.
    To minimize and get the most likely alpha using this, use the following example:
    
    from scipy.optimize import minimize_scalar
    bestfit_alpha = minimize_scalar(dd.loglikelihood,bounds=(low_bound,upper_bound),method='bounded',
                args=(list_of_n,list_of_n0))
    '''
    
    assert len(n)==len(n0), 'Lengths do not match!'
    # Sum over likelihoods
    likelihood = 0
    for i,j in zip(n,n0):
        prob = dd_prob(j,alpha)[i]
        assert prob!=0, 'Requires tighter bounds on allowed alpha!'
        likelihood -= np.log(prob)
    return likelihood

def contours(alpha,pc,nmax):
    '''
    Will generate contour data for a given alpha and percentage.
    pc is 1-the percent of the contour desired, ie. for 95% contours pc = 0.05
    nmax is the maximum abundance we want the contour for.
    Return the range used and the intervals.
    '''
    interval = []
    # Use a unique logrange to speed up sampling. Uses base 10 and 50 points by default.
    logrange = np.unique(np.logspace(0,np.log10(nmax),dtype=int))
    for n in logrange:
        dd_cdf = np.cumsum(dd_prob(n,alpha))
        interval.append((np.searchsorted(dd_cdf,pc/2)/n,np.searchsorted(dd_cdf,1-pc/2)/n))
    return logrange,interval
