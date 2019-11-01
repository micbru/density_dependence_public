# density_dependence_public
Code supplementing "Density dependence and a colonization rule predict spatial patterning"

dd_functions.py is the main python file with the following functions:

    dd_prob(n0,alpha) - Gives the Pi function numerically for a given n0 and alpha.
    bisect(df,xmax,ymax,level=1,xkey='gx',ykey='gy',skey='sp') - Can bisect an explicitly spatial dataset at any level.
    create_f(df,thresh=0) - Creates fractions from bisection data.
    loglikelihood(alpha,n,n0) - Calculates the loglikelihood for a given alpha. 
                                Can be used in combination with minimize_scalar to find most likely alpha.
    contours(alpha,pc,nmax) - For plotting contour data, as in Fig. 2 in the paper.

The other files included here are Python notebooks to replicate the figures in the paper, and to demonstrate usage of dd_functions.py. The required data are available at: https://doi.org/10.6078/D1MQ2V (Serpentine) and https://doi.org/10.15146/5xcp-0d46 (BCI).
    
    PlotPiFunction_Fig1.ipynb - Examples for directly plotting the Pi function.
    BCI.ipynb - Replicates analysis for the BCI dataset, generating the BCI part of Figs. 2, 4, and 5.
    Serpentine.ipynb - As BCI, but for the serpentine dataset.
    SerpentineBCIScaling_Fig3.ipynb - Replicates the aggregate scaling we showed in Fig. 3.

Please contact me if you have any question about how to use this code. If you do use it, please cite the Density dependence paper (To be published) or my github page for the time being.
