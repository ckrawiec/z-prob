# z-prob

Calculate the Bayesian posterior probability that target galaxies belong to a given range in redshift.

runzprob.py: reads a config file and calls zprobability.py, writes outputs to fits file

zprobability.py: describes Target and Template classes which use various implementations of the Gaussian likelihood function to calculate probabilities
