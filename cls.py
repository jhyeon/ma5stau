def cls(NumObserved, ExpectedBG, BGError, SigHypothesis, NumToyExperiments):
    import scipy.stats
    # generate a set of expected-number-of-background-events, one for each toy
    # experiment, distributed according to a Gaussian with the specified mean
    # and uncertainty
    ExpectedBGs = scipy.stats.norm.rvs(loc=ExpectedBG, scale=BGError, size=NumToyExperiments)

    # Ignore values in the tail of the Gaussian extending to negative numbers
    ExpectedBGs = [value for value in ExpectedBGs if value > 0]

    # For each toy experiment, get the actual number of background events by
    # taking one value from a Poisson distribution created using the expected
    # number of events.
    ToyBGs = scipy.stats.poisson.rvs(ExpectedBGs)
    ToyBGs = map(float, ToyBGs)

    # The probability for the background alone to fluctutate as LOW as
    # observed = the fraction of the toy experiments with backgrounds as low as
    # observed = p_b.
    # NB (1 - this p_b) corresponds to what is usually called p_b for CLs.
    p_b = scipy.stats.percentileofscore(ToyBGs, NumObserved, kind='weak')*.01
    #print 1. - p_b

    # Toy MC for background+signal
    ExpectedBGandS = [expectedbg + SigHypothesis for expectedbg in ExpectedBGs]
    ExpectedBGandS = [x for x in ExpectedBGandS if x > 0]
    if len(ExpectedBGandS)==0:
        return 0.
    ToyBplusS = scipy.stats.poisson.rvs(ExpectedBGandS)
    ToyBplusS = map(float, ToyBplusS)

    # Calculate the fraction of these that are >= the number observed,
    # giving p_(S+B). Divide by (1 - p_b) a la the CLs prescription.
    p_SplusB = scipy.stats.percentileofscore(ToyBplusS, NumObserved, kind='weak')*.01

    if p_SplusB>p_b:
        print 0.
    else:
        print 1.-(p_SplusB / p_b) # 1 - CLs


import math
#rejected = rejected(low) and/or rejected(high)
#accepted = accepted(low) and accepted(high)

print '####',100,1
cll = cls(10, 6.0, 1.7, 7.32, 100000) #low
clh = cls(7, 10.2, 3.3, 4.44, 100000) #high

print '####',120,1
cll = cls(10, 6.0, 1.7, 14.65, 100000) #low
clh = cls(7, 10.2, 3.3, 8.24, 100000) #high

print '####',120,20
cll = cls(10, 6.0, 1.7, 9.82, 100000) #low
clh = cls(7, 10.2, 3.3, 5.24, 100000) #high

print '####',120,30
cll = cls(10, 6.0, 1.7, 7.55, 100000) #low
clh = cls(7, 10.2, 3.3, 4.98, 100000) #high

print '####',120,50
cll = cls(10, 6.0, 1.7, 2.98, 100000) #low
clh = cls(7, 10.2, 3.3, 1.66, 100000) #high

print '####',140,1
cll = cls(10, 6.0, 1.7, 16.82, 100000) #low
clh = cls(7, 10.2, 3.3, 11.54, 100000) #high

print '####',140,30
cll = cls(10, 6.0, 1.7, 11.95, 100000) #low
clh = cls(7, 10.2, 3.3, 7.81, 100000) #high

print '####',140,50
cll = cls(10, 6.0, 1.7, 7.23, 100000) #low
clh = cls(7, 10.2, 3.3, 4.32, 100000) #high

print '####',160,1
cll = cls(10, 6.0, 1.7, 17.67, 100000) #low
clh = cls(7, 10.2, 3.3, 12.87, 100000) #high

print '####',180,60
cll = cls(10, 6.0, 1.7, 10.88, 100000) #low
clh = cls(7, 10.2, 3.3, 8.87, 100000) #high

print '####',180,80
cll = cls(10, 6.0, 1.7, 7.41, 100000) #low
clh = cls(7, 10.2, 3.3, 5.45, 100000) #high

print '####',200,50
cll = cls(10, 6.0, 1.7, 12.58, 100000) #low
clh = cls(7, 10.2, 3.3, 13.16, 100000) #high

print '####',200,70
cll = cls(10, 6.0, 1.7, 10.54, 100000) #low
clh = cls(7, 10.2, 3.3, 10.96, 100000) #high

print '####',200,80
cll = cls(10, 6.0, 1.7, 8.95, 100000) #low
clh = cls(7, 10.2, 3.3, 8.1, 100000) #high

print '####',200,100
cll = cls(10, 6.0, 1.7, 6.06, 100000) #low
clh = cls(7, 10.2, 3.3, 4.94, 100000) #high

print '####',220,80
cll = cls(10, 6.0, 1.7, 8.64, 100000) #low
clh = cls(7, 10.2, 3.3, 9.56, 100000) #high

print '####',220,100
cll = cls(10, 6.0, 1.7, 7.4, 100000) #low
clh = cls(7, 10.2, 3.3, 6.92, 100000) #high

print '####',250,70
cll = cls(10, 6.0, 1.7, 7.39, 100000) #low
clh = cls(7, 10.2, 3.3, 12.86, 100000) #high

print '####',250,100
cll = cls(10, 6.0, 1.7, 6.59, 100000) #low
clh = cls(7, 10.2, 3.3, 9.77, 100000) #high

print '####',250,120
cll = cls(10, 6.0, 1.7, 5.32, 100000) #low
clh = cls(7, 10.2, 3.3, 6.75, 100000) #high

print '####',250,150
cll = cls(10, 6.0, 1.7, 2.96, 100000) #low
clh = cls(7, 10.2, 3.3, 3.15, 100000) #high

print '####',280,1
cll = cls(10, 6.0, 1.7, 5.64, 100000) #low
clh = cls(7, 10.2, 3.3, 13.37, 100000) #high

print '####',300,70
cll = cls(10, 6.0, 1.7, 4.47, 100000) #low
clh = cls(7, 10.2, 3.3, 11.51, 100000) #high

print '####',300,100
cll = cls(10, 6.0, 1.7, 4.35, 100000) #low
clh = cls(7, 10.2, 3.3, 10.01, 100000) #high

print '####',300,120
cll = cls(10, 6.0, 1.7, 4.1, 100000) #low
clh = cls(7, 10.2, 3.3, 8.27, 100000) #high

print '####',300,140
cll = cls(10, 6.0, 1.7, 3.83, 100000) #low
clh = cls(7, 10.2, 3.3, 6.95, 100000) #high

print '####',320,1
cll = cls(10, 6.0, 1.7, 3.48, 100000) #low
clh = cls(7, 10.2, 3.3, 11.85, 100000) #high

print '####',320,50
cll = cls(10, 6.0, 1.7, 3.59, 100000) #low
clh = cls(7, 10.2, 3.3, 10.99, 100000) #high

print '####',340,100
cll = cls(10, 6.0, 1.7, 2.64, 100000) #low
clh = cls(7, 10.2, 3.3, 8.51, 100000) #high

print '####',340,120
cll = cls(10, 6.0, 1.7, 2.83, 100000) #low
clh = cls(7, 10.2, 3.3, 7.63, 100000) #high

print '####',360,30
cll = cls(10, 6.0, 1.7, 2.23, 100000) #low
clh = cls(7, 10.2, 3.3, 8.84, 100000) #high

print '####',360,50
cll = cls(10, 6.0, 1.7, 2.16, 100000) #low
clh = cls(7, 10.2, 3.3, 8.54, 100000) #high

print '####',360,60
cll = cls(10, 6.0, 1.7, 2.18, 100000) #low
clh = cls(7, 10.2, 3.3, 8.62, 100000) #high

print '####',360,80
cll = cls(10, 6.0, 1.7, 2.20, 100000) #low
clh = cls(7, 10.2, 3.3, 8.20, 100000) #high

print '####',360,100
cll = cls(10, 6.0, 1.7, 2.13, 100000) #low
clh = cls(7, 10.2, 3.3, 7.92, 100000) #high

print '####',360,110
cll = cls(10, 6.0, 1.7, 2.15, 100000) #low
clh = cls(7, 10.2, 3.3, 7.36, 100000) #high

print '####',380,30
cll = cls(10, 6.0, 1.7, 1.71, 100000) #low
clh = cls(7, 10.2, 3.3, 7.62, 100000) #high

print '####',380,50
cll = cls(10, 6.0, 1.7, 1.75, 100000) #low
clh = cls(7, 10.2, 3.3, 7.5, 100000) #high

print '####',380,80
cll = cls(10, 6.0, 1.7, 1.72, 100000) #low
clh = cls(7, 10.2, 3.3, 7.08, 100000) #high
