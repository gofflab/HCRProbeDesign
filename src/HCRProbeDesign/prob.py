#!/usr/bin/env python
import math,operator,random,sys
import numpy as np

#######
#Probability Tools for DNA sequence analysis
#######
def snr(observed,expected):
    '''
    This function calculates the signal-to-noise ratio of an observed value compared to an expected
    value
    
    :param observed: the actual number of observed counts
    :param expected: the expected value of the distribution
    :return: The SNR value
    '''
    return observed/expected

def zscore(observed,expected):
    '''
    The z-score is the number of standard deviations from the mean a data point is
    
    :param observed: the actual observed count of the test statistic
    :param expected: the expected number of counts under the null hypothesis
    :return: The z-score for the observed value.
    '''
    return (observed-expected)/math.sqrt(expected)

def which_bin(bins, x, safe=0):
    """
    # if we're interested in binning x with boundaries
    # 0, 5, 10, 15
    # then it will return which boundary it belongs in.
    # if x<0: -1
    # if 0<=x<5: 0
    # if 5<=x<10: 1
    # if 10<=x<15: 2
    # if x>=15: 3
    """
    if x<bins[0]: return -1
    for i in range(1,len(bins)):
        if x<bins[i]: return i-1
    if safe and x==bins[-1]: return len(bins)
    return len(i)+1

def cumulative_sum(quality):
    '''
    Given a list of numbers, return a list of the cumulative sum of the numbers in the list
    
    :param quality: a list of integers representing the quality of the video
    :return: The cumulative sum of the quality scores.
    '''
    if not quality: return quality
    sum_q = quality[:]
    for i in range(1,len(quality)):
        sum_q[i] = sum_q[i-1]+quality[i]
    return sum_q

def frequency_dic(seq):
    """Generates dictionary of k,v='nucleotide':'frequency' from seq"""
    dic = {}
    bases = ['A','C','G','T']
    seq=seq.upper()
    for b in bases:
        dic[b]=seq.count(b)/float(len(seq))
    return dic

def pick_one(dic):
    '''
    Given a dictionary of key:value pairs, 
    where the values are cumulative sums of the keys, 
    pick a key with probability proportional to its value
    
    :param dic: the dictionary to pick from
    :return: The character with the probability of occurrence specified by the dictionary.
    '''
    # {'A': .18, 'C': .32, 'G': .32, 'T': .18}
    # will generate A with probability .18 and so on
    items = dic.items()
    cums = cumulative_sum(cget(items,1))
    if 1: #debug:
        #print cums
        x = random.uniform(0,cums[-1])
        bin = which_bin(cums, x, safe=1)
        #print "%s is in bin %s and char %s. items=%s"%(
        #    x,bin,items[bin][0],items)
        return items[bin+1][0]
    else:
        return items[which_bin(cums, random.uniform(0,cums[-1]), safe=1)][0]

def pick_many(dic, n):
    '''
    Given a dictionary of items and their probabilities, 
    pick n items according to those probabilities
    
    :param dic: a dictionary of the form {'A': .18, 'C': .32, 'G': .32, 'T': .18}
    :param n: number of random choices
    :return: A list of n characters.
    '''
    # {'A': .18, 'C': .32, 'G': .32, 'T': .18}
    # will generate A with probability .18 and so on
    items = dic.items()
    cums = cumulative_sum(cget(items,1))
    choices = []
    for i in range(0,n):
        x = random.uniform(0,cums[-1])
        bin = which_bin(cums, x, safe=1)
        choices.append(items[bin+1][0])
    return choices

def gaussian(x,mu,sigma):
    """
    Evaluate N(mu,sigma) at x.
    where N(mu,sigma) is a gaussian of mean mu and stdev sigma
    """

    return ( (1.0/math.sqrt(2*math.pi*sigma)) * (math.e**(-((x-mu)**2)/(2*sigma**2))))

def make_gaussian(mu,sigma):
    """
    usage:
    N2_3 = make_gaussian(2,3)
    N2_3(4) -> gaussianN(2,3) evaluated at 4
    """
    return lambda x,mu=mu,sigma=sigma: ( (1.0/math.sqrt(2*math.pi*sigma)) * (math.e**(-((x-mu)**2)/(2*sigma**2))))

def make_adder(n):
    """
    usage:
    Add2=make_adder(2)
    Add2(3) -> 5
    """
    return lambda x,n=n: x+n

#############
#Math Primitives
#############
loge_2 = math.log(2)

def avg(l,precise=0):
    '''
    Return the average of a list of numbers
    
    :param l: The list of numbers to be averaged
    :param precise: If True, the average will be rounded to the nearest integer, defaults to 0
    (optional)
    :return: The average of the list of numbers.
    '''
    if not l: return 0
    if precise:
        return reduce(operator.add,l,0)/float(len(l))
    else:
        return reduce(operator.add,l,0)/len(l)

def movavg(s, n):
    ''' returns an n period moving average for the time series s

        s is a list ordered from oldest (index 0) to most recent (index -1)
        n is an integer

        returns a numeric array of the moving average
    '''
    s = np.array(s)
    c = np.cumsum(s)
    return (c[n-1:] - c[:-n+1]) / float(n)


def median(l):
    if not l: return None
    l = my_sort(l)
    if len(l)%2: return my_sort(l)[len(l)/2]
    else: return (l[len(l)/2]+l[len(l)/2-1])/2.0

def stdev(l, failfast=1):
    return math.sqrt(variance(l,failfast=failfast))

def variance(l,failfast=1):
    if (not l) or len(l)==1:
        if failfast: raise "tools.variance: Not enough samples.  Need >= 2, got %s"%len(l)
        else: return 0#'N/A'
    m = avg(l,1)
    s = 0
    for i in l:
        s = s + (i-m)*(i-m)
    return s / (len(l)-1)

def log2(x):
    #converting bases: log_a(b) = log_c(b)/log_c(a)
    #i.e. log_2(x) = log_e(2)/log_e(x) = log_10(2)/log_10(x)
    return math.log(x)/float(loge_2)

def log_k(x,k):
    return math.log(x)/math.log(k)

def prob2score(prob):
    '''
    Given a probability, return the corresponding score
    
    :param prob: the probability of the event
    :return: The score for each SNP.
    '''
    #1/100 -> 20
    try:
        return -10*float(math.log10(float(prob)))
    except:
        return -1

def p2bits(p):
    """Takes p-value and returns negative log2"""
    return -log2(p)

def factorial(n):
    result = 1
    for i in range(n,0,-1):
        #print i
        result = result * i
    return result

###########
#Poisson
###########
def poisson_expected(rate):
    '''
    The function poisson_expected() takes a rate and prints out the probability of a number of events
    occuring in a given time period
    
    :param rate: the average number of events per interval
    '''
    for x in range(1,50,1):
        p = poisson(rate,x)
        print(f"{x}\t{p}\t{12000000*p}")

def poisson(rate, x):
    """Returns the probability of observing a count of x"""
    return math.exp(-rate)*(rate**x)/factorial(x)

######################
#Binomial Distribution
#######################
def binomial_likelihood_ratio(ps,k,n):
    '''
    The likelihood ratio is the log of the likelihood of the hypothesis being tested divided by the
    likelihood of the null hypothesis
    
    :param ps: the two hypotheses
    :param k: number of successes
    :param n: number of trials
    :return: The log-likelihood ratio.
    '''
    # p[0] is the null hypothesis
    # p[1] is the hypothesis being tested
    assert(len(ps)==2)
    likelihoods = []
    for p in ps:
        likelihoods.append(binomial(p,k,n))
    #i = argmax(likelihoods)
    #p = likelihoods[i] / sum(likelihoods)
    #return p
    if likelihoods[0]: return np.log(likelihoods[1]) / likelihoods[0]
    else:
        print(f"Warning: likelihood ratio set to sys.maxint.  p(H1)={p[1]}, p(H0)=0")
        return sys.maxint

def binomial_log_likelihood_ratio(ps,k,n):
    return log_binomial(ps[1],k,n) - log_binomial(ps[0],k,n)

def log_binomial(p,k,n):
    # the log probability of seeing exactly k successes in n trials
    # given the probability of success is p
    return log_n_choose_k(n,k)+math.log(p)*k+math.log(1-p)*(n-k)

def binomial(p,k,n):
    # probability of seeing exactly k successes in n trials, given
    # the probability of success is p
    #return n_choose_k(n,k)*(p**k)*((1-p)**(n-k))
    return n_choose_k(n,k)*(p**k)*((1-p)**(n-k))

def cumBinomial(p,k,n):
    #Returns the cumulative probability from the binomaial distribution
    Pval = 0.0
    for j in range(0,k+1):
        Pval+=binomial(p,j,n)
    return Pval

def n_choose_k(n,k):
    # (n k) = n! / (k! (n-k)!)
    #
    #         n*(n-1)*(n-2)*....*(n-k+1)
    #       = --------------------------
    #              k*(k-1)*...*1
    assert(k<=n)
    k = min(k, n-k)
    nominator   = range(n,n-k,-1)
    denominator = range(k,0,-1)

    result = 1.0
    for nom, den in map(None, nominator, denominator):
        result = (result * nom) / den
        #result = result*nom
        #print result
        #result = result/den
        #print result

    return result

def log_n_choose_k(n,k):
    '''
    Given n and k, return the log of n choose k
    
    :param n: number of items in the list
    :param k: the number of successes
    :return: The log of the number of ways to choose k objects from n objects.
    '''
    # (n k) = n! / (k! (n-k)!)
    #
    #         n*(n-1)*(n-2)*....*(n-k+1)
    #       = --------------------------
    #              k*(k-1)*...*1
    assert(k<=n)
    k = min(k, n-k)
    nominator   = range(n,n-k,-1)
    denominator = range(k,0,-1)

    result = 0
    for nom, den in map(None, nominator, denominator):
        result = (result + math.log(nom)) - math.log(den)
    return result

#################
#Dictionary Tools
#################
def cget(diclist, key, strict=1):
    # cross_get was: gather(diclist,key)
    # gathers the same key from a list of dictionaries
    # can also be used in lists

    # input: a list of dictionaries all of which contains key
    # output: a list of elements d[key] for each d in diclist
    if strict:
        # return map(lambda d,key=key: d[key], diclist)
        result = [None]*len(diclist)
        for i in range(0,len(diclist)):
            result[i] = diclist[i][key]
        return result
    else:
        results = []
        for dic in diclist:
            if dic and generic_has_key(dic,key):
                results.append(dic[key])
        return results
