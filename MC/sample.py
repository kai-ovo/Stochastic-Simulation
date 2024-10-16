from venv import create
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import standard_cauchy as npcauchy
from numpy.random import normal, uniform 
from utils import *
import pylab as pl

        
class sampler(object):

    def draw(self):
        """
        draw n samples 
        """
        pass

class cauchy(object):

    def __init__(self, gamma) -> None:
        super().__init__()
        self.gamma = gamma
    
    def draw(self,n):
        y = np.random.rand(n)
        x = self.gamma * np.tan(np.pi*(y-0.5))
        return x

    def sample_mean_seq(self, n, S):
        """
        n: #iid samples to compute sample mean
        S: #sample mean realizations
        """
        return [np.mean(self.draw(n)) for _ in range(S)]

class AR(object):
    """
    Simple rejection sampler 
    Biasing distribution is always standard normal
    """
    def __init__(self, p, C) -> None:
        """
        p: p(x) returns the unnormalized density value at x
            biasing density is standard normal
        C: scaling factor
        """
        super().__init__()
        self.g = normal 
        self.p = p
        self.C = C
    
    def draw(self, n):
        """
        draw n valid samples from rejection sampler

        output:
            samples: RS samples, length n
            t: totals samples been drawn
        """
        samples = []
        phi = lambda x : 1/np.sqrt(2*np.pi) * np.exp(-x**2/2) # standard normal density
        t = 0
        while True:
            t += 1
            y = self.g()
            w = self.p(y) / (self.C * phi(y))
            u = uniform()
            if u <= w:
                samples.append(y)
            if len(samples) == n:
                break
        
        return samples, t

    def draw_until(self, t):
        """
        draw a total (valid + invalid) of t samples
        output:
            samples: valid samples
            n: #valid samples
        """
        samples = []
        phi = lambda x : 1/np.sqrt(2*np.pi) * np.exp(-x**2/2) # standard normal density
        for _ in range(t):
            y = self.g()
            w = self.p(y) / (self.C * phi(y))
            u = uniform()
            if u <= w:
                samples.append(y)
        n = len(samples)
        return samples, n

    def target_estimate(self, t, N):
        """
        Estimate I = \int x^2*p(x) dx for N times; each time draws a total of t samples 

        t: stopping time for one draw of samples
        N: number of estimates 
        output: sample estimates, length N
        """
        estimates = []
        for _ in range(N):
            samples, _ = self.draw_until(t)
            samples = np.asarray(samples)
            estimate = np.mean(np.power(samples, 2))
            estimates.append(estimate)
        
        return estimates



class SNIS(object):
    """
    Simple self-normalized importance sampler
    Biasing distribution is always standard normal
    """
    def __init__(self) -> None:
        super().__init__()
        self.g = normal

    def draw(self):
        pass

    def target_estimate(self, t, N):
        pass


def sample_variance(samples):
    """
    compute sample variance with samples
    """
    samples = np.asarray(samples)
    mu = np.mean(samples)
    var = 0
    for s in samples:
        var += (s - mu)**2
    var *= 1/(len(samples)-1)    
    return var

def sample_mean_hist(sampler, n, S, figname):
    xs = sampler.sample_mean_seq(n,S)
    lbl = f'n={n}'
    plt.hist(xs)
    print(xs)
    print(max(xs))
    plt.title(f'n = {n}', fontsize=16)
    plt.ylabel('Probability')
    plt.xlabel('Sample Mean')
    plt.savefig(figname)

figpath = '../figs/'
create_path(figpath)

gamma = 1
n = 100
s = 10

sampler = cauchy(gamma)
samples = sampler.draw(n)
plt.figure()
plt.hist(samples, bins=50, range=(min(samples), max(samples)))
figname = figpath + 'cauchy-samples.png'
plt.savefig(figname)


np_samples = npcauchy(100)
plt.figure()
plt.hist(np_samples, bins=50, range=(min(samples), max(samples)))
figname = figpath + 'npcauchy-samples.png'
plt.savefig(figname)



# figname = f'cauchy-gamma{gamma}-n{n}-s{s}.png'
# figname = figpath + figname
# sample_mean_hist(sampler, n, s, figname)