import numpy as np
from scipy.fft import fft, ifft

BDF1 = 0
BDF2 = 1
TR = 2
TTR = 3

def dlt_lms(z,lms):
    if lms == BDF1:
        return 1-z
    elif lms == BDF2:
        return 1-z+(1/2)*(1-z)**2
    elif lms == TR:
        return 2*(1-z)/(1+z)
    elif lms == TTR:
        cs = [0.893817850529318,   0.684154908023834,   0.629642997466429]
        return 1-z+(1/2)*(1-z)**2+(1/4)*cs[0]*(1-z)**3+(1/8)*cs[1]*(1-z)**4+(1/16)*cs[2]*(1-z)**5
    else:
        return "No such LMS method"

def evalCQ(g, Ks, N, T, lms = TR):
    dt = T/N
    lam = 10**(-14/(2*N+1))
    lams = lam**np.arange(N+1)
    ts = np.linspace(0,T,N+1)
    zs = lam*np.exp(-2*np.pi*1j*np.arange(N+1)/(N+1))
    u = (1/lams)*ifft(Ks(dlt_lms(zs,lms)/dt)*fft(lams*g(ts)))
    
    return u


