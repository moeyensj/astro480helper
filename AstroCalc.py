import numpy as np
import matplotlib.pyplot as plt

GAIN = 2.3 # e-/ADU

def snrCalc(N_star, npix, N_s, N_d, N_r, t):
    return N_star*t / np.sqrt(N_star*t + npix*(N_s*t + N_d*t + N_r**2))

def expCalc(snr, N_star, npix, N_s, N_d, N_r):
    A = N_star**2
    B = (-(snr)**2)*(N_star + npix*(N_s + N_d))
    C = (-(snr)**2)*(npix*N_r**2)
    return (-B + (B**2 - 4*A*C)**0.5) / (2*A)

def pixelCalc(r):
    return int(np.pi*r**2)

def darkCurrentCalc(T, A, I, E):
    k = 8.6175 * 10**-5
    return (2.5*10**15)*(A*I*T**5)*np.exp(-E/(2*k*T))

def cToK(temp_C):
    return temp_C + 273.15

def kToC(temp_K):
    return temp_K - 273.15

def eToADU(e):
    return e / GAIN

def aduToE(adu):
    return adu * GAIN

def fluxRatioCalc(m1, m2):
    return 10.0**((m1 - m2)/-2.5)

def apparentMagCalc(m0,f1,f0):
    return -2.5*np.log10(f1/f0) + m0

def linearityPlot(exposure,counts):

    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(10,6)

    ax.scatter(exposure,counts)
    ax.set_xlabel(r'Exposure Time, $s$',size=14)
    ax.set_ylabel(r'Counts, $ADUs$',size=14)
    ax.set_title(r'CCD Linearity for ST-10XME',size=14)
    #ax.set_xlim(0,125)
    #ax.set_ylim(0,60000)
    #ax.axvline(90,linestyle='--')
    #ax.text(30,47000,'Linear Region',size=15)
    #ax.text(92,18000,'Non-Linear Region',size=15)
    #ax.grid()

    #plt.savefig('CCDlinearity.jpg')

    return