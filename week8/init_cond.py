import numpy as np







def scale(m1,m2,a,ecc,fase,G):
    Gcgs  = 6.667e-8   #G in cgs
    pc = 3.086e18   #pc in cgs
    msun = 1.989e33 #msun in cgs
    rsun = 6.95e10  #rsun in cgs
    AU = 1.49e13    #AU in cgs
    L_scale = 1./AU                          #1cm/1pc
    M_scale = 1./msun                        #1g/1M_sun
    T_scale = (Gcgs * L_scale**3/M_scale)**0.5  #tscale in nbody units
    V_scale = (L_scale/T_scale)              # vel in nbody units
    return (L_scale,M_scale,T_scale,V_scale)

def gen_pos(m1,m2,a,ecc,fase,G):
    mtot=m1+m2
    x=np.zeros([2,3],float)

    x[0,0] = -m2/mtot * a * (1.-ecc**2) * np.cos(fase)/(1.+(ecc * np.cos(fase)))
    x[0,1] = -m2/mtot * a * (1.-ecc**2) * np.sin(fase)/(1.+(ecc * np.cos(fase)))
    x[0,2] = 0.0

    x[1,0] = m1/mtot * a * (1.-ecc**2) * np.cos(fase)/(1.+(ecc * np.cos(fase)))
    x[1,1] = m1/mtot * a * (1.-ecc**2) * np.sin(fase)/(1.+(ecc * np.cos(fase)))
    x[1,2] = 0.0
    return x

def gen_vel(m1,m2,a,ecc,fase,G):
    mtot=m1+m2
    v=np.zeros([2,3],float)
    
    v[0,0] = -m2/mtot * (ecc * np.cos(fase) /(1+ ecc * np.cos(fase)) - 1.) * np.sin(fase) * (1+ ecc * np.cos(fase)) * (G * mtot / (a*(1- ecc*ecc)))**0.5
    v[0,1] = -m2/mtot * (ecc * (np.sin(fase))**2 /(1+ ecc * np.cos(fase)) + np.cos(fase)) * (1+ ecc * np.cos(fase)) * (G * mtot / (a*(1- ecc*ecc)))**0.5
    v[0,2] = 0.0

    v[1,0] = m1/mtot * (ecc * np.cos(fase) /(1+ ecc * np.cos(fase)) - 1.) * np.sin(fase) * (1+ ecc * np.cos(fase)) * (G * mtot / (a*(1- ecc*ecc)))**0.5   
    v[1,1] = m1/mtot * (ecc * (np.sin(fase))**2 /(1+ ecc * np.cos(fase)) + np.cos(fase)) * (1+ ecc * np.cos(fase)) * (G * mtot / (a*(1- ecc*ecc)))**0.5
    v[1,2] = 0.0
    return v
