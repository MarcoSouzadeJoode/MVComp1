import numpy as np

def advect_conservative_donorcell(xi,rho,ui,dt):
    rhonew             = np.zeros_like(rho)
    rhonew[0]          = rho[0]
    rhonew[-1]         = rho[-1]
    F                  = np.zeros(len(ui))
    mask               = ui[1:-1]>0
    F[1:-1][mask]      = rho[:-1][mask] * ui[1:-1][mask]
    mask               = ui[1:-1]<0
    F[1:-1][mask]      = rho[1:][mask]  * ui[1:-1][mask]
    rhonew[1:-1]       = rho[1:-1] - dt * ( F[2:-1] - F[1:-2] ) / ( xi[2:-1] - xi[1:-2] )
    return rhonew

def advect_translation_upwind(xc,rho,uc,dt):
    rhonew             = np.zeros_like(rho)
    rhonew[0]          = rho[0]
    rhonew[-1]         = rho[-1]
    mask               = uc[1:-1]>0
    eps                = uc[1:-1]*dt/(xc[1:-1]-xc[:-2])
    rhonew[1:-1][mask] = eps[mask] * rho[:-2][mask] + (1-eps[mask]) * rho[1:-1][mask]
    mask               = uc[1:-1]<0
    eps                = uc[1:-1]*dt/(xc[1:-1]-xc[2:])
    rhonew[1:-1][mask] = eps[mask] * rho[2:][mask] + (1-eps[mask]) * rho[1:-1][mask]
    return rhonew
