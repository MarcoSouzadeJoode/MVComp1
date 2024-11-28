##########################################
## Toy code by M.Mapelli, Nov 21st 2024 ##
##########################################

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from numpy.fft import fft2,ifft2,rfft2,irfft2,fftfreq
import cmath



                   
""" auxiliary functions for Plummer model"""
def set_phi(x): #generate phi
    phi=x*2.*np.pi
    return phi 
def set_theta(x): #generate theta
    theta=1.0-x*2.
    theta=np.arccos(theta)
    return theta
def set_radius(a,x): #generate radius
    radius=a/(x**(-2./3.)-1.)**0.5
    return radius
def gocart(r,theta,phi): #change from spherical to cartesian coord.
    x=r*np.sin(theta)*np.cos(phi)
    y=r*np.sin(theta)*np.sin(phi)
    z=r*np.cos(theta)
    return x,y
""" end auxiliary functions for Plummer model (wait for Monte Carlo)"""



def generate_particles(N,maxi):
    """ Generate particles """
 
    if(plummer):
        np.random.seed(42)
        x=np.zeros(N,float)
        y=np.zeros(N,float)
        
        x1=np.random.rand(N)
        x2=np.random.rand(N)
        x3=np.random.rand(N)
        phi=set_phi(x1)
        theta=set_theta(x2)
        a=0.1 #plummer scale radius
        radius=set_radius(a,x3)
        x,y=gocart(radius,theta,phi)

        print(min(x),max(x))
        print(min(y),max(y))
        print(len(x),len(y))

        # remove particles outside boundaries
        # (alternative: use periodic conditions)
        d=np.where(abs(x)<maxi)
        x,y=x[d[0]],y[d[0]]
        d=np.where(abs(y)<maxi)
        x,y=x[d[0]],y[d[0]]
        print(len(x))

        N=len(x)              
    else:
        x,y=(2.*np.random.random(N))-1.,(2.*np.random.random(N))-1.
    
    m=np.ones(len(x))
    mtot=sum(m)
    
    return(x,y,m,mtot,N)






def density_NGP(x,y,m,Np,h,mini):
    """
    Compute the density on a grid using the NGP algorithm.

    Parameters:
    - x, y: x, y positions of particles.
    - m: masses of particles.
    - Np: int, the number of grid points along each axis (NpxNp grid).
    - h: length of a cell
    - mini: minimum value of the box size
    
    Returns:
    - density: 2D numpy array (Np x Np) containing the density values.
    """
    
    rhop=np.zeros([Np,Np],float)
    
    for l in range(len(x)):
        # Calculate the grid indices for particles based on their positions
        i = ((x[l] - mini) / h).astype(int)
        j = ((y[l] - mini) / h).astype(int)

        # Increment the density at the grid point
        rhop[i, j] += m[l]
    

    rhop=rhop/(h*h)

    print("Density field computed")
    plt.title("Density Field (NGP)")
    plt.imshow(rhop,norm="log",vmin=1e4,vmax=3e6,extent=(mini,maxi,mini,maxi))
    plt.xlabel("x")
    plt.ylabel("y")
    cbar=plt.colorbar()
    cbar.set_label("Density")
    plt.savefig("density.pdf")
    plt.close()
    return rhop





def density_CIC(x,y, m, Np, box_size, mini):
    """
    Compute the density on a grid using the Cloud-In-Cell (CIC) algorithm.

    Parameters:
    - x,y: Each row contains the (x, y) position of a particle.
    - m: mass of each particle.
    - Np: int, the number of grid points along each axis (NpxNp grid).
    - box_size: float, the size of the simulation box (assumed square).
    - mini: minimum value of the box size
    
    Returns:
    - density: 2D numpy array (Np x Np) containing the density values.
    """
    rhop = np.zeros([Np, Np],float)
    h = box_size / Np  # Grid spacing
    print("cell size=",h)


    for l in range(len(x)):
        """
        TO BE FILLED IN
        ADD YOUR OWN CODE HERE TO CALCULATE
        DENSITY ON THE MESH WITH CIC METHOD
        """
        print("TO BE FILLED IN")

        
    return rhop

    

def calc_pot(rhop,Np,maxi,mini):
    """ calculate potential with FFT """
    
    h = (maxi - mini) / Np

    # write the frequencies
    kx = np.fft.fftfreq(Np, d=h) * 2 * np.pi
    ky = np.fft.fftfreq(Np, d=h) * 2 * np.pi
    KX, KY = np.meshgrid(kx, ky)
    k2 = KX**2 + KY**2
    ### NOTE: it is very important that the frequencies
    ### are centered around zero.
    ## This is ensured by fft.fftfreq
    """
    TO BE FILLED IN
    
    ADD YOUR OWN CODE HERE TO CALCULATE
    POTENTIAL ON THE MESH WITH FOURIER TRANSFORM
    IN THE CURRENT VERSION THE CODE READS
    A PRE-COMPUTED POTENTIAL FROM THE FILE "fft.dat"
    NOTE: THE EXAMPLE WORKS ONLY WITH A MESH OF 100x100!
    HINT: USE FFT2 and THEN IFFT2 TO GO BACK TO REAL SPACE
    HINT2: BE CAREFUL NOT TO DIVIDE BY ZERO FOR kx=0, ky=0
    """
    phi=np.zeros([Np,Np])

    
    """
    TO BE FILLED IN:
    PLOT POTENTIAL
    """

    return phi





if __name__ == "__main__":

    plummer=True
    # if True calculates plummer model, else uniform density field



    N=100000   #number of particles 

    Np=100    #number of grid points

    mini=-1.  # min of the box
    maxi=1.   # max of the box

    h=(maxi-mini)/Np #cell length

    
    xp=np.linspace(mini,maxi,int(N)) 
    yp=np.copy(xp)


    ## generate particles
    x,y,m,mtot,N=generate_particles(N,maxi)
    print(N)


    ## choose how to calculate the density on the mesh
    mesh="NGP"

    ## calculate density on the mesh
    if(mesh=="NGP"):
        rhop=density_NGP(x,y,m,Np,h,mini)
    elif(mesh=="CIC"):
        rhop=density_CIC(x,y, m, Np, (maxi-mini), mini)
    else:
        print("Error, no way to produce a density field")
        exit()

    ## calculate potential solving Poisson's eq. via Fourier method
    phi=calc_pot(rhop,Np,maxi,mini)


    

