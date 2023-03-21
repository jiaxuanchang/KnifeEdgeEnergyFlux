from matplotlib import pylab


def knifeedge(delta,nmodes):
    """" an,Fn=knifeedge(delta,nmodes)
    
    Compute the amplitude of vertical modes and energy flux of the internal wave solution on a knife-edge ridge.
    
    Parameters:
    ----------- 
    delta : is relative height (obstacle height/water depth)
    nmodes : number of modes to return.
         
    Returns:
    --------
    an : amplitude of vertical modes
    Fn : =Fknife/F0=n^{-1}an^2, nondimensional modes energy flux. eq(21) in St. Laurent et al., 2003
             
    Notes: 
    ------

             
    JX. Chang  Mar 20, 2023
    """
    
    import numpy as np
    
    
    n=np.arange(1,nmodes+1)
    n=np.array([n]).T
    m=np.arange(nmodes)[np.newaxis]

    d = np.pi*(1 - delta)
    Amn = (n * np.sin(n*d) * np.cos(m*d) -
           m * np.cos(n*d) * np.sin(m*d)) + n
    Amn += (-n* np.cos(n*d) * np.cos(m*d) -
            m * np.sin(n*d) * np.sin(m*d)) 
    Amn *= 1/(m**2-n**2) #(17)

    for ni in range(1, nmodes):
        Amn[ni-1, ni] = (ni*np.pi*delta - np.sin(ni*d) * np.cos(ni*d) - np.sin(ni*d)**2) / 2 / ni  #(19)

    cm=np.sin((m*np.pi)*(1-delta))/m #(18)
    cm=cm.T 
    cm[0]=-np.pi*delta #(20)
    A_inv= np.linalg.inv(Amn.T)
    an=A_inv.dot(cm) #an=A-1mn cm 

    Fn=an**2/n #(21)

    return an,Fn
