import matplotlib.pyplot as plt
import numpy as np
import map_equ
from ipfnpytools.closest import closest
import ipywidgets as widgets


class Equilibrium:
    
    def __init__(self, equilibrium, countour, axes):
        self._eq = equilibrium 
        self._cs = countour
        self._ax = axes
        self._contourlevels = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.02, 1.04, 1.06])
        
    def __call__(self, time):
        index = closest(self._eq.t_eq, time)
        for col in self._cs.collections:
            col.remove()
            
        self._cs = self._ax.contour(
            self._eq.Rmesh, self._eq.Zmesh, 
            np.sqrt((self._eq.pfm[:,:,:] - self._eq.psi0)/(self._eq.psix-self._eq.psi0))[:,:,index].T,
            levels=self._contourlevels, linewidths=1.2, cmap='viridis_r'
        ) 

def plot_equilibrium(shot, time, equilibrium='EQH', axes=None):

    eq = map_equ.equ_map()
    eq.Open(shot, equilibrium)
    #Populate the fields
    eq.read_pfm()
    eq.read_scalars()
    #Get the index closest to the requested time
    index = closest(eq.t_eq, time)

    if axes is None:
        ax = plt.gca()
    else:
        ax = axes
        
    fig = ax.get_figure()
    
    contourlevels = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.02, 1.04, 1.06])

    cs = ax.contour(eq.Rmesh, eq.Zmesh, np.sqrt((eq.pfm[:,:,:] - eq.psi0)/(eq.psix-eq.psi0))[:,:,index].T,
                     levels=contourlevels, linewidths=1.2, cmap='viridis_r') 
    
    equil = Equilibrium(eq, cs, ax)
    
    def update(i):
        equil(i)
        
    slider = widgets.FloatSlider(
        value=eq.t_eq[int(len(eq.t_eq)/2)],
        min=eq.t_eq[0],
        max=eq.t_eq[-1],
        step=eq.t_eq[1]-eq.t_eq[0],
        description='Time',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='.6f',
    )

    widgets.interact(update, i=slider)
    
    return equil