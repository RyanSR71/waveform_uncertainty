import numpy as np
import bilby
import logging
import tqdm
from ..wfu.utils import ProgressBar
import matplotlib.pyplot as plt
from pesummary.gw.file.strain import StrainData
from pesummary.io import read



def match(signal,data,duration,**kwargs):

    PSDs = kwargs.get('PSDs',None)

    if PSDs is None:
        PSDs = np.ones(len(signal))

    if len(signal) != len(data):
        raise Exception('Signal and Data do not have the same shape!')
    
    signal_match = np.sqrt(bilby.gw.utils.matched_filter_snr(signal,signal,PSDs,4))
    data_match = np.sqrt(bilby.gw.utils.matched_filter_snr(data,data,PSDs,4))
    normalized_match = np.abs(bilby.gw.utils.matched_filter_snr(signal,data,PSDs,4)/(signal_match*data_match))
    
    return normalized_match



def match_plot(ppE_waveform_generator,GR_waveform_generator,injection,beta_tildes,delta_epsilon_tildes,**kwargs):
    PSDs = kwargs.get('PSDs',None)
    save = kwargs.get('save',False)
    path = kwargs.get('path',None)
    levels = kwargs.get('levels',np.linspace(0,1,21))
    
    GR_waveform = GR_waveform_generator.frequency_domain_strain(parameters=injection)['plus']
    
    if len(beta_tildes) != len(delta_epsilon_tildes):
        raise Exception(f'beta_tilde array and delta_epsilon_tilde array do not have the same length! {len(beta_tildes)}, {len(delta_epsilon_tildes)}')
        
    resolution = len(beta_tildes)
    match_grid = np.zeros([resolution,resolution])
    
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    log.addHandler(ProgressBar())
    with tqdm.tqdm(total=resolution**2, desc="Generating Contour Plot",position=0,leave=True) as pbar:
        for i in range(resolution):
            for j in range(resolution):

                injection['beta_tilde'] = beta_tildes[i]
                injection['delta_epsilon_tilde'] = delta_epsilon_tildes[j]

                ppE_waveform = ppE_waveform_generator.frequency_domain_strain(parameters=injection)['plus']

                match = match(ppE_waveform,GR_waveform,ppE_waveform_generator.duration,PSDs=PSDs)

                match_grid[j,i] = match

                pbar.update(1)
                
    X,Y = np.meshgrid(beta_tildes,delta_epsilon_tildes)

    plt.contourf(X,Y,match_grid,levels=levels)
    plt.colorbar(label=r'$\mathrm{match}(\tilde{h}_\mathrm{ppE},\tilde{h}_\mathrm{GR})$')
    
    plt.grid(False)
    plt.xlim(beta_tildes[0],beta_tildes[-1])
    plt.ylim(delta_epsilon_tildes[0],delta_epsilon_tildes[-1])
    plt.xlabel(r'$\tilde\beta$')
    plt.ylabel(r'$\delta\tilde\epsilon$')
    
    M = np.round(bilby.gw.conversion.generate_mass_parameters(injection)['total_mass'])
    plt.title(r'$M=val1\ [\mathrm{M}_\odot],\ b=val2$'.replace('val1',str(M)).replace('val2',str(injection['b'])))
    
    if save is True:
        if path is None:
            plt.savefig(f"match_plot_b={injection['b']}_M={M}.png")
        else:
            plt.savefig(f"{path}/match_plot_b={injection['b']}_M={M}.png")
    
    plt.show()
