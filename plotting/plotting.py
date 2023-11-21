import os
import numpy as np
import json
from skimage.io import imread
from skimage.io import imsave
import pandas
import warnings
import time

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

# Optional, uncomment to use color maps from Fabio Crameri:
# https://doi.org/10.5281/zenodo.1243862
#from cmcrameri import cm




def plot_scatter(h_ax,
                 v_ax,
                 save_name,
                 data_base_path,
                 data_in,
                 sampling,
                 v_ax_lim=[],
                 h_ax_lims=[],
                 save_plot=True,
                 h_label='',
                 v_label='',
                 plot_title=''):
    
    """
    Scatter plot of bead sizes (FWHM in one dimension) vs. 
    bead position (spatial coordinate in one dimension).

    Parameters
    ----------
    h_ax : string matching a column in the input dataframe
        data on the horizontal axis - use a bead position column, e.g. 'x_center'
        
    v_ax : string matching a column in the input dataframe
        data on the horizontal axis - use a bead size column, e.g. 'FWHM_x'

    save_name : string, file name for saved plot

    data_in : input pandas dataframe of bead positions and 
       measurements (generated below)

    save_plot : (optional) bool, determines whether plot is saved

    h_label : (optional) string, custom label for the horizontal axis.
        
    v_label : (optional) string, custom label for the vertical axis.
    
    plot_title : (optional) string, custom title displayed on plot.
    """  
    
    # plot data
    fig, ax = plt.subplots()
    plt.scatter(data_in[h_ax]*sampling,
                data_in[v_ax],
                marker='.')
    
    # update labels
    if h_label == '':
        h_label = h_ax + ' (μm)'
    plt.xlabel(h_label)

    if v_label == '':
        v_label = v_ax + ' (μm)'
    plt.ylabel(v_label)
    
    plt.title(plot_title)

    plt.ylim(bottom=0)
    if v_ax_lim != []:
        plt.ylim(top=v_ax_lim)
    if h_ax_lims != []:
        plt.xlim(h_ax_lims[0],h_ax_lims[1])
    
    
    # save plot
    if save_plot:
        print('Saving plot')
        plt.savefig(os.path.join(data_base_path,save_name))
    else:
        print('Plot NOT saved!')
        
    return fig, ax


def plot_histo(h_ax,
               h_max_lim,
               save_name,
               data_base_path,
               data_in,
               anno='',
               save_plot=True,
               h_label='',
               v_label='',
               percentile_lim=[0.1,99.9]):
    
    """
    Histogram of bead sizes (FWHM in one dimension)

    Parameters
    ----------
    h_ax : string matching a column in the input dataframe
        data on the horizontal axis - use a bead size column, e.g. 'FWHM_x'
        
    h_max_lim : float, upper limit of values on the horizontal axis
        
    save_name : string, file name for saved plot

    data_in : input pandas dataframe of bead positions and 
       measurements (generated below)

    save_plot : (optional) bool, determines whether plot is saved

    h_label : (optional) string, custom label for the horizontal axis.
        
    v_label : (optional) string, custom label for the vertical axis.

    percentile_lim : (optional) two element float list, indicates
        the percentile range of data to plot
    """  
    
    
    # drop data outside percentile limits. To disable, set percentile_lim=[0,100]
    min_cutoff = np.percentile(data_in[h_ax],percentile_lim[0])
    max_cutoff = np.percentile(data_in[h_ax],percentile_lim[1])
    data_in = data_in.drop(data_in[(data_in[h_ax] < min_cutoff) | (data_in[h_ax] > max_cutoff)].index)

    # plot data
    fig, ax = plt.subplots()
    bins = np.linspace(0, h_max_lim, 31)
    plt.hist(data_in[h_ax], bins = bins, edgecolor = 'black',linewidth=0.5)
    
    # update labels
    if h_label == '':
        h_label = h_ax + ' (μm)'
    plt.xlabel(h_label)
    
    if v_label == '':
        v_label = 'N'
    plt.ylabel(v_label)
    
    data_med = round(np.median(data_in[h_ax]),3)
    data_std = round(np.std(data_in[h_ax]),3)
    
    plt.title(r'Median = ' + str(data_med) + r' μm   Stand. Dev. = ' + 
              str(data_std) + ' μm')
    
    # add annotation
    plt.annotate(anno, xy=(0.0*plt.xlim()[1], 0.9*plt.ylim()[1]))
    
    # save plot
    if save_plot:
        print('Saving plot')
        plt.savefig(os.path.join(data_base_path,save_name))
    else:
        print('Plot NOT saved!')
    
    return fig, ax

def plot_scatter_heat(h_ax,
                      v_ax,
                      color_ax,
                      save_name,
                      data_base_path,
                      data_in,
                      sampling,
                      save_plot=True,
                      h_label='',
                      v_label='',
                      c_label='',
                      plot_title='',
                      clim=[],
                      hlim=[],
                      vlim=[],
                      hticks=[],
                      vticks=[],
                      percentile_lim=[0.1,99.9],
                      base_map='',
                      equal=False
                     ):
        
    """
    Scatter plot heat map. Color shows bead sizes (FWHM in one dimension),
    horizontal and vertical position shows bead position (spatial coordinate
    in two dimensions).

    Parameters
    ----------
    h_ax : string matching a column in the input dataframe
        data on the horizontal axis - use a bead position column, e.g. 'x_center'
        
    v_ax : string matching a column in the input dataframe
        data on the vertical axis - use a position column, e.g. 'y_center'

    save_name : string, file name for saved plot

    data_in : input pandas dataframe of bead positions and 
       measurements (generated below)

    save_plot : (optional) bool, determines whether plot is saved.

    h_label : (optional) string, custom label for the horizontal axis.
        
    v_label : (optional) string, custom label for the vertical axis.
    
    c_label : (optional) string, custom label for colorbar.
    
    plot_title : (optional) string, custom title displayed on plot.
    
    clim : (optional) two element float list, indicates colormap range.
        Data above/below the range is still displayed but will all be
        the same color as the max/min of the range (use percentile_lim
        to truncate data). Omit to autoscale the colormap. Autoscaling
        is meant for use with a diverging colormap. The median data value
        will be mapped to the center of the colormap (i.e. neutral color)
        and the most extreme value (max or min) will be mapped to one end
        of the colormap. The rest of the data is mapped linearly based
        on these values.
        
    hlim : (optional) two element float list, indicates horizontal axis range.
    
    vlim : (optional) two element float list, indicates vertical axis range.
    
    hticks : (optional) float list, custom tick marks for horizontal axis.
    
    vticks : (optional) float list, custom tick marks for vertical axis.
    
    percentile_lim : (optional) two element float list, indicates
        the percentile range of data to plot
        
    base_map : (optional) custom colormap
        For matplotlib colormaps, use e.g. <base_map=mpl.colormaps['coolwarm']>
        For cmcrameri colormaps, use e.g. <base_map=cm.batlow>
        
    equal : (optional) bool, use the scale for horizontal and vertical axes
        
    """  
    
    if base_map == '':
        # to use built-in colormaps, uncomment line below
        base_map = mpl.colormaps['coolwarm']
        # base_map = cm.batlow
    
    # drop data outside percentile limits. To disable, set percentile_lim=[0,100]
    min_cutoff = np.percentile(data_in[color_ax],percentile_lim[0])
    max_cutoff = np.percentile(data_in[color_ax],percentile_lim[1])
    data_in = data_in.drop(data_in[(data_in[color_ax] < min_cutoff) | (data_in[color_ax] > max_cutoff)].index)

    # use autoscaled colormap (recommended)
    if clim == []:

        # find colormap range to center map at median data value
        med_val = np.median(data_in[color_ax])
        max_val = np.max(data_in[color_ax])
        min_val = np.min(data_in[color_ax])        

        reach = np.max([max_val - med_val, med_val - min_val])

        cmap_lower = 0.5 - 0.5*(med_val-min_val)/reach
        cmap_upper = 0.5 + 0.5*(max_val-med_val)/reach
        cmap_nsteps = int((cmap_upper - cmap_lower)*256)

        if(cmap_nsteps < 128):
            warnings.warn('Autoscaled colormap is significantly truncated - specify explicit clim')
            
        # truncate colormap to new range
        cmap_full = base_map
        cmp = ListedColormap(cmap_full(np.linspace(cmap_lower, cmap_upper, cmap_nsteps)))
    
    # use explicitly definied colormap range
    else:
        cmp = base_map
            
    # plot data
    fig,ax = plt.subplots()
    if equal:
        ax.set_aspect('equal')

    #fig = plt.figure()
    #ax = fig.add_axes([0.1,0.1,0.8,0.8])
    plt.scatter(data_in[h_ax]*sampling,
                data_in[v_ax]*sampling,
                marker='.',
                c=data_in[color_ax],
                cmap=cmp)
    cbar = plt.colorbar(format='%0.1f')
    
    # use explicitly definied colormap range
    if not clim == []:
        plt.clim(clim)
        
    # update labels
    if c_label == '':
        c_label = color_ax + ' (μm)'
    cbar.set_label(c_label)
    
    if h_label == '':
        h_label = h_ax + ' (μm)'
    plt.xlabel(h_label)

    if v_label == '':
        v_label = v_ax + ' (μm)'
    plt.ylabel(v_label)
        
    plt.title(plot_title)
    
    # set axis bounds
    if hlim != []:
        plt.xlim(hlim)
    if vlim != []:
        plt.ylim(vlim)
    if hticks != []:
        ax.xaxis.set_ticks(hticks)
    if vticks != []:
        ax.yaxis.set_ticks(vticks)

    # save plot
    if save_plot:
        print('Saving plot')
        plt.savefig(os.path.join(data_base_path,save_name))
    else:
        print('Plot NOT saved!')
    
    return fig, ax

def plot_avg_fit(h_ax,
                 v_ax,
                 save_name,
                 data_base_path,
                 data_in,
                 save_plot=True,
                 grid_dim=300):
    
    """
    2D plot of average Gaussian fit (i.e. a 2D slide of the 3D Gaussian function).

    Parameters
    ----------
    h_ax : int, horizontal plot axis (0: x, 1: y, 2: z)
    
    v_ax : int, vertical plot axis (0: x, 1: y, 2: z) - must be different from h_ax
    
    save_name : string, file name for saved plot

    data_in : input pandas dataframe of bead positions and 
       measurements (generated below)

    save_plot : (optional) bool, determines whether plot is saved.
    
    grid_dim : (optional) int, grid size
    """
    
    if h_ax not in [0,1,2] or v_ax not in [0,1,2]:
        raise Exception('h_ax and v_ax must be 0, 1, or 2')
    if h_ax==v_ax:
        raise Exception('h_ax and v_ax must be different')
        
    flat_dim=[0,1,2]
    flat_dim.remove(h_ax)
    flat_dim.remove(v_ax)
    flat_dim=flat_dim[0]
    
    offset = 0
    amplitude = 1
    xo = 0
    yo = 0
    zo = 0

    mean_FWHM = [0, 0, 0]
    mean_FWHM[0] = data_in.loc[:, 'FWHM_x'].mean()
    mean_FWHM[1] = data_in.loc[:, 'FWHM_y'].mean()
    mean_FWHM[2] = data_in.loc[:, 'FWHM_z'].mean()
    mean_sigma = [x / 2.355 for x in mean_FWHM]
    mean_sigma

    sigma_x = mean_sigma[0]
    sigma_y = mean_sigma[1]
    sigma_z = mean_sigma[2]
    
    grid_max = 1.5*max([mean_FWHM[h_ax], mean_FWHM[v_ax]])
    h_domain = np.linspace(-grid_max, grid_max, grid_dim)
    v_domain = np.linspace(-grid_max, grid_max, grid_dim)
    h_cord, v_cord = np.meshgrid(h_domain, v_domain)
    flat_cord = np.zeros((grid_dim, grid_dim))

    if h_ax==0:
        x_cord = h_cord
    elif h_ax==1:
        y_cord = h_cord
    elif h_ax==2:
        z_cord = h_cord

    if v_ax==0:
        x_cord = v_cord
    elif v_ax==1:
        y_cord = v_cord
    elif v_ax==2:
        z_cord = v_cord

    if flat_dim==0:
        x_cord = flat_cord
    elif flat_dim==1:
        y_cord = flat_cord
    elif flat_dim==2:
        z_cord = flat_cord
        
    g = offset + amplitude*np.exp(-(((x_cord-xo)**2)/(2*sigma_x**2) + 
                            ((y_cord-yo)**2)/(2*sigma_y**2) + 
                            ((z_cord-zo)**2)/(2*sigma_z**2)))
    
    fig,ax = plt.subplots()
    ax.imshow(g, cmap='gray');
    ax.axis('off');
    
    # save plot
    if save_plot:
        print('Saving plot')
        save_path=os.path.join(data_base_path,save_name)
        print(save_path)
        imsave(save_path,(g*255).astype('uint8'))
        print('DONE')
    else:
        print('Plot NOT saved!')
    
    return fig, ax