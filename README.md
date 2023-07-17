# Bead_PSF_computation

This software can be used to compute the PSF of an optical sectioning microscope (e.g. light-sheet) from a 3D image volume of diffraction limited beads. It locates well-separated beads and then fits a 2D Gaussian function to the maximum projection of the bead. The location and PSF dimenisons of each bead are saved as a .csv file. Various visualization functions allow the PSF to be characterized with respect to imaging depth, field position, etc.

Note that computing the PSFs from the 3D image volume can take a long time (several hours for a few thousand beads). Once you compute the PSFs once, you can load in the results from the .csv again later to more quickly generate new visualizations. 

This software is based on a package by [Nicholas Sofroniew](https://github.com/sofroniewn/psf) (Chan Zuckerberg Initiative), and edited by [Adam Glaser](https://github.com/adamkglaser) and [Kevin Bishop](https://github.com/kevinwbishop) (University of Washington).


### System requirements
This package is written in Python and run as a Jupyter Notebook. Users should install Anaconda (tested on Conda 4.3.30, Python 3.9.16, Windows 10).

In addition, the following packages are required, many of which come with Python/Anaconda. This code has been tested with the version numbers indicated, though other versions may also work.

Numpy (1.23.5)
Scikit-image (0.19.3)
Matplotlib (3.7.1)
Pandas (1.5.3)
Cmcrameri (1.4) - optional, for use with Scientific Colormaps from Fabio Crameri: http://doi.org/10.5281/zenodo.1243862

No specific computing hardware is required. However PSF calculation is computationally intensiv, and will in general run faster on a more powerful machine.

### Installation
After installing the required packages above, run the following command in Anaconda to clone this repository:
```bash
git clone https://github.com/kevinwbishop/Bead_PSF_computation
```

Typical installation time (after setting up Anaconda environment): less than 5 minutes.

### Demo
A sample input dataset (Tiff stack of a bead phantom imaged on a light-sheet microscope: fused_tp_0_ch_0.tif) and sample PSF results (PSF_data.csv) are provided in the demo folder. To run the demo, open the PSF_notebook Jupyter notebook and run, following the instructions. The expected outputs are provided in the notebook. The provided dataset is very small (14 beads) for demonstration purposes and should run in ~1 minute on a standard desktop PC. Computing PSFs for a typical dataset of several thousand beads may take a few hours.

### Instructions for use
To run, open PSF_notebook.ipynb as a Jupyter notebook and follow the instructions.
