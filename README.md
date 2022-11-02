# Bead_PSF_computation

This software can be used to compute the PSF of an optical sectioning microscope (e.g. light-sheet) from a 3D image volume of diffraction limited beads. It locates well-separated beads and then fits a 2D Gaussian function to the maximum projection of the bead. The location and PSF dimenisons of each bead are saved as a .csv file. Various visualization functions allow the PSF to be characterized with respect to imaging depth, field position, etc.

Note that computing the PSFs from the 3D image volume can take a long time (several hours for a few thousand beads). Once you compute the PSFs once, you can load in the results from the .csv again later for to more quickly generate new visualizations. 

To run, open PSF_notebook.ipynb as a Jupyter notebook and follow the instructions.

This package was originally written by [Nicholas Sofroniew](https://github.com/sofroniewn/psf) (Chan Zuckerberg Initiative), and edited by [Adam Glaser](https://github.com/adamkglaser) and [Kevin Bishop](https://github.com/kevinwbishop) (University of Washington).


### installation

```bash
git clone https://github.com/kevinwbishop/Bead_PSF_computation
```
