# psf

Originally written by Nicholas Sofroniew (CZI: https://github.com/sofroniewn/psf), edited by Adam Glaser and Kevin Bishop (UW).

This is currently an internal version for Liu Lab at UW and should NOT be shared publically. To run, open PSF_notebook.ipynb as a Jupyter notebook and follow the instructions

----------------
> Compute the point spread function from a 3d image of beads

This package computes the point spread function (psf) from a 3d image of beads. It first finds the centers of each well separated bead and then fits 2d Gaussians to max projections of the bead. It returns the resulting psfs as a table.


### installation

```bash
git clone https://github.com/kevinwbishop/Bead_PSF_computation
```
