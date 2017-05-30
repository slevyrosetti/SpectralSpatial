# SpectralSpatial
Files for producing spectral and spectral-spatial prewinding RF pulses for MRI
Based on S. N. Williams, J-F. Nielsen, J.A. Fessler, and D.C. Noll, "Design of spectral-spatial phase prewinding pulses and their use in small-tip fast recovery steady-state imaging", Magnetic Resonance in Medicine, 2017 Early View

main files:
-example.m              runs example design/plots on human data
-PrewindingPulse.mat    human data for pulse design and simulation
-prewinding_pulse.m     sets up spectral and spectral-spatial RF pulse design
-spectralRF.m           spectral prewinding design
-spectralspatialRF.m    spectral prewinding design
... additional file details listed in example.m

needs the following external packages:
-Michigan IRT           Jeff Fessler's Image Recon Toolbox, http://web.eecs.umich.edu/~fessler/code/index.html
-Var. Dens. Spiral      Brian Hargreaves' VDS Functions, http://www-mrsrl.stanford.edu/~brian/vdspiral/
-PTx Bloch Simulator    Hao Sun's Mex Bloch Simulator, http://www-personal.umich.edu/~sunhao/

disclaimer: first Github repo, hoping to learn more soon. Please contact me at sydneynw_at_umich_dot_edu for questions or comments, thanks!
