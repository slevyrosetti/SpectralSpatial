# SpectralSpatial
Files for producing spectral and spectral-spatial prewinding RF pulses for MRI <br />
Based on S. N. Williams, J-F. Nielsen, J.A. Fessler, and D.C. Noll, "Design of spectral-spatial phase prewinding pulses and their use in small-tip fast recovery steady-state imaging", Magnetic Resonance in Medicine, 2017 http://onlinelibrary.wiley.com/doi/10.1002/mrm.26794/abstract <br />
<br />
Files found at: [https://github.com/sydneynw/SpectralSpatial/](https://github.com/sydneynw/SpectralSpatial/)<br />
<br />
main files:<br />
-example.m              runs example design/plots on human data<br />
-PrewindingPulse.mat    human data for pulse design and simulation<br />
-prewinding_pulse.m     sets up spectral and spectral-spatial RF pulse design<br />
-spectralRF.m           spectral prewinding design<br />
-spectralspatialRF.m    spectral prewinding design<br />
... additional file details listed in example.m<br />
<br />
needs the following external packages:<br />
-Michigan IRT           Jeff Fessler's Image Recon Toolbox, http://web.eecs.umich.edu/~fessler/code/index.html<br />
-Var. Dens. Spiral      Brian Hargreaves' VDS Functions, http://www-mrsrl.stanford.edu/~brian/vdspiral/<br />
-PTx Bloch Simulator    Hao Sun's Mex Bloch Simulator, http://www-personal.umich.edu/~sunhao/<br />
<br />
disclaimer: first Github repo, hoping to learn more soon. Please contact me at sydneynw_at_umich_dot_edu for questions or comments, thanks!<br />
