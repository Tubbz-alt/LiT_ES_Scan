% load wake
global wake;
wake = load('slac.dat');

% load params
global PARAM;
param_5mm_R56;

% load beamline
FACET_NOTCH_bl;

% set beam init params
init_beam = 1;

% run LiT
tic;
out = LiT(beamline,init_beam,init_param,PARAM.SIMU.BIN,0);
display(['LiTrack took ' num2str(toc) ' seconds to run.']);

% plot phase space
plot_ps(out,PARAM.SIMU.BIN,1,1,0,0);