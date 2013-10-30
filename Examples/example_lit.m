addpath(genpath('../../LiT_ES_Scan/'));

% load wake
global wake;
wake = load('slac.dat');

% load params
global PARAM;
param_test1;

% tweak params here
% for a full list of parameters see 'Parameters/param_test1.m'
PARAM.INIT.NPART = 2.09E10; % # beam eleectrons
PARAM.NRTL.AMPL  = 0.0397;  % Compressor Phase
PARAM.LONE.PHAS  = -21.0;   % 2-10 phase
PARAM.LI20.NLO   = -0.006;  % notch low energy
PARAM.LI20.NHI   = 0.008;   % notch high energy
PARAM.LI20.EHI   = 0.035;   % S20 low energy
PARAM.LI20.ELO   = -0.025;  % S20 high energy

% load beamline
FACET_NOTCH_bl;

% set beam init params
init_beam = 1;

% run LiT
tic;
out = LiT(beamline,init_beam,init_param,PARAM.SIMU.BIN,0);
display(['LiTrack took ' num2str(toc) ' seconds to run.']);

% plot phase space
plot_ps(out.BEAM,PARAM.SIMU.BIN,1,0);