% load wake
global wake;
wake = load('slac.dat');

% load params
global PARAM;
param_5mm_R56;
init_beam = 1;

% scan parameter (see file Parameters/param_5mm_R56.m 
% for a full list of available parameters)
param_name = {'NRTL AMPL','LONE DECK'};
low_limit  = [0.0390,-22.0];
high_limit = [0.0415,-19.5];
n_steps    = [10, 10];
nPar       = numel(n_steps);

% create scan values
nrtl_vals = linspace(low_limit(1),high_limit(1),n_steps(1));
phas_vals = linspace(low_limit(2),high_limit(2),n_steps(2));

% allocate space for results
SIGMA_Z = zeros(n_steps);
RMS_Z   = zeros(n_steps);
I_PEAK  = zeros(n_steps);
FWHM_D  = zeros(n_steps);
RMS_D   = zeros(n_steps);

for i = 1:numel(nrtl_vals)
    for j = 1:numel(phas_vals)
    
        % Update parameters by calling SetPars and then executing beamline file
        SetPars([nrtl_vals(i) phas_vals(j)], param_name, nPar);
        FACET_NOTCH_bl;
    
        % Run LiTrack
        LiT_struct = LiT(beamline,init_beam,init_param,PARAM.SIMU.BIN,0);
    
        % Generate phase spaces for plots
        [profile,spectrum,phase_space,beam_zd] = gen_ps(LiT_struct,PARAM.SIMU.BIN);
        
        % Calculate beam moments (outside of LiTrack)
        [SIGMA_Z(i,j),RMS_Z(i,j),I_PEAK(i,j),FWHM_D(i,j),RMS_D(i,j)] = calculate_moments(profile,spectrum);
        
        % display progress
        if i == 1 && j == 1; del = 0; else del = 1; end;
        disp_prog(j+numel(nrtl_vals)*(i-1),numel(nrtl_vals)*numel(phas_vals),del);
    end
end
fprintf('\n');

% Contour plot of peak current
[c,h] = contourf(1000*nrtl_vals,phas_vals,I_PEAK');
clabel(c);
xlabel('NRTL Compressor Amplitude [MV]','fontsize',14);
ylabel('LI02-LI10 Phase [deg]','fontsize',14);
title('Peak Current [kA]','fontsize',14);