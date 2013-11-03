% load wake
global wake;
wake = load('slac.dat');

% load params
global PARAM;
param_5mm_R56;
init_beam = 1;

% scan parameter (see file Parameters/param_5mm_R56.m 
% for a full list of available parameters)
param_name = {'INIT NPART','NRTL AMPL','LONE DECK'};
low_limit  = [1.95E10, 0.0390,-22.0];
high_limit = [2.05E10, 0.0410,-20.0];
n_steps    = [5, 5, 5];
nPar       = numel(n_steps);
nSims      = prod(n_steps);

% create scan values 
%(this is the most brilliant function ever written, you're welcome universe)
[scan_inds, scan_vals] = ScanSpace(low_limit,high_limit,n_steps);

% allocate space for results
SIGMA_Z = zeros(n_steps);
RMS_Z   = zeros(n_steps);
I_PEAK  = zeros(n_steps);
FWHM_D  = zeros(n_steps);
RMS_D   = zeros(n_steps);

% estimate amount of time and data needed
nVec = 4;
time_per_sim = 0.25;
guess_time_and_data(nVec,PARAM.SIMU.BIN,time_per_sim,nSims);

% allocate space for projections
spectra = zeros([PARAM.SIMU.BIN n_steps]);
spec_ax = zeros([PARAM.SIMU.BIN n_steps]);
profils = zeros([PARAM.SIMU.BIN n_steps]);
prof_ax = zeros([PARAM.SIMU.BIN n_steps]);

display('Starting scan');
for i = 1:nSims
    
    % Update parameters by calling SetPars and then executing beamline file
    SetPars(scan_vals(:,i), param_name, nPar);
    FACET_NOTCH_bl;
    
    % Run LiTrack
    LiT_struct = LiT(beamline,init_beam,init_param,PARAM.SIMU.BIN,0);
    
    % Generate phase spaces for plots
    [profile,spectrum,phase_space,beam_zd] = gen_ps(LiT_struct,PARAM.SIMU.BIN);
    spectra(:,i) = spectrum(:,2);
    spec_ax(:,i) = spectrum(:,1);
    profils(:,i) = profile(:,2);
    prof_ax(:,i) = profile(:,1);
    
    % Calculate beam moments (outside of LiTrack)
    [SIGMA_Z(i),RMS_Z(i),I_PEAK(i),FWHM_D(i),RMS_D(i)] = calculate_moments(profile,spectrum);
    
    % display progress
    if i == 1; del = 0; else del = 1; end;
    disp_prog(i,nSims,del);
    
end
fprintf('\n');
display('Finished!');