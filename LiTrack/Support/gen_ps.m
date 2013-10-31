function [z_proj,e_proj,phase_space,beam_zd] = gen_ps(LiT_struct,Nbin)
% [z_proj,e_proj,beam_zd] = gen_ps(beam,Nbin)
%
% This function generates phase spaces and projections from LiTrack output.
% 
% INPUTS: 
%
%   beam:       z = beam(:,1) [m]
%               e = beam(:,2) [GeV]
%   Nbin:       # of bins for histogram
%
% OUTPUTS:
%
%   z_proj:     Current profile projection of phase space
%               z_axis = z_proj(:,1) [microns]
%               z_prof = z_proj(:,2) [kA] 
%   e_proj:     Energy spectrum projection of phase space
%               e_axis = e_proj(:,1) [%]
%               e_spec = e_proj(:,2) [N/bin] 
%   beam_zd:    z = beam_zd(:,1) [microns]
%               d = beam_zd(:,2) [%]

% Convert LiTrack beam from [m,GeV] to [um,%]
lit_z = 1e6*LiT_struct.BEAM(:,1);
lit_d = 100*(LiT_struct.BEAM(:,2)-mean(LiT_struct.BEAM(:,2)))/mean(LiT_struct.BEAM(:,2));

% Histogram beam to determine, max and position
[nz,zb] = hist(lit_z,Nbin);
[ne,eb] = hist(lit_d,Nbin);
phase_space = hist2(lit_z,lit_d,zb,eb);


% load SI_constants
SI_consts;

% Convert z-distribution to current
dz = 1e-6*(zb(2)-zb(1)); % bin spacing [m]
qb = LiT_struct.QP*nz*SI_e; % charge/bin [C]
I  = qb*SI_c/(1000*dz); % current [kA]
z_proj = [zb' I'];

% Convert e-distribution to N electron/bin
qb = LiT_struct.QP*ne;
e_proj = [eb' qb'];

beam_zd = [lit_z lit_d];