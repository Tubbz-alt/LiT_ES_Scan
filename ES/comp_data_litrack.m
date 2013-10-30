function [beam_zd, z_ax, lit_z_prof, e_ax, lit_e_spec, yag_spec, tcav_prof] = comp_data_litrack(beam,Nbin,yag_axis,yag_spec,tcav_axis,tcav_prof)
% [beam_zd,z_ax,lit_z_prof,e_ax,lit_e_prof,yag_spec,tcav_prof] = comp_data_litrack(beam,Nbin,yag_axis,yag_spec,tcav_axis,tcav_spec)
%
% This function compares data with LiTrack simulations.
% 
% INPUTS: 
%
%   beam:       z = beam(:,1) [m]
%               e = beam(:,2) [GeV]
%   Nbin:       # of bins for histogram
%   yag_axis:   data Delta axis [%], must be same length as Nbin
%   yag_spec:   data Delta spectrum, must be same length as Nbin
%   tcav_axis:  (optional) data Z axis [microns], must be same length as Nbin
%   tcav_prof:  (optional) data Z current profile, must be same length as Nbin
%
% OUTPUTS:
%
%   beam_zd:    z = beam_zd(:,1) [microns]
%               d = beam_zd(:,2) [%]
%   z_ax:       data and LiTrack Z axis [microns]
%   lit_z_prof: LiTrack current profile binned on z_ax
%   e_ax:       data and LiTrack delta axis [%]
%   lit_e_spec: LiTrack Delta spectrum binned on e_ax

% Convert LiTrack beam from [m,GeV] to [um,%]
lit_z = 1e6*beam(:,1);
lit_d = 100*(beam(:,2)-mean(beam(:,2)))/mean(beam(:,2));

% Histogram beam to determine, max and position
[nz,zb] = hist(lit_z,Nbin);
[ne,eb] = hist(lit_d,Nbin);

if nargin == 6 % compare syag and tcav with LiTrack
    
    % Shift LiTrack max to zero
    [~,lit_ind] = max(nz);
    lit_z = lit_z - zb(lit_ind);

    % Shift TCAV max to zero
    [~,dat_ind] = max(tcav_prof);
    tcav_axis = tcav_axis - tcav_axis(dat_ind);
    
    % Histogram on shared axis, zero end bins, scale
    [nz,zb] = hist(lit_z,tcav_axis);
    nz(1) = 0;
    nz(end) = 0;
    lit_zmax = max(nz);
    dat_zmax = max(tcav_prof);
    tcav_prof = (lit_zmax/dat_zmax)*tcav_prof;

    % Histogram on shared axis, zero end bins, scale
    [ne,eb] = hist(lit_d,yag_axis);
    ne(1) = 0;
    ne(end) = 0;
    lit_emax = max(ne);
    dat_emax = max(yag_spec);
    yag_spec = (lit_emax/dat_emax)*yag_spec;

elseif nargin == 4 % compare syag with LiTrack

    % Histogram on shared axis, zero end bins, scale
    [nz,zb] = hist(lit_z,Nbin);
    [ne,eb] = hist(lit_d,yag_axis);
    ne(1) = 0;
    ne(end) = 0;
    lit_emax = max(ne);
    dat_emax = max(yag_spec);
    yag_spec = (lit_emax/dat_emax)*yag_spec;

end

% Copy results
beam_zd = [lit_z lit_d];
z_ax = zb;
e_ax = eb;
lit_z_prof = nz;
lit_e_spec = ne;