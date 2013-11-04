function [profile, spectrum, yag_spec, residual] = compare_spectra(beam,Nbin,yag_axis,yag_spec,gaussFilter)
% [z_ax,lit_z_prof,e_ax,lit_e_prof,yag_spec] = comp_data_litrack(beam,Nbin,yag_axis,yag_spec)
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
%   smear:      (optional) a struct of smear parameters
%
% OUTPUTS:
%
%   z_ax:       data and LiTrack Z axis [microns]
%   lit_z_prof: LiTrack current profile binned on z_ax
%   e_ax:       data and LiTrack delta axis [%]
%   lit_e_spec: LiTrack Delta spectrum binned on e_ax
%   yag_spec:   Normalized YAG spectrum
%   residual:   Difference between LiTrack and data

% Convert LiTrack beam from [m,GeV] to [um,%]
lit_z = 1e6*beam(:,1);
lit_d = 100*(beam(:,2)-mean(beam(:,2)))/mean(beam(:,2));

% Histogram on shared axis, zero end bins, scale
[nz,zb] = hist(lit_z,Nbin);
[ne,eb] = hist(lit_d,yag_axis);
ne(1) = 0;
ne(end) = 0;

if nargin > 4
    
    ne = conv(ne, gaussFilter, 'same');
end

lit_emax = max(ne);
dat_emax = max(yag_spec);
yag_spec = (lit_emax/dat_emax)*yag_spec;

% Copy results
profile = [zb' nz'];
spectrum = [eb' ne'];

residual = sum((yag_spec - ne).^2);