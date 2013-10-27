function [dE,zc] = calc_wake(z,L,Qp,Nbin)
%===============================================================================
%  [dE,zc] = calc_wake(z,L,Ne,Nbin);
%
%  Function to return the wakefield induced energy profile vs. z for
%  a set of given axial coordinates "z".
%
%  INPUTS:	z:		The internal axial coordinates, within the bunch, of
%					each electron with respect to any fixed point [m]
%			L:		The length of the linac [m]
%			Qp:		Particle macro charge
%			Nbin:   The number of bins to use
%			
%
%  OUTPUTS:	dE:		The energy loss per binned bunch slice [MeV]
%			zc:		The sample points along z where dE is calculated [m]
%					[e.g. plot(zc,dE)]
%			
%===============================================================================
SI_consts;
global wake;

% Histogram particles into Nbins with zeros at ends for interpolation
[N,zc] = hist(z,Nbin-2);
% Bin spacing (m)
dzc = zc(2)-zc(1);
% Zero padded axis
zc = [zc(1)-dzc zc zc(Nbin-2)+dzc];
zc = zc';
% Zero padded histogram
N = [0 N 0];
N = N';

% delta vector
dE = zeros(Nbin,1);

% scale variable (mC/V ?)
scl = -L*Qp*SI_e*1E-6;

% Bin separation vector
dzi = dzc*((1:Nbin)' - 1);

% Interpolated wake vector
Wf = interp1(wake(:,1),wake(:,2),dzi);

% Self wake bin
Wf(1) = Wf(1)/2;

% Sum delta due to wake from leading bins
for j =1:Nbin
    dE(j) = sum(scl*N(j:-1:1).*Wf(1:j));
end
