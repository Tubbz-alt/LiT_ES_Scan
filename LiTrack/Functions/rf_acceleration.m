function [beam_out, Nb] = rf_acceleration(beam_in,Nb,Qp,Nbin,params)

SI_consts;

[Z,ind] = sort(beam_in(:,1));
E = beam_in(ind,2);
E0 = mean(E);

E_acc = params(2);
phi   = params(3)*pi/180;
rf    = params(4);
fb    = params(5);
Lacc  = params(6);



[dE_wake,zc_wake] = calc_wake(Z,Lacc,Qp,Nbin);  % convolve wake function with beam in Nbin histogram
int_wake = interp1(zc_wake,dE_wake,Z,'linear');	% inerpolate wake response onto particles
wakeloss = 1E-3*mean(int_wake);                 % mean wake loss [GeV]

if fb
    Egain = (E_acc - E0 - wakeloss)/cos(phi);
else
    Egain = E_acc;
end

Erf = Egain*cos(phi + 2*pi*Z/(SI_c/SI_sband));

E = E + Erf + 1E-3*int_wake;

beam_out = [Z E];

Nb = length(beam_out);