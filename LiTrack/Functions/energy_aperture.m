function [beam_out, Nb] = energy_aperture(beam_in,Nb,params)

[E,ind] = sort(beam_in(:,2));
Z = beam_in(ind,1);
%E0 = mean(E);

d_lo = params(2);	% minimum dE/E to allow through [ ]
d_hi = params(3);	% maximum dE/E "   "     "      [ ]
E0   = params(4);

ind_lo = (E-E0)/E0 < d_hi;
ind_hi = (E-E0)/E0 > d_lo;

E = E(ind_lo & ind_hi);
Z = Z(ind_lo & ind_hi);
    
beam_out = [Z E];
Nb = numel(E);