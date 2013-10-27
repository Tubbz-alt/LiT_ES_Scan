function [beam_out, Nb] = notch_collimator(beam_in,Nb,params)

[E,ind] = sort(beam_in(:,2));
Z = beam_in(ind,1);
E0 = mean(E);

n_lo = params(2);	% minimum dE/E to cut
n_hi = params(3);	% maximum dE/E to cut

ind_lo = (E-E0)/E0 > n_lo;
ind_hi = (E-E0)/E0 < n_hi;

E = E(~(ind_lo & ind_hi));
Z = Z(~(ind_lo & ind_hi));
    
beam_out = [Z E];
Nb = numel(E);