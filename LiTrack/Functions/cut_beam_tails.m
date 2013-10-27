function [beam_out, Nb] = cut_beam_tails(beam_in,Nb,params)

[Z,ind] = sort(beam_in(:,1));
E = beam_in(ind,2);

Z0 = mean(Z);
sigz = std(Z);

z_sig = params(2);	% sigma_z cut

ind = abs(Z-Z0)/sigz < z_sig;

Z = Z(ind);
E = E(ind);
    
beam_out = [Z E];
Nb = numel(Z);