function [beam_out, Nb] = ISR_energySpread(beam_in,Nb,params)

E = beam_in(:,2);

id_rms = params(2);
d = id_rms*randn(Nb,1);

E = E.*(1+d);

beam_out = [beam_in(:,1) E];

Nb = length(beam_out);