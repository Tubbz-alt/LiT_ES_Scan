function [beam_out, Nb] = bunch_compressor(beam_in,Nb,params)

z = beam_in(:,1);
E = beam_in(:,2);

R56   = params(2);       % R56 value [m]
T566  = params(3);       % T566 value [m] (=-3*R56/2 for non-quad system)
U5666 = params(4);       % U5666 value [m] (=2*R56 for non-quad system)
E56   = params(5);       % Nominal energy of compressor [GeV]

if params(6)
    
    T566  = -1.5*R56;            % T566 value [m]
    U5666 =  2.0*R56;            % U5666 value [m]
    
end

% relative energy error w.r.t. nominal compressor energy
dd = (E-E56)/E56;			

% compress or anti-compress bunch [m]
z = R56*dd + T566*dd.^2 + U5666*dd.^3 + z;		

beam_out = [z E];
Nb = length(beam_out);