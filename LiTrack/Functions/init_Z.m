function x = init_Z(N,xsig,asym,cut)

%	x = asym_gaussian(N,xsig,xbar,asym,cut,tail,halo,halo_pop);
%
%	Function to generate an asymmetric gaussian truncated at
%	+/- cut*xsig and with added tails described by "tail".
%
%    INPTS:
%			N:			Total number of particles, even after cuts.
%			xsig:		rms of distribution after asymmetry and cuts
%			asym:		asymmetry parameter (-1 < asym < 1) - [asym
%						of +-0.2 is weak and +-0.8 is strong (0 is none)
%						[asym>0 gives slope>0; i.e. larger sigma @ x<0]
%			cut:		Number of "xsig" to truncate dist. at
%						(0.5 <= cut < inf) - still get length(x)=N.
%
%
%    OUTPUTS:
%			x:			Array of random numbers with length N

%===================================================================

f1 = 1 + asym;
f2 = 1 - asym;

N1 = round(N*f1/4/erfn(cut));		% boost N1 to accomodate cuts
N2 = round(N*f2/4/erfn(cut));		% boost N2 to accomodate cuts

x1      = f1*xsig*randn(2*N1,1);	% generate a gaussian distribution + offset
cut_ind = x1<0 & x1>(-cut*xsig*f1);	% eliminate positive values and cuts
x       = x1(cut_ind);

x2      = f2*xsig*randn(2*N2,1);    % generate a gaussian distribution + offset
cut_ind = x2>0 & x2<(cut*xsig*f2);	% eliminate positive values
x       = [x; x2(cut_ind)];

% add or subtract particles to match number asked for
NN = length(x);
dN = N - NN;
if dN > 0
  x1 = min([f1 f2])*xsig*(rand(dN,1)-0.5);	% generate a gaussian distribution + offset
  x = [x; x1];
elseif dN < 0
  adN = abs(dN);
  i = round((NN/adN)*(1:adN));
  x(i) = [];
end

x = x - mean(x);
x = x*xsig/std(x);