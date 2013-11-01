function rms = calc_rms(x,y,h)

if nargin==2; h=0.2; end

max_y = max(y);
ind_y = find(y > h*max_y);

mean_x = sum(y(ind_y).*x(ind_y))/sum(y(ind_y));

rms = sqrt(sum(y(ind_y).*(x(ind_y)-mean_x).^2)/sum(y(ind_y)));

