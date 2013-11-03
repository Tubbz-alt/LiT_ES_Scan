function [scan_inds, scan_vals] = ScanSpace(low_limit,high_limit,n_steps)

n_vals = prod(n_steps);
idx = 1:n_vals;


npars = numel(n_steps);
INDS      = cell(1,npars);
[INDS{:}] = ind2sub(n_steps,idx);

scan_inds = [];
scan_vals = [];

for i = 1:npars
    
    par_vals  = linspace(low_limit(i),high_limit(i),n_steps(i));
    scan_inds = [scan_inds; INDS{i}];
    cur_vals  = par_vals(INDS{i});
    scan_vals = [scan_vals; cur_vals];
    
end