function gain = set_gain(omega,estimated_cost,scale_factor)
%set_gain: This function determines the gain based on the estimated cost

w = mean(omega);
gain = w/(scale_factor*estimated_cost);