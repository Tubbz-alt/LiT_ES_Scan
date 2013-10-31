function gain = set_gain(omega,estimated_cost)
%set_gain: This function determines the gain based on the estimated cost

w = mean(omega);
ratio = 20;
gain = w/(ratio*estimated_cost);