function f = fCA (v,cA)

c_inf  = 0.28; % CA-units

f = v/(1+(cA/c_inf)^2);
