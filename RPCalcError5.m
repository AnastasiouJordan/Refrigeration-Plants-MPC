function [E5, h_RP, v] = RPCalcError5(pUAvec, u, p, t, N)

% Simulate the outputs of the system based on parameter estimates.
% First convert p_vector of inital parameter guesses to a structure:
p = V2S(pUAvec, p.regressedparameterfields, p);


h_RP0 = p.C_p * (u.T_RP(0) - p.T_0) + p.h_0;

[~, h_RP] = ode45(@(t, h_RP) FridgePlantsODEs(p, h_RP, u, t), t(1:N), h_RP0);
v = RPIntermediates(h_RP', u, p, t(1:N));

% Finally, we calculate the error based on this output of the simulation
% using the estimated parameters

E5 = [v.T_RP(:,5) - u.T_RP.Values(1:N,5)];
