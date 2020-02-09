function [N2,zDomain] = StratificationProfileWithMixedLayer(z_p,dRhoML)
%StratificationProfileWithMixedLayer
%   z_p is the depth of the mixed layer
%   dRhoML is the density jump across the transition layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stratification parameters
%
N0 = 5e-4; % bouyancy frequency at the surface
delta_p = 15; % width of the transition layer

Lgm = 1300; % e-fold scale of stratification below the mixed layer
Nb = (6.5e-3)*exp(-4000/Lgm); % buoyancy freqency at the bottom. This is tuned so that the equivalent depth is approx 80 cm.
D = 4000; % depth of the ocean

N2_ml = dRhoML*(9.81)/(2*1025*delta_p);

% Stratification profile---deep exponential, then rapid transition to mixed
% layer at the surface.
N2 = @(z) N0*N0*exp( ((z+D)/Lgm+log(Nb)-log(N0)) .* (1-tanh((z-z_p)/delta_p)) ) + N2_ml*sech( (z-z_p)/delta_p ).^2;

zDomain = [-D, 0];
end

