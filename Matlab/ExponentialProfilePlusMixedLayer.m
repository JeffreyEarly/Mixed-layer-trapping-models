%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stratification parameters
%
N0 = 5e-4; % bouyancy frequency at the surface
z_p = -100; % mixed layer depth
delta_p = 15; % width of the transition layer

Lgm = 1300; % e-fold scale of stratification below the mixed layer
Nb = (6.5e-3)*exp(-4000/Lgm); % buoyancy freqency at the bottom. This is tuned so that the equivalent depth is approx 80 cm.
D = 4000; % depth of the ocean

dRhoML = 0; % change in density across the mixed layer.. certainly 2 kg/m^3 is reasonable in some places.
N2_ml = dRhoML*(9.81)/(2*1025*delta_p);

% Stratification profile---deep exponential, then rapid transition to mixed
% layer at the surface.
N2 = @(z) N0*N0*exp( ((z+D)/Lgm+log(Nb)-log(N0)) .* (1-tanh((z-z_p)/delta_p)) ) + N2_ml*sech( (z-z_p)/delta_p ).^2;


z = linspace(-D,0,1000)';
im = InternalModes(N2,[-D 0],z,31,'rho0',1025,'nEVP', 513, 'N2',1);
im.upperBoundary = UpperBoundary.freeSurface;
im.ShowLowestModesAtFrequency(0);

[F,G] = im.ModesAtFrequency(0.0);
Fmode= -F(1,2)*F(:,1) + F(:,2);
Gmode= -F(1,2)*G(:,1) + G(:,2);

Fmode2= -F(1,3)*F(:,1) + F(:,3);
Gmode2= -F(1,3)*G(:,1) + G(:,3);

% Normalize it to be 1 at the surface
% Fmode = Fmode/Fmode(end);
% Gmode = Gmode/Fmode(end);
% 
% Fmode2 = Fmode2/Fmode2(end);
% Gmode2 = Gmode2/Fmode2(end);

figure
subplot(1,4,1)
plot(im.rho,z, 'LineWidth', 2)
title('density (kg/m^3)')
ylabel('depth (meters)')

subplot(1,4,2)
plot(im.N2,z, 'LineWidth', 2)
title('N^2 (1/s^2)')
yticks([]);

subplot(1,4,3)
plot([Fmode, Fmode2],z, 'LineWidth', 2)
title('(u,v)-mode')
yticks([]);

subplot(1,4,4)
plot([Gmode, Gmode2],z, 'LineWidth', 2)
title('density-mode');
yticks([]);

return
EddyRadius = 50e3;
EddyAmplitude = 0.1; % cm
z0 = -400;
eddy = @(r,zd) exp( -r.^2/EddyRadius/EddyRadius/2) .* interp1(Fmode,z,zd);



