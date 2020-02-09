%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stratification parameters
%
z_p = -100; % mixed layer depth
dRhoML = 3; % density jump across the transition layer
[N2,zDomain] = StratificationProfileWithMixedLayer(z_p,dRhoML);

z = linspace(min(zDomain),max(zDomain),1000)';
im = InternalModes(N2,zDomain,z,31,'rho0',1025,'nEVP', 513, 'N2',1);
im.upperBoundary = UpperBoundary.freeSurface;

[F,G] = im.ModesAtFrequency(0.0);
Fmode= -F(1,2)*F(:,1) + F(:,2);
Gmode= -F(1,2)*G(:,1) + G(:,2);

Fmode2= -F(1,3)*F(:,1) + F(:,3);
Gmode2= -F(1,3)*G(:,1) + G(:,3);

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

packfig(1,4)

print('-depsc',sprintf('../ModesWithDeltaRho%d.eps',round(dRhoML)))
