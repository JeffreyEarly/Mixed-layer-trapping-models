%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stratification parameters
%
z_p = -100; % mixed layer depth
dRhoML = 0; % density jump across the transition layer
[N2,zDomain] = StratificationProfileWithMixedLayer(z_p,dRhoML);

z = linspace(min(zDomain),max(zDomain),1000)';
im = InternalModes(N2,zDomain,z,31,'rho0',1025,'nEVP', 513, 'N2',1);
im.upperBoundary = UpperBoundary.freeSurface;

[F,G] = im.ModesAtFrequency(0.0);
Fmode= -F(1,2)*F(:,1) + F(:,2);
Gmode= -F(1,2)*G(:,1) + G(:,2);


dRhoML = 3; % density jump across the transition layer
[N2,zDomain] = StratificationProfileWithMixedLayer(z_p,dRhoML);

z = linspace(min(zDomain),max(zDomain),1000)';
im = InternalModes(N2,zDomain,z,31,'rho0',1025,'nEVP', 513, 'N2',1);
im.upperBoundary = UpperBoundary.freeSurface;

[F,G] = im.ModesAtFrequency(0.0);

Fmode2= -F(1,2)*F(:,1) + F(:,2);
Gmode2= -F(1,2)*G(:,1) + G(:,2);

if 1==0
    % match interior density anomaly
    rescale = max(abs(Gmode))/max(abs(Gmode2));
    Gmode2 = Gmode2*rescale;
    Fmode2 = Fmode2*rescale;
else
    % match SSH
    Fmode = Fmode/Fmode(end);
    Gmode = Gmode/Fmode(end);
    
    Fmode2 = Fmode2/Fmode2(end);
    Gmode2 = Gmode2/Fmode2(end);
end

figure
subplot(1,2,1)
plot([Fmode, Fmode2],z, 'LineWidth', 2)
title('(u,v)-mode')
ylabel('depth (meters)')

subplot(1,2,2)
plot([Gmode, Gmode2],z, 'LineWidth', 2)
title('density-mode');
yticks([]);
legend('\Delta \rho=0', '\Delta \rho=3','Location','east')

packfig(1,2)

print('-depsc','../ModesWithDeltaRhoComparisonSSHMatched.eps')
