%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stratification parameters
%
z_p = -100; % mixed layer depth
dRhoML = 2; % density jump across the transition layer
[N2,zDomain] = StratificationProfileWithMixedLayer(z_p,dRhoML);
latitude = 31;

z = linspace(min(zDomain),max(zDomain),1000)';
im = InternalModes(N2,zDomain,z,latitude,'rho0',1025,'nEVP', 513, 'N2',1);
im.upperBoundary = UpperBoundary.freeSurface;

% Let's use the barotopic mode to create a combined mode with no bottom
% velocity.
[F,G] = im.ModesAtFrequency(0.0);
Fmode= -F(1,2)*F(:,1) + F(:,2);
Gmode= -F(1,2)*G(:,1) + G(:,2);

% Normalize it to be 1 at the surface
Fmode = Fmode/Fmode(end);
Gmode = Gmode/Fmode(end);

% Stream function for an eddy
EddyRadius = 50e3;
EddyAmplitude = 0.1; % cm
psi_eddy = @(r,zd) (EddyAmplitude*9.81/im.f0)*exp( -r.^2/EddyRadius/EddyRadius/2) .* interp1(Fmode,z,zd);
u_eddy = @(x,y,zd) (y/EddyRadius/EddyRadius).*(EddyAmplitude*9.81/im.f0).*exp( -(x.^2 + y.^2)/EddyRadius/EddyRadius/2) .* interp1(im.z,Fmode,zd);
v_eddy = @(x,y,zd) - (x/EddyRadius/EddyRadius).*(EddyAmplitude*9.81/im.f0).*exp( -(x.^2 + y.^2)/EddyRadius/EddyRadius/2) .* interp1(im.z,Fmode,zd);


Lx = 350e3;
Ly = Lx;
D = 400;
N = 128;
x = linspace(-Lx/2,Lx/2,N).';
y = linspace(-Ly/2,Ly/2,N).';
z = linspace(-D,0,N).';
[X,Y,Z] = ndgrid(x,y,z);

U = u_eddy(X,Y,Z);
V = v_eddy(X,Y,Z);

figure, pcolor(squeeze(X(:,N/2,:)),squeeze(Z(:,N/2,:)),squeeze(U(:,N/2,:))), shading interp, xlabel('x'), title('u')
figure, pcolor(squeeze(X(:,N/2,:)),squeeze(Z(:,N/2,:)),squeeze(V(:,N/2,:))), shading interp, xlabel('x'), title('v')
figure, quiver(X(:,:,1),Y(:,:,1),U(:,:,end),V(:,:,end))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now fetch a wind record
startDate = datenum('2007-06-08 04:00:00');
endDate = datenum('2008-05-29 00:00:00');
[t_wind, u_wind, v_wind] = PapaWindRecord(startDate, endDate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use a damped slab model
slab_depth = abs(z_p);
slab_damp = 6; % damping time, in days
[t, u, v] = OBLModel_DampedSlab( t_wind, u_wind, v_wind, slab_depth, latitude, slab_damp );
figure, plot(t,[u,v])

% So now you can simply add this u and v time series, over the mixed layer
% region of the eddy.


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Try a completely different model---this one has depth dependence.
% K0 = 100e-5;
% [t_model1, u_model1, v_model1] = OBLModel_InfiniteLayerConstantK( t_wind, u_wind, v_wind, z, latitude, K0 );