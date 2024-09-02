

addpath(genpath('.\utils'))
addpath(genpath('.\SeisLab_10.0301'))

%% Possible improvements:
% Use of sgs and rockphysics to create a more realistic facies for the upper and lower zones 
% Use SGS to disturbe net2gross to make it more continuous/realistic
% Improve sgy inline and xline definicion. Currently it is the same of x y 

% Sampling intervals 2
dx = 25;
dy = 25;
dz = 1;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STATIC PETROPHYSICS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Transform grids
[phi_seismic, mask, X, Y, Z] = upscale_geo2seis(G,petro_phi,dx,dy,dz);
[vsh_seismic] = upscale_geo2seis(G,petro_vshale,dx,dy,dz);
[p2013_seismic] = upscale_geo2seis(G,petro_p1,dx,dy,dz);
[p2024_seismic] = upscale_geo2seis(G,petro_p2,dx,dy,dz);
[pe2013_seismic] = upscale_geo2seis(G,petro_pe1,dx,dy,dz);
[pe2024_seismic] = upscale_geo2seis(G,petro_pe2,dx,dy,dz);
[sw2013_seismic] = upscale_geo2seis(G,petro_sw1,dx,dy,dz);
[sw2024_seismic] = upscale_geo2seis(G,petro_sw2,dx,dy,dz);

% Permute to have time in the first index, sgy libs work in this way
phi_seismic = permute(phi_seismic, [3, 2, 1]); 
vsh_seismic = permute(vsh_seismic, [3, 2, 1]); 
p2013_seismic = permute(p2013_seismic, [3, 2, 1]); 
p2024_seismic = permute(p2024_seismic, [3, 2, 1]); 
pe2013_seismic = permute(pe2013_seismic, [3, 2, 1]); 
pe2024_seismic = permute(pe2024_seismic, [3, 2, 1]); 
sw2013_seismic = permute(sw2013_seismic, [3, 2, 1]); 
sw2024_seismic = permute(sw2024_seismic, [3, 2, 1]); 
Y = permute(Y, [3, 2, 1]); 
X = permute(X, [3, 2, 1]); 
Z = permute(Z, [3, 2, 1]); 

%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME-LAPSE 2013 %%%
%%%%%%%%%%%%%%%%%%%%%%%

%%  Transform grids
[Vp2013_seismic] = upscale_geo2seis(G,Vp1,dx,dy,dz);
[Vs2013_seismic ] = upscale_geo2seis(G,Vs1,dx,dy,dz);
[Rho2013_seismic ] = upscale_geo2seis(G,Rho1,dx,dy,dz);

% Permute to have time in the first index, sgy libs work in this way
Vp2013_seismic = permute(Vp2013_seismic, [3, 2, 1]); 
Vs2013_seismic = permute(Vs2013_seismic, [3, 2, 1]); 
Rho2013_seismic = permute(Rho2013_seismic, [3, 2, 1]); 

% Define a new facies for the upper and lower zones 
vp_outer = 4000;
vs_outer = 2000;
rho_outer = 2.5;
Vp2013_seismic(isnan(Vp2013_seismic)) = vp_outer; 
Vs2013_seismic(isnan(Vs2013_seismic)) = vs_outer; 
Rho2013_seismic(isnan(Rho2013_seismic)) = rho_outer; 

%% Convert do time
t0 = 2000; % velocity of the upper layer not included in the grid: min(Z(:))/(t0/2) = 2.8688km/s for t0=2000
dt_fine = 0.5;
dt_seis = 2; % in meters, it is approximatly mean(Vp1(:))*(dt_seis/1000)/2 = 3.66 for dt = 2
[Vp2013_seismic_time, time2013] = convert2time(Vp2013_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);
Vs2013_seismic_time = convert2time(Vs2013_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);
Rho20130_seismic_time = convert2time(Rho2013_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);

%treat nans
Vp2013_seismic_time(isnan(Vp2013_seismic_time)) = vp_outer;
Vs2013_seismic_time(isnan(Vs2013_seismic_time)) = vs_outer;
Rho20130_seismic_time(isnan(Rho20130_seismic_time)) = rho_outer;


figure
subplot(211)
imagesc(squeeze(Vp2013_seismic(:,150,:)))
subplot(212)
imagesc(squeeze(Vp2013_seismic_time(:,150,:)))

%% Make SGY files and compute seismic
load('C:\Users\leandro\Documents\GIT\Grids - Felipe\Data\wavelet.mat')
% Convert petrophysics to time:
phi_seismic_time = convert2time(phi_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);
vsh_seismic_time = convert2time(vsh_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);
p2013_seismic_time = convert2time(p2013_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);
p2024_seismic_time = convert2time(p2024_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);
pe2013_seismic_time = convert2time(pe2013_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);
pe2024_seismic_time = convert2time(pe2024_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);
sw2013_seismic_time = convert2time(sw2013_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);
sw2024_seismic_time = convert2time(sw2024_seismic,dz,Vp2013_seismic,t0,dt_fine,dt_seis);

Seismic2013 = compute_seismic(4, Rho20130_seismic_time,Vp2013_seismic_time, Vs2013_seismic_time, [10, 20, 30, 40],wavelet);

%write_segy(dt_seis, Vp2013_seismic_time, X, Y, t0, 'Pvelocity2013.sgy')
%write_segy(dt_seis, phi_seismic_time, X, Y, t0, 'phi.sgy')
%write_segy(dt_seis, sw2024_seismic_time, X, Y, t0, 'sw2024.sgy')
%write_segy(dt_seis, Seismic2013(:,:,:,1), X, Y, t0, 'seismic2013_10deg_NoiseFree.sgy')
%write_segy(dt_seis, Seismic2013(:,:,:,2), X, Y, t0, 'seismic2013_20deg_NoiseFree.sgy')
%write_segy(dt_seis, Seismic2013(:,:,:,3), X, Y, t0, 'seismic2013_30deg_NoiseFree.sgy')
%write_segy(dt_seis, Seismic2013(:,:,:,4), X, Y, t0, 'seismic2013_40deg_NoiseFree.sgy')



%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME-LAPSE 2024 %%%
%%%%%%%%%%%%%%%%%%%%%%%

[Vp2024_seismic ] = upscale_geo2seis(G,Vp2,dx,dy,dz);
[Vs2024_seismic ] = upscale_geo2seis(G,Vs2,dx,dy,dz);
[Rho2024_seismic ] = upscale_geo2seis(G,Rho2,dx,dy,dz);

% Permute to have time in the first index, sgy libs work in this way
Vp2024_seismic = permute(Vp2024_seismic, [3, 2, 1]); 
Vs2024_seismic = permute(Vs2024_seismic, [3, 2, 1]); 
Rho2024_seismic = permute(Rho2024_seismic, [3, 2, 1]); 

% Define a new facies for the upper and lower zones 
Vp2024_seismic(isnan(Vp2024_seismic)) = vp_outer; 
Vs2024_seismic(isnan(Vs2024_seismic)) = vs_outer; 
Rho2024_seismic(isnan(Rho2024_seismic)) = rho_outer; 

%% Convert do time
t0 = 2000;
dt_fine = 0.5;
dt_seis = 2;
[Vp2024_seismic_time, time2024]= convert2time(Vp2024_seismic,dz,Vp2024_seismic,t0,dt_fine,dt_seis);
Vs2024_seismic_time = convert2time(Vs2024_seismic,dz,Vp2024_seismic,t0,dt_fine,dt_seis);
Rho2024_seismic_time = convert2time(Rho2024_seismic,dz,Vp2024_seismic,t0,dt_fine,dt_seis);


Seismic2024 = compute_seismic(4, Rho2024_seismic_time,Vp2024_seismic_time, Vs2024_seismic_time, [10, 20, 30, 40],wavelet);


figure
subplot(211)
imagesc(squeeze(Vp2024_seismic(:,150,:)))
subplot(212)
imagesc(squeeze(Vp2024_seismic_time(:,150,:)))


%% PLOTS
figure
subplot(211)
h1 = pcolor([1:size(Seismic2013,3)],time2013,squeeze(Seismic2013(:,150,:,end)));
caxis([-3 3])
shading interp
set(gca,'Ydir','reverse')
subplot(212)
h2 = pcolor([1:size(Seismic2024,3)],time2024,squeeze(Seismic2024(:,150,:,end)));
colormap(seismic_simple)
caxis([-3 3])
shading interp
set(gca,'Ydir','reverse')
ax1 = gca;
ax2 = subplot(211);

linkaxes([ax1, ax2], 'xy'); % Link x and y axes

