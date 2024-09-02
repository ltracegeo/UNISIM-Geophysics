
addpath(genpath('.\SeReM'))

%%   COMPUTE TIME-LAPSE 2013 USING SOFT SAND MODEL

load('.\Data\petrophysics.mat')

temp = 80;
salinity = 0.06;
GOR = interp1([193.36 213.26],[105.42 115.01],210); % Reference value from the CMG case assuming PB=210
gas_gravity = 0.75;
api=19;

coordnumber = 9;
criticalporo = 0.4;


%% Unit conversion, and treating extreme values
petro_phi(petro_phi>0.9)=0.01; % shale has porosity = 1 in the model
%petro_phi(petro_vshale==1) = 0.0; % spurious values of pure shale with high porosity in the model
petro_phi(petro_vshale==1) = petro_phi(petro_vshale==1)*0.12; % spurious values of pure shale with high porosity in the model
petro_phi(petro_phi>=criticalporo) = criticalporo - 0.01;
petro_phi(petro_phi<=0) = 0.01;
petro_vshale(petro_vshale>=0.99) = 0.99;
petro_vshale(petro_vshale<=0) = 0.02;
%petro_phi = petro_phi .* (1 - petro_vshale);
petro_sw1(petro_sw1>=0.99) = 0.99;
petro_sw1(petro_sw1==0) = 0.99; % shale has porosity = 1 in the model
petro_pe1(petro_pe1==0) = median(petro_pe1(:));
petro_p1(petro_p1==0) = median(petro_p1(:));
petro_pe1 = petro_pe1/1e6; % convert to GB
petro_p1 = petro_p1/1e6; % convert to GB

petro_sw2(petro_sw2>=0.99) = 0.99;
petro_sw2(petro_sw2==0) = 0.99;
petro_pe2(petro_pe2==0) = median(petro_pe2(:));
petro_p2(petro_p2==0) = median(petro_p2(:));
petro_pe2 = petro_pe2/1e6; % convert to GB
petro_p2 = petro_p2/1e6; % convert to GB


%% Apply rock-physics/petroelastic model using SeReM
%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME-LAPSE 2013 %%%
%%%%%%%%%%%%%%%%%%%%%%%
for pto=1:numel(petro_phi)
    
    Kminc = [37 21];
    Gminc = [44 7];
    Rhominc = [2.65 2.58];
    Volminc = [ 1-petro_vshale(pto) petro_vshale(pto) ];
    
    [K_bri(pto), rho_bri(pto)] = BatzleWangBrine(temp, petro_p1(pto), salinity); % ok
    [K_oil(pto), rho_oil(pto)] = BatzleWangOil(temp, petro_p1(pto), GOR, api, gas_gravity);
    
    Kflc = [K_bri(pto) K_oil(pto)];
    Rhoflc = [rho_bri(pto) rho_oil(pto)];
    Sflc = [petro_sw1(pto) 1-petro_sw1(pto)];
    
    patchy = 0;
    
    [Kmat(pto), Gmat(pto), Rhomat(pto), Kfl(pto), Rhofl(pto)] = MatrixFluidModel (Kminc, Gminc, Rhominc, Volminc, Kflc, Rhoflc, Sflc, patchy);
    
    Rho1(pto) = DensityModel(petro_phi(pto), Rhomat(pto), Rhofl(pto));
    
    [Vp1(pto), Vs1(pto)] = SoftsandModel(petro_phi(pto), Rho1(pto), Kmat(pto), Gmat(pto), Kfl(pto), criticalporo, coordnumber, petro_pe1(pto));
    
end

% Convert to m/s
Vp1 = Vp1*1000;
Vs1 = Vs1*1000;

%% Cross-plots
sample_size = 50000; % Ajuste conforme necessário
sample_indices = randsample(length(Vp1), sample_size)';

figure
subplot(221)
scatter(petro_phi(sample_indices),petro_vshale(sample_indices),10,petro_sw1(sample_indices),'filled')
grid
xlabel('Phi')
ylabel('Vshale')
cb = colorbar;
cb.Label.String = 'Sw1';
subplot(222)
scatter(Vp1(sample_indices).*Rho1(sample_indices),Vp1(sample_indices)./Vs1(sample_indices),10,petro_vshale(sample_indices),'filled')
grid
xlabel('AI')
ylabel('Vp/Vs')
cb = colorbar;
cb.Label.String = 'Vshale';
subplot(223)
scatter(Vp1(sample_indices).*Rho1(sample_indices),Vp1(sample_indices)./Vs1(sample_indices),10,petro_sw1(sample_indices),'filled')
grid
xlabel('AI')
ylabel('Vp/Vs')
cb = colorbar;
cb.Label.String = 'Sw1';
subplot(224)
scatter(Vp1(sample_indices).*Rho1(sample_indices),Vp1(sample_indices)./Vs1(sample_indices),10,petro_phi(sample_indices),'filled')
grid
xlabel('AI')
ylabel('Vp/Vs')
cb = colorbar;
cb.Label.String = 'Phi';

%% Apply rock-physics/petroelastic model using SeReM
%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME-LAPSE 2024 %%%
%%%%%%%%%%%%%%%%%%%%%%%
for pto=1:numel(petro_phi)
    
    Kminc = [37 21];
    Gminc = [44 7];
    Rhominc = [2.65 2.58];
    Volminc = [ 1-petro_vshale(pto) petro_vshale(pto) ];
    
    [K_bri(pto), rho_bri(pto)] = BatzleWangBrine(temp, petro_p2(pto), salinity); % ok
    [K_oil(pto), rho_oil(pto)] = BatzleWangOil(temp, petro_p2(pto), GOR, api, gas_gravity);
    
    Kflc = [K_bri(pto) K_oil(pto)];
    Rhoflc = [rho_bri(pto) rho_oil(pto)];
    Sflc = [petro_sw2(pto) 1-petro_sw2(pto)];
    
    patchy = 0;
    
    [Kmat(pto), Gmat(pto), Rhomat(pto), Kfl(pto), Rhofl(pto)] = MatrixFluidModel (Kminc, Gminc, Rhominc, Volminc, Kflc, Rhoflc, Sflc, patchy);
    
    Rho2(pto) = DensityModel(petro_phi(pto), Rhomat(pto), Rhofl(pto));
    
    [Vp2(pto), Vs2(pto)] = SoftsandModel(petro_phi(pto), Rho2(pto), Kmat(pto), Gmat(pto), Kfl(pto), criticalporo, coordnumber, petro_pe2(pto));
    
end

% Convert to m/s
Vp2 = Vp2*1000;
Vs2 = Vs2*1000;

%% Cross-plots
sample_size = 50000; % Ajuste conforme necessário
sample_indices = randsample(length(Vp2), sample_size)';

figure
subplot(221)
scatter(petro_phi(sample_indices),petro_vshale(sample_indices),10,petro_sw2(sample_indices),'filled')
grid
xlabel('Phi')
ylabel('Vshale')
cb = colorbar;
cb.Label.String = 'Sw1';
subplot(222)
scatter(Vp2(sample_indices).*Rho2(sample_indices),Vp2(sample_indices)./Vs2(sample_indices),10,petro_vshale(sample_indices),'filled')
grid
xlabel('AI')
ylabel('Vp/Vs')
cb = colorbar;
cb.Label.String = 'Vshale';
subplot(223)
scatter(Vp2(sample_indices).*Rho2(sample_indices),Vp2(sample_indices)./Vs2(sample_indices),10,petro_sw2(sample_indices),'filled')
grid
xlabel('AI')
ylabel('Vp/Vs')
cb = colorbar;
cb.Label.String = 'Sw1';
subplot(224)
scatter(Vp2(sample_indices).*Rho2(sample_indices),Vp2(sample_indices)./Vs2(sample_indices),10,petro_phi(sample_indices),'filled')
grid
xlabel('AI')
ylabel('Vp/Vs')
cb = colorbar;
cb.Label.String = 'Phi';







