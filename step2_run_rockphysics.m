
addpath(genpath('.\SeReM'))
addpath(genpath('.\SeisLab_10.0301'))
addpath(genpath('.\utils'))

%% Rocl-physics parameters: 
Kminc = [37 21];
Gminc = [44 7];
Rhominc = [2.65 2.58];

temp = 80;
salinity = 0.06;
GOR = interp1([193.36 213.26],[105.42 115.01],210); % Reference value from the CMG case assuming PB=210
gas_gravity = 0.75;
api=19;

coordnumber = 9;
criticalporo = 0.4;


%% Load well, run rock-physics and export wells
% 
% load('.\Data\WellNames.mat')
% 
% well_logs_table = [];
% for w = 1:size(WellNames,1)
% 
%     %%% Load wells
%     file_name = WellNames.Well(w) + '.las';
%     file_path = char('.\Data\well_las\' + file_name);    
%     
%     %%% Create elastic logs in the well structure
%     wells{w} = read_las_file(file_path);    
%     wells{w}.curve_info{4,1} = 'SW_2013'; wells{w}.curve_info{4,3} = 'SW_2013';
%     wells{w}.curve_info{5,1} = 'SW_2024'; wells{w}.curve_info{5,3} = 'SW_2024';
%     wells{w}.curve_info{6,1} = 'Press_2013'; wells{w}.curve_info{6,3} = 'Press_2013';
%     wells{w}.curve_info{7,1} = 'Press_2024'; wells{w}.curve_info{7,3} = 'Press_2024';
%     wells{w}.curve_info{8,1} = 'EffPress_2013'; wells{w}.curve_info{8,3} = 'EffPress_2013';
%     wells{w}.curve_info{9,1} = 'EffPress_2024'; wells{w}.curve_info{9,3} = 'EffPress_2024';    
%     wells{w}.curve_info{end+1,1} = 'vshale'; wells{w}.curve_info{end,2} = 'm3/m3'; wells{w}.curve_info{end,3} = 'Vshale';
%     wells{w}.curve_info{end+1,1} = 'Vp2013'; wells{w}.curve_info{end,2} = 'm/s'; wells{w}.curve_info{end,3} = 'Vp2013';
%     wells{w}.curve_info{end+1,1} = 'Vs2013'; wells{w}.curve_info{end,2} = 'm/s'; wells{w}.curve_info{end,3} = 'Vs2013';
%     wells{w}.curve_info{end+1,1} = 'Rho2013'; wells{w}.curve_info{end,2} = 'g/cm3'; wells{w}.curve_info{end,3} = 'Rho2013';
%     wells{w}.curve_info{end+1,1} = 'Ip2013'; wells{w}.curve_info{end,2} = 'kg/(m2.s)'; wells{w}.curve_info{end,3} = 'Ip2013';
%     wells{w}.curve_info{end+1,1} = 'VpVs2013'; wells{w}.curve_info{end,2} = '_'; wells{w}.curve_info{end,3} = 'VpVs2013';    
%     wells{w}.curve_info{end+1,1} = 'Vp2024'; wells{w}.curve_info{end,2} = 'm/s'; wells{w}.curve_info{end,3} = 'Vp2024';
%     wells{w}.curve_info{end+1,1} = 'Vs2024'; wells{w}.curve_info{end,2} = 'm/s'; wells{w}.curve_info{end,3} = 'Vs2024';
%     wells{w}.curve_info{end+1,1} = 'Rho2024'; wells{w}.curve_info{end,2} = 'g/cm3'; wells{w}.curve_info{end,3} = 'Rho2024';
%     wells{w}.curve_info{end+1,1} = 'Ip2024'; wells{w}.curve_info{end,2} = 'kg/(m2.s)'; wells{w}.curve_info{end,3} = 'Ip2024';
%     wells{w}.curve_info{end+1,1} = 'VpVs2024'; wells{w}.curve_info{end,2} = '_'; wells{w}.curve_info{end,3} = 'VpVs2024';
%     wells{w}.curve_info{end+1,1} = 'Time'; wells{w}.curve_info{end,2} = 'ms'; wells{w}.curve_info{end,3} = 'OWT';
%         
%     %%% Interp NANs:
%     for curve=2:size(wells{w}.curves,2)
%         wells{w}.curves(:,curve) = naninterp(wells{w}.curves(:,curve));
%     end
%     
%     %%% Treat
%     petro_phi_log = wells{w}.curves(:,2);
%     petro_vshale_log = wells{w}.curves(:,3);
%     petro_vshale_log = 1 - petro_vshale_log;   
%     petro_sw1_log = wells{w}.curves(:,4);
%     petro_sw2_log = wells{w}.curves(:,5);
%     petro_p1_log = wells{w}.curves(:,6);
%     petro_p2_log = wells{w}.curves(:,7);
%     petro_pe1_log = wells{w}.curves(:,8);
%     petro_pe2_log = wells{w}.curves(:,9);
%     
%     petro_phi_log(petro_phi_log >0.9)=0.01; % shale has porosity = 1 in the model
%     %petro_phi_log(petro_vshale_log==1) = 0.0; % spurious values of pure shale with high porosity in the model
%     petro_phi_log(petro_vshale_log==1) = petro_phi_log(petro_vshale_log==1)*0.12; % spurious values of pure shale with high porosity in the model
%     petro_phi_log(petro_phi_log>=criticalporo) = criticalporo - 0.01;
%     petro_phi_log(petro_phi_log<=0) = 0.01;
%     petro_vshale_log(petro_vshale_log>=0.99) = 0.99;
%     petro_vshale_log(petro_vshale_log<=0) = 0.02;
%     %petro_vshale_log = petro_vshale_log .* (1 - petro_vshale_log);
%     petro_sw1_log(petro_sw1_log>=0.99) = 0.99;
%     petro_sw1_log(petro_sw1_log==0) = 0.99; % shale has porosity = 1 in the model
%     petro_pe1_log(petro_pe1_log==0) = median(petro_pe1_log(:));
%     petro_p1_log(petro_p1_log==0) = median(petro_p1_log(:));
%     petro_pe1_log = petro_pe1_log/1e6; % convert to GB
%     petro_p1_log = petro_p1_log/1e6; % convert to GB
% 
%     petro_sw2_log(petro_sw2_log>=0.99) = 0.99;
%     petro_sw2_log(petro_sw2_log==0) = 0.99;
%     petro_pe2_log(petro_pe2_log==0) = median(petro_pe2_log(:));
%     petro_p2_log(petro_p2_log==0) = median(petro_p2_log(:));
%     petro_pe2_log = petro_pe2_log/1e6; % convert to GB
%     petro_p2_log = petro_p2_log/1e6; % convert to GB    
%     
%     wells{w}.curves(:,2) = petro_phi_log;
%     wells{w}.curves(:,10) = petro_vshale_log;
%     wells{w}.curves(:,4) = petro_sw1_log;
%     wells{w}.curves(:,5) = petro_sw2_log;
%     wells{w}.curves(:,6) = petro_p1_log;
%     wells{w}.curves(:,7) = petro_p2_log;
%     wells{w}.curves(:,8) = petro_pe1_log;
%     wells{w}.curves(:,9) = petro_pe2_log;     
%     
%     %%% Well markers (top and base):
%     markers(w,1) = wells{w}.curves(1,1);
%     markers(w,2) = wells{w}.curves(end,1);    
%                
%     %%% Run rock-physics Time-lapse 1:       
%     for pto = 1:size(wells{w}.curves,1)
%         Volminc = [ 1 - petro_vshale_log(pto) petro_vshale_log(pto) ]; 
%         [K_bri, rho_bri] = BatzleWangBrine(temp, petro_p1_log(pto), salinity);
%         [K_oil, rho_oil] = BatzleWangOil(temp, petro_p1_log(pto), GOR, api, gas_gravity);
%         Kflc = [K_bri K_oil];
%         Rhoflc = [rho_bri rho_oil];
%         Sflc = [petro_sw1_log(pto) 1-petro_sw1_log(pto)];
%         patchy = 0;
%         [Kmat, Gmat, Rhomat, Kfl, Rhofl] = MatrixFluidModel (Kminc, Gminc, Rhominc, Volminc, Kflc, Rhoflc, Sflc, patchy);
%         Rho1(pto) = DensityModel(petro_phi_log(pto), Rhomat, Rhofl);
%         [Vp1(pto), Vs1(pto)] = SoftsandModel(petro_phi_log(pto), Rho1(pto), Kmat, Gmat, Kfl, criticalporo, coordnumber, petro_pe1_log(pto));        
%     end
%     Vp1 = Vp1*1000; Vs1 = Vs1*1000;
%     Ip1 = Vp1.*Rho1;
%     VpVs1 = Vp1./Vs1;
%     
%     %%% Run rock-physics Time-lapse 2:
%     for pto = 1:size(wells{w}.curves,1)
%         Volminc = [ 1 - petro_vshale_log(pto) petro_vshale_log(pto) ];
%         [K_bri, rho_bri] = BatzleWangBrine(temp, petro_p2_log(pto), salinity); % ok
%         [K_oil, rho_oil] = BatzleWangOil(temp, petro_p2_log(pto), GOR, api, gas_gravity);
%         Kflc = [K_bri K_oil];
%         Rhoflc = [rho_bri rho_oil];
%         Sflc = [petro_sw2_log(pto) 1-petro_sw2_log(pto)];
%         patchy = 0;
%         [Kmat, Gmat, Rhomat, Kfl, Rhofl] = MatrixFluidModel (Kminc, Gminc, Rhominc, Volminc, Kflc, Rhoflc, Sflc, patchy);
%         Rho2(pto) = DensityModel(petro_phi_log(pto), Rhomat, Rhofl);
%         [Vp2(pto), Vs2(pto)] = SoftsandModel(petro_phi_log(pto), Rho2(pto), Kmat, Gmat, Kfl, criticalporo, coordnumber, petro_pe2_log(pto));        
%     end
%     Vp2 = Vp2*1000; Vs2 = Vs2*1000;
%     Ip2 = Vp2.*Rho2;
%     VpVs2 = Vp2./Vs2;        
% 
%     %%% Compute time:
%     depth = wells{w}.curves(:,1);
%     v0 = 2868.8; %m/s      
%     v0_outer = 4000; %m/s %v0_outer = 4005.11; %m/s    
%     
%     ways = 1; % if is OneWayT or TwoWayT
%     t0 = ways*1000;
%     t0_outer = ways*1000*(depth(1) - v0/(t0/(ways*1000)))/v0_outer;    
%     t0 = t0 + t0_outer;
%     dz = [0; diff(depth)];
%     time = t0 + ways*1000*cumsum(dz./Vp1',1);
%     time = time - (time(2)-time(1));
%     
%     %%% Save elastic properties in struct:
%     wells{w}.curves(:,11) = Vp1;
%     wells{w}.curves(:,12) = Vs1;
%     wells{w}.curves(:,13) = Rho1;
%     wells{w}.curves(:,14) = Ip1;
%     wells{w}.curves(:,15) = VpVs1;
%     wells{w}.curves(:,16) = Vp2;
%     wells{w}.curves(:,17) = Vs2;
%     wells{w}.curves(:,18) = Rho2;
%     wells{w}.curves(:,19) = Ip2;
%     wells{w}.curves(:,20) = VpVs2;            
%     wells{w}.curves(:,21) = time;            
%     
%     clear Vp1 Vs1 Rho1 Vp2 Vs2 Rho2
% 
%     file_path = char('.\Export\Wells\' + file_name);
%     write_las_file(wells{w},file_path); 
%     
% end
% 

%%   COMPUTE TIME-LAPSE 2013 USING SOFT SAND MODEL

load('.\Data\petrophysics.mat')


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







