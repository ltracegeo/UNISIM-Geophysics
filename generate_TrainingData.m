    
clear

addpath(genpath('.\SeReM'))
addpath(genpath('.\SeisLab_10.0301'))
addpath(genpath('.\utils'))

%% MONTE CARLO SAMPLING

% mid porosity
n_sim = 50000;
Phi_train = 0.1 + 0.05*randn(n_sim,1);
vshale = 0.05 + 0.05*randn(n_sim,1);
sw_train13 = rand(n_sim,1);
sw_train24 = rand(n_sim,1);
facies_train = ones(n_sim,1);


% high porosity
n_sim = 50000;
Phi_train = [ Phi_train ; 0.25 + 0.05*randn(n_sim,1) ];
vshale = [vshale; 0.25 + 0.075*randn(n_sim,1)];
sw_train13 = [sw_train13; rand(n_sim,1)];
sw_train24 = [sw_train24; rand(n_sim,1)];
facies_train = [facies_train; 2*ones(n_sim,1)];

% Shale
n_sim = 50000;
Phi_train = [ Phi_train ; 0.02 + 0.005*randn(n_sim,1) ];
vshale = [vshale; 0.95 + 0.002*randn(n_sim,1)];
sw_train13 = [ sw_train13 ; 0.9 + 0.02*randn(n_sim,1) ];
sw_train24 = [ sw_train24 ; 0.9 + 0.02*randn(n_sim,1) ];
facies_train = [facies_train; 3*ones(n_sim,1)];

% truncate values
Phi_train(Phi_train<0) = 0.001;
Phi_train(Phi_train>=0.4) = 0.39;
sw_train13(sw_train13>1) = 0.999;
sw_train24(sw_train24>1) = 0.999;
vshale(vshale>=0.99) = 0.99;
vshale(vshale<=0) = 0.02;

% pressure for all facies
n_sim_total = length(Phi_train);
petro_p13 = 0.03248 + 0.0003*randn(n_sim_total,1);
petro_p24 = 0.03 + 0.0025*randn(n_sim_total,1);
petro_pe13 = 0.0377 + 0.0008*randn(n_sim_total,1);
petro_pe24 = 0.0399 + 0.0016*randn(n_sim_total,1);



%%  Simulate observed data (elastic properties) with noise
criticalporo = 0.4;
std_vp = 75;
std_vs = 37.5;
std_rho = 0.075;

[Vp, Vs, Rho] = RPM_unisim(Phi_train, sw_train13, vshale, petro_p13, petro_pe13 );
Vp = Vp + std_vp*randn(size(Vp));
Vs = Vs + std_vs*randn(size(Vs));
Rho = Rho + std_rho*randn(size(Rho));
Ip_train13 = Vp.*Rho;
VPVS_train13 = Vp./Vs;

[Vp, Vs, Rho] = RPM_unisim(Phi_train, sw_train24, vshale, petro_p24, petro_pe24 );
Vp = Vp + std_vp*randn(size(Vp));
Vs = Vs + std_vs*randn(size(Vs));
Rho = Rho + std_rho*randn(size(Rho));
Ip_train24 = Vp.*Rho;
VPVS_train24 = Vp./Vs;

%% EXPORT TO LAS FILE

well = read_las_file('.\Data\well_las\INJ003.las');

well.curve_info{3,1} = 'VSH'; well.curve_info{3,2} = '_'; well.curve_info{3,3} = 'VSH';
well.curve_info{4,1} = 'SW_2013'; well.curve_info{3,2} = '_'; well.curve_info{3,3} = 'SW_2013';
well.curve_info{5,1} = 'SW_2024'; well.curve_info{4,2} = '_'; well.curve_info{4,3} = 'SW_2024';
well.curve_info{6,1} = 'Press_2013'; well.curve_info{5,2} = 'kPa'; well.curve_info{5,3} = 'Press_2013';
well.curve_info{7,1} = 'EffPress_2013'; well.curve_info{6,2} = 'kPa'; well.curve_info{6,3} = 'EffPress_2013';
well.curve_info{8,1} = 'Press_2024'; well.curve_info{7,2} = 'kPa'; well.curve_info{7,3} = 'Press_2024';
well.curve_info{9,1} = 'EffPress_2024'; well.curve_info{8,2} = 'kPa'; well.curve_info{8,3} = 'EffPress_2024';
well.curve_info{10,1} = 'Ip2013'; well.curve_info{9,2} = 'kg/(m2.s)'; well.curve_info{9,3} = 'Ip2013';
well.curve_info{11,1} = 'VpVs2013'; well.curve_info{10,2} = '_'; well.curve_info{10,3} = 'VpVs2013';
well.curve_info{12,1} = 'Ip2024'; well.curve_info{11,2} = 'kg/(m2.s)'; well.curve_info{11,3} = 'Ip2024';
well.curve_info{13,1} = 'VpVs2024'; well.curve_info{12,2} = '_'; well.curve_info{12,3} = 'VpVs2024';


%depth = linspace(depth(1),depth(end),size(TrainingData,1));
depth = 1:length(Phi_train);
well.step = depth(2)-depth(1);
well.first = depth(1);
well.last = depth(end);

well.curves = [];
well.curves(:,1) = depth;
well.curves(:,2) = Phi_train;
well.curves(:,3) = vshale;
well.curves(:,4) = sw_train13;
well.curves(:,5) = sw_train24;
well.curves(:,6) = petro_p13;
well.curves(:,7) = petro_pe13;
well.curves(:,8) = petro_p24;
well.curves(:,9) = petro_pe24;
well.curves(:,10) = Ip_train13;
well.curves(:,11) = VPVS_train13;
well.curves(:,12) = Ip_train24;
well.curves(:,13) = VPVS_train13;

write_las_file(well,'.\Export\Wells\TraininData_well.las'); 























