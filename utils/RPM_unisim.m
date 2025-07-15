function [Vp, Vs, Rho] = RPM_unisim(Phi, sw, vshale, pressure, eff_pressure )



Kminc = [37 21];
Gminc = [44 7];
Rhominc = [2.65 2.58];

temp = 80;
salinity = 0.06;
GOR = interp1([193.36 213.26],[105.42 115.01],210); % Reference value from the CMG case assuming PB=210
gas_gravity = 0.75;
api=19;
patchy = 0;
coordnumber = 9;
criticalporo = 0.4;


for pto = 1:length(Phi)
    Volminc = [ 1 - vshale(pto) vshale(pto)];
    [K_bri, rho_bri] = BatzleWangBrine(temp, pressure(pto), salinity);
    [K_oil, rho_oil] = BatzleWangOil(temp, pressure(pto), GOR, api, gas_gravity);
    Kflc = [K_bri K_oil];
    Rhoflc = [rho_bri rho_oil];
    Sflc = [sw(pto) 1-sw(pto)];
    [Kmat, Gmat, Rhomat, Kfl, Rhofl] = MatrixFluidModel (Kminc, Gminc, Rhominc, Volminc, Kflc, Rhoflc, Sflc, patchy);
    Rho(pto) = DensityModel(Phi(pto), Rhomat, Rhofl);
    [Vp(pto), Vs(pto)] = SoftsandModel(Phi(pto), Rho(pto), Kmat, Gmat, Kfl, criticalporo, coordnumber, eff_pressure(pto));
end
Vp = Vp*1000; Vs = Vs*1000;

%%%%%%%%%%%%%%%%%%%%  ANTIGO  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Temperature = 80;
% 
% SIZE = size(Phi);
% 
% Phi(Phi<=0) = 0.001;
% Phi(Phi>=0.40) = 0.399;
% 
% sw(sw<=0) = 0.001;
% sw(sw>=1.0) = 0.999;
% 
% Phi = Phi(:);
% sw = sw(:);
% 
% 
% coordnumber=9;
% pressure=0.060;
% 
% Kminc = [37];
% Gminc = [44];
% Rhominc = [2.65];
% Volminc = ones(size(sw));
% 
% Kflc = [3.14 0.53];
% Rhoflc = [1.06 0.52];
% Sflc = [sw 1-sw];
% 
% [Kmat, Gmat, Rhomat, Kfl, Rhofl] = MatrixFluidModel (Kminc, Gminc, Rhominc, Volminc, Kflc, Rhoflc, Sflc, 0);
% 
% Rho = DensityModel(Phi, Rhomat, Rhofl);
% 
% [Vp, Vs] = SoftsandModel(Phi, Rho, Kmat, Gmat, Kfl, criticalporo, coordnumber, pressure);
% 
% Vp = Vp * 1000;
% Vs = Vs * 1000;
% 
% Vp  = reshape(Vp,SIZE);
% Vs  = reshape(Vs,SIZE);
% Rho = reshape(Rho,SIZE);
% 
% %Ip = Vp .* Rho;
% %VPVS = Vp./Vs;
% %Ip = reshape(Ip,SIZE);
% %VPVS = reshape(VPVS,SIZE);
