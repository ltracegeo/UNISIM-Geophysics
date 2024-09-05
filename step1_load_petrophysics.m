%
%% initializing
close all
addpath(genpath('.\MRST'))
mrstModule add deckformat

%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME-LAPSE 2013 %%%
%%%%%%%%%%%%%%%%%%%%%%%

%% Load, process Porosity
% Process
grdecl = fullfile('.\Data\porosity\data.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);
%Processa a grid e obtem os centroides
G = processGRDECL(grdecl,'checkgrid', true);
G(:,2) = []; %Remove todas as geometrias incluidas no arquivo menos a primeira;
G = computeGeometry(G);

petro_phi = grdecl.PORO(G.cells.indexMap);

%% Load, process netgross/vshale
% Process
grdecl = fullfile('.\Data\ntg\data.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

petro_vshale = grdecl.PORO(G.cells.indexMap);
petro_vshale = 1 - petro_vshale;

%% Load, process sw1
% Process
grdecl = fullfile('.\Data\sw1\data.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

petro_sw1 = grdecl.PORO(G.cells.indexMap);

%% Load, process pe1
% Process
grdecl = fullfile('.\Data\pe1\data.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

petro_pe1 = grdecl.PORO(G.cells.indexMap);

%% Load, process p1
% Process
grdecl = fullfile('.\Data\p1\data.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

petro_p1 = grdecl.PORO(G.cells.indexMap);

%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME-LAPSE 2024 %%%
%%%%%%%%%%%%%%%%%%%%%%%

%% Load, process sw2
% Process
grdecl = fullfile('.\Data\sw2\data.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

petro_sw2 = grdecl.PORO(G.cells.indexMap);

%% Load, process pe2
% Process
grdecl = fullfile('.\Data\pe2\data.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

petro_pe2 = grdecl.PORO(G.cells.indexMap);

%% Load, process p2
% Process
grdecl = fullfile('.\Data\p2\data.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

petro_p2 = grdecl.PORO(G.cells.indexMap);


%% Save .mat
%save('.\Data\petrophysics.mat');





