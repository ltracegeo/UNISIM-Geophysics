function [property_seismic,mask_f,x_seismic,y_seismic,z_seismic] = upscale_geo2seis_krig(G,property_grid,dx,dy,dz)

x_grid = G.cells.centroids(:,1);
y_grid = G.cells.centroids(:,2);
z_grid = G.cells.centroids(:,3);

x_seismic = [-dx/2+min(x_grid):dx:max(x_grid)+dx/2]';
y_seismic = [-dy/2+min(y_grid):dy:max(y_grid)+dy/2]';
z_seismic = [-2*dz/2+min(z_grid):dz:max(z_grid)+2*dz]';

[X,Y,Z] = meshgrid(x_seismic,y_seismic,z_seismic);

%% Transforming the geo grid coords x y z to seismic i j k 
A = (size(X,1) - 1) / (max(x_seismic(:)) - min(x_seismic(:)));
B = 1 - A*min(x_seismic(:));
i_grid = round(A*x_grid(:) + B);

A = (size(X,2) - 1) / (max(y_seismic(:)) - min(y_seismic(:)));
B = 1 - A*min(y_seismic(:));
j_grid = round(A*y_grid(:) + B);

A = (size(X,3) - 1) / (max(z_seismic(:)) - min(z_seismic(:)));
B = 1 - A*min(z_seismic(:));
k_grid = round(A*z_grid(:) + B);

%% Using kriging methods for better accuracy
% TODO:
% Levar Z em consideração
% Diminuir L
depth = size(Z,3);
krig_cube = zeros(size(Z));
% for d=1:depth - 1
for d=250
    %Copies the grids for editing
    x_slice = x_grid;
    y_slice = y_grid;
    z_slice = z_grid;
    dvalues = property_grid';
    
    %Uses only coordenates from current slice
    Xslice = X(:,:,d);
    Yslice = Y(:,:,d);
    Zslice = Z(:,:,d);
    xcoords = [Yslice(:) Xslice(:) Zslice(:)];
    
    %Filters out all values which are too far away from current slice
    filter = Z(1,1,d);
    upper = filter + 1;
    lower = filter - 1;
    outside = find(z_grid > upper | z_grid < lower); 
    x_slice(outside) = [];
    y_slice(outside) = [];
    z_slice(outside) = [];
    dvalues(outside) = [];
    
    dcoords = [y_slice, x_slice, z_slice];
    xvar = var(dvalues,0,'all');
    l = 25;
    type = 'exp';
    krig = 1;
    angles = [0, 0, 0];
    dim = size(Z);
    dimX = dim(1);
    dimY = dim(2);
    
    [xok, ~] = Kriging_options(xcoords, dcoords, dvalues', xvar, l, type, krig, angles);
    slice = reshape(xok, [dimX, dimY]);
    krig_cube(:,:,d) = slice;
end

%{
firstX = X(:,:,250);
firstY = Y(:,:,250);
xcoords = [firstX(:) firstY(:)];
dcoords = [x_grid y_grid];
dvalues = property_grid';

xvar = var(dvalues,0,'all');
l = 1;
type = 'exp';
krig = 1;
angles = [0, 0, 0];

[xok, vok] = Kriging_options(xcoords, dcoords, dvalues, xvar, l, type, krig, angles);
%}
%% Accumulating geo grid matrix indexes into seismic i j k index
positions = [i_grid j_grid k_grid];
index2upscale = cell(size(X));
for pos = 1:length(positions)
    index2upscale{positions(pos,1),positions(pos,2),positions(pos,3)}(end+1) = pos;
end

%% Performing upscaling by computing the median of the acumulated geo grid values within a seismic index/voxel i j k
property_seismic = zeros(size(X)) + nan;
parfor i=1:numel(X)
    i_neig = i;
    % if there is not geo grid values within seismic i j k, add neighbor to compute median
    if isempty(index2upscale{i})
        [ii,jj,kk] = ind2sub(size(X),i);
        
        neighbors = [
            -1, 0, 0;  % left neighbor
            1, 0, 0;   % right neighbor
            0, -1, 0;  % front neighbor
            0, 1, 0;   % back neighbor
            0, 0, -1;  % bottom neighbor
            0, 0, 1;   % top neighbor
            ];
        %     neighbors = [
        %         -1, -1, -1; -1, -1, 0; -1, -1, 1;
        %         -1, 0, -1; -1, 0, 0; -1, 0, 1;
        %         -1, 1, -1; -1, 1, 0; -1, 1, 1;
        %         0, -1, -1; 0, -1, 0; 0, -1, 1;
        %         0, 0, -1; 0, 0, 1;
        %         0, 0, 0;
        %         0, 1, -1; 0, 1, 0; 0, 1, 1;
        %         1, -1, -1; 1, -1, 0; 1, -1, 1;
        %         1, 0, -1; 1, 0, 0; 1, 0, 1;
        %         1, 1, -1; 1, 1, 0; 1, 1, 1;
        %     ];
        
        iijjkk = [ii,jj,kk] + neighbors;        
        iijjkk(any(iijjkk<1,2),:)=[]; % Removing lower than 1
        iijjkk(any(iijjkk>size(X),2),:)=[]; % Removing higher than matrix size
        
        i_neig = sub2ind(size(X),iijjkk(:,1),iijjkk(:,2),iijjkk(:,3));
    end
    if ~isempty(cat(2,index2upscale{i_neig}))
        property_seismic(i) = median(property_grid(cat(2,index2upscale{i_neig})));
    end
     
end

%% Define a mask with the valid points in the seismis domain
mask = zeros(size(X));
mask(~isnan(property_seismic)) = 10;
mask = (mask);
mask_f = imgaussfilt3(mask,1);
mask_f(mask_f<4) = 0;
mask_f(mask_f>=4) = 1;

%% Interpolate few positions within the mask and out of the property_seismic
indexes2interpolate = find( mask_f==1 & isnan(property_seismic));
indexes_already_upscaled = find( mask_f==1 &  ~isnan(property_seismic));

F = scatteredInterpolant(X(indexes_already_upscaled),Y(indexes_already_upscaled),Z(indexes_already_upscaled),property_seismic(indexes_already_upscaled));
Vq = F(X(indexes2interpolate),Y(indexes2interpolate),Z(indexes2interpolate));

property_seismic (indexes2interpolate) = Vq ;

krig_cube(property_seismic==0) = 0;
property_seismic = krig_cube;
