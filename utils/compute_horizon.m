function [horizons] = compute_horizon(phi)

[~,J,K] = size(phi);

for j = 1:size(phi,2)
    for k = 1:size(phi,3)        
        trace = phi(:,j,k);
        horizon_pos_volume = abs(diff(isnan(trace)));
        positions = find(horizon_pos_volume>0);
        if length(positions)>0
            horizons(:,j,k) = [positions(1); positions(end)];
        end
    end
end


horizonte1 = squeeze(horizons(1,:,:));
valid_ind = find(horizonte1>0);
tointerp_ind = find(horizonte1==0);

[X,Y] = meshgrid([1:K],[1:J]);

l=20;
type = 'sph';
krig = 0;
xcoords = [X(tointerp_ind),Y(tointerp_ind)];
dcoords = [X(valid_ind),Y(valid_ind)];
dz = horizonte1(valid_ind);
zvar = var(dz);
[v_grid] = Kriging_options(xcoords, dcoords, dz, zvar, l, type, krig, [0 0 0]);

horizonte1(tointerp_ind) = v_grid;

figure
imagesc(horizonte1)


horizonte2 = squeeze(horizons(2,:,:));
l=5;
type = 'sph';
krig = 0;
xcoords = [X(tointerp_ind),Y(tointerp_ind)];
dcoords = [X(valid_ind),Y(valid_ind)];
dz = horizonte2(valid_ind);
zvar = var(dz);
[v_grid] = Kriging_options(xcoords, dcoords, dz, zvar, l, type, krig, [0 0 0]);

horizonte2(tointerp_ind) = v_grid;


horizons(1,:,:) = horizonte1;
horizons(2,:,:) = horizonte2;
