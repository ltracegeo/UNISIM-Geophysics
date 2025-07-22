function [property_depth,depth_regular] = convert2depth(property,dt,Vp,z0,dz_fine,dz_coarse)
% NOT TESTED
disp('NOT TESTED')

dt = dt/2000;

[~,J,K] = size(property);

time = z0 + cumsum(dt.*Vp,1);
time_regular_fine = z0:dz_fine:max(time(:));
depth_regular = z0+dz_coarse/2:dz_coarse:time_regular_fine(end)-dz_coarse/2;
winMA = ones(dz_coarse/dz_fine,1);
winMA = winMA/sum(winMA);
property_depth = zeros(numel(depth_regular),J,K);
for j = 1:size(property,2)
    for k = 1:size(property,3)
        property_trace = interp1(time(:,j,k),property(:,j,k),time_regular_fine)';
        property_trace = conv(property_trace,winMA,'same');
        property_trace = interp1(time_regular_fine,property_trace,depth_regular);
        property_depth(:,j,k) = property_trace;
    end
end

end
