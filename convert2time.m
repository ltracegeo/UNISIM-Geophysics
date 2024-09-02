function [property_seis_time,time_regular_seis] = convert2time(property,dz,Vp,t0,dt_fine,dt_seis)
%Convert do time. First we interpolate in a fine regular time grid (sampling rate dt_fine) and
% then we upscale it to a coarse grid of seismic (sampling rate dt_seis)


[~,J,K] = size(property);

time = t0 + 2*1000*cumsum(dz./Vp,1);
time_regular_fine = t0:dt_fine:max(time(:));
time_regular_seis = t0+dt_seis/2:dt_seis:time_regular_fine(end)-dt_seis/2;
winMA = ones(dt_seis/dt_fine,1);
winMA = winMA/sum(winMA);
property_seis_time = zeros(numel(time_regular_seis),J,K);
for j = 1:size(property,2)
    for k = 1:size(property,3)        
        property_trace = interp1(time(:,j,k),property(:,j,k),time_regular_fine)';
        property_trace = conv(property_trace,winMA,'same');
        property_trace = interp1(time_regular_fine,property_trace,time_regular_seis);
        property_seis_time(:,j,k) = property_trace;
    end
end

end
