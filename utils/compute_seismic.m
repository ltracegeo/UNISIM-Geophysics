function [Seismic, Rpp] = compute_seismic(type, rho_cube,vp_cube, vs_cube, angles, wavelets)

n_angles = numel(angles);
Rpp = zeros( size(vp_cube,1), size(vp_cube,2), size(vp_cube,3), n_angles );
Seismic = zeros( size(vp_cube,1), size(vp_cube,2), size(vp_cube,3), n_angles );

for j = 1:size(vp_cube,2)
    for k = 1:size(vp_cube,3)
        
        Rpp_trace = reflec_coef(type, rho_cube(:,j,k), vp_cube(:,j,k), vs_cube(:,j,k), angles);
        Rpp_trace (isnan(Rpp_trace )) = 0;
        
        for ang=1:n_angles
            seismic_trace(:,ang) = conv(Rpp_trace(:,ang),wavelets(:,ang),'same');
        end        
        
        Rpp(:,j,k,:) = reshape(Rpp_trace,size(vp_cube,1),1,1,[]);
        Seismic(:,j,k,:) = reshape(seismic_trace,size(vp_cube,1),1,1,[]);
        
    end
end