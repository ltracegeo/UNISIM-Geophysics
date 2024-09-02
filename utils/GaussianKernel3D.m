function G = GaussianKernel3D(I, J, K, std_i, std_j, std_k)
    % Create coordinate grids
    [x, y, z] = ndgrid(1:I, 1:J, 1:K);
    
    % Center coordinates
    cx = (I + 1) / 2;
    cy = (J + 1) / 2;
    cz = (K + 1) / 2;
    
    % Compute the Gaussian function
    G = exp(-(((x - cx).^2) / (2 * std_i^2) + ((y - cy).^2) / (2 * std_j^2) + ((z - cz).^2) / (2 * std_k^2)));
    
    % Normalize the Gaussian to ensure the sum of all elements is 1 (optional)
    G = G / sum(G(:));
end