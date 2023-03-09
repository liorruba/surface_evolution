function Z = linear_slope_diffusion(Z, slope_of_repose, res)
[u,v] = gradient(Z, res);
slope_map = sqrt(u.^2 + v.^2);

% Choose k based on the Courant condition
K = 0.25 * res;
K = K  ./ res.^2;
kk = 1;
while any(slope_map(:) > slope_of_repose)
    [u,v] = gradient(Z, res);
    slope_map = sqrt(u.^2 + v.^2);
    for ii = 2:(length(Z)-1)
        for jj = 2:(length(Z)-1)
            buffii = (Z(ii + 1,jj) - 2*Z(ii,jj) + Z(ii - 1,jj));
            buffjj = (Z(ii,jj + 1) - 2*Z(ii,jj) + Z(ii,jj - 1));
            Z(ii,jj) = Z(ii,jj) + K * (buffii + buffjj);
        end
    end
    plot((Z(250,:)))
    drawnow
    kk = kk + 1;
end
end
