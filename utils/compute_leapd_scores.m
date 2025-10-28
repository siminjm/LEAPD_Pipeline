function s = compute_leapd_scores(V, P0,m0, P1,m1, dstar, is_norm)
% COMPUTE_LEAPD_SCORES - distance-to-hyperplane ratio (in [0,1])
if nargin<6 || isempty(dstar), dstar = size(P0,2); end
if nargin<7, is_norm = 0; end

if isvector(V), V = V(:)'; end
Y0 = bsxfun(@minus, V, m0);
Y1 = bsxfun(@minus, V, m1);
proj0 = (Y0*P0) * P0';
proj1 = (Y1*P1) * P1';
if is_norm==0
    d1 = vecnorm(Y0 - proj0, 2, 2);
    d2 = vecnorm(Y1 - proj1, 2, 2);
else
    d1 = vecnorm(proj0,2,2)./max(vecnorm(Y0,2,2),eps);
    d2 = vecnorm(proj1,2,2)./max(vecnorm(Y1,2,2),eps);
end
s = d2 ./ max(d1 + d2, eps);
end
