function [P, m] = build_hyperplanes(V, d)
% BUILD_HYPERPLANES - returns basis P (Kxd) and mean m (1xK)
% V: NxK
m = mean(V,1);
Vc = V - m;
[~,~,Vr] = svd(Vc/ sqrt(size(V,1)-1), 'econ');
d = min(d, size(Vr,2));
P = Vr(:,1:d);
end
