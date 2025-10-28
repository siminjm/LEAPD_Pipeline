function LPC = compute_yw(Xf, order)
% COMPUTE_YW - AR (Burg) coefficients per subject
n = numel(Xf); LPC = zeros(n, order-1);
for i=1:n
    x = Xf{i};
    a = arburg(x, order); % returns [1 a2 a3 ...]
    LPC(i,:) = a(2:end);
end
end
