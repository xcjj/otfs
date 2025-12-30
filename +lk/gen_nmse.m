function nmse = gen_nmse(estimate, truth)
%gen_nmse Compute normalized mean-square error.
%
% Description:
%   nmse = sum(|estimate-truth|^2) / sum(|truth|^2)
%
% Developer: Jia. Institution: PML. Date: 2024/06/02

num = sum(abs(estimate(:) - truth(:)).^2);
den = sum(abs(truth(:)).^2);
nmse = num / max(den, eps);
end
