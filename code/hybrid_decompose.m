function [F_RF,F_BB,W] = hybrid_decompose(Wd,Nrf,maxiter,F_RF_init)

if nargin < 4
    F_RF_init = [];
end
if nargin < 3 || isempty(maxiter)
    maxiter = 20;
end

[M,~] = size(Wd);
if isempty(F_RF_init) || size(F_RF_init,2) ~= Nrf || size(F_RF_init,1) ~= M
    F_RF = exp(1i*angle(randn(M,Nrf)+1i*randn(M,Nrf)))/sqrt(M);
else
    F_RF = exp(1i*angle(F_RF_init))/sqrt(M);
end

for iter = 1:1:maxiter
    F_BB = (F_RF'*F_RF)\(F_RF'*Wd);
    F_RF = exp(1i*angle(Wd*F_BB'))/sqrt(M);
end

W = F_RF*F_BB;

end
