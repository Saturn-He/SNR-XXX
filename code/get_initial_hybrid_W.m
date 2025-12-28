% Obtain an initial hybrid beamformer by decomposing a digital solution.
% Inputs: Wd: full-digital beamformer; Nrf: number of RF chains
% Outputs: F_RF: analog beamformer; F_BB: digital beamformer; W: hybrid beamformer

function [F_RF,F_BB,W] = get_initial_hybrid_W(Wd,Nrf)

maxiter = 20;
[F_RF,F_BB,W] = hybrid_decompose(Wd,Nrf,maxiter,[]);

end
