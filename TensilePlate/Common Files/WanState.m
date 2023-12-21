function Ts = WanState(wf, De, iK, ih, jh, Z, Z_)
% WanState
% Refer to 
% Wan, J., Chen, Z., Chu, X. et al. 
% Improved method for zero-energy mode suppression in peridynamic correspondence model. 
% Acta Mech. Sin. 35, 1021–1032 (2019). 
% https://doi.org/10.1007/s10409-019-00873-y
CK = (De*iK')'; %
Ts = matmul(CK(ih,:),Z,2,2,1) - matmul(CK(jh,:),Z_,2,2,1); % 稳定附加力态
Ts = Ts.*[wf, wf];
end