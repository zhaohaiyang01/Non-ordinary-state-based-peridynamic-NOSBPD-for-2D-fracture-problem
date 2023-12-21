function Ts = LiState(wf, X, Y, x, bc, Z, Z_, pv, ih, jh, iK,dof)
% Li State
% Refer to
% Li, P., Hao, Z.M., Zhen, W.Q.: 
% A stabilized non-ordinary state-based peridynamic model. 
% Comput. Methods Appl. Mech. Eng. 339, 262–280 (2018)
wf0 = bc*1/2*wf./x.^3;
C = matmul(X.*[wf0,wf0],X,dof,1,dof); % 0.5*w*bc
KN = matmul(C,matmul(Z,Y-X,dof,1,dof),dof,dof,dof);
SS = zeros(size(pv,1),dof^2);
for i = 1:1:size(KN,2)
    SS(:,i) = accumarray(ih,KN(:,i).*pv(jh));
end
    SSiK = matmul(SS, iK, dof,dof,dof);
    ST = SSiK(ih,:) + SSiK(jh,:);
    Ts = matmul(C, Z-Z_, dof, dof, 1) + matmul(ST,X.*[wf,wf],dof,dof,1);
end