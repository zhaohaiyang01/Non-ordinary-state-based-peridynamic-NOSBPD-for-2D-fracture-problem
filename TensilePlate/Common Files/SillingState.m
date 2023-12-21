function Ts = SillingState(C, G, Z, Z_, wf, pv, ih, jh,dof)
% Silling State
% Refer to
% Silling, S.A.: 
% Stability of peridynamic correspondence material models and their particle discretizations. 
% Comput. Methods Appl. Mech. Eng. 322, 42–57 (2017)
OMG0 = accumarray(ih, wf.*pv(jh));
CGOI = C*G*OMG0(jh)./OMG0(ih);
CGOJ = C*G*OMG0(ih)./OMG0(jh);
Ts = CGOI*ones(1,dof).*Z - CGOJ*ones(1,dof).*Z_;
end