clear
clc
close all
format long
addpath('Common Files\')
% Input information
[Geome,Mater] = AAInput();
dt = Geome.dt;
dx = Geome.dx;
dof = Geome.dof;
nt = Geome.nt;
sc = Mater.sc;
emod = Mater.emod;
delta = Geome.delta;

% Generate material points
[coord,totint] = Gen_nodes(Geome);

% Find displacement and force boundary
[coord,totnode,fnodes,inrf,dire,inbu] = AAFindBoundary(Geome,coord,totint);

% Visualize nodes and boundaries
plot_node(coord,totint,Geome,fnodes)

[allPointMember,numNeiPoint,pointfamtot,ih,fail] = two_neighboornodes(coord,Geome);
pv = dx^2.*ones(max(ih),1);

% Computing initial relative position state and shape tensor 
X = coord(allPointMember,:) - coord(ih,:);        % Relative position state at bond
x = vecnorm(X,2,2);
wf = exp(-x.^2/delta^2);                          % Gauss weight function
wf = wf.*fail;

K = innerPro(wf,X,X,pv(allPointMember),ih,totnode,dof);   % Initial shape tensor of bond
K_matrix = reshape(K',dof,dof,[]);          % Matrix form of shape tensor
K_inv = pageinv(K_matrix);
K_inv = pagetranspose(K_inv);
K_inv = reshape(K_inv,dof^2,[])';           % Inverse of shape tensor
[Dmat] = D_elastic(Mater,Geome);

% Initialization
bc = 12.0d0 * emod / (pi * (delta^4));
massvec = 0.25d0 * dt * dt * ((4.0d0/3.0d0)*pi*(delta)^3) * bc / dx*8; % ADR parameters
bforce = zeros(size(coord));
% bforce(fnodes(1:inrf),dire) = 10e6/dx;
% bforce(fnodes(inrf+1:end),dire) = -10e6/dx;
velhalfold = zeros(size(coord));
pforceold = zeros(size(coord));
disp = zeros(size(coord));
dmg = zeros(totnode,1);
withZeroEngyMode = 0; % 0-without zero energy mode,1-with zero energy mode

for tt = 1:nt
    fprintf("%d/%d\n",tt,nt);
    wf = wf.*fail;

    % Application of boundary conditions
    if isfield(Geome,'uinc') && (sum(Geome.uinc) ~= 0)
        disp(totint+1:totint+inbu,:) = -ones(inbu,1)*Geome.uinc(1:dof).*tt*dt;
        disp(totint+inbu+1:totnode,:) = ones(totnode-inbu-totint,1)*Geome.uinc(1:dof).*tt*dt;
        vel(totint+1:totint+inbu,:) = -ones(inbu,1)*Geome.uinc(1:dof);
        vel(totint+inbu+1:totnode,:) = ones(totnode-inbu-totint,1)*Geome.uinc(1:dof);
    end

    K = innerPro(wf,X,X,pv(allPointMember),ih,totnode,dof);   % Initial shape tensor of bond
    K_matrix = reshape(K',dof,dof,[]);          % Matrix form of shape tensor
    K_inv = pageinv(K_matrix);
    K_inv = pagetranspose(K_inv);
    K_inv = reshape(K_inv,dof^2,[])';           % Inverse of shape tensor

    condK = abs(K(:,1))+abs(K(:,2))+abs(K(:,3))+abs(K(:,4));
    condKinv = abs(K_inv(:,1))+abs(K_inv(:,2))+abs(K_inv(:,3))+abs(K_inv(:,4));
    rcondK = 1./(condK.*condKinv);
    rcondK(isnan(rcondK)) = 0;
    index1 = find(rcondK<1e-7);
    K_inv(index1,:) = [1,0,0,1].*ones(size(index1,1),1);  

    acoord = coord+disp;
    Y = acoord(allPointMember,:)-acoord(ih,:); % Deformation state
    K = innerPro(wf,Y,X,pv(allPointMember),ih,totnode,dof);            % Shape tensor after deformation
    F = matmul(K, K_inv, dof,dof,dof);            % Deformation gradient tensor

    condF = abs(F(:,1))+abs(F(:,2))+abs(F(:,3))+abs(F(:,4));
    detF = F(:,1).*F(:,4)-F(:,2).*F(:,3);
    Finv = [F(:,4),-F(:,2),-F(:,3),F(:,1)]./detF;
    condFinv = abs(Finv(:,1))+abs(Finv(:,2))+abs(Finv(:,3))+abs(Finv(:,4));
    rcondF = 1./(condF.*condFinv);
    rcondF(isnan(rcondF)) = 0;
    index1 = find(rcondF<1e-7);
    F(index1,:) = [1,0,0,1].*ones(size(index1,1),1);  

    epsilon = GreenStrainTensor(F,dof);
    S = (Dmat*epsilon')'; % Cauchy stress tensor

    F_matrix = reshape(F',dof,dof,[]);          % Matrix form of shape tensor
    FT = pagetranspose(F_matrix);
    FT_inv = pageinv(FT);
    detF = F(:,1).*F(:,4)-F(:,2).*F(:,3);
    S_matrix = reshape(S',dof,dof,[]);
    PK1_matrix = pagemtimes(S_matrix,FT_inv);
    PK1 = pagetranspose(PK1_matrix);
    S1 = detF.*reshape(PK1,dof^2,[])';

    SiK = matmul(S1,K_inv,dof,dof,dof);
    ST = SiK(ih,:) + SiK(allPointMember,:);% Sij*(xj-xi)-Sji*(xi-xj)=(Sij+Sji)*(xj-xi)=ST*X
    Tbfore = matmul(ST,X.*(wf.*ones(1,dof)),dof,dof,1); % 原先力态
    if(withZeroEngyMode==1) % 如果需要有零能模式存在，则不附加力态
        T =  Tbfore;
    else
        Z = Y - matmul(F(ih,:),X,dof,dof,1);          % Non-uniform deformation gradient tensor Z
        Z_ = -Y - matmul(F(allPointMember,:),-X,dof,dof,1);    % Non-uniform deformation gradient tensor Z'
        %Ts = SillingState(10e3*emod/dx, 1, Z, Z_, wf, pv, ih, allPointMember,dof); % Silling
        %Ts = LiState(wf, X, Y, x,18*emod/2/(1-Mater.pratio-2*Mater.pratio^2)/pi/delta^4, Z, Z_, pv,ih, allPointMember, K_inv,dof); %Li Pan
        Ts  = WanState(wf, Dmat, K_inv, ih, allPointMember, Z, Z_);
        T =  Tbfore+Ts;
    end
    pforce = zeros(size(coord));
    for m = 1:dof
        pforce(:,m) =  accumarray( ih, T(:,m).*pv(allPointMember,1));
    end

    if Geome.ADR == 1
        [vel,disp,velhalfold,pforceold] = ADR(velhalfold,disp,pforce,massvec,pforceold,dt,tt,bforce);
    else
        acc = (pforce + bforce)./dens;
        vel = vel+acc.*dt;
        disp = disp+vel.*dt;
    end

    if isfield(Geome,'uinc') && (sum(Geome.uinc) ~= 0)
        disp(totint+1:totint+inbu,:) = -ones(inbu,1)*Geome.uinc(1:dof).*tt*dt;
        disp(totint+inbu+1:totnode,:) = ones(totnode-inbu-totint,1)*Geome.uinc(1:dof).*tt*dt;
        vel(totint+1:totint+inbu,:) = -ones(inbu,1)*Geome.uinc(1:dof);
        vel(totint+inbu+1:totnode,:) = ones(totnode-inbu-totint,1)*Geome.uinc(1:dof);
    end
    if anynan(disp) || max(disp(:,1)) > 0.000011
         c = 1;
    end
    deformX = coord(allPointMember,:) - coord(ih,:)+disp(allPointMember,:) - disp(ih,:);
    deformx = vecnorm(deformX,2,2);
    s = (abs(deformx)-abs(x))./abs(x);
    index = s > Mater.sc;
    fail(index,1) = 0.0001;
    conbond = accumarray(ih,fail);
    a = fail>0;
    a = fail(a);
    if max(index) > 0 && size(a,1) < 1137584
        b = 1;
    end
    dmg = 1-accumarray(ih,fail.*pv(allPointMember,1))./accumarray(ih,pv(allPointMember,1));

if tt == 500
    d = 1;
end
    if mod(tt,20) == 0 %|| tt >22%|| tt < 101
        Processing(Geome,Mater,totint,coord,disp,dmg,tt)
        drawnow
    end
end





