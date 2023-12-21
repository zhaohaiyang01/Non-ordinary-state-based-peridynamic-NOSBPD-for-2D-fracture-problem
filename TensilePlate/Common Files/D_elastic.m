function [Dmat]=D_elastic(Mater,Geome)
emod = Mater.emod;            %  Young's modulus
pratio = Mater.pratio;          % Poison ratio
dof = Geome.dof;
if dof == 2
    state = Mater.state;
    switch state
        case 'Plane stress'
            Dmat = [1         0              0       pratio; ...
                    0   (1-pratio)/2   (1-pratio)/2     0;...
                    0   (1-pratio)/2   (1-pratio)/2     0;
                    pratio       0              0          1]*emod/(1-pratio^2);
        case 'Plane strain'
            fac = emod*(1-pratio)/(1+pratio)/(1-2*pratio);
            a = pratio/(1-pratio);
            b = (1-2*pratio)/2/(1-pratio);
            Dmat = [1 0 0 a;
                    0 b b 0;
                    0 b b 0;
                    a 0 0 1]*fac;
    end
elseif dof == 3
    miu = emod/2/(1+pratio);
    namda = pratio*emod/((1+pratio)*(1-2*pratio));
    a = namda+2*miu;
    miu2 = 2*miu;
    Dmat = [ a,   0,   0,   0,  namda,  0,  0,   0, namda;
             0,  miu2, 0,   0,    0,    0,  0,   0,   0;
             0,   0,  miu2, 0,    0,    0,  0,   0,   0;
             0,   0,   0,  miu2,  0,    0,  0,   0,   0;
           namda, 0,   0,   0,    a,    0,  0,   0, namda;
             0,   0,   0,   0,    0,  miu2, 0,   0,   0;
             0,   0,   0,   0,    0,    0, miu2, 0,   0;
             0,   0,   0,   0,    0,    0,  0,  miu2, 0;
           namda, 0,   0,   0,  namda,  0,  0,   0,   a];
end