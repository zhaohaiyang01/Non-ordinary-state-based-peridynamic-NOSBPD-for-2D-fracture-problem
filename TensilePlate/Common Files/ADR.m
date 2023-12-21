function [vel,disp,velhalfold,pforceold] = ADR(velhalfold,disp,pforce,massvec,pforceold,dt,tt,bforce)
    % Adaptive dynamic relaxation ⬇⬇⬇
    q = abs(velhalfold)>1e-12;
    cn1 = - sum(sum(disp(q) .* disp(q) .* (pforce(q) / massvec - pforceold(q) / massvec) ./ (dt * velhalfold(q))));
    cn2 = sum(sum(disp.*disp));
    if(abs(cn2)>1e-12)
        if ((cn1 / cn2) > 0.0)
            cn = 2.0 * sqrt(cn1 / cn2);
        else
            cn = 0.0;
        end
    else
        cn = 0.0;
    end
    if (cn > 2.0d0)
        cn = 1.9d0;
    end
    
    if (tt==1)
        velhalf = 1.0d0 * dt / massvec * (pforce + bforce) / 2.0;
    else
        velhalf = ((2.0d0 - cn * dt) * velhalfold + 2.0 * dt / massvec * (pforce + bforce)) / (2.0 + cn * dt);
    end
    vel = 0.5d0 * (velhalfold + velhalf);
    disp = disp + velhalf * dt;
    velhalfold = velhalf;
    pforceold = pforce;
    % Adaptive dynamic relaxation ⬆⬆⬆

