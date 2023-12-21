function Processing(Geome,Mater,totint,coord,disp,dmg,tt,vel)
% Analytical Solution

t = Geome.dt:Geome.dt:Geome.nt*Geome.dt;

colormap jet;
subplot(222)
scale = 150;
scatter(coord(:, 1)+scale*disp(:, 1), coord(:, 2)+scale*disp(:, 2), [], dmg, "filled")
colorbar
axis equal

subplot(223)
scatter(coord(:, 1)+scale*disp(:, 1), coord(:, 2)+scale*disp(:, 2), [], disp(:,1), "filled")
colorbar
axis equal

