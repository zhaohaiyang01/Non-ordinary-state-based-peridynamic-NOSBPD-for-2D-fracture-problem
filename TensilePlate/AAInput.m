function [Geome,Mater] = AAInput()
%% Geometry information
Geome.length = 0.05;
Geome.width = 0.05;
Geome.thick = 1;
Geome.ndivx = 150;
Geome.ndivy = 150;
Geome.ndivz = 1;
Geome.dof = 2;
Geome.nbnd = 3; 
Geome.nt = 2000;                    
% nt: Total number of time step
Geome.dt = 1;
Geome.ADR = 1;
Geome.uinc = [0,1e-8,0];                  
% Displacement BC,uinc:velcoty
% Geome.P = 10e6;                       
% P:Force
Geome.prec = 1;    
Geome.TipL = [-0.00501,0];
Geome.TipR = [0.00501,0];
% Prec:Whether there is a pre-crack 1 - exist,0 - non-exist
% *** The precrack process applies only to 2D problems ***
Geome.dx = Geome.length/Geome.ndivx;
Geome.radij = Geome.dx/2;
Geome.area = Geome.dx^2;
Geome.delta = 2.515 * Geome.dx;
%Geome.delta = 0.0015075;
Geome.deltak = 3.015*Geome.dx;
Geome.ADR = 1; % 1:Using ADR,0:Forward Euler
%% Material information
Mater.dens = 2400;%8000.0;                   
% dens: Density
Mater.emod = 30e9;%192.0e9;                  
% emod: Elastic modulus
Mater.pratio = 0.1;%1/3;              
% pratio: Poisson's ratio

if Geome.dof ==  2
    Geome.vol = Geome.area*Geome.thick;
    Mater.state = 'Plane stress';
end

Mater.sc = 0.0012;%0.001
