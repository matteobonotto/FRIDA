%% MATLAB 2 GMSH x RFX MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Write ITER geometrical data write in a suitable format for GMSH mesh
%   generation.
%   
%   Matteo Bonotto 01/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
clear;  clc;  close all;
restoredefaultpath

% % addpath ../../../source/routines_FRIDA/AutoMESH_VI_ver_1.0/
dir_FRIDA = '../../../source/routines_FRIDA/';
addpath(genpath(dir_FRIDA))

input_Data = '../data/geo/';


%%

%% geometry per evol (VI)

%%% Plasma domain (PDE solution)

%RFX-mod2
rfw = 0.49; 
r_PSS_int = 0.5115;
rbc = .5*(r_PSS_int + rfw);

%RFX-mod
rfw = 0.459; 
r_VV_int = 0.475;
rbc = .5*(r_VV_int + rfw);



nn_FW = 120;

FW = [1.995+rfw*cos(linspace(0,2*pi-2*pi/nn_FW,nn_FW))' ...
    rfw*sin(linspace(0,2*pi-2*pi/nn_FW,nn_FW))'];

bc=[1.995+rbc*cos(linspace(0,2*pi-2*pi/nn_FW,nn_FW))' ...
    rbc*sin(linspace(0,2*pi-2*pi/nn_FW,nn_FW))'];

n_FW = [];
n_bc = nn_FW;
factor_bc = [];
N_order = 2;
MEX_OPT = true;

spacing = zeros(2,1);
spacing(1) = .03;
spacing(2) = .03;

dir_AutoMESH = [dir_FRIDA, '/AutoMESH_ver_2/'];
tic
[meshData_pla] = fun_buildmesh_pla_FRIDA_evol(...
    FW,n_FW,bc,n_bc,factor_bc,spacing,N_order,dir_AutoMESH,MEX_OPT);
toc

%% Define a coupling surface (here we are taking the BCs curve as CS)
tri = meshData_pla.t;
ind_n_B = meshData_pla.ind_n_bc;
[ind_t_B] = fun_indTriOnEdge(tri,ind_n_B);
ind_surf = zeros(size(ind_t_B));

for ii = 1 : length(ind_t_B)
    tri_ii = meshData_pla.t(ind_t_B(ii),:);
    tmp = intersect(tri_ii,meshData_pla.ind_n_bc);
    ind_surf(ii) = tmp(end);
end

C_surf = meshData_pla.n(ind_surf,:);
meshData_pla.C_surf = C_surf;

plot(meshData_pla.C_surf(:,1), meshData_pla.C_surf(:,2), '*b')
legend('vac + plasma mesh', ...
    'plasma mesh', ...
    'limiter curve', ...
    'BCs curve', ...
    'normal unit vec',...
    'coupling surface')



%% Save mesh
save([input_Data, 'RFXmod_mesh_pla.mat'], 'meshData_pla')







