

clc; clearvars; close all;


restoredefaultpath

dir_FRIDA = '../../FRIDA/routines_FRIDA_ver_3.2';
% % dir_FRIDA = '../../EQUILIBRIA/routines_FRIDA_ver_3.1';
% % dir_VI = '../../Integral_MQS_axi/routines_VI_axi_ver_1.02_beta';

% % addpath ../..\Inductance_Coefficient\fun_Inductance_ver_1.01

addpath(genpath(dir_FRIDA))


set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% geometry per evol

% % filename_geo = 'RFX_evol_2';
% % % % prepara_geo_FRIDA_RFX_evol
% % 
% % load(['INPUT_FRIDA_geo_' filename_geo],'meshData', 'meshData_loc')


%% Load scenario data

qq = load('./RFX/36922/data_RFX_shot_36922.mat');


[qq.t_p,ind] = unique(qq.t_p);

qq.im = qq.im(ind,:);
qq.ifs = qq.ifs(ind,:);

figure
subplot(2,1,1)
plot(qq.t_p,qq.im)
xlim(qq.t_p([1 end]))
subplot(2,1,2)
plot(qq.t_p,qq.ifs)
xlim(qq.t_p([1 end]))

run ./CARMA0NL/CarMa0NL_2_AXI_only_shell_MB_beta_drop_50ms/for_FRIDA_evol.m
BETAP0_TIME = repmat(BETAP0_TIME(1),size(BETAP0_TIME));

% % load('./CARMA0NL/BETAP0_TIME_36922_raw');
% % BETAP0_TIME = smoothdata(BETAP0_TIME,'SmoothingFactor',.1);
% % BETAP0_TIME = BETAP0_TIME - BETAP0_TIME(1) + .3;

n_time_new = 200;
t_start_new = 1.02	;
t_end_new = 1.1134;
% % n_time_new = 25;
% % t_start_new = 1.0	;
% % t_end_new = 1.04;


% % n_time_new = 500;
% % t_start_new = .6	;
% % % % dt = 1.202405e-03;
% % dt = 5.01e-4;
% % t_end_new = t_start_new + n_time_new*dt;
% % t_start_new = .6	;
% % t_end_new = 1.2;
% % t_end_new = .605;


time_sim_new = linspace(t_start_new,t_end_new,n_time_new);

imm = double(interp1(qq.t_imm_at,qq.imm_at,time_sim_new));
ifs = double(interp1(qq.t_ifs_at,qq.ifs_at,time_sim_new));
iSC = double(interp1(qq.t_SC_ref,qq.SCCurrentRef_eda1,time_sim_new));
ipla = double(interp1(qq.t_ipla_at,qq.ipla_at,time_sim_new));
% %  = repmat(58240.090000,1,n_time_new);

% % BETAP0_TIME = interp1(time_sim,BETAP0_TIME,time_sim_new);

% % imm = repmat(imm(1,:).',1,n_time_new).';
% % ifs = repmat(ifs(1,:).',1,n_time_new).';
% % iSC = repmat(iSC(1,:).',1,n_time_new).';


% % imm = repmat(imm(1,:).',1,n_time_new).';
% % ifs = repmat(ifs(1,:).',1,n_time_new).';
% % iSC = repmat(iSC(1,:).',1,n_time_new).';
% % ipla = [58240.0952723379 ipla];
% % BETAP0_TIME = [.3 BETAP0_TIME];

% % n_fac = 10;
% % fac_red = .3;
% % fac = [1 1 linspace(1,fac_red,n_fac) linspace(fac_red,fac_red,n_time_new-n_fac-2)];


% % qqqq = ifs(:,3).*fac.';
% % 
% % figure
% % plot(time_sim,ifs(:,3),'k')
% % hold on
% % plot(time_sim, ifs(:,3) + ifs(:,3).*(1-fac.'),'r')
% % 
% % ifs(:,3) = ifs(:,3) + ifs(:,3).*(1-fac.');


% % n_time_new = length(ipla);
% % 
% % time_sim_new = linspace(t_start_new,t_end_new,n_time_new);

time_sim = time_sim_new;
n_time = length(time_sim);

figure
subplot(3,1,1)
plot(time_sim,imm); hold on
plot(time_sim,imm,'o')
xlim(time_sim([1 end]))
subplot(3,1,2)
plot(time_sim,ifs); hold on
plot(time_sim,ifs,'o')
xlim(time_sim([1 end]))
subplot(3,1,3)
plot(qq.t_ipla_at,qq.ipla_at); hold on
plot(time_sim,ipla,'o')
xlim(time_sim([1 end]))



% % if 0
% %     ipla = double(interp1(qq.t_ipla_at,qq.ipla_at,time_sim_new));
% %     
% %     figure
% %     plot(time_sim_new,ipla)
% %     
% %     fac_curr_quench = ipla/max(ipla);
% %     
% %     % % IP_TIME = [IP_TIME(1); IP_TIME(1)*fac_currr_quench.'];
% %     IP_TIME = [IP_TIME(1)*fac_curr_quench.'];
% %     IP_TIME(IP_TIME<0) = 0;
% %     
% %     figure
% %     plot(time_sim_new,IP_TIME)
% %     
% %     t_start_new = 1.053	;
% %     t_end_new = 1.125;
% %     
% %     dt = 2*(t_end_new - t_start_new)/n_time
% %     dt = (t_end_new - t_start_new)/n_time
% %     
% %     fprintf('%f,\n',IP_TIME)
% %     
% %     fac_beta_drop = fac_curr_quench + .5.*(1 - fac_curr_quench);
% %     BETAP0_TIME = .3*fac_beta_drop;
% %     dt = .4*(t_end_new - t_start_new)/n_time
% %     
% %     BETAP0_TIME = [BETAP0_TIME repmat(BETAP0_TIME(end),1,250)]
% %     
% %     fprintf('%f,\n',BETAP0_TIME)
    
    load data_beta_drop 
    n_fac = 25;
    fac = [linspace(1,0.8,n_fac) linspace(0.8,0.8,n_time-n_fac)];
    figure
    plot(time_sim_new,fac)
% % end
% % RADD = 1000;
% % qq_print = (ifs(:,3) - ifs(1,3))*RADD
% % qq_print2 = (ifs(:,3) - ifs(1,3))/0.654528941789845     ;
% % 
% % figure
% % plot(time_sim,ifs(:,3)); hold on
% % plot(time_sim,ifs(1,3) + qq_print/RADD,'--')
% % 
% % 
% % fprintf('%f,',qq_print)


%% Plot for article

% % 
% % figure
% % subplot(3,1,1)
% % plot(qq.t_imm_at,qq.imm_at,'linewidth',2); hold on
% % plot(time_sim,imm,'o')
% % xlim([0 1.2])
% % ylabel('$I$ [A]')
% % title('$I$ OH')
% % set(gca,'fontsize',10)
% % subplot(3,1,2)
% % plot(qq.t_ifs_at,qq.ifs_at,'linewidth',2); hold on
% % plot(time_sim,ifs,'o')
% % xlim([0 1.2])
% % title('$I$ FS')
% % ylabel('$I$ [A]')
% % set(gca,'fontsize',10)
% % subplot(3,1,3)
% % plot(qq.t_ipla_at,qq.ipla_at,'linewidth',2); hold on
% % plot(time_sim,ipla,'o')
% % xlim([0 1.2])
% % title('$I$ plasma')
% % ylabel('$I$ [A]')
% % xlabel('$t$ [s]')
% % set(gca,'fontsize',10)
% % 
% % savefig current_waveforms_rampdown
% % print('current_waveforms_rampdown', '-depsc', '-r500')


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

Conductors_geo = [];
Conductors_index = [];

spacing = zeros(size(Conductors_geo,2)+2,1);
spacing(1) = .03;
spacing(2) = .03;

tic
[meshData_pla] = fun_buildmesh_pla_FRIDA_evol(...
    FW,n_FW,bc,n_bc,factor_bc,spacing,N_order,dir_FRIDA,MEX_OPT);
toc

%%% Coupling surface
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

% % plot(meshData_pla.C_surf(:,1),meshData_pla.C_surf(:,2),'k*')

%%% external sources (active + passive conductors)
load('./RFX/Data/RFXmod_mesh_VI_linear.mat')
% % load('./RFX/Data/RFXmod_noTSS_asCerma0nl_mesh_VI_linear.mat')
% % load('./RFX/Data/RFXmod_noTSS_no_saddle_mesh_VI_linear.mat')
meshData_ext = meshData;

% % warning('neglecting TSS (putting fictitious high resistivity)') 
% % meshData_ext.rho_pas(3) = 10;

% % warning('neglecting VV (putting fictitious high resistivity)') 
% % meshData_ext.rho_pas(1) = 10;

% % warning('neglecting gaps on toroidal direction for PSS') 
% % meshData_ext.ind_passive_cut = 3;

figure
hold on;    grid on

mesh_vac = triplot(meshData_ext.t(:,1:3),meshData_ext.n(:,1),meshData_ext.n(:,2));
set(mesh_vac,'color',0*[1 1 1]);

mesh_vac = triplot(meshData_pla.t(:,1:3),meshData_pla.n(:,1),meshData_pla.n(:,2));
set(mesh_vac,'color','r'); axis equal

addpath('.\CARMA0NL\CARIDDI_2')





%% Magnetic sensors

%%%%% flux-loops
theta_fp =  [22.5 67.5 112.5 157.5 202.5 247.5 292.5 337.5];
% % theta_fp =  linspace(0,360,100); theta_fp = theta_fp(1:end-1);
thetarad = theta_fp*pi/180;
r_fluxloop = 0.5065;
PSISENS_R = [1.995 + r_fluxloop*cos(thetarad)]';
PSISENS_Z = [r_fluxloop*sin(thetarad)]';

flux_loops = [PSISENS_R PSISENS_Z];

SETTINGS.flux_loops = flux_loops;



%%%%% pick-up
theta = [27,72,117,162,207,252,297,342];
% % theta =  linspace(0,360,100); theta = theta(1:end-1);
thetarad=theta*pi/180;
r_pickup = 0.5085;
BSENS_R = [1.995 + r_pickup*cos(thetarad)]';
BSENS_Z = [r_pickup*sin(thetarad)]';
BSENS_T = [-sin(thetarad)' cos(thetarad)'];

pickup = [BSENS_R BSENS_Z BSENS_T];

SETTINGS.pickup = pickup;

phi_0 = 70*pi/180;
phi = linspace(0,2*pi,25) - phi_0;
phi = phi(1:end-1);

XFLD_mat = BSENS_R*cos(phi);
YFLD_mat = BSENS_R*sin(phi);
ZFLD_mat = repmat(BSENS_Z,1,length(phi));

XFLD = XFLD_mat(:);
YFLD = YFLD_mat(:);
ZFLD = ZFLD_mat(:);

figure
plot_mesh_CARIDDI('RFX_active_passive_renumber.msh',3,[.5 .5 .5])
% % plot_mesh_CARIDDI('RFX_active_passive_renumber.msh')
hold on
plot3(XFLD(1:8),YFLD(1:8),ZFLD(1:8),'bo', 'LineWidth',2)

fprintf('%f, ',XFLD.')

meshData_ext.n(:,1)
set(mesh_vac,'color',0*[1 1 1]);


%% Currents

load('./RFX/Data/Nturns_RFX_64')
Currents_evol = zeros(64,n_time);
KONNAX_ACT = meshData_ext.KONNAX_ACT;

Currents_evol = zeros(64,n_time);
for ii=1:n_time
    
    Currents = double([imm(ii,:) ifs(ii,:) iSC(ii,[1 2])]');
% %     warning('imposing zero current on saddle coils')
    % Active Currents
    Currents_all = KONNAX_ACT'*Currents;
    
    Currents_all(48) = 0;
    Currents_all(50) = 0;
    
    % %     Currents_evol(:,ii) = Currents_all(1:56).*Nturns_RFX(1:56);
    Currents_evol(:,ii) = Currents_all.*Nturns_RFX;
end

% % load data_curr_check_carma0nl.mat
% % 
% % I_act = [I_act repmat(I_act(:,end),1,length(time_sim)-size(I_act,2))];
% % 
% % Currents_evol = Currents_evol;
% % Currents_evol([45 46],:) = Currents_evol([45 46],1) - I_act([23 51],:);

Conductors.Nconductors = size(Currents_all,1);
Conductors.Currents = Currents_evol(:,1);
Conductors.Nturns = ones(size(Nturns_RFX));

% % openfig aaihwevf.fig
% % load Gmat.mat
% % dCurrent = Currents_evol - Currents_evol(:,1);
% % for ii = 1:8
% %    subplot(2,4,ii)
% %    plot(time_sim,Gmat(ii,:)*dCurrent,'c--') 
% % end


%% IE_evol entries

ref_IC = load('Ipasma_36922_VI_new.mat');

J_c_0 = double(interp1(ref_IC.time_sim, ref_IC.J_VI_rr.', time_sim(1))).';
v_gap_t0 = double(interp1(ref_IC.time_sim, ref_IC.xx_VI(end-1:end,:).', time_sim(1))).';

% %     IE_1 = load(['./INPUT_FRIDA_equil_RFX_36922']);
% %     Iplasma_ref = load('Ipasma_36922_VI_RFXmod2.mat');

ALFAMF_TIME= 1.15.*1.2846158822057*ones(1,n_time);
ALFANF_TIME= 0.998531120938779*ones(1,n_time);
BETAP0_TIME= .9*0.3000*ones(1,n_time).*fac;
IE_evol.plaparameter.Ipla  = 1.0625.*ipla;

% % load('prova.mat')
% % ALFAMF_TIME= v_bestDes(1)*ones(1,n_time);
% % ALFANF_TIME= v_bestDes(2)*ones(1,n_time);
% % BETAP0_TIME= v_bestDes(3)*ones(1,n_time).*fac;
% % IE_evol.plaparameter.Ipla  = v_bestDes(4)*ipla/ipla(1);

% % ALFAMF_TIME= 1.2846158822057*ones(1,n_time);
% % ALFANF_TIME= 0.998531120938779*ones(1,n_time);
% % BETAP0_TIME= .0*ones(1,n_time);
% % IE_evol.plaparameter.Ipla  = double(interp1(ref_IC.time_sim, ref_IC.Ipla_tot_rec.', time_sim)).';


IE_evol.plaparameter.Centroid = vec_linspace([2 0]',[2 0]',n_time);
IE_evol.plaparameter.R_0      = linspace(1.995,1.995,n_time);

IE_evol.plaparameter.beta_0   = BETAP0_TIME;
IE_evol.plaparameter.alpha_M  = ALFAMF_TIME;
IE_evol.plaparameter.alpha_N  = ALFANF_TIME;



% %     BETAP0_TIME = 0.300000000000000*ones(1,n_time);

% % load('./CARMA0NL/BETAP0_TIME_36922_raw');
% % 
% % % % BETAP0_TIME = BETAP0_TIME_smooth;
% % % % BETAP0_TIME = smoothdata(BETAP0_TIME,'loess',20);
% % BETAP0_TIME2 = smoothdata(BETAP0_TIME,'SmoothingFactor',.1);
% % BETAP0_TIME2 = BETAP0_TIME2 - BETAP0_TIME2(1) + .3;
% % 
% % % % figure
% % % % plot(linspace(.6,1.2,900),BETAP0_TIME); hold on
% % % % plot(linspace(.6,1.2,900),BETAP0_TIME2);
% % 
% % BETAP0_TIME = interp1(linspace(.6,1.2,900),BETAP0_TIME2,time_sim);
% % 
% % figure
% % plot(time_sim,BETAP0_TIME,'k'); hold on
% % % %     plot(time_sim,fac_VDE.*BETAP0_TIME,'ro'); hold on
% % % % plot(time_sim,smoothdata(BETAP0_TIME,'loess',10),'r', 'LineWidth',2)
% % 
% % % %     BETAP0_TIME = fac_VDE.*BETAP0_TIME;
% % 
% % IE_evol.plaparameter.beta_0 = BETAP0_TIME.';
% % IE_evol.plaparameter.alpha_M = ALFAMF_TIME.';
% % IE_evol.plaparameter.alpha_N = ALFANF_TIME.';
% % 
% % 
% % fac_Ipla = 1;
% % IE_evol.plaparameter.Ipla = double(interp1(Iplasma_ref.time_sim, Iplasma_ref.Ipla_tot_rec, time_sim));
% % fac_Ipla = 5.6628e+4/IE_evol.plaparameter.Ipla(1);
% % IE_evol.plaparameter.Ipla = fac_Ipla*IE_evol.plaparameter.Ipla;
% % IE_evol.plaparameter.Ipla = fac_VDE.*IE_evol.plaparameter.Ipla;



% % IE_evol.plaparameter.Ipla = IP_TIME.';



IE_evol.Conductors.Nconductors  = Conductors.Nconductors;
IE_evol.Conductors.Nturns       = Conductors.Nturns;
IE_evol.Conductors.Currents     = Currents_evol;
IE_evol.time_sim                = time_sim;


%%

% % figure
% % plot(time_sim,IE_evol.Conductors.Currents); hold on
% % plot(time_sim(322),IE_evol.Conductors.Currents(:,322),'o', 'LineWidth',2)
% % 
% % figure
% % plot(time_sim,IE_evol.plaparameter.Ipla); hold on
% % plot(time_sim(322),IE_evol.plaparameter.Ipla(322),'o', 'LineWidth',2)


%% run FRIDA evol
filename_out = 'temp_out';
filename_preproc = 'RFXmod_36922_rampdown';
filename_save = 'temp_save';

filename_geo = 'RFXmod_36922_rampdown';
save(['INPUT_FRIDA_evol_VI_geo_' filename_geo],'meshData_pla', 'meshData_ext')

filename_evol = 'RFXmod_36922_rampdown';
save(['INPUT_FRIDA_evol_VI_' filename_evol],'IE_evol')


SETTINGS.filename_geo = filename_geo;
SETTINGS.filename_evol = filename_evol;
SETTINGS.filename_preproc = filename_preproc;

SETTINGS.filename_save = filename_save;
SETTINGS.filename_out = filename_out;

SETTINGS.FIRST_SOLUTION = false;
SETTINGS.SAVE_OUTPUT = false;
SETTINGS.SOLVER = 'NR';
SETTINGS.RUN_MEX_ROUTINE = true;

SETTINGS.VAC_METHOD = 2;

SETTINGS.plaparameter.Centroid = [2 -.05];
SETTINGS.FIGURES_DEBUG = 0;

SETTINGS.J_c_t0 = J_c_0;
SETTINGS.v_gap_t0 = v_gap_t0;

SETTINGS.QuasiNewton = true;
SETTINGS.QuasiNewton_factor = 0.5;

SETTINGS.ii_START = 1   ;
SETTINGS.PREPROC = true;
SETTINGS.FIGURES = true;
SETTINGS.IS_EVOL = 1;
close all
[OUT_FRIDA_TD] = main_run_FRIDA_evol_VI(SETTINGS);

save('OUT_FRIDA_TD_RFXmod_36922_rampdown_11.mat', ...
    'OUT_FRIDA_TD' ,...
    'IE_evol', ...
    'meshData_pla', ...
    'meshData_ext');
 

%%

if 0
    
    load data_debug.mat

SETTINGS.SOLVER_RELAX_TRESHOLD_1 = 20
SETTINGS.ii_freeB_max = 25;
    SETTINGS.FIGURES = true;
    [OUT_th] = run_FRIDA_TD_CN(meshData_pla,... % CN = Crank Nicolson
        meshData_ext,...
        plaparameter_th,...
        Conductors_th,...
        solk_old, ...
        SETTINGS);
    
end


%%
[RR_grid,ZZ_grid] = meshgrid(linspace(1.4,2.6,100),linspace(-.6,.6,100)); 

n_time_max = numel(find(OUT_FRIDA_TD.Ax_r_TD));

close all
nist = 1:10:n_time_max;
for ii_time =nist
    
    figure(ii_time)
    axis equal; hold on; colormap jet,
    % %     contour(RR_grid,ZZ_grid,griddata(meshData_pla.n(:,1),meshData_pla.n(:,2), ...
    % %         OUT_FRIDA_TD.Psi_TD(:,ii_time)/pi/2,RR_grid,ZZ_grid),OUT_FRIDA_TD.Psi_B_TD(ii_time)*[1 1]/2/pi,'k-', 'LineWidth',3)
    contour(RR_grid,ZZ_grid,griddata(meshData_pla.n(:,1),meshData_pla.n(:,2), ...
        OUT_FRIDA_TD.Psi_TD(:,ii_time)/pi/2,RR_grid,ZZ_grid), ...
        linspace(OUT_FRIDA_TD.Psi_B_TD(ii_time),OUT_FRIDA_TD.Psi_Ax_TD(ii_time), 25)/pi/2)
    colormap jet
    title(['t = ' num2str(time_sim(ii_time)) 's'])
% %     pause
end

showfilm(nist,0.1);


figure
plot(OUT_FRIDA_TD.Separatrix_r_TD(:,1:10:n_time_max),OUT_FRIDA_TD.Separatrix_z_TD(:,1:10:n_time_max))
hold on; axis equal
print('LCFS_evol', '-dpng', '-r500')

% % print('LCFS_evol_volt_kick', '-dpng', '-r500')







