function [solk1]=fun_NR_FRIDA_TD_CN(INPUT_solver)

%%
mu0=4*pi*1.e-7;

%%

solk         = INPUT_solver.solk;
meshData     = INPUT_solver.meshData;
meshData_loc = INPUT_solver.meshData_loc;
plaparameter = INPUT_solver.plaparameter;
ii_freeB     = INPUT_solver.ii_freeB;
SETTINGS     = INPUT_solver.SETTINGS;   
Conductors   = INPUT_solver.Conductors;

if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    Centroid = plaparameter.Centroid;
    Ipla =     plaparameter.Ipla;
    psibar =   plaparameter.psibar;
    FdF =      plaparameter.FdF;
    dP =      plaparameter.dP;
    
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    R_0 =      plaparameter.R_0;
    beta_0 =   plaparameter.beta_0;
    alpha_M =  plaparameter.alpha_M;
    alpha_N =  plaparameter.alpha_N;
    Ipla =     plaparameter.Ipla;
end


%%
nn              = meshData_loc.nn;
nn_D            = length(meshData_loc.ind_D);
nn_B            = length(meshData_loc.ind_B);
n_pas           = meshData.n_pas;
nt              = meshData_loc.nt;
N_order         = meshData_loc.shape_functions_N_order;
ww_Gauss        = meshData_loc.ww_Gauss_pla;
P_Gauss         = meshData_loc.nodes_Gauss_pla;
shape_functions = meshData_loc.shape_functions;
Iphi            = zeros(nn,1);
ind_D           = meshData_loc.ind_D;
ind_B           = meshData_loc.ind_B;
KK_nobc         = meshData_loc.KK_nobc;
KK_D            = KK_nobc(ind_D,ind_D);
KK_bc           = KK_nobc(ind_D,ind_B);
GG_BEM          = meshData_loc.GG_BEM;
nodes_Gauss_pla = meshData_loc.nodes_Gauss_pla;
ww_Gauss_pla    = meshData_loc.ww_Gauss_pla;
rr              = meshData_loc.n(:,1);
rr_Gauss        = meshData_loc.nodes_Gauss_pla(:,1);
ind_n_INpla     = meshData_loc.ind_n_INpla;
% % meshData_loc.KK_D_tilde = -(GG_BEM.psi/mu0)*KK_D;
% % meshData_loc.KK_bc_tilde = (speye(nn_B)-(GG_BEM.psi/mu0)*KK_bc);


%%
[~,ind] = max(solk.Psi(solk.ind_n_axis));
solk.Psi_axis_node = solk.Psi(solk.ind_n_axis(ind));


%% Fixed boundary step to converge on current density

% % if solk.residuo_norm >= SETTINGS.THRESHOLD_FIXB && ii_freeB <= SETTINGS.NFREEB_FIXB
% %     
% %     fprintf('Fixed-Boundary step \n')
% %     
% %     input_FixB.plaparameter = plaparameter;
% %     input_FixB.meshData_loc = meshData_loc;
% %     input_FixB.solk         = solk;
% %     input_FixB.solk.xx_k    = solk.xx_k(1:end-n_pas);
% %     input_FixB.SETTINGS     = SETTINGS;
% %     
% %     tic
% %     [solk_fix] = fun_FixB_FRIDA(input_FixB);
% %     toc
% %     
% %     solk.Psi = solk_fix.Psi;
% %     solk.Psi_axis = solk_fix.Psi_axis;
% %     solk.ind_n_axis = solk_fix.ind_n_axis;
% %     solk.delta_a_weights = solk_fix.delta_a_weights;
% %     solk.lambda = solk_fix.lambda;
% %     solk.lambda_star = mu0*solk_fix.lambda;
% %     
% %     [~,ind] = max(solk.Psi(solk.ind_n_axis));
% %     % % solk.ind_n_axis = solk.ind_n_axis(ind);
% %     solk.Psi_axis_node = solk.Psi(solk.ind_n_axis(ind));
% %     % % solk.delta_a_weights = 1;
% %     
% %     clear solk_fix
% %     
% % else
    
% %     fprintf('Skipping Fixed-Boundary step \n')
    
% % end

    Rplot=[.9*min(meshData_loc.n(:,1)) 1.1*max(meshData_loc.n(:,1))];
    Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
    rgrid=linspace(Rplot(1),Rplot(2),200);
    zgrid=linspace(Zplot(1),Zplot(2),200);
    [RR,ZZ] = meshgrid(rgrid,zgrid);
    PSI_2 = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk.Psi,RR,ZZ);
    % %
    
    if SETTINGS.FIGURES_DEBUG == 1
        figure
        edges=pdemesh(meshData.n',meshData.e_interface',[]);
        set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
        axis([Rplot Zplot]),  hold on,
        contour(RR,ZZ,PSI_2,100);
        contour(RR,ZZ,PSI_2,[solk.Psi_Bpla solk.Psi_Bpla],'r','LineWidth',2);
        colormap('cool'); colorbar vert;
        plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k.')
    end
    
lambda_k = solk.lambda;
lambda_star_k = solk.lambda_star;


%%


%%
tri             = meshData_loc.t;
N_order         = meshData_loc.shape_functions_N_order;
P_Gauss         = meshData_loc.nodes_Gauss_pla;
n_Gauss         = meshData_loc.n_Gauss;
shape_functions = meshData_loc.shape_functions;

Psi_k=solk.Psi;
Psi_axis=solk.Psi_axis_node;
Psi_Bpla=sum(solk.Psi_Bpla)./numel(solk.Psi_Bpla);
psi_bar=(Psi_k-Psi_axis)/(Psi_Bpla-Psi_axis); %normalized poloidal magnetic flux

if SETTINGS.RUN_MEX_ROUTINE == false
    [psi_bar] = fun_evaluateFluxGaussPoints_v2(tri,psi_bar,N_order,P_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [psi_bar] = fun_evaluateFluxGaussPoints_v2_mex(tri,psi_bar,N_order,P_Gauss,n_Gauss,shape_functions);
end
solk.psi_bar = psi_bar;

if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    
% %     psi_bar(psi_bar>1) = NaN;
    psi_bar(psi_bar<0) = 0;

    fdfn=interp1(psibar,FdF,psi_bar,'linear');   
    dpn=interp1(psibar,dP,psi_bar,'linear');   
    gg_tot=(rr_Gauss.*dpn+fdfn./(mu0*rr_Gauss)); %computation of nodes plasma current density
    gg_tot=gg_tot/2/pi;
    
    qq = isnan(gg_tot);
    gg_tot(isnan(gg_tot)) = 0;

    gg_tot_in = gg_tot;
    gg_tot_in(meshData_loc.ind_G_out) = 0;
    
    nn = meshData_loc.nn;
    ind_t_selected = meshData_loc.ind_t_INpla;
    ww_Gauss = meshData_loc.ww_Gauss_pla;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Source] = fun_calcSourceFEM_v2(tri,nn,gg_tot_in,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
        [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg_tot_in,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    end
    
    % %     figure;  plot3(meshData_loc.n(:,1),meshData_loc.n(:,2),Source,'o')
    
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    
    tmp_gg = ((1-psi_bar.^alpha_M).^alpha_N);
    tmp_gg(imag(tmp_gg) ~= 0) = abs(tmp_gg(imag(tmp_gg) ~= 0)); % occhio!!!
    aa_r_tot = (rr_Gauss*beta_0/R_0+R_0*(1-beta_0)./rr_Gauss);
    gg_tot = tmp_gg.*aa_r_tot;
    
    gg_tot_in = gg_tot;
    gg_tot_in(meshData_loc.ind_G_out) = 0;

    nn = meshData_loc.nn;
    ind_t_selected = meshData_loc.ind_t_INpla;
    ww_Gauss = meshData_loc.ww_Gauss_pla;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Source] = fun_calcSourceFEM_v2(tri,nn,gg_tot_in,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
        [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg_tot_in,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    end
    
end

Source_D = Source(ind_D);

% % figure
% % triplot(tri(:,1:3),meshData_loc.n(:,1),meshData_loc.n(:,2),'k')
% % hold on; axis equal
% % triplot(tri(ind_t_selected,1:3),meshData_loc.n(:,1),meshData_loc.n(:,2),'r')
% % plot(P_Gauss(meshData_loc.ind_G_in,1),P_Gauss(meshData_loc.ind_G_in,2),'o')


%% Sigma rows of stiffness matrix

sigma_a_tmp = zeros(1,nn);
sigma_b_tmp = zeros(1,nn);

sigma_a_tmp(solk.ind_n_axis) = solk.delta_a_weights;
sigma_b_tmp(solk.ind_n_XP)   = solk.delta_b_weights;

if SETTINGS.FIGURES_DEBUG == 1
    plot(meshData_loc.n(solk.ind_n_axis,1),meshData_loc.n(solk.ind_n_axis,2),'r*','LineWidth',2)
    plot(meshData_loc.n(solk.ind_n_XP,1),meshData_loc.n(solk.ind_n_XP,2),'b*','LineWidth',2)
end

sigma_a = [sigma_a_tmp(ind_D) sigma_a_tmp(ind_B)];
sigma_b = [sigma_b_tmp(ind_D) sigma_b_tmp(ind_B)];


%% Eddy current equation

A_CN_1_inv = meshData.A_CN_1_inv;
A_CN_2 = meshData.A_CN_2;

h_step = meshData.h_step;

% % KONNAX_passive = meshData.KONNAX_passive;

% % A = meshData.A;
C = meshData.C;
B_pla = meshData.B_pla;
B_act = meshData.B_act;
D_pla = meshData.D_pla;
D_act = meshData.D_act;

P_tilde = meshData.P_tilde;
P_tilde_bc = meshData.P_tilde_bc;

Q_tilde = meshData.Q_tilde;
Q_tilde_bc = meshData.Q_tilde_bc;

I_act_old = Conductors.Currents_old;
Iphi_old  = solk.Iphi_old(ind_D);
I_act     = Conductors.Currents;


b_x = A_CN_1_inv*A_CN_2*solk.xx_DAE_old + ...  
    .5*h_step*A_CN_1_inv*(B_act*I_act_old + B_pla*Iphi_old + B_act*I_act);

% % A_CN = A_CN_1_inv*A_CN_2;
% % B_CN_act = .5*h_step*A_CN_1_inv*B_act;
% % B_CN_pla = .5*h_step*A_CN_1_inv*B_pla;
% % 
% % xx_prova = (eye(size(A_CN)) - A_CN)\(2*B_pla*Iphi_old + 2*B_act*I_act);
% % 
% % qq = load('a_pas_t0');
% % sum(qq.lambda_Ipla_th - Iphi_old)
% % sum(I_act_old - I_act)
% % sum(qq.lambda_Ipla_th - Iphi_old)
% % 
% % xx_prova = -meshData.A\(B_act*I_act_old + B_pla*qq.lambda_Ipla_th);
% % norm(xx_prova - solk.xx_DAE_old)/norm(solk.xx_DAE_old)
% % 
% % xx_prova2 = (eye(size(A_CN)) - A_CN)\(2*B_pla*Iphi_old + 2*B_act*I_act);
% % 
% % 
% % 
% % xx_t1 = A_CN*solk.xx_DAE_old + B_CN_act*(I_act_old + I_act) + B_CN_pla*(qq.lambda_Ipla_th + Iphi_old);
% % norm(xx_t1 - solk.xx_DAE_old)/norm(xx_t1)
% % 
% % qq2 = load('data_test_DAE.mat');
% % 
% % xx_t11 = qq2.A_CN*qq2.xx_t + qq2.B_CN*(qq2.uu_t + qq2.uu_t1);
% % norm(xx_t11 - qq2.xx_t)/norm(qq2.xx_t)
% % xx_prova = (eye(size(qq2.A_CN)) - qq2.A_CN)\(2*qq2.B_CN*qq2.uu_t);
% % norm(xx_prova - qq2.xx_t)/norm(qq2.xx_t)


% % qq = [A_CN_2*solk.xx_DAE_old .5*h_step*(B_act*I_act_old + B_pla*Iphi_old + B_act*I_act)];

GG_tilde_bc_pas = meshData.GG_tilde_bc_pas;
b_psic_bc = GG_tilde_bc_pas*(C*b_x + D_act*I_act);

% % [C*b_x D_act*I_act]

GG_BEM_D = zeros(nn_B,nn);
GG_BEM_D(:,ind_D) = GG_BEM.psi; % (ALREADY DIVIDED BY 2*PI)
GG_BEM_D(:,meshData_loc.ind_n_NOTpla) = 0;
GG_BEM_D = GG_BEM_D(:,ind_D);

KK_D_tilde = -(GG_BEM_D/mu0)*KK_D;
KK_bc_tilde = sparse(speye(nn_B)-(GG_BEM_D/mu0)*KK_bc);

S_tilde = KK_D_tilde + P_tilde;
S_tilde_bc = KK_bc_tilde + P_tilde_bc;



%% Build stiffness matrix and RHS


AA_psi    = [KK_D KK_bc sparse(nn_D,2) -Source_D];
AA_psibc  = [S_tilde S_tilde_bc sparse(nn_B,3)];
AA_psia   = [sigma_a -1 0 0];
AA_psib   = [sigma_b 0 -1 0];
AA_lambda = [sparse(1,nn) 0 0 sum(Source_D)];

AA = [AA_psi; ...
    AA_psibc; ...
    AA_psia; ...
    AA_psib; ...
    AA_lambda];

% % solk.Psi

b_psi_bc_act = SETTINGS.Psi_MS_bc_act + b_psic_bc;


RHS = [zeros(nn_D,1); b_psi_bc_act; 0; 0; mu0*Ipla];

xx_k = solk.xx_k;

RR_psi_k = AA*xx_k - RHS;

% % xx_k = [solk.Psi(ind_D); ...
% %     solk.Psi(ind_B); ...
% %     solk.Psi_axis; ...
% %     solk.Psi_B; ...
% %     solk.lambda_star; ...
% %     solk.psi_bc_pas_old];

% % xx_new = AA\RHS;
% % 
% % xx_new(i_lambda) = xx_new(i_lambda)/fac_scale_lambda;
% % 
% % Ipla_check = xx_new(i_lambda)*sum(Source_D)/mu0;

ind_plot = floor(1000*rand(1));
if SETTINGS.FIGURES_DEBUG == 1
    figure(ind_plot);  hold on;
% %     plot(plaparameter.Psi_CREATE_loc_B,'k') ; hold on;
    plot(solk.Psi(ind_B),'r--','LineWidth',2)
end


%%

i_nn_D = (1:nn_D)';
i_nn_B = (i_nn_D(end)+1:nn_D+nn_B)';
i_psia = i_nn_B(end)+1;
i_psib = i_psia(end)+1;
i_lambda = i_psib(end)+1;


%%

% % % JAC_0
% % DR_psi_0   = [KK_D KK_bc sparse(nn_D,3)];
% % DRpsibc_0  = [S_tilde S_tilde_bc sparse(nn_B,3)];
% % DRpsia_0   = [sparse(1,nn) -1 0 0];
% % DRpsib_0   = [sparse(1,nn) 0 -1 0];
% % DRlambda_0 = [sparse(1,nn+3)];
% % 
% % JAC_0 = [DR_psi_0; ...
% %     DRpsibc_0; ...
% %     DRpsia_0; ...
% %     DRpsib_0; ...
% %     DRlambda_0];
% % 
% % % JAC_k
% % DR_psi_k   = [Jac_source sparse(nn_D,nn_B)  DRpsi_Dpsia DRpsi_Dpsib DRpsi_Dlambda];
% % DRpsibc_k  = [sparse(nn_B,nn+3)];
% % DRpsia_k   = [sigma_a 0 0 0];
% % DRpsib_k   = [sigma_b 0 0 0];
% % DRlambda_k = [DRlambda_Dpsi DRlambda_Dpsibc DRlambda_Dpsia DRlambda_Dpsib DRlambda_Dlambda];
% % 
% % JAC_k = [DR_psi_k; ...
% %     DRpsibc_k; ...
% %     DRpsia_k; ...
% %     DRpsib_k; ...
% %     DRlambda_k];

%% Solution 

% % fprintf('Computing solution (Direct - sparse LU)\n')
time_solver_start = tic;

if SETTINGS.QuasiNewton == false
    
    calc_Jacobian_FRIDA

    if SETTINGS.RUN_MEX_ROUTINE == false
        h_k = mldivide(JAC,RR_psi_k); % usually faster for purely FEM problems (purely sparse Matrix)
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        h_k = cs_lusol (JAC,RR_psi_k,1,[]); % usually faster for FEM-BEM problems
        % %     opts.ordering = 'amd';
        % %     h_k = umfpack (JAC, '\', RR_psi_k,[],opts); % usually faster than cs_lusol
    end
    
else
    
    if isfield(solk,'RR_psi_k')
        
        % %         if norm(RR_psi_k)/norm(xx_k) < 1e-6
        if norm(RR_psi_k)/norm(solk.RR_psi_k) < SETTINGS.QuasiNewton_factor
            fprintf('Running Quasi-Newton\n')
            SETTINGS.QuasiNewton_fullNR = false;
            h_k = klu (solk.LU, '\', RR_psi_k);
            solk1.LU = solk.LU;
        else 
            SETTINGS.QuasiNewton_fullNR = true;            
        end
        
    else
        SETTINGS.QuasiNewton_fullNR = true;
        
    end
    
    
    if SETTINGS.QuasiNewton_fullNR == true
        fprintf('Running Full-Newton\n')
        calc_Jacobian_FRIDA
        if SETTINGS.RUN_MEX_ROUTINE == false
            [L,U,P,Q,R] = lu(JAC);
            h_k = Q*(U\(L\(P*(R\RR_psi_k))));
        elseif SETTINGS.RUN_MEX_ROUTINE == true
            [LU] = klu (JAC) ; %
            solk1.LU = LU;
            h_k = klu (LU, '\', RR_psi_k) ;
        end
    end
    
end
% % toc(time_solver_start);

time_solver = toc(time_solver_start);

if SETTINGS.VERBOSE
    fprintf('Time per NR step is %5.1s seconds\n', time_solver)
end

%% Post-Processing

% % fprintf('\n')
% % fprintf('Residue psi(end) = %1.3e \n',h_k(nn))
% % fprintf('Residue psi_a    = %1.3e \n',h_k(nn+1))
% % fprintf('Residue psi_b    = %1.3e \n',h_k(nn+2))
% % fprintf('Residue lambda   = %1.3e \n',h_k(nn+3))
% % fprintf('\n')

% % eigs(JAC,1,'sa')
% % eigs(JAC,5,0)/eigs(JAC,1,'la')



xx_k1 = xx_k - h_k;


alfa = 1;


%%% 
% % alfa = linspace(.1,10,100);
% % 
% % for ii = 1:length(alfa)
% %     alfa_k = alfa(ii);
% %     [Residue_new,JAC_new] = fun_Residue_FRIDA_TD_CN(meshData_loc,meshData,solk,Conductors,plaparameter,(xx_k + alfa_k*h_k),SETTINGS);
% %     
% %     c1 = 1e-4;
% %     c2 = .9;
% %     
% %     wolfe_1 = Residue_new  <= (RR_psi_k + (c1*alfa_k*h_k.'*JAC).');
% %     wolfe_2 = -(h_k.'*JAC_new).' <= -(c2*h_k.'*JAC).';
% %     
% %     [numel(wolfe_1) numel(find(wolfe_1)) numel(find(wolfe_2))]
% %     
% % end


%%%

xx_k1(i_lambda) = xx_k1(i_lambda);

Ipla_check = xx_k1(i_lambda)*sum(Source_D)/mu0;
error_Ipla = 100*abs((Ipla_check - plaparameter.Ipla))/plaparameter.Ipla;

if error_Ipla > 1
% %     fprintf('Preliminary check on lambda*Iphi:          error=%6.4f\n',error_Ipla)
% % else
    fprintf('Preliminary check on lambda*Iphi: error=%6.4f\n',error_Ipla)
end

if SETTINGS.SOLVER_RELAXED == true

    treshold_1 = SETTINGS.SOLVER_RELAX_TRESHOLD_1;
    treshold_2 = SETTINGS.SOLVER_RELAX_TRESHOLD_2;
    treshold_3 = SETTINGS.SOLVER_RELAX_TRESHOLD_3;
    treshold_4 = SETTINGS.SOLVER_RELAX_TRESHOLD_4;

    alfa_1 = SETTINGS.SOLVER_RELAX_ALPHA_1;
    alfa_2 = SETTINGS.SOLVER_RELAX_ALPHA_2;
    alfa_3 = SETTINGS.SOLVER_RELAX_ALPHA_3;
    alfa_4 = SETTINGS.SOLVER_RELAX_ALPHA_4;
    alfa_5 = SETTINGS.SOLVER_RELAX_ALPHA_5;

    if error_Ipla < treshold_1
        alfa = alfa_1;
        
    elseif error_Ipla >= treshold_1 && error_Ipla < treshold_2
        alfa = alfa_2;
        
    elseif error_Ipla >= treshold_2 && error_Ipla < treshold_3
        alfa = alfa_3;
        
    elseif error_Ipla >= treshold_3 && error_Ipla < treshold_4
        alfa = alfa_4;
        
    elseif error_Ipla >= treshold_4
        alfa = alfa_5;
        
    end
        
end

xx_k1 = xx_k - alfa*h_k;
    

% % norm(h_k)

Psi_D = xx_k1(1:nn_D);
Psi_B = xx_k1(nn_D+1:nn);

Psi_tot = zeros(meshData_loc.nn,1);
Psi_tot(ind_D) = Psi_D;
Psi_tot(ind_B) = Psi_B;

lambda_star = xx_k1(i_lambda);
lambda = lambda_star/mu0;

Psi_a = xx_k1(i_psia);
Psi_b = xx_k1(i_psib);

Ipla_check = lambda*sum(Source_D);

error_Ipla = 100*abs((Ipla_check - plaparameter.Ipla))/plaparameter.Ipla;

if alfa ~= 1
    fprintf('Check on lambda*Iphi after regularization: error=%6.4f\n',error_Ipla)
end

Iphi = Source*lambda;

Jphi = gg_tot*lambda;


%%
xx_DAE = Q_tilde*Psi_D + Q_tilde_bc*Psi_B + b_x;

% % solk.xx_DAE_old = xx_DAE;

a_pas_new = xx_DAE(1:meshData.n_pas);
v_gap_new = xx_DAE(meshData.n_pas+1:end);

J_c_new = C*xx_DAE + D_act*Conductors.Currents + D_pla*Iphi(ind_D);

% % psi_c_BC_new = GG_tilde_bc_pas*J_c_new;


% % load meshData_pas;
% % figure
% % 
% % ind_n_pas = unique(meshData.t(find(sum(meshData.type.' == meshData.ind_pas.')),:));
% % figure
% % plot3(meshData_pas.n(:,1),meshData_pas.n(:,2),J_c_new,'o'); hold on;
% % plot3(meshData.n(ind_n_pas,1),meshData.n(ind_n_pas,2),J_c_new,'*r')



%%
if SETTINGS.FIGURES_DEBUG == 1
    figure(ind_plot)
    plot((Psi_tot(ind_B)),'b:','LineWidth',2)
end
      
if SETTINGS.FIGURES_DEBUG == 1
    figure
    edges=pdemesh(meshData.n',meshData.e_interface',[]);
    set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
    Rplot=[.9*min(meshData_loc.n(:,1)) 1.08*max(meshData_loc.n(:,1))];
    Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
    rgrid=linspace(Rplot(1),Rplot(2),200);
    zgrid=linspace(Zplot(1),Zplot(2),200);
    [RR,ZZ] = meshgrid(rgrid,zgrid);
    axis([Rplot Zplot]),  hold on,
    PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),Psi_tot,RR,ZZ);
    contour(RR,ZZ,PSI,50);
    colormap('cool'); colorbar vert;
    contour(RR,ZZ,PSI,[Psi_b Psi_b],'k');
    plot(meshData_loc.n(solk.ind_n_axis,1),meshData_loc.n(solk.ind_n_axis,2),'r*','LineWidth',2)
    plot(meshData_loc.n(solk.ind_n_XP,1),meshData_loc.n(solk.ind_n_XP,2),'b*','LineWidth',2)
    plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k.')
end



%% Output

solk1.Psi_MS_bc      = solk.Psi_MS_bc;
% % solk1.Psi_MS_bc_act  = solk.Psi_MS_bc_act;
% % solk1.Psi_MS_bc_pas  = solk.Psi_MS_bc_pas;
solk1.psi_bar        = psi_bar;
solk1.Ipla_check     = Ipla_check;
solk1.xx_k           = xx_k1;
solk1.h_k_psi        = (Psi_tot - solk.Psi);
solk1.h_k            = h_k;
solk1.Psi            = Psi_tot;
solk1.RR_psi_k       = RR_psi_k;
solk1.lambda         = lambda;
solk1.lambda_star    = lambda_star;
solk1.error_Ipla     = error_Ipla;
solk1.Iphi           = Iphi;
solk1.Jphi           = Jphi;
solk1.Psi_a          = Psi_a;
solk1.Psi_b          = Psi_b;
solk1.time_solver    = time_solver;

solk1.xx_DAE        = xx_DAE;
solk1.a_pas         = a_pas_new;
solk1.v_gap         = v_gap_new;
solk1.J_c           = J_c_new;
% % solk1.Iphi_old       = Iphi;



