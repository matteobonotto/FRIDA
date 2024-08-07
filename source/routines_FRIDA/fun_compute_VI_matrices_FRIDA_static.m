function [meshData_pas,...
    G_BCs_VI_act, ...
    G_BCs_VI_pas, ...
    G_flux_VI_pas, ...
    G_B_pickup_VI_pas, ...
    G_flux_VI_act, ...
    G_B_pickup_VI_act, ...
    vec_A] = fun_compute_VI_matrices_FRIDA_static(meshData_ext,meshData_pla,SETTINGS)

%% Select only passive elements

ind_pas = meshData_ext.ind_pas;
if numel(ind_pas) == 1
    ind_t_pas =  find(meshData_ext.type == ind_pas);
else
    ind_t_pas = fun_vec_find(meshData_ext.type,ind_pas);
    % %     ind_t_pas =  find(sum(meshData_ext.type.' == (ind_pas)))';
end

ind_t_comp = ind_t_pas;

temp = meshData_ext.t(ind_t_comp,:)';
ind_n_comp=unique(temp(:),'stable');



node_new=meshData_ext.n(ind_n_comp,:);
N_order= meshData_ext.N_order;

if N_order == 1
    tri_new=zeros(length(ind_t_comp),3);
elseif N_order == 2
    tri_new=zeros(length(ind_t_comp),6);
end

tri_chosen = meshData_ext.t(ind_t_comp,:);
for ii=1:length(ind_n_comp)
    [a,b]=find(tri_chosen==ind_n_comp(ii));
    for jj=1:length(a)
        tri_new(a(jj),b(jj))=ii;
    end
end

if N_order == 1
    tri_new=tri_new(:,1:3);
end


% sort triangles
if size(ind_pas,1)>size(ind_pas,2)
    ind_t_map = repmat(ind_pas,1,2);
else
    ind_t_map = repmat(ind_pas.',1,2);
end

type_old_loc = meshData_ext.type(ind_t_comp);
type_new_loc = zeros(size(tri_new,1),1);

for ii = 1:size(ind_t_map,1)
    if ind_t_map(ii,1) == ind_t_map(ii,2)
        ind_sel = find((type_old_loc == (ind_t_map(ii,1):ind_t_map(ii,2))));
    else
        ind_sel = find(sum(type_old_loc.' == (ind_t_map(ii,1):ind_t_map(ii,2)).'))';
    end
    type_new_loc(ind_sel) = ii;
end

temp_node_new = zeros(size(node_new));
temp_tri_new  = zeros(size(tri_new));

kk = 1;
for ii = 1:size(ind_t_map,1)
    ind_nodes_ii = unique(tri_new(find(type_new_loc == ii),:));
    for jj = 1:numel(ind_nodes_ii)
        temp_node_new(kk,:) = node_new(ind_nodes_ii(jj),:);
        [ind_row,ind_col,~] = find(tri_new == ind_nodes_ii(jj));
        for hh = 1:numel(ind_row)
            temp_tri_new(ind_row(hh),ind_col(hh)) = kk;
        end
        kk = kk+1;
    end
end

tri_new = temp_tri_new;
node_new = temp_node_new;

meshData_pas.t = tri_new;
meshData_pas.n = node_new;
meshData_pas.type = meshData_ext.type(ind_t_comp);
meshData_pas.nn = size(meshData_pas.n,1);
meshData_pas.nt = size(meshData_pas.t,1);
meshData_pas.N_order = N_order;


if SETTINGS.FIGURES == true
    pdemesh(node_new',[],[tri_new(:,1:3)'; ones(1,size(tri_new,1))],'EdgeColor','k'); hold on
end

%%
tri = meshData_pas.t;
nodes = meshData_pas.n;
type_pas = meshData_pas.type;

nodes_matrix = [nodes(tri(:,1),:) ...
    nodes(tri(:,2),:) ...
    nodes(tri(:,3),:)];

P1 = nodes_matrix(:,[1 2]);
P2 = nodes_matrix(:,[3 4]);
P3 = nodes_matrix(:,[5 6]);

degree_G_target = 8;
degree_G_source = 12;
[w_G_target,P_G_target,n_G_target] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_target);
[w_G_source,P_G_source,n_G_source] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_source);

if SETTINGS.RUN_MEX_ROUTINE
    [MatInd_nodes_tri] = fun_MatInd_nodes_tri_fast_mex(meshData_pas.N_order,meshData_pas.t,meshData_pas.nn);
else
    [MatInd_nodes_tri] = fun_MatInd_nodes_tri(meshData_pas.N_order,meshData_pas.t,meshData_pas.nn);
end

tri_act = meshData_ext.t;
nodes_act = meshData_ext.n;
keyreg_act = meshData_ext.type;
ind_act = meshData_ext.ind_act;
n_act = length(ind_act);

degree_G_source_ext = degree_G_target;


%%
eta = zeros(meshData_pas.nt,1);

for ii = 1:length(meshData_ext.ind_pas)
    eta(meshData_pas.type == meshData_ext.ind_pas(ii)) = meshData_ext.rho_pas(ii);
end

%%
% % tic
% % L_VI = fun_L_VI_stable_symmetric_fast(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     degree_G_source, ...
% %     n_G_source, ...
% %     P_G_source,...
% %     degree_G_target, ...
% %     n_G_target, ....
% %     P_G_target);
% % toc

% % tic
% % R_VI = fun_R_VI_stable(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     degree_G_source, ...
% %     n_G_source, ...
% %     P_G_source,...
% %     eta);
% % toc

% % tic
% % U_VI = fun_U_VI_stable_fast(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     degree_G_source, ...
% %     n_G_source, ...
% %     P_G_source);
% % toc

% % tic
% % V_VI = fun_V_VI_stable_fast(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     degree_G_source, ...
% %     n_G_source, ...
% %     P_G_source);
% % toc


% % tic
% % D_VI = fun_D_VI(tri,...
% %     type_pas, ...
% %     ind_t_map);
% % toc


% % tic
% % [M_VI,~] = fun_M_act_VI_stable_fast(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     degree_G_source_ext, ...
% %     degree_G_target, ...
% %     n_G_target, ....
% %     P_G_target,...
% %     n_act, ...
% %     tri_act, ...
% %     nodes_act, ...
% %     ind_act, ...
% %     keyreg_act);
% % toc


%% Green Mat for BCs

R_BCs = meshData_pla.n(meshData_pla.ind_n_bc,1);
Z_BCs = meshData_pla.n(meshData_pla.ind_n_bc,2);

% Passive contribution
tic
G_BCs_VI_pas = fun_G_flux_passive_VI_stable_fast(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    R_BCs, ...
    Z_BCs);
toc

%%% Active contribution
tic
G_BCs_VI_act =fun_G_flux_active_VI_stable_fast(n_act, ...
    tri_act, ...
    nodes_act, ...
    ind_act, ...
    keyreg_act, ...
    degree_G_source, ...
    R_BCs, ...
    Z_BCs);
toc


%% Green Mat for Coupling surface

% % R_CS = meshData_pla.C_surf(:,1);
% % Z_CS = meshData_pla.C_surf(:,2);
% % 
% % % Passive contribution
% % tic
% % G_CS_VI_pas = fun_G_flux_passive_VI_stable_fast(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     degree_G_source, ...
% %     n_G_source, ...
% %     P_G_source,...
% %     R_CS, ...
% %     Z_CS);
% % toc
% % G_CS_VI_pas = G_CS_VI_pas.';

% % fil_source = meshData_pla.C_surf;
% % [G_CS_VI_pas2] = fun_M_fil_VI_stable_fast(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     degree_G_target, ...
% %     n_G_target, ....
% %     P_G_target,...
% %     fil_source);
% % toc

%%% Active contribution
% % tic
% % G_CS_VI_act =fun_G_flux_active_VI_stable_fast(n_act, ...
% %     tri_act, ...
% %     nodes_act, ...
% %     ind_act, ...
% %     keyreg_act, ...
% %     degree_G_source, ...
% %     R_CS, ...
% %     Z_CS);
% % toc



%% Green Mat for Probes
%%% Flux loops
if any(strcmp('flux_loops',fieldnames(SETTINGS)))
    
    R_fluxloop = SETTINGS.flux_loops(:,1);
    Z_fluxloop = SETTINGS.flux_loops(:,2);
   
    %%% Passive contribution
    tic
    G_flux_VI_pas = fun_G_flux_passive_VI_stable_fast(tri, ...
        nodes, ...
        N_order, ...
        degree_G_source, ...
        n_G_source, ...
        P_G_source,...
        R_fluxloop, ...
        Z_fluxloop);
    toc
    
    %%% Active contribution
    tic
    G_flux_VI_act =fun_G_flux_active_VI_stable_fast(n_act, ...
        tri_act, ...
        nodes_act, ...
        ind_act, ...
        keyreg_act, ...
        degree_G_source, ...
        R_fluxloop, ...
        Z_fluxloop);
    toc
    
else
    
    G_flux_VI_pas = [];
    G_flux_VI_act = [];
    
end

%%% Pick-up
if any(strcmp('pickup',fieldnames(SETTINGS)))
        
    R_pickup = SETTINGS.pickup(:,1);
    Z_pickup = SETTINGS.pickup(:,2);
    Vers     = SETTINGS.pickup(:,3:4);
    
    %%% Passive contribution
    tic
    [G_Br_VI_pas,G_Bz_VI_pas] = fun_G_BrBz_passive_VI_stable_fast(tri, ...
        nodes, ...
        N_order, ...
        degree_G_source, ...
        n_G_source, ...
        P_G_source,...
        R_pickup, ...
        Z_pickup);
    toc
    
    G_B_pickup_VI_pas = G_Br_VI_pas.*Vers(:,1) + G_Bz_VI_pas.*Vers(:,2);
    
    %%% Active contribution
    tic
    [G_Br_VI_act,G_Bz_VI_act] =fun_G_BrBz_active_VI_stable_fast(n_act, ...
        tri_act, ...
        nodes_act, ...
        ind_act, ...
        keyreg_act, ...
        degree_G_source, ...
        R_pickup, ...
        Z_pickup);
    toc
    
    G_B_pickup_VI_act = G_Br_VI_act.*Vers(:,1) + G_Bz_VI_act.*Vers(:,2);
    
    else
    
    G_B_pickup_VI_pas = [];
    G_B_pickup_VI_act = [];
    
end



%%

vec_A = fun_vec_Area_element_VI_stable(tri, ...
    nodes, ...
    N_order, ...
    MatInd_nodes_tri,...
    degree_G_source, ...
    n_G_source, ...
    ind_pas, ...
    meshData_pas.type);




