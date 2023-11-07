%% MATLAB 2 GMSH x ITER MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Write ITER geometrical data write in a suitable format for GMSH mesh
%   generation.
%   
%   Matteo Bonotto 01/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
clear;  clc;  close all;

%% Input Geometry
load('ITER_mesh_CREATE');
load('ITER_geo_act');

%%
nt=size(meshData.t,1);
P1=[meshData.n(meshData.t(:,1),1) meshData.n(meshData.t(:,1),2)];
P2=[meshData.n(meshData.t(:,2),1) meshData.n(meshData.t(:,2),2)];
P3=[meshData.n(meshData.t(:,3),1) meshData.n(meshData.t(:,3),2)];
centro_t=(P1+P2+P3)/3;

keyreg = meshData.type;


figure
    hold on;    grid on
% %     pdemesh(p,[],t)
    mesh=triplot(meshData.t,meshData.n(:,1),meshData.n(:,2),'k');
    axis equal,  hold on, xlabel('r [m]'), ylabel('z [m]')
% %     axis([0.8 4.5 -1.5 1.5])
% Conduttori
ind_conduttori=find(keyreg>0 & keyreg<500);
    triplot(meshData.t(ind_conduttori,:),meshData.n(:,1),meshData.n(:,2),'r')
% % % shell 
% % ind_shell=find(keyreg>=117 & keyreg<=176);
% %     triplot(meshData.t(ind_shell,:),meshData.n(:,1),meshData.n(:,2),'g')
% % % vessel
% % ind_vessel=find(keyreg>=57 & keyreg<=116);
% %     triplot(meshData.t(ind_vessel,:),meshData.n(:,1),meshData.n(:,2),'y')
% % % TSS
% % ind_TSS=find(keyreg>=177 & keyreg<=236);
% %     triplot(meshData.t(ind_TSS,:),meshData.n(:,1),meshData.n(:,2),'b')
% % % %     plot(centro_t(ind_TSS,1),centro_t(ind_TSS,2),'bo')
% % % saddle
% % ind_saddle=find(keyreg>=237 & keyreg<=244);
% %     triplot(meshData.t(ind_saddle,:),meshData.n(:,1),meshData.n(:,2),'c')
% % % %     plot(centro_t(ind_saddle,1),centro_t(ind_saddle,2),'co')
% % % Aactive coils (256 conduttori collegati in 194 circuiti)
% % ind_active=find(keyreg>=1 & keyreg<=56);
% %     triplot(meshData.t(ind_active,:),meshData.n(:,1),meshData.n(:,2),'r')
%     
%% VV, Shell,
load('ITER_fw');
FW = fw;

n_FW = 800;  n_bc = ceil(n_FW/1.8); mesh_name = 'ITER_meshData_FRIDA_evol_HHR';
n_FW = 600;  n_bc = ceil(n_FW/1.8); mesh_name = 'ITER_meshData_FRIDA_evol_HR';
n_FW = 400;  n_bc = ceil(n_FW/1.8); mesh_name = 'ITER_meshData_FRIDA_evol_MR';
n_FW = 260;  n_bc = ceil(n_FW/2); mesh_name = 'ITER_meshData_FRIDA_evol_LR';
n_FW = 200;  n_bc = ceil(n_FW/1.6); mesh_name = 'ITER_meshData_FRIDA_evol_LLR';
% % n_bc = 150;

[FW_r,FW_z]=equispaced(FW(:,1),FW(:,2),n_FW);
FW = [FW_r(1:end) FW_z(1:end)];


load('ITER_VV');
VV_in = .5*(shell_int(1:2:end-1,:) + shell_int(2:2:end,:));

[VV_in_r,VV_in_z]=equispaced(VV_in(:,1),VV_in(:,2),n_bc);
VV_in = [VV_in_r(1:end) VV_in_z(1:end)];


% % figure
% % plot(shell_int(:,1),shell_int(:,2),'o-')
% % axis equal; hold on
% % plot(VV_in(:,1),VV_in(:,2),'*-r')

figure
plot(VV_in(:,1),VV_in(:,2),'o-')
axis equal; hold on
plot(FW(:,1),FW(:,2),'*-r')

VV = VV_in;



% boundary
ind = boundary(FW(:,1),FW(:,2),1);
VV = FW(ind,:);

VV = unique(VV,'rows');



[VV,~] = fun_ordinapunti(VV);
[RR_EQUI,ZZ_EQUI]=equispaced(VV(:,1),VV(:,2),200);
VV = [RR_EQUI ZZ_EQUI];

U_t = zeros(size(VV));
U_n = zeros(size(VV));
nfil=length(VV(:,1));
rstart=VV(:,1);
zstart=VV(:,2);
for kpoint=1:nfil
    kprev=kpoint-1; if(kprev<1) kprev=nfil; end
    kfol=kpoint+1; if(kfol>nfil) kfol=1; end
    t=[rstart(kfol)-rstart(kprev); zstart(kfol)-zstart(kprev)];
    t=t/norm(t);
    U_t(kpoint,:)=t';
    
    U_n(kpoint,:)=[t(2) -t(1)]; 

end

factor = -.2;

VV = VV + factor*U_n;

plot(VV(:,1),VV(:,2),'or-')


Ne = floor(.15*n_FW);
r0=6.1;
z0=0.4;
a=.65;
b=1.05;
theta=linspace(0,2*pi,Ne)';
ell.RR=r0+a*cos(theta(1:end-1));
ell.ZZ=z0+b*sin(theta(1:end-1));
ellisse1 = [ell.RR ell.ZZ];
[ellisse1r,ellisse1z]=equispaced(ellisse1(:,1),ellisse1(:,2),Ne);
ellisse1 = [ellisse1r(1:end) ellisse1z(1:end)];


Ne = floor(.47*n_FW);
r0=6.1;
z0=0.4;
a=.75*2;
b=1.25*2;
theta=linspace(0,2*pi,Ne)';
ell.RR=r0+a*cos(theta(1:end-1));
ell.ZZ=z0+b*sin(theta(1:end-1));
ellisse2 = [ell.RR ell.ZZ];
[ellisse2r,ellisse2z]=equispaced(ellisse2(:,1),ellisse2(:,2),Ne);
ellisse2 = [ellisse2r(1:end) ellisse2z(1:end)];


plot(ellisse1(:,1),ellisse1(:,2),'-k')
plot(ellisse2(:,1),ellisse2(:,2),'-r')

%% Geometria Saddle coils
for j=1:14
% %     pause
    conduttori{j}=[source_act.RR(j,:)' source_act.ZZ(j,:)'] ;
end

%%
addpath ./AutoMESH_ver_2;

% % % mesh_name = 'ITER_meshData_FRIDA_evol_HR';

GEOMETRY.shape{1} = ellisse1;
GEOMETRY.shape{2} = ellisse2;
GEOMETRY.shape{3} = FW;
GEOMETRY.shape{4} = VV;
for ii=1:numel(conduttori)
    GEOMETRY.shape{ii+4}=conduttori{ii};
end

GEOMETRY.spacing{1} = 1;
GEOMETRY.spacing{2} = 1;
GEOMETRY.spacing{3} = 1;
GEOMETRY.spacing{4} = 1;
for ii=1:numel(conduttori)
    GEOMETRY.spacing{ii+4}=.2;
end

GEOMETRY.surface{1} = 1;
GEOMETRY.surface{2} = [2 1];
GEOMETRY.surface{3} = [3 2];
GEOMETRY.surface{4} = [4 3];
for ii=1:numel(conduttori)
    GEOMETRY.surface{ii+4}=ii+4;
end

[p,e,t]=autoMESH_complex(mesh_name,GEOMETRY,2);

figure
triplot(t(:,1:3),p(:,1),p(:,2),'k')
axis equal; hold on;


%%

two_d_nodes = p;
two_d_elements = t;
N_order = 2;


%%
figure
c=[0.5 0.5 0.5];
pdemesh(two_d_nodes',[],[two_d_elements(:,1:3)'; ones(1,size(two_d_elements,1))])
hold on;
% % plot(two_d_nodes(:,1),two_d_nodes(:,2),'ro')
% % plot(two_d_nodes(:,1),two_d_nodes(:,2),'or')
% % 
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('GMsh to MATLAB import','fontsize',14)


% % plot(two_d_nodes(two_d_elements(1:280,1),1),...
% %     two_d_nodes(two_d_elements(1:280,1),2),'r*');
% % hold on

% % for ii=1:size(cond,2)
% %     conduttore=cond{ii};
% %     plot(conduttore(:,1),conduttore(:,2),'*c')
% % end

%baricentro triangoli
if N_order == 1
    P1=[two_d_nodes(two_d_elements(:,1),1), ...
        two_d_nodes(two_d_elements(:,1),2)];
    P2=[two_d_nodes(two_d_elements(:,2),1), ...
        two_d_nodes(two_d_elements(:,2),2)];
    P3=[two_d_nodes(two_d_elements(:,3),1), ...
        two_d_nodes(two_d_elements(:,3),2)];
    
    centro=[sum(P1(:,1)+P2(:,1)+P3(:,1),2)/3, ...
        sum(P1(:,2)+P2(:,2)+P3(:,2),2)/3];
elseif N_order == 2
    P1=[two_d_nodes(two_d_elements(:,1),1), ...
        two_d_nodes(two_d_elements(:,1),2)];
    P2=[two_d_nodes(two_d_elements(:,2),1), ...
        two_d_nodes(two_d_elements(:,2),2)];
    P3=[two_d_nodes(two_d_elements(:,3),1), ...
        two_d_nodes(two_d_elements(:,3),2)];
    P4=[two_d_nodes(two_d_elements(:,4),1), ...
        two_d_nodes(two_d_elements(:,4),2)];
    P5=[two_d_nodes(two_d_elements(:,5),1), ...
        two_d_nodes(two_d_elements(:,5),2)];
    P6=[two_d_nodes(two_d_elements(:,6),1), ...
        two_d_nodes(two_d_elements(:,6),2)];
    
    centro=[sum(P1(:,1)+P2(:,1)+P3(:,1)+P4(:,1)+P5(:,1)+P6(:,1),2)/6, ...
        sum(P1(:,2)+P2(:,2)+P3(:,2)+P4(:,2)+P5(:,2)+P6(:,2),2)/6];
elseif N_order == 3
    P1=[two_d_nodes(two_d_elements(:,1),1), ...
        two_d_nodes(two_d_elements(:,1),2)];
    P2=[two_d_nodes(two_d_elements(:,2),1), ...
        two_d_nodes(two_d_elements(:,2),2)];
    P3=[two_d_nodes(two_d_elements(:,3),1), ...
        two_d_nodes(two_d_elements(:,3),2)];
    P4=[two_d_nodes(two_d_elements(:,4),1), ...
        two_d_nodes(two_d_elements(:,4),2)];
    P5=[two_d_nodes(two_d_elements(:,5),1), ...
        two_d_nodes(two_d_elements(:,5),2)];
    P6=[two_d_nodes(two_d_elements(:,6),1), ...
        two_d_nodes(two_d_elements(:,6),2)];
    P7=[two_d_nodes(two_d_elements(:,7),1), ...
        two_d_nodes(two_d_elements(:,7),2)];
    P8=[two_d_nodes(two_d_elements(:,8),1), ...
        two_d_nodes(two_d_elements(:,8),2)];
    P9=[two_d_nodes(two_d_elements(:,9),1), ...
        two_d_nodes(two_d_elements(:,9),2)];
    
    centro=[sum(P1(:,1)+P2(:,1)+P3(:,1)+P4(:,1)+P5(:,1)+P6(:,1)+P7(:,1)+P8(:,1)+P9(:,1),2)/9, ...
        sum(P1(:,2)+P2(:,2)+P3(:,2)+P4(:,2)+P5(:,2)+P6(:,2)+P7(:,2)+P8(:,2)+P9(:,2),2)/9];
end



% % plot(centro(:,1),centro(:,2),'.r


% % n_points=250;
% % VV_in=[1.995+0.475*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
% %     0.475*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];
% % VV_ex=[1.995+0.505*sin(linspace(0,2*pi-2*pi/n_points,n_points))' ...
% %     0.505*cos(linspace(0,2*pi-2*pi/n_points,n_points))'];
% % n_points=200;
% % Shel_in=[1.995+0.5115*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
% %     0.5115*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];
% % Shel_ex=[1.995+0.5145*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
% %     0.5145*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];
% % TSS_in=[2+0.553*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
% %     0.553*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];
% % TSS_ex=[2+0.6*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
% %     0.6*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];

type=zeros(size(two_d_elements,1),1);
% elementi interni al vessel
% % [ind_t_pla,~]=inpolygon(centro(:,1),centro(:,2),...
% %     VV_in(:,1),VV_in(:,2));
% % ind_t_pla=find(ind_t_pla);


% % ind_conduttori=find(keyreg>0);
%     triplot(meshData.t(ind_conduttori,:),meshData.n(:,1),meshData.n(:,2),'r')
% plasma 
[ind_t_pla,~]=inpolygon(centro(:,1),centro(:,2),...
    FW(:,1),FW(:,2));
ind_t_pla=find(ind_t_pla);


% % plot(centro(ind_t_pla,1),centro(ind_t_pla,2),'o')


% % mesh = pdemesh(two_d_nodes',[],[two_d_elements(ind_t_pla,:)'; ones(size(ind_t_pla))']','EdgeColor','r');
triplot(two_d_elements(ind_t_pla,1:3),two_d_nodes(:,1),two_d_nodes(:,2),'r');



for j=1:14
% %     pause
    cond{j}=[source_act.RR(j,:)' source_act.ZZ(j,:)'] ;
end

keyreg_create=meshData.type;
meshData_create=meshData;

nt=size(two_d_elements,1);
keyreg=zeros(1,nt);
point.RR=meshData_create.n(:,1);
point.ZZ=meshData_create.n(:,2);
centro_coils=zeros(size(cond,2),2);
for ii=1:size(cond,2)
    centro_coils(ii,:)=sum(cond{ii})/size(cond{ii},1);
end


ind_conduttori=[];
for ii=1:size(cond,2)
    Points=cond{ii};
    ind=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
    ind_conduttori=ind;
    c=rand(3,1);
    triplot(two_d_elements(ind,1:3),two_d_nodes(:,1),two_d_nodes(:,2),'c');
    keyreg(ind_conduttori)=ii;
end


% Plasma
keyreg(ind_t_pla)=-1;

%%
% nodi FW
ind_n_in=reshape(two_d_elements(ind_t_pla,:),numel(two_d_elements(ind_t_pla,:)),1);
ind_n_in=unique(ind_n_in);
ind_n=boundary(two_d_nodes(ind_n_in,1),two_d_nodes(ind_n_in,2),1);
ind_n=ind_n_in(ind_n);
ind_n_FW=(ind_n(1:end-1));
        
plot(two_d_nodes(ind_n_FW,1),two_d_nodes(ind_n_FW,2),'ob','LineWidth',2) 

% nodi boundary
ind_t_vac = find(keyreg == 0);
ind_n_in=reshape(two_d_elements(ind_t_vac,:),numel(two_d_elements(ind_t_vac,:)),1);
ind_n_in=unique(ind_n_in);
ind_n=boundary(two_d_nodes(ind_n_in,1),two_d_nodes(ind_n_in,2),.9);
ind_n=ind_n_in(ind_n);
ind_n_bc=(ind_n(1:end-1));
        
plot(two_d_nodes(ind_n_bc,1),two_d_nodes(ind_n_bc,2),'og','LineWidth',2) 

%%
ind_bound=boundary(two_d_nodes(:,1),two_d_nodes(:,2),0);
ind_bound=ind_bound(1:end-1);

%% Lati di bordo (tutto tranne il vuoto)
ind_vacuum=find(keyreg>=0);
ind_not_vacuum=1:numel(keyreg);
figure
pdemesh(two_d_nodes',[],[two_d_elements(ind_vacuum,1:3)'; ones(1,length(ind_vacuum))],'EdgeColor','k');
% % mesh = pdemesh(two_d_nodes',[],two_d_elements(ind_vacuum,1:4)','EdgeColor','k');
LATIBORDO=1;
meshData.t_choosen=two_d_elements(ind_vacuum,:);
[~,lati_bordo]=fun_lati(meshData,LATIBORDO);
figure
edges=pdemesh(two_d_nodes',lati_bordo',[]);
    set(edges,'color','k')
    set(edges,'linewidth',1)  

%% Midpoints (quadratic interpolation)
% % % Lati (tutto)
% % meshData.t_choosen=two_d_elements;
% % [lati,lati_bordo]=fun_lati(meshData,0);
% %     
% % midpoints=zeros(size(lati,1),2);
% % for ii=1:size(lati,1)
% %     midpoints(ii,:)=sum(two_d_nodes(lati(ii,:),:))*0.5;
% % end
% % two_d_nodes=[two_d_nodes; midpoints];
% % 
% % jj=size(two_d_nodes,1)+1;
% % midpoints=zeros(size(two_d_nodes));
% % midpoints_index=zeros(size(two_d_elements));
% % tic
% % for ii=1:size(two_d_elements,1)
% % % %     ii
% %     triangle=two_d_elements(ii,:);
% %     P=two_d_nodes(triangle,:);
% %     M=[0.5*(P(1,:)+P(2,:))
% %         0.5*(P(2,:)+P(3,:));...
% %         0.5*(P(3,:)+P(1,:))];
% %     a=ismember(M(1,:),midpoints,'rows');
% %     b=ismember(M(2,:),midpoints,'rows');
% %     c=ismember(M(3,:),midpoints,'rows');
% %     if ~a
% %         midpoints(ii,:)=M(1,:);
% %         mid(1)=jj;
% %         jj=jj+1;
% %     else
% %         mid(1)=size(two_d_nodes,1)+find(ismember(midpoints,M(1,:),'rows'));
% %     end
% %     if ~b
% %         midpoints(ii,:)=M(2,:);
% %         mid(2)=jj;
% %         jj=jj+1;
% %     else
% %         mid(2)=size(two_d_nodes,1)+find(ismember(midpoints,M(2,:),'rows'));
% %     end
% %     if ~c
% %         midpoints(ii,:)=M(3,:);
% %         mid(3)=jj;
% %         jj=jj+1;
% %     else
% %         mid(3)=size(two_d_nodes,1)+find(ismember(midpoints,M(3,:),'rows'));
% %     end
% %     midpoints_index(ii,:)=mid;
% % end
% % toc
% % two_d_elements=[two_d_elements midpoints_index];
% % two_d_nodes=[two_d_nodes; midpoints];
% % 
% % figure
% % pdemesh(two_d_nodes',[],[two_d_elements ones(size(two_d_elements,1),1)]');

%% Salvataggio dati (file meshData)
clear meshData
meshData.t=two_d_elements;
meshData.n=two_d_nodes;
meshData.e=lati_bordo;

%% Save Data
meshData.type=keyreg;
meshData.b=ind_bound;
meshData.vess=ind_n_FW;
meshData.ind_n_FW=ind_n_FW;
meshData.ind_n_bc=ind_n_bc;
meshData.type_explanation{1,1}='keyreg -1: plasma region';
meshData.type_explanation{2,1}='keyreg > 0: coils';
% % meshData.type_explanation{3,1}='keyreg = 117: shell';
% % meshData.type_explanation{4,1}='keyreg = 57: vessel';
% % meshData.type_explanation{5,1}='keyreg = 237: saddle coils';
% % meshData.type_explanation{6,1}='keyreg = 177: TSS';


save([mesh_name '.mat'], 'meshData');




