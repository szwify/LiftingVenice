%% Problem INPUT
Htot = 30;
L = 2*Htot;
H1 = 15;
H2 = 20;
H3 = 22;
H4 = Htot;

% Mesh PARAMETERS
Hs = [0,2;
    -H1,1;
    -H2,0.5; 
    -H3,0.5;
    -H4,1;
    ];
% Mesh GENERATION
[vert,etri, tria,tnum] = stratimesher(L, Hs);
the_coor=vert; 
connect=tria;

% Mesh PLOT
    figure(1);
    patch('faces',tria(tnum==1,1:3),'vertices',vert, ...
        'facecolor',[0.9,0.9,0.9], ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; %axis image off;
    patch('faces',tria(tnum==2,1:3),'vertices',vert, ...
        'facecolor',[0.7,0.7,0.7], ...
        'edgecolor',[.2,.2,.2]) ;
    patch('faces',tria(tnum==3,1:3),'vertices',vert, ...
        'facecolor',[0.4,0.4,0.4], ...
        'edgecolor',[.2,.2,.2]) ;
    patch('faces',tria(tnum==4,1:3),'vertices',vert, ...
        'facecolor',[0.7,0.7,0.7], ...
        'edgecolor',[.2,.2,.2]) ;    
%     patch('faces',tria(tnum==5,1:3),'vertices',vert, ...
%         'facecolor',[.5,.5,.5], ...
%         'edgecolor',[.2,.2,.2]) ;    
%     patch('faces',edg(:,1:2),'vertices',node, ...
%         'facecolor','w', ...
%         'edgecolor',[.1,.1,.1], ...
%         'linewidth',1.5) ;
    title('Lifting Venice - MESH') ;
%     axis ij;
    ylabel('Depth [m]');
    xlabel('Length [m]');
    legend({'Layer XYZ','Clay','Injection','Clay'},'Location','best')


mesh = FEmesh(the_coor,connect);


%% 2D axisymmetry
Config='Axis';
%
% MATERIAL PROPERTIES
% all stiffness in MPa, ohio Sandstone

k=8.4e3; % elastic drained bulk modulus
g=6.8e3;  % shear modulus
b=0.707692; % Biot coefficient 
%phi=0.19;  % porosity 
%k_f=2.2e3;  % fluid bulk modulus
M=9.18478e3; % Biot Modulus
%overN=1./M-phi/k_f;N=1/overN  % Biot intrinsic Modulus (fluid independent)
k_u=k+b^2*M     % undrained bulk modulus
perm =0.137549e-3; % permeability value adjusted such that c= k /(mu S)=1 , where S is storage coefficient

mu_f=1; rho=1;
kappa = perm / mu_f  ; 

[E_u,nu_u]=Enu_from_kg(k_u,g);

[E,nu]=Enu_from_kg(k,g);

% undrained elastic stiffness tensor
L_u=Elastic_Isotropic_Stiffness(k_u,g,'Axis');

% drained elastic stiffness tensor
L_d=Elastic_Isotropic_Stiffness(k,g,'Axis');


%% Boundary Condition
klt_down = find(mesh.XY(:,2)==-H4);
Imp_displacement=[klt_down  ones(length(klt_down),1)  zeros(length(klt_down),1) ];

klt_left = find(mesh.XY(:,1)==0);
Imp_displacement=[Imp_displacement;
    klt_left  zeros(length(klt_left),1)  ones(length(klt_left),1) ];

klt_right = find(mesh.XY(:,1)==L);
Imp_displacement=[Imp_displacement;
    klt_right  zeros(length(klt_right),1)  ones(length(klt_right),1) ];

 
hold on;
plot(mesh.XY(klt_down,1),mesh.XY(klt_down,2),'og','HandleVisibility','off')
hold on;
plot(mesh.XY(klt_left,1),mesh.XY(klt_left,2),'og','HandleVisibility','off')
hold on;
plot(mesh.XY(klt_right,1),mesh.XY(klt_right,2),'og','HandleVisibility','off')





%%%% create boundary loads matrix ....format is dof node1 intensity1 node2
%%%% intensity 2   (where node 1 node 2 define a segment)
%%%% 
deltap=1;
ktl_flux=find(mesh.XY(:,1)==0 & mesh.XY(:,2)<-H2 & mesh.XY(:,2)>-H3);
flux_load =[];
hold on;
plot(mesh.XY(ktl_flux,1),mesh.XY(ktl_flux,2),'dr')
hold on;
%    switch from fr,0=  to----- fx fy    - to be recheck
for seg=1:length(ktl_flux)-1
    
    flux = 1;
    
    flux_load =[flux_load ;...
                   1 ktl_flux(seg) flux ktl_flux(seg+1) flux;
                   2 ktl_flux(seg) 0 ktl_flux(seg+1) 0]; 
end

%Boundary_loads =[];
%
% no Initial stress field 
mySig_o= zeros(mesh.Nelts,4);  % 4 components in 2D axi-symmetry (srr, szz, srz, stt
% mySig_o(:,1)=1; 
% mySig_o(:,2)=1;
% mySig_o(:,4)=1;

%% different matrices

% elasticity
proplist={L_d};
[K,ID_array]=AssembleMatrix(mesh,'Axis','Elasticity',proplist,3);

[Fflux]=AssembleVectorBoundaryTerm(mesh,'Axis','BoundaryLoads',flux_load,ID_array,2);
[Fbody]=AssembleVectorVolumeTerm(mesh,'Axis','InitialStress',mySig_o,ID_array,3);

Fload=Fbody*0;

[eq_free_u,fix_nonZero_u,eq_fix_u]=PrepareDirichletBC(Imp_displacement,ID_array);

% mass matrix term
proplist={1./M};
[S,Id_p]=AssembleMatrix(mesh,'Axis','Mass',proplist,3);
% pore-pressure fixed to zero on the outer radius.
eq_fix_p=[Id_p(ktl_rad)  ];
eq_free_p=setdiff(Id_p(:),eq_fix_p(:));

% laplacian term
proplist={kappa};
[D,Id_p]=AssembleMatrix(mesh,'Axis','Laplacian',proplist,3);

% Coupling term 
proplist={b};
[C,Id_u,Id_p]=AssembleCouplingMatrix(mesh,mesh,'Axis',proplist,3);

ntot_u=length(Id_u(:));
ntot_p=length(Id_p(:));
%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%
% implicit time - integration

%% Solving directly for the undrained solution at t=0 ! 
%  We dont fix the pore pressure to zero on the outer radius for the undrained response.
dt= 0.;
AA=S+dt*D;

ntot=ntot_u+ntot_p;

TotMat=sparse(ntot,ntot);
TotMat(1:ntot_u,1:ntot_u)=K;
TotMat(ntot_u+1:ntot,ntot_u+1:ntot)=-AA;
TotMat(ntot_u+1:ntot,1:ntot_u)=-C';
TotMat(1:ntot_u,ntot_u+1:ntot)=-C;
%
Ftot=sparse(ntot,1);
Ftot(1:ntot_u)=Fbody;

eq_free=[eq_free_u;ntot_u   + Id_p(:)];%+eq_free_p];
Undrained_sol=sparse(ntot,1);
Undrained_sol(eq_free) = TotMat(eq_free,eq_free)\Ftot(eq_free);

pp_o=Undrained_sol(ntot_u   + Id_p(:));
U_o=sparse(ntot_u,1);
U_o(eq_free_u) = Undrained_sol(eq_free_u);
%  
AnalyticSolUr= @(r) -(r*(1-2* nu_u)/(2.* g*(1.+ nu_u)));

% radial displacement at z=0 and at at r=0 -> should be the same....
%
figure(2)
plot(mesh.XY(klt_z,1), U_o(ID_array(klt_z,1)),'ok' )
hold on
plot(mesh.XY(klt_r,2), U_o(ID_array(klt_r,2)),'*k' )
hold on
plot(mesh.XY(klt_r,2),AnalyticSolUr(mesh.XY(klt_r,2)),'-r');

% reshape Usol_u
udisp = [U_o(ID_array(:,1)) U_o(ID_array(:,2))];

% plot deformed mesh
figure(3)
plotmesh(mesh.XY,connect,[.2 .2 .2],'w')
hold on;
plotmesh(mesh.XY+udisp*1e3,connect,[.8 .2 .2],'none')

% t=0+ solution is the undrained response
%pp_o=ones(mesh_fem.FEnode.n_tot,1)*Pp_u_u;
pp_o(eq_fix_p)=0.   ; % NOW we fix the pore pressure to zero. 

figure(4) 
trisurf(mesh.conn,mesh.XY(:,1),mesh.XY(:,2),full(pp_o))
 
%%

%%%%%%
dt=0.002;
%----- matrix for fluid flow for ct time step
AA=S+dt*D;

ntot=ntot_u+ntot_p;

TotMat=sparse(ntot,ntot);
TotMat(1:ntot_u,1:ntot_u)=K;
TotMat(ntot_u+1:ntot,ntot_u+1:ntot)=-AA;
TotMat(ntot_u+1:ntot,1:ntot_u)=-C';
TotMat(1:ntot_u,ntot_u+1:ntot)=-C;


eq_free=[eq_free_u;ntot_u+eq_free_p];

% time
t_k(1)=0.;
% history of pp at the sphere center
hist_pp_1(1)=pp_o(1);

% Algorithm ct
n_step=501;

% time loop -- 
for i=2:n_step;
    
    t_k(i)=t_k(i-1)+dt;
    
    Dp=0.*pp_o;   % at nodes for fluid flow
    D_U=0.*U_o;
     
    D_x=[D_U; Dp];
    
    po_lhs=dt*D(eq_free_p,:)*pp_o(:);
    
    P_s=sparse(ntot,1);
    
    P_s(eq_free_p+(ntot_u))=po_lhs;
     
    D_x(eq_free)=TotMat(eq_free,eq_free)\P_s(eq_free);
    
    Dp(eq_free_p)=D_x(eq_free_p+(ntot_u));
    D_U(eq_free_u)=D_x(eq_free_u);
    
    pp_o=pp_o+Dp  ;
    U_o=U_o+D_U;
    
    hist_pp_1(i)=pp_o(1);
    
 %   hist_ur_10(i)=U_o(10);
    
end

%%
 %trisurf(Ien,mesh_fem.FEnode.coor(:,1),mesh_fem.FEnode.coor(:,2), pp_o)
ResA = csvread("PP-Sphere-ohio.csv"); % load the analytical solution vector t,ppat center (dt 0.002 time from 0 to 1

figure(5)
title(' Pore pressure evolution at the sphere center'); hold on;
 plot(t_k,hist_pp_1,'k') ; hold on;
plot(ResA(:,1),ResA(:,2),'r') 
xlabel('time');
ylabel(' pore pressure / applied load');

figure(6)
rel_error = abs(ResA(1:n_step,2)-hist_pp_1')./ResA(1:n_step,2)
plot(t_k,rel_error);
xlabel('time');
ylabel(' relative error on pressure '); 

figure(7) 
 % reshape Usol_u
udisp = [U_o(ID_array(:,1)) U_o(ID_array(:,2))];

plotmesh(mesh.XY,connect,[.2 .2 .2],'w')
hold on;
plotmesh(mesh.XY+udisp*1e3,connect,[.8 .2 .2],'none')

figure(8) 
trisurf(mesh.conn,mesh.XY(:,1),mesh.XY(:,2),full(pp_o))
 