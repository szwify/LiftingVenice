%
%     Poroelastic SPHERE - sudden isotropic loading at t=0, drained BC at
%     r=1
%
%      
%%%% geometry and mesh
Radius = 1.;
nrad = 80;
thetas=linspace(0.,pi/2.,nrad);
node = [ 0., 0. ;
    Radius*cos(thetas)' , Radius*sin(thetas)'
    ];

edge = [];
for e=1:nrad,
    edge = [ edge ; e e+1 ];
end
edge = [ edge ; nrad+1 1 ];

%-- call mesh-gen.
[vert,etri, tria,tnum] = refine2(node,edge,[],[],0.08) ; % do not touch
the_coor=vert;
connect=tria;

plotmesh(the_coor,connect,[.2 .2 .2],'w')
% boundary
patch('faces',edge(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;


mesh = FEmesh(the_coor,connect);

% 2D axisymmetry
Config='Axis';
%
% MATERIAL PROPERTIES
% all stiffness in MPa, ohio Sandstone

k=8.4e3; % elastic drained bulk modulus
g=6.8e3;  % shear modulus
b=0.707692; % Biot coefficient 
phi=0.19;  % porosity 
k_f=2.2e3;  % fluid bulk modulus
M=9.18478e3; % Biot Modulus
overN=1./M-phi/k_f;N=1/overN  % Biot intrinsic Modulus (fluid independent)
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

%  SURFACES FOR BOUNDARY CONDITIONS

% impose z displacement bottom  (dof 2)
% impose r displacement left  (dof 1)

klt_z = find(mesh.XY(:,2)==0.);
Imp_displacement=[klt_z  2*ones(length(klt_z),1)  0.*ones(length(klt_z),1) ];

klt_r = find(round(mesh.XY(:,1),4) ==0.);
Imp_displacement=[Imp_displacement;
    klt_r  1*ones(length(klt_r),1)  0.*ones(length(klt_r),1) ];

%%% fixed pressure on the outside
% nodes on radius
ktl_rad=find(round(sqrt(mesh.XY(:,1).^2.+mesh.XY(:,2).^2),4)==Radius )

hold on;
plot(mesh.XY(klt_z,1),mesh.XY(klt_z,2),'or')
hold on;
plot(mesh.XY(klt_r,1),mesh.XY(klt_r,2),'or')
hold on;
plot(mesh.XY(ktl_rad,1),mesh.XY(ktl_rad,2),'*b')

% Initial stress and external tractions
%%%%  
Boundary_loads =[];
%
% Initial stress field,
mySig_o= zeros(mesh.Nelts,4);  % 4 components in 2D axi-symmetry (srr, szz, srz, stt
mySig_o(:,1)=-1;
mySig_o(:,2)=-1;
mySig_o(:,4)=-1;
 

%% different matrices

% elasticity
proplist={L_d};
[K,ID_array]=AssembleMatrix(mesh,'Axis','Elasticity',proplist,3);

%[Fload]=AssembleVectorBoundaryTerm(mesh,'Axis','BoundaryLoads',Boundary_loads,ID_array,3);
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
 