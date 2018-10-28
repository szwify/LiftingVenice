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

%%
ne_t=length(connect(:,1))

objN=FEnode(the_coor);  % object FEnode

Ien = connect;          % connectivity

% add a FEM On it and interpolation
mat_zones=ones(length(Ien(:,1))); %     single material
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones); % linear FE
% 2D axisymmetry
Config='Axis';
%%
% MATERIAL PROPERTIES
% all stiffness in MPa, ohio Sandstone

k=8.4e3; % elastic drained bulk modulus
g=6.8e3;  % shear modulus
b=0.707692; % Biot coefficient 
phi=0.19;  % porosity 
k_f=2.2e3;  % fluid bulk modulus
M=9.18478e3; % Biot Modulus
overN=1./M-phi/k_f;N=1/overN  % Biot intrinsic Modulus (fluid independent)
ku=k+b^2*M     % undrained bulk modulus
perm =0.137549e-3; % permeability value adjusted such that c= k /(mu S)=1 , where S is storage coefficient
mu_f=1; rho=1;

propObject=Properties_PoroElastic_Isotropic(k,g,b,N,phi,perm,mu_f,k_f,rho);

%%
% elasticity dof
dof_h_e=DOF_handle(mesh_fem,2,'Matrix');  % the dof_handle
% pressure dof
dof_h_p=DOF_handle(mesh_fem,1,'Matrix');  % the dof_handle for each problems is the same in axis

% impose z displacement bottom  (dof 2)
% impose r displacement left  (dof 1)

klt_z = find(mesh_fem.FEnode.coor(:,2)==0.);
Imp_displacement=[klt_z  2*ones(length(klt_z),1)  0.*ones(length(klt_z),1) ];

klt_r = find(round(mesh_fem.FEnode.coor(:,1),4) ==0.);
Imp_displacement=[Imp_displacement;
    klt_r  1*ones(length(klt_r),1)  0.*ones(length(klt_r),1) ];

%%% fixed pressure on the outside
% nodes on radius
ktl_rad=find(round(sqrt(mesh_fem.FEnode.coor(:,1).^2.+mesh_fem.FEnode.coor(:,2).^2),4)==Radius )

eq_fix_p=[dof_h_p.ID_array(ktl_rad)  ];
eq_free_p=setdiff(dof_h_p.ID_array(:),eq_fix_p(:));

hold on;
plot(mesh_fem.FEnode.coor(klt_z,1),mesh_fem.FEnode.coor(klt_z,2),'or')
hold on;
plot(mesh_fem.FEnode.coor(klt_r,1),mesh_fem.FEnode.coor(klt_r,2),'or')
hold on;
plot(mesh_fem.FEnode.coor(ktl_rad,1),mesh_fem.FEnode.coor(ktl_rad,2),'*b')

%%
%%%%  Preparation of the different FE matrix
Boundary_loads =[];
%
% Initial stress field,
mySig_o= zeros(mesh_fem.Nelts,4);  % 4 components in 2D axi-symmetry (srr, szz, srz, stt
mySig_o(:,1)=1;
mySig_o(:,2)=1;
mySig_o(:,4)=1;
 
%
%%%%%%%%%%%%%%%%%%%%%%%%% Staggered Solution of the time evolution problem
%%%%%%%%%%%%%%%%%%%%%%%%%
% implicit time - integration

% prepare elastic sub-problem matrices
obj_elas=Elasticity_Block('Axis',mesh_fem,propObject,Imp_displacement,Boundary_loads,mySig_o);
[K,dof_aux]=BuildStiffness(obj_elas);   %% Assemble stiffness matrix...
[Fbody]=BuildBodyLoad(obj_elas) ;

[eq_free_u,fix_nonZero,eq_fix_u]=PrepareDisplacementBC(obj_elas)


% prepare diffusion sub-problem matrices
% Storage matrix...
obj_mass=PDE_part_matrix('Mass','Axis',mesh_fem,propObject,dof_h_p);
tic;
S=(1/(propObject.M))*Assembly(obj_mass,3);  % with unit coefficient
toc;
% permeability matrix L   
obj_perm=PDE_part_matrix('Laplacian','Axis',mesh_fem,Extract_HydroProp(propObject),dof_h_p);
L=(Assembly(obj_perm,3));

%%%%% coupling matrix   coupling coef. is propObject.b (isotropy only)
obj_coupl=PDE_part_coupling_matrix('Axis',mesh_fem,mesh_fem,propObject,dof_h_e,dof_h_p);
C= (Assembly(obj_coupl,3));

%% Solving directly for the undrained solution at t=0 ! 
dt= 0.;
AA=S+dt*L;

ntot=dof_h_e.nrow+dof_h_p.nrow;

TotMat=sparse(ntot,ntot);
TotMat(1:dof_h_e.nrow,1:dof_h_e.nrow)=K;
TotMat(dof_h_e.nrow+1:ntot,dof_h_e.nrow+1:ntot)=-AA;
TotMat(dof_h_e.nrow+1:ntot,1:dof_h_e.nrow)=-C';
TotMat(1:dof_h_e.nrow,dof_h_e.nrow+1:ntot)=-C;

Ftot=sparse(ntot,1);
Ftot(1:dof_h_e.nrow)=Fbody;

eq_free=[eq_free_u;(dof_h_e.nrow)   + dof_h_p.ID_array(:)];%+eq_free_p]; We dont fix pore pressure for the undrained response.

Undrained_sol=sparse(ntot,1);
Undrained_sol(eq_free) = TotMat(eq_free,eq_free)\Ftot(eq_free);

pp_o=sparse(dof_h_p.nrow,1);
pp_o=Undrained_sol((dof_h_e.nrow)   + dof_h_p.ID_array(:));
eq_u =sort([dof_h_e.ID_array(:,1),dof_h_e.ID_array(:,2)]);

U_o=sparse(dof_h_e.nrow,1);
U_o(eq_u) = Undrained_sol(eq_u);

AnalyticSolUr= @(r) -(r*(1-2*propObject.nu_u)/(2.*propObject.g*(1.+propObject.nu_u)));

% radial displacement at z=0 and at at r=0 -> should be the same....
%
figure(2)
plot(mesh_fem.FEnode.coor(klt_z,1), U_o(obj_elas.DOF_handle.ID_array(klt_z,1)),'ok' )
hold on
plot(mesh_fem.FEnode.coor(klt_r,2), U_o(obj_elas.DOF_handle.ID_array(klt_r,2)),'*k' )
hold on
plot(mesh_fem.FEnode.coor(klt_r,2),AnalyticSolUr(mesh_fem.FEnode.coor(klt_r,2)),'-r');

% reshape Usol_u
udisp = [U_o(obj_elas.DOF_handle.ID_array(:,1)) U_o(obj_elas.DOF_handle.ID_array(:,2))];

% plot deformed mesh
figure(3)
plotmesh(the_coor,connect,[.2 .2 .2],'w')
hold on;
plotmesh(the_coor+udisp*1e3,connect,[.8 .2 .2],'none')

% t=0+ solution is the undrained response
%pp_o=ones(mesh_fem.FEnode.n_tot,1)*Pp_u_u;
pp_o(eq_fix_p)=0.   ; % NOW we fix the pore pressure to zero. 

figure(4) 
trisurf(Ien,mesh_fem.FEnode.coor(:,1),mesh_fem.FEnode.coor(:,2),full(pp_o))
 
%%

%%%%%%
dt=0.002;
%----- matrix for fluid flow for ct time step
AA=S+dt*L;

ntot=dof_h_e.nrow+dof_h_p.nrow;
TotMat=sparse(ntot,ntot);
TotMat(1:dof_h_e.nrow,1:dof_h_e.nrow)=K;
TotMat(dof_h_e.nrow+1:ntot,dof_h_e.nrow+1:ntot)=-AA;
TotMat(dof_h_e.nrow+1:ntot,1:dof_h_e.nrow)=-C';
TotMat(1:dof_h_e.nrow,dof_h_e.nrow+1:ntot)=-C;

eq_free=[eq_free_u;(dof_h_e.nrow)+eq_free_p];

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
    
    po_lhs=dt*L(eq_free_p,:)*pp_o(:);
    
    P_s=sparse(ntot,1);
    
    P_s(eq_free_p+(dof_h_e.nrow))=po_lhs;
     
    D_x(eq_free)=TotMat(eq_free,eq_free)\P_s(eq_free);
    
    Dp(eq_free_p)=D_x(eq_free_p+(dof_h_e.nrow));
    D_U(eq_free_u)=D_x(eq_free_u);
    
    pp_o=pp_o+Dp  ;
    U_o=U_o+D_U;
    
    hist_pp_1(i)=pp_o(1);
    
 %   hist_ur_10(i)=U_o(10);
    
end

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



%scatter3(AvgCoor(:,1),AvgCoor(:,2),NodeToElt*pp_o)
 

 figure(7) 
 % reshape Usol_u
udisp = [U_o(obj_elas.DOF_handle.ID_array(:,1)) U_o(obj_elas.DOF_handle.ID_array(:,2))];

plotmesh(the_coor,connect,[.2 .2 .2],'w')
hold on;
plotmesh(the_coor+udisp*1e3,connect,[.8 .2 .2],'none')



