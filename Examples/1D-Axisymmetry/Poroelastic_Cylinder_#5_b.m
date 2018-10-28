% TEST:: POROELASTICITY  - Staggered Algorithm
% CYLINDRICAL CORE  retrieved
%  MODE-1 Loading
%  at t=0+=eps : initial condition of the numerical model set-up from the
%  analytical solution.
%  s_rr(r=1)=-p*   , and pore-pressure and displacement profile from the
%  inverse laplace transform of the analytical solution p(r=1)=0
%  See Detournay & Cheng 1993 for the Analytical solution
%%
% simple 1d mesh
R_b=1.;
coor=[0.:0.05:R_b]';
Ien=[];
for i=1:length(coor)-1;
    Ien(i,:)=[i i+1];
end;

% Mesh object
objN=FEnode(coor);
mat_zones=ones(length(Ien(:,1)));
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones);

Config='Axis';
dof_h=DOF_handle(mesh_fem,1,'Matrix');  % the dof_handle for each problems is the same in axis

%%% MATERIAL PROPERTIES
% all stiffness in MPa,

% Poroelastic parameters
k=5e3;       % Drained Bulk modulus
g=4.565e3;   % shear modulus
b=0.8;       % Biot's coefficient
K_hyd=1e-3; % hydraulic conductivity
perm=1e-8; % intrinsic perm
k_f=2.2e3;     % fluid bulk moduli
nu_f=perm/K_hyd; % fluid visco
phi=0.2;

M=8.58e3;     % Biot Modulus
N=1/(1/M-phi/k_f);

%propObject=Properties_PoroElastic_Isotropic(k,g,b,M,K_hyd,1.); % poroelastici mat. with rho=1

propObject=Properties_PoroElastic_Isotropic(k,g,b,N,phi,perm,nu_f,k_f,1.); % poroelastici mat. with rho=1

%%% Construction of the different Finite element matrix on the mesh
% Storage (i.e. Mass type)  for the term in d p/ dt  in the fluid cont. eq.
% Permeability   for the term in K_hyd M \nabla^2 p  in the fluid cont. eq.
%
% Stiffness matrix (drained coefficient for the elasticity eq.

% Storage matrix...
obj_mass=PDE_part_matrix('Mass',Config,mesh_fem,propObject,dof_h);
tic;
S=Assembly(obj_mass,4);  % unit coefficient
toc;
% permeability matrix P   (multiply by M)
obj_perm=PDE_part_matrix('Laplacian',Config,mesh_fem,Extract_HydroProp(propObject),dof_h);
P=propObject.M*(Assembly(obj_perm,4));

% dof_handle  the  for source term in fluid cont. equation
dof_h_l=DOF_handle(mesh_fem,1,'Load');
% the boundary conditions for the fluid cont eq. are : zero pp at r=1
% (cylinder surface).
Imp_pp=[dof_h.nrow 1 0];

%%% Elasticity problem.....
ProblemType='Elasticity';
Imp_displacement=[1 1 0  ]; % zero displacement @ cylinder center
Boundary_loads=[dof_h.nrow 1 -1]; %

%%% Initial stress field (compression negative)
mySig_o=-ones(length(coor)-1,3);
%%% The drained problem stiffness matrix
obj_elas=Elasticity_Block('Axis',mesh_fem,propObject,Imp_displacement,Boundary_loads,mySig_o*0);
[K,dof_aux]=BuildStiffness(obj_elas);   %% Assemble stiffness matrix...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      INITIAL CONDITIONS
%      UN_DRAINED PROBLEM
% AT t=0+, the initial compressive stress field  is instanneously released
% we first solve an undrained problem to get the pore-pressure field and
% displacemnt at t=0+

[Stress_u,Strain_u,AvgCoor]=Stress_And_Strain(obj_elas,zeros(dof_h.nrow,1));
%%%% get the corresponding pore pressure profile via Skempton coefficient
% variation of mean stress
 
% We impose the analytical solution at a time t=eps=1e-4 as initial
% condition of the numericals scheme.
%%%%% COMPUTE ANALYTICAL SOLUTION USING NUMERICAL INVERSE LAPLACE TRANSFORM
t_o=1e-4;
for i=1:length(coor)
    p_ini(i)=gavsteh(@(s)Cylinder_Poro_LaplaceSolution_Mode1_p(propObject,s,coor(i)),t_o,12);

    ur_ini(i)=gavsteh(@(s)Cylinder_Poro_LaplaceSolution_Mode1_ur(propObject,s,coor(i)),t_o,12);
end


for i=1:mesh_fem.Nelts
    pp_o_e(i)=gavsteh(@(s)Cylinder_Poro_LaplaceSolution_Mode1_p(propObject,s,AvgCoor(i)),t_o,12);
end
pp_o_e=pp_o_e';

% these are the initial conditions for the further evolution
U_o=ur_ini';
U_o(1)=0;


Stress_o=Stress_u; % we don't care about stresses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOME matrices to switch fromt elt to node values
% matrix from elt to node values...:: only work for constant elt size
Elt_to_Node=zeros(dof_h.nrow,mesh_fem.Nelts);
Elt_to_Node(1,1)=1;
for i=2:mesh_fem.Nelts,
    Elt_to_Node(i,i)=0.5;
    Elt_to_Node(i,i-1)=0.5;
end
Elt_to_Node(dof_h.nrow,mesh_fem.Nelts)=1;

% Matrice to switch from  Node to Elt value (simple averaging), ct elt size
Node_to_Elt=zeros(mesh_fem.Nelts,dof_h.nrow);
for i=1:dof_h.nrow-1,
    Node_to_Elt(i,i)=0.5;
    Node_to_Elt(i,i+1)=0.5;
end

pp_o=p_ini'; %Elt_to_Node*pp_o_e; % initial pp field on nodes

pp_o(dof_h.nrow)=0; % put last node to zero for further evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% IMPLICIT (staggered algorithm)....
dt=0.0005;  % ct time step
%----- matrix for fluid flow for ct time step
AA=S+dt*P;

% history of pp at the cylinder center
hist_pp_1(1)=pp_o(1);
% history of dispalcement at the node 10
hist_ur_10(1)=U_o(10);
% time
t_k(1)=0.;

% Algorithm ct
beta=0.8;   % relaxation coefficient
eps_p=1e-4; % cvgence test on pp increment
eps_u=1e-4; % cvgence test on U increment
iter_max=30;

for i=2:300;
    
    t_k(i)=t_k(i-1)+dt;
    
    Dp_k_e=0.*pp_o_e; %-- at the middle of the element (for elasticity)
    Dp_k=0.*pp_o;   % at nodes for fluid flow
    
    D_U=0.*U_o;
    
    % START loop for convergence
    %   Fixed point iteration and staggered coupling between fluid cont. and
    %   drained elastic system
    
    cvgence=0;
    k=0;
    
    while ( (k<iter_max) && (cvgence~=1) )
        
        k=k+1;
        
        % Increment of pore pressure from cont. equation
        % compute source term due to increment of volumetric strain
        [Ds,Dstrain,AvgCoor]=Stress_And_Strain(obj_elas,D_U);
        D_vol_Strain=3*mean(Dstrain(:,1:3)')';
        
        Source_k=-(propObject.b)*(propObject.M)*D_vol_Strain; % source is -b M Deps_v
        obj_source=PDE_part_load_body('Laplacian',Config,mesh_fem,propObject,Source_k,dof_h_l);
        [Fs]=Assembly(obj_source,4);
        
        % solution of pore pressure at previous time is a load    -dt*P*p_o
        % except nodes for r=1 where p=0
        po_load=-dt*P(1:dof_h.nrow-1,1:dof_h.nrow-1)*pp_o(1:dof_h.nrow-1);
        %        if (i>2)
        sol=AA(1:dof_h.nrow-1,1:dof_h.nrow-1)\(Fs(1:dof_h.nrow-1)+po_load);
        Dp_k_new=[sol ;0.];
        
        % relaxation
        
        Dp_k=beta*Dp_k_new+(1-beta)*Dp_k;
        Dp_k_e=Node_to_Elt*Dp_k; % back to elt  value.
        
        % SOLVE for Elasticity increment
        % drained problem with current pp field
        
        Sig_k=-(propObject.b)*[Dp_k_e Dp_k_e Dp_k_e];
        obj_elas=Elasticity_Block('Axis',mesh_fem,propObject,Imp_displacement,Boundary_loads,Sig_k);
        [D_U_k]=Solve(obj_elas);
        
        % relaxation
        D_U=beta*D_U_k+(1-beta)*D_U;
        
        if ( (max(abs((Dp_k_new-Dp_k)./Dp_k_new)) < eps_p)  && (max(abs((D_U-D_U_k)./D_U)) < eps_u ) )
            cvgence=1;
              disp('convergence AT its ');disp(k);
            %  disp(norm(Dp_k_new-Dp_k));
        end
        
    end;
    
    pp_o=pp_o+Dp_k_new;
    pp_k_e=pp_o_e+Dp_k_e;
    U_o=U_o+D_U;
    
    hist_pp_1(i)=pp_o(1);
    hist_ur_10(i)=U_o(10);
    
   %plot(coor,pp_o,'.-'); hold on;
    
end;

%%

%%%%% COMPUTE ANALYTICAL SOLUTION USING NUMERICAL INVERSE LAPLACE TRANSFORM
for i=1:length(t_k)
    p_A_1(i)=gavsteh(@(s)Cylinder_Poro_LaplaceSolution_Mode1_p(propObject,s,0.),t_k(i),12);
    ur_A_10(i)=gavsteh(@(s)Cylinder_Poro_LaplaceSolution_Mode1_ur(propObject,s,coor(10)),t_k(i),12);
end


%%%% 
figure(1)
plot(t_k*propObject.c_dif,p_A_1 ); hold on;
plot(t_k*propObject.c_dif,hist_pp_1,'r.')
title('Pore Pressure evolution at cylinder center');

figure(2)
plot(t_k*propObject.c_dif,ur_A_10); hold on;
plot(t_k*propObject.c_dif,hist_ur_10,'r.')
title('Radial evolution evolution at nodes 10');

% limits

