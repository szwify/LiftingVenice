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

ne_t=length(connect(:,1))

objN=FEnode(the_coor); % obj FEnode

Ien = connect;

% add a FEM On it and interpolation
mat_zones=ones(length(Ien(:,1))); %     single material
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones); % linear FE
% 2D plane-strain axisymmetry
Config='Axis';

% MATERIAL PROPERTIES
% all stiffness in MPa,

k=4.2e3;
g=3.1e3;
b=0.8;
phi=0.2;
k_f=2.2e3;
M=5.e3;
overN=1./M-phi/k_f;N=1/overN
ku=k+b^2*M
perm =1/3.6127e+03;mu_f=1; rho=1;

propObject=Properties_PoroElastic_Isotropic(k,g,b,N,phi,perm,mu_f,k_f,rho);

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


%%% fixed pressure.
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


%%%% SOLUTION OF THE UNDRAINED RESPONSE.
% solve the undrained response from imposed initial stress field.
Boundary_loads =[];
%
% Initial stress field,
mySig_o= zeros(mesh_fem.Nelts,4);  % 4 components in 2D axi-symmetry (srr, szz, srz, stt
mySig_o(:,1)=1;
mySig_o(:,2)=1;
mySig_o(:,4)=1;

% Elastic Undrained Problem at t=0;
propObject_Undrained = Properties_Elastic_Isotropic(propObject.k_u,propObject.g,propObject.rho);
% create elasticity Block
obj_elas_u=Elasticity_Block(Config,mesh_fem,propObject_Undrained,Imp_displacement,...
    Boundary_loads,mySig_o);
Usol_u=full(Solve(obj_elas_u));

%
AnalyticSolUr= @(r) -(r*(1-2*propObject.nu_u)/(2.*propObject.g*(1.+propObject.nu_u)));

% radial displacement at z=0 and at at r=0 -> should be the same....
%
figure(2)
plot(mesh_fem.FEnode.coor(klt_z,1), Usol_u(obj_elas_u.DOF_handle.ID_array(klt_z,1)),'ok' )
hold on
plot(mesh_fem.FEnode.coor(klt_r,2), Usol_u(obj_elas_u.DOF_handle.ID_array(klt_r,2)),'*k' )
hold on
plot(mesh_fem.FEnode.coor(klt_r,2),AnalyticSolUr(mesh_fem.FEnode.coor(klt_r,2)),'-r');

% reshape Usol_u
udisp = [Usol_u(obj_elas_u.DOF_handle.ID_array(:,1)) Usol_u(obj_elas_u.DOF_handle.ID_array(:,2))];

% plot deformed mesh
figure(3)
plotmesh(the_coor,connect,[.2 .2 .2],'w')
hold on;
plotmesh(the_coor+udisp*1e3,connect,[.8 .2 .2],'none')

% stress
[StressG,StrainG,AvgCoor]=Stress_And_Strain(obj_elas_u,Usol_u,'Gauss');

[Stress,Strain,AvgCoor]=Stress_And_Strain(obj_elas_u,Usol_u);

% stress -should be all relaxed to zero !

% undrained Porre pressure response
Pp_u=-mean((Stress(:,[1 2 4])-mySig_o(:,[1 2 4])),2)*propObject.B

Pp_u_u=mean(Pp_u)
%%%%%%%%%%%%%%%%%%%%%%%%% Staggered Solution of the time evolution problem
%%%%%%%%%%%%%%%%%%%%%%%%%
% implicit time - integration

% prepare elastic sub-problem matrices
obj_elas=Elasticity_Block('Axis',mesh_fem,propObject,Imp_displacement,Boundary_loads,mySig_o*0);
[K,dof_aux]=BuildStiffness(obj_elas);   %% Assemble stiffness matrix...

[eq_free_u,fix_nonZero,eq_fix_u]=PrepareDisplacementBC(obj_elas)
 % prepare diffusion sub-problem matrices

% Storage matrix...
obj_mass=PDE_part_matrix('Mass',Config,mesh_fem,propObject,dof_h_p);
tic;
S=Assembly(obj_mass,3);  % unit coefficient
toc;
% permeability matrix L   (multiplied by M)
obj_perm=PDE_part_matrix('Laplacian',Config,mesh_fem,Extract_HydroProp(propObject),dof_h_p);
L=propObject.M*(Assembly(obj_perm,3));
% dof_handle  the  for source term in fluid cont. equation
dof_h_l=DOF_handle(mesh_fem,1,'Load');

%%
%%%%% coupling matrix   coef. is propObject.b
obj_coupl=PDE_part_coupling_matrix(Config,mesh_fem,mesh_fem,propObject,dof_h_e,dof_h_p);
C= (Assembly(obj_coupl,3));

 
%%
%%%%%%%
% t=0+ solution is the undrained response
pp_o=ones(mesh_fem.FEnode.n_tot,1)*Pp_u_u;
pp_o(eq_fix_p)=0.
U_o = Usol_u;

% 
% figure(3)
% trisurf(Ien,mesh_fem.FEnode.coor(:,1),mesh_fem.FEnode.coor(:,2),pp_o)


%%   Staggered solution ----- 
%%%%%%
dt=0.001;
%----- matrix for fluid flow for ct time step
AA=S+dt*L;
% time
t_k(1)=0.;
% history of pp at the sphere center
hist_pp_1(1)=pp_o(1);

% Algorithm ct
beta=.8;   % relaxation coefficient
eps_p=1e-3; % cvgence test on pp increment
eps_u=1e-3; % cvgence test on U increment
iter_max=40;

% time loop
for i=2:50;
    
    t_k(i)=t_k(i-1)+dt;
    
    Dp_k_e=0.*pp_o_e; %-- at the middle of the element (for elasticity)
    Dp_k=0.*pp_o;   % at nodes for fluid flow
    Dp_k_new=Dp_k;
    
    D_U=0.*U_o;
    D_U_k=D_U;
    % staggered coupling loop
    cvgence=0;
    k=0;
    
    % Solution of pore pressure at previous time is a load    -dt*P*p_o
    po_lhs=-dt*L(eq_free_p,:)*pp_o(:);
    
    while ( (k<iter_max) && (cvgence~=1) )
        
        k=k+1;
        % Increment of pore pressure from cont. equation
        % compute source term due to increment of volumetric strain
         
        Fs=-(propObject.M)*(C'*D_U);
        % 
        Dp_k_new(eq_free_p)=(AA(eq_free_p,eq_free_p)\(Fs(eq_free_p)+po_lhs));
        
        % relaxation
        Dp_k=beta*Dp_k_new+(1-beta)*Dp_k;
        
        Dp_k_e=NodeToElt*Dp_k; %   to elt  value.
        % SOLVE for Elasticity increment
        % drained problem with current Delta pp field
 
         Fc=C*Dp_k;
         D_U_k(eq_free_u)=K(eq_free_u,eq_free_u)\Fc(eq_free_u);%full(Solve(obj_elas,K));
         % relaxation 
         D_U=beta*D_U_k+(1-beta)*D_U;

        if ( (max(abs((Dp_k_new(eq_free_p)-Dp_k(eq_free_p))./Dp_k_new(eq_free_p))) < eps_p)  && (max(abs((D_U(eq_free_u)-D_U_k(eq_free_u))./D_U(eq_free_u))) < eps_u ) )
            cvgence=1;
            disp(['convergence AT its : ' num2str(k)]);
              disp(norm((Dp_k_new-Dp_k) ));
        end
%           disp(['its ',k]);
%           disp(max(abs((Dp_k_new-Dp_k)./Dp_k_new))) ;
%           disp((max(abs((D_U-D_U_k)./D_U))));
    end
    
    pp_o=pp_o+Dp_k_new;
    pp_k_e=pp_o_e+Dp_k_e;
    U_o=U_o+D_U;
    
    hist_pp_1(i)=pp_o(1);
    % hist_ur_10(i)=U_o(10);
    
end

 %trisurf(Ien,mesh_fem.FEnode.coor(:,1),mesh_fem.FEnode.coor(:,2), pp_o)
 figure(4)
 plot(t_k,hist_pp_1)

%scatter3(AvgCoor(:,1),AvgCoor(:,2),NodeToElt*pp_o)
 