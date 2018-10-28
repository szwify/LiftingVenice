% RESERVOIR LAYER IN AXISYMMETRY SURROUNDED BY a material with 100 times
% lower permeability. 

%%%% geometry and mesh
Depth = 10.;
H_res = 1;
Ly = 2*Depth+H_res;
Lr = 2*Ly;

% mes reservoir
nr_y =3; 
hy=H_res/nr_y;

node =[];
for i=0:nr_y
    node=[node; 0. Depth+i*hy];
end

nr_x=80;
hx=Lr/nr_x;
for i=1:nr_x
    node=[node; i*hx Depth+H_res];
end

for i=1:nr_y
    node=[node; Lr Depth+H_res-i*hy];
end

for i=1:nr_x-1
    node=[node; Lr-i*hx Depth];
end

nn_res = length(node);

edge = [];
for e=1:nn_res-1,
   edge = [edge; e e+1 ]
end
edge= [edge; e+1 1 ];
nr_edge = length(edge);

node=[node; 0 0 ;
    Lr 0.; 
     0. Ly ;
     Lr Ly;]

edge =[edge; 1 nn_res+1;
    nn_res+1 nn_res+2; 
    nn_res+2 nr_y+nr_x+nr_y+1;
     nr_y+1 length(node)-1;
     length(node)-1 length(node);
     length(node) nr_y+nr_x+1;
     ]

Part={[1:nr_edge],
    [nr_edge+1 nr_edge+2 nr_edge+3 nr_y+nr_x+nr_y+1:nr_edge]
        ,
 [ length(edge)-2 length(edge)-1 length(edge) (nr_y+nr_x):-1:nr_y+1]  
 };


%-- call mesh-gen.
[vert,etri, tria,tnum] = refine2(node,edge,Part,[],2.5) ; % do not touch
the_coor=vert;
connect=tria;
mat_id=tnum;  % 1 is the reservoir

plotmesh(the_coor,connect,[.2 .2 .2],'w')

% boundary
patch('faces',edge(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;


%%
ne_t=length(connect(:,1))

objN=FEnode(the_coor);  % object FEnode


% add a FEM On it and interpolation
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,connect,mat_id); % linear FE

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

propObject_res=Properties_PoroElastic_Isotropic(k,g,b,N,phi,perm,mu_f,k_f,rho);


% surrounding material.
propObject_sur=Properties_PoroElastic_Isotropic(k,g,b,N,phi,perm*0.01,mu_f,k_f,rho);

%%

% elasticity dof
dof_h_e=DOF_handle(mesh_fem,2,'Matrix');  % the dof_handle
% pressure dof
dof_h_p=DOF_handle(mesh_fem,1,'Matrix');  % the dof_handle for each problems is the same in axis

% impose r at r=0
% impose z at z=0

%  imp_disp [ node number,  dof number, value of imp. displacement] 
% dof number 1:: u_x, 2:: u_y
klt_z = find(mesh_fem.FEnode.coor(:,2)==0.);
Imp_displacement=[klt_z  2*ones(length(klt_z),1)  0.*ones(length(klt_z),1) ];

%%
klt_r = find(round(mesh_fem.FEnode.coor(:,1),4) ==0.);
Imp_displacement=[Imp_displacement;
    klt_r  1*ones(length(klt_r),1)  0.*ones(length(klt_r),1) ];

%%%%  Preparation of the different FE matrix
Boundary_loads =[ ];
%

% Initial stress field,   % 3 in plane strain (sxx,syy,sxy)
mySig_o= zeros(mesh_fem.Nelts,4);  % 4 components in 2D axi-symmetry (srr, szz, srz, stt
mySig_o(:,1)=1;
mySig_o(:,2)=1;
mySig_o(:,4)=1;
 

%%
listProp=[propObject_res,propObject_sur,propObject_sur];

obj_elas=Elasticity_Block(Config,mesh_fem,listProp,Imp_displacement,Boundary_loads,mySig_o);
[K,dof_aux]=BuildStiffness(obj_elas);   %% Assemble stiffness matrix...


[eq_free_u,fix_nonZero,eq_fix_u]=PrepareDisplacementBC(obj_elas)


% prepare diffusion sub-problem matrices
% Storage matrix...
obj_mass=PDE_part_matrix('Mass','Axis',mesh_fem,listProp,dof_h_p);
tic;
S=(1/(propObject_sur.M))*Assembly(obj_mass,3);  % with unit coefficient
toc;

%%
prop_hydro=[Extract_HydroProp(propObject_res),Extract_HydroProp(propObject_sur),Extract_HydroProp(propObject_sur)];

% permeability matrix L   
obj_perm=PDE_part_matrix('Laplacian','Axis',mesh_fem,prop_hydro,dof_h_p);
L=(Assembly(obj_perm,3));

%%%%% coupling matrix   coupling coef. is propObject.b (isotropy only)
obj_coupl=PDE_part_coupling_matrix('Axis',mesh_fem,mesh_fem,listProp,dof_h_e,dof_h_p);
C= (Assembly(obj_coupl,3));
%%


