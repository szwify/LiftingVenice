%GEOMETRY
 
%constants
Radius = 0.1;
Lplug = 0.01;
Lend = 0.6;
Lfrac = Lend-Lplug;
Hfrac = -0.003;
nrad = 15;
 
%nodes
thetas=linspace(0.,pi/2.,nrad);
node = [   
      Radius*cos(thetas)' , Radius*sin(thetas)' ;
      0. Lend;
      Lend Lend;
      Lend 0.;    
      Radius+Lplug, 0;
       Radius, Hfrac;
      Radius+Lplug, Hfrac;
      Lend, Hfrac;
      ];
 
%edges
edge = [];
for e=1:nrad+3,
edge = [ edge ; e e+1 ];
end
edge = [ edge ; e+1 1 ];
n1=length(edge);  % we have finish material 1.

nn_res = length(node);
 
edge= [edge; 1 20 ; 20 21; 21 19;]; % add segment for the boundary of 2
edge = [edge; 21 22;22 18 ] ;   %  add segment for boundary 3
%edge= [edge; 18 19; 19 20; 20 21; 21 18];
 
nr_edge = length(edge);
 
%% 
for i=1:length(Part{1}),
    plot(node(edge(i,1:2),1),node(edge(i,1:2),2),'-.'); hold on;
end,
%%
%

 Part={[1:n1],[n1 n1+1 n1+2 n1+3],[n1-1 n1+3 n1+4 n1+5] };
 
%%
% Part={[1:nr_edge];[nr_edge+1 nr_edge+2 nr_edge+3 nr_y+nr_x+nr_y+1:nr_edge];
%         ,
%  [ length(edge)-2 length(edge)-1 length(edge) (nr_y+nr_x):-1:nr_y+1];
%  
%  []; };
 %%
%Part={[1:nr_edge];[1:19];[20:23];[21:27]; };
 
%-- call mesh-gen.
[vert,etri, tria,tnum] = refine2(node,edge,Part,[],2.5) ; % do not touch
the_coor=vert;
connect=tria;
mat_id=tnum;  % 1 is the reservoir
 
%plotmesh(the_coor,connect,mat_id,[0.2 0.2 0.2],['y' ,'r', 'g'])
 

elt_mat_1 = find(mat_id(:)==1);
elt_mat_2 = find(mat_id(:)==2);
elt_mat_3 = find(mat_id(:)==3);


patch('faces',connect(elt_mat_1,:),'vertices',the_coor, ...
                            'facecolor','r', ...
                        'edgecolor','k') ;
    
                        hold on; 

patch('faces',connect(elt_mat_2,:),'vertices',the_coor, ...
                            'facecolor','g', ...
                        'edgecolor','k') ;
    
                        hold on; 
                        
patch('faces',connect(elt_mat_3,:),'vertices',the_coor, ...
                            'facecolor','b', ...
                        'edgecolor','k') ;
    
                        hold on; 
                        
                        
                        axis image off;

%PlainStrain
ne_t=length(connect(:,1))
 
objN=FEnode(the_coor);  % object FEnode
Ien = connect;          % connectivity
 
% add a FEM On it and interpolation
mat_zones=mat_id; %     single material
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones); % linear FE



% Plane strain model
Config='PlaneStrain';
 
 


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
 
%rock reservoir
propObject_res=Properties_PoroElastic_Isotropic(k,g,b,N,phi,perm*0.0001,mu_f,k_f,rho);
 
% plug
propObject_plug=Properties_PoroElastic_Isotropic(k,g,b,N,phi,perm,mu_f,k_f,rho);
 
% fracture
propObject_frac=Properties_PoroElastic_Isotropic(k,g,b,N,phi,perm*100000,mu_f,k_f,rho);
 
 
%BOUNDARY CONDITIONS
 
% elasticity dof
dof_h_e=DOF_handle(mesh_fem,2,'Matrix');  % the dof_handle (2 DOF for displacement)
% pressure dof
dof_h_p=DOF_handle(mesh_fem,1,'Matrix');  % the dof_handle for each problems is the same in axis
 
% impose x at x=0
% impose y at y=Hfrac
 
% find 
klt_y = find(mesh_fem.FEnode.coor(:,2)==Hfrac);
Imp_displacement=[klt_y  2*ones(length(klt_y),1)  0.*ones(length(klt_y),1) ];
 
 
klt_x = find(round(mesh_fem.FEnode.coor(:,1),4) ==0.);
Imp_displacement=[Imp_displacement; 
    klt_x  1*ones(length(klt_x),1)  0.*ones(length(klt_x),1)];

% Boundary loads.

%%%% create boundary loads matrix ....format is dof node1 intensity1 node2
%%%% intensity 2   (where node 1 node 2 define a segment)
%%%% 
deltap=1;
ktl_rad=find(round(sqrt(mesh_fem.FEnode.coor(:,1).^2.+mesh_fem.FEnode.coor(:,2).^2),4)==Radius )
Boundary_loads =[];
hold on;
plot(mesh_fem.FEnode.coor(ktl_rad,1),mesh_fem.FEnode.coor(ktl_rad,2),'og')
hold on;
%    switch from fr,0=  to----- fx fy    - to be recheck
for seg=1:length(ktl_rad)-2
    
    theta_1 = atan2(mesh_fem.FEnode.coor(ktl_rad(seg),2),mesh_fem.FEnode.coor(ktl_rad(seg),1));
    theta_2 = atan2(mesh_fem.FEnode.coor(ktl_rad(seg+1),2),mesh_fem.FEnode.coor(ktl_rad(seg+1),1));
    
    Boundary_loads =[Boundary_loads ;...
                   1 ktl_rad(seg) deltap*cos(theta_1) ktl_rad(seg+1) deltap*cos(theta_2); ...
                   2 ktl_rad(seg) deltap*sin(theta_1) ktl_rad(seg+1) deltap*sin(theta_2); 
                  ]; 
end

 Boundary_loads =[Boundary_loads ; 
                 1 1 deltap 20 deltap];
             
            
%ktl_rad=find(round(sqrt(mesh_fem.FEnode.coor(:,1).^2.+mesh_fem.FEnode.coor(:,2).^2),4)==Radius )

eq_fix_p=[dof_h_p.ID_array(ktl_rad)  ];
eq_free_p=setdiff(dof_h_p.ID_array(:),eq_fix_p(:)); 

%%
% no Initial stress field 
mySig_o= zeros(mesh_fem.Nelts,3);  % 4 components in 2D axi-symmetry (srr, szz, srz, stt
% mySig_o(:,1)=1; 
% mySig_o(:,2)=1;
 

ListProp=[propObject_res,propObject_plug, propObject_frac];

%%

% elasticity dof
dof_h_e=DOF_handle(mesh_fem,2,'Matrix');  % the dof_handle
% pressure dof
dof_h_p=DOF_handle(mesh_fem,1,'Matrix');  % the dof_handle for each problems is the same in axis

% create elasticity Block
obj_elas=Elasticity_Block(Config,mesh_fem,ListProp,Imp_displacement,...
    Boundary_loads,mySig_o);
%%
%%elastic problem to check that everything is fine 
 Usol=full(Solve(obj_elas));
% Reshape UsolF
udisp = [Usol(obj_elas.DOF_handle.ID_array(:,1)) Usol(obj_elas.DOF_handle.ID_array(:,2))
];

% plot deformed mesh
figure(3)
plotmesh(the_coor,connect,[.2 .2 .2],'w')
hold on;
plotmesh(the_coor+udisp*1e3,connect,[.8 .2 .2],'none') % deformesh
   
%% 

[K,dof_aux]=BuildStiffness(obj_elas);   %% Assemble stiffness matrix...
[Fbody]=BuildBodyLoad(obj_elas);
[eq_free_u,fix_nonZero,eq_fix_u]=PrepareDisplacementBC(obj_elas)
[Fload]=BuildBoundaryLoad(obj_elas);
% prepare diffusion sub-problem matrices
% Storage matrix...
obj_mass=PDE_part_matrix('Mass',Config,mesh_fem,ListProp,dof_h_p);
tic;
S=(1/(propObject_res.M))*Assembly(obj_mass,3);  % with unit coefficient
toc;
 

prop_hydro=[Extract_HydroProp(propObject_res),Extract_HydroProp(propObject_plug),Extract_HydroProp(propObject_frac)];

% permeability matrix L   
obj_perm=PDE_part_matrix('Laplacian',Config,mesh_fem,prop_hydro,dof_h_p);
L=(Assembly(obj_perm,3));

%%%%% coupling matrix   coupling coef. is propObject.b (isotropy only)
obj_coupl=PDE_part_coupling_matrix(Config,mesh_fem,mesh_fem,ListProp,dof_h_e,dof_h_p);
C= (Assembly(obj_coupl,3));


%% 
dt= 0.;
AA=S+dt*L;

ntot=dof_h_e.nrow+dof_h_p.nrow;

TotMat=sparse(ntot,ntot);
TotMat(1:dof_h_e.nrow,1:dof_h_e.nrow)=K;
TotMat(dof_h_e.nrow+1:ntot,dof_h_e.nrow+1:ntot)=-AA;
TotMat(dof_h_e.nrow+1:ntot,1:dof_h_e.nrow)=-C';
TotMat(1:dof_h_e.nrow,dof_h_e.nrow+1:ntot)=-C;

Ftot=sparse(ntot,1);
Ftot(1:dof_h_e.nrow)=Fbody+Fload;
 

eq_free=[eq_free_u;(dof_h_e.nrow) + dof_h_p.ID_array(:)];   %+eq_free_p]; We dont fix pore pressure for the undrained response.

Undrained_sol=sparse(ntot,1);
Undrained_sol(eq_free) = TotMat(eq_free,eq_free)\Ftot(eq_free);

pp_o=sparse(dof_h_p.nrow,1);
pp_o=Undrained_sol((dof_h_e.nrow) + dof_h_p.ID_array(:));
eq_u =sort([dof_h_e.ID_array(:,1),dof_h_e.ID_array(:,2)]);

U_o=sparse(dof_h_e.nrow,1);
U_o(eq_u) = Undrained_sol(eq_u);
          
% Reshape UsolF
udisp = [Usol(obj_elas.DOF_handle.ID_array(:,1)) Usol(obj_elas.DOF_handle.ID_array(:,2))
];

% plot deformed mesh
figure(3)
plotmesh(the_coor,connect,[.2 .2 .2],'w')
hold on;
plotmesh(the_coor+udisp*1e3,connect,[.8 .2 .2],'none') % deformesh


% t=0+ solution is the undrained response
%pp_o=ones(mesh_fem.FEnode.n_tot,1)*Pp_u_u;
eq_fix_p
%pp_o(eq_fix_p)=0.   ; % NOW we fix the pore pressure to zero. 

figure(4) 
trisurf(Ien,mesh_fem.FEnode.coor(:,1),mesh_fem.FEnode.coor(:,2),full(pp_o))


%%
eq_fix_p=
