%    Plane-strain Hole in a plate - to compare with Kirsh solution 
%     here just checking displacement - internally pressurized hole with
%     unit pressure
% ----   seems to work - needs to check stress, 
% ----      as well as case with deviatoric far field loading....

Radius = 1.;
Lend = 50.;
nrad = 20;

thetas=linspace(0.,pi/2.,nrad); 
node = [   Radius*cos(thetas)' , Radius*sin(thetas)' ;
      0. Lend;
      Lend Lend;
      Lend 0.;
];

edge = [];
for e=1:nrad+2,
edge = [ edge ; e e+1 ];
end
edge = [ edge ; e+1 1 ];


%------------------------------------------- call mesh-gen.
   [vert,etri, tria,tnum] = refine2(node,edge,[],[],5) ; % do not touch
  
  the_coor=vert; 
 vert=the_coor

    figure(1);
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
hold on; 


%% 
connect=tria;
 
ne_t=length(connect(:,1)) 

objN=FEnode(the_coor); % obj FEnode

Ien = connect;

% add a FEM On it and interpolation
mat_zones=ones(length(Ien(:,1))); %single material
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones); % linear FE

% MATERIAL PROPERTIES
% all stiffness in MPa, 

k=4.2e3;
g=3.1e3; 

propObject=Properties_Elastic_Isotropic(k,g,1.);

% Elasticity problem  1D plane-strain axisymmetry 
% ProblemType='Elasticity';
Config='PlaneStrain';

dof_h=DOF_handle(mesh_fem,2,'Matrix');  % the dof_handle

% impose z displacement bottom  (dof 2) 
klt_z = find(mesh_fem.FEnode.coor(:,2)==0.);
Imp_displacement=[klt_z  2*ones(length(klt_z),1)  0.*ones(length(klt_z),1) ];

% impose r displacement left  (dof 1)
klt_r = find(round(mesh_fem.FEnode.coor(:,1),4)==0.);
Imp_displacement=[Imp_displacement;
    klt_r  1*ones(length(klt_r),1)  0.*ones(length(klt_r),1) ];
 
hold on;
plot(mesh_fem.FEnode.coor(klt_z,1),mesh_fem.FEnode.coor(klt_z,2),'or')
hold on;
plot(mesh_fem.FEnode.coor(klt_r,1),mesh_fem.FEnode.coor(klt_r,2),'or')

%Imp_displacement=[Imp_displacement; klt_r(end) 2 0; klt_z(end) 1 0];

%%%% create boundary loads matrix ....format is dof node1 intensity1 node2
%%%% intensity 2   (where node 1 node 2 define a segment)
%%%% 

% ktFy = find(mesh_fem.FEnode.coor(:,2)==Lend);
% Boundary_loads =[];
% %    
% for seg=1:length(ktFy)-1
%     Boundary_loads = [Boundary_loads ;...
%                    2 ktFy(seg) 1 ktFy(seg+1) 1; 
%                    ]; 
% end
% ktFx = find(mesh_fem.FEnode.coor(:,1)==Lend);
% %    
% for seg=1:length(ktFx)-1
%     Boundary_loads = [Boundary_loads ;...
%                    1 ktFx(seg) 1 ktFx(seg+1) 1; 
%                    ]; 
% end

ktl_rad=find(round(sqrt(mesh_fem.FEnode.coor(:,1).^2.+mesh_fem.FEnode.coor(:,2).^2),4)==Radius )
Boundary_loads =[];

%    switch from fr,0=  to----- fx fy    - to be recheck
for seg=1:length(ktl_rad)-1
    
    theta_1 = atan2(mesh_fem.FEnode.coor(ktl_rad(seg),2),mesh_fem.FEnode.coor(ktl_rad(seg),1));
    theta_2 = atan2(mesh_fem.FEnode.coor(ktl_rad(seg+1),2),mesh_fem.FEnode.coor(ktl_rad(seg+1),1));
    
    Boundary_loads = [Boundary_loads ;...
                   1 ktl_rad(seg) cos(theta_1) ktl_rad(seg+1) cos(theta_2); ...
                   2 ktl_rad(seg) sin(theta_1) ktl_rad(seg+1) sin(theta_2); 
                  ]; 
end
%
% s Initial stress field 
mySig_o= zeros(mesh_fem.Nelts,3);  % 4 components in 2D axi-symmetry, 3 in plane-strain
  mySig_o(:,1)=0.3;
  mySig_o(:,2)=0.4;


% create elasticity Block
obj_elas=Elasticity_Block(Config,mesh_fem,propObject,Imp_displacement,...
    Boundary_loads,mySig_o);

[K,dof_aux]=BuildStiffness(obj_elas); 

[F]=BuildBoundaryLoad(obj_elas); 

Usol=full(Solve(obj_elas)); 

% radial displacement at z=0 and at at r=0 -> should be the same....
%

figure(2)
plot(mesh_fem.FEnode.coor(klt_z,1), Usol(obj_elas.DOF_handle.ID_array(klt_z,1)),'ok' )
hold on
plot(mesh_fem.FEnode.coor(klt_r,2), Usol(obj_elas.DOF_handle.ID_array(klt_r,2)),'*k' )
hold on
% analytical solution for the case of wellbore pressure...
myr=linspace(Radius,Lend,40);
ur =1.*Radius./myr/(2*propObject.g);
plot(myr,ur,'-r');

% Need to script things for stress plot 
[StressG,StrainG,AvgCoor]=Stress_And_Strain(obj_elas,Usol,'Gauss');
 
[Stress,Strain,AvgCoor]=Stress_And_Strain(obj_elas,Usol);

Strain 

