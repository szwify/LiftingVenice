% implementation of TRI3 element for axi-symmetry elasticity
%    
%     SPHERE
% 
% ---- there is still something weird ... check load integration in
% axisymmetry
%  - works ok with inital stress state. 
%  we apply an unit external pressure here
% it looks like there is a slight inbalance in the applied load (on z=0
% axis)

Lc=1.;
node=[ 0. 0. ;Lc 0.; 20*Lc  0. ; 20*Lc  20*Lc ; 0.  20*Lc ];
nrad=4;
edge = [];
for e=1:nrad,
edge = [ edge ; e e+1 ];
end
edge = [ edge ; nrad+1 1 ];
 
%------------------------------------------- call mesh-gen.
[vert,etri, tria,tnum] = refine2(node,edge,[],[],.5) ; % do not touch
  
  the_coor=vert; 
  connect=tria;
 
  figure(1);
  plotmesh(the_coor,connect,[.2 .2 .2],'w')
 %%
%    % boundary 
%    patch('faces',edge(:,1:2),'vertices',node, ...
%         'facecolor','w', ...
%         'edgecolor',[.1,.1,.1], ...
%         'linewidth',1.5) ;

 
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

%%
% Elasticity problem  1D plane-strain axisymmetry 
% ProblemType='Elasticity';
Config='PlaneStrain';

dof_h=DOF_handle(mesh_fem,2,'Matrix');  % the dof_handle

% impose y displacement bottom  (dof 2) 
% impose x displacement left  (dof 1)

klt_y = find(mesh_fem.FEnode.coor(:,2)==0. && mesh_fem.FEnode.coor(:,1)>=Lc);
Imp_displacement=[klt_z  2*ones(length(klt_z),1)  0.*ones(length(klt_z),1) ];

klt_x = find(round(mesh_fem.FEnode.coor(:,1),4) ==0.);

Imp_displacement=[Imp_displacement; 
    klt_r  1*ones(length(klt_r),1)  0.*ones(length(klt_r),1) ];
 
%%
hold on;
plot(mesh_fem.FEnode.coor(klt_z,1),mesh_fem.FEnode.coor(klt_z,2),'or')
hold on;
plot(mesh_fem.FEnode.coor(klt_r,1),mesh_fem.FEnode.coor(klt_r,2),'or')

%%
%%%% create boundary loads matrix ....format is dof node1 intensity1 node2
%%%% intensity 2   (where node 1 node 2 define a segment)
%%%% 

ktl_rad=find(round(sqrt(mesh_fem.FEnode.coor(:,1).^2.+mesh_fem.FEnode.coor(:,2).^2),4)==Radius )
Boundary_loads =[];
hold on;
plot(mesh_fem.FEnode.coor(ktl_rad,1),mesh_fem.FEnode.coor(ktl_rad,2),'og')
hold on;
%    switch from fr,0=  to----- fx fy    - to be recheck
for seg=1:length(ktl_rad)-1
    
    theta_1 = atan2(mesh_fem.FEnode.coor(ktl_rad(seg),2),mesh_fem.FEnode.coor(ktl_rad(seg),1));
    theta_2 = atan2(mesh_fem.FEnode.coor(ktl_rad(seg+1),2),mesh_fem.FEnode.coor(ktl_rad(seg+1),1));
    
    Boundary_loads =[Boundary_loads ;...
                   1 ktl_rad(seg) cos(theta_1) ktl_rad(seg+1) cos(theta_2); ...
                   2 ktl_rad(seg) sin(theta_1) ktl_rad(seg+1) sin(theta_2); 
                  ]; 
end

%Boundary_loads =[];
%
% no Initial stress field 
mySig_o= zeros(mesh_fem.Nelts,4);  % 4 components in 2D axi-symmetry (srr, szz, srz, stt
% mySig_o(:,1)=1; 
% mySig_o(:,2)=1;
% mySig_o(:,4)=1;

% create elasticity Block
obj_elas=Elasticity_Block(Config,mesh_fem,propObject,Imp_displacement,...
    Boundary_loads,mySig_o);

% [K,dof_aux]=BuildStiffness(obj_elas); 
% 
% [F]=BuildBoundaryLoad(obj_elas); 

Usol=full(Solve(obj_elas)); 


%
AnalyticSolUr= @(r) (r*(1-2*propObject.nu)/(2.*propObject.g*(1.+propObject.nu)));

% radial displacement at z=0 and at at r=0 -> should be the same....
% 
figure(2)
plot(mesh_fem.FEnode.coor(klt_z,1), Usol(obj_elas.DOF_handle.ID_array(klt_z,1)),'ok' )
hold on
plot(mesh_fem.FEnode.coor(klt_r,2), Usol(obj_elas.DOF_handle.ID_array(klt_r,2)),'*k' )
hold on
plot(mesh_fem.FEnode.coor(klt_r,2),AnalyticSolUr(mesh_fem.FEnode.coor(klt_r,2)),'-r');

% Reshape UsolF
udisp = [Usol(obj_elas.DOF_handle.ID_array(:,1)) Usol(obj_elas.DOF_handle.ID_array(:,2))
];

% plot deformed mesh
figure(3)
plotmesh(the_coor,connect,[.2 .2 .2],'w')
hold on;
plotmesh(the_coor+udisp*1e3,connect,[.8 .2 .2],'none')

% Stress
[StressG,StrainG,AvgCoor]=Stress_And_Strain(obj_elas,Usol,'Gauss');
 
[Stress,Strain,AvgCoor]=Stress_And_Strain(obj_elas,Usol);

% stress should be uniform & equal to applied pressure (srr,szz,stt) and
% srz=0

Stress
