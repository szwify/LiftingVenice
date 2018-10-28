% implementation of TRI3 element for plane-strain elasticity

% IMPOSED STress PATCH TEST (e_xx =1, e_yy=1 exy=0 )
% create a  Tri3 grid  of 8 elements for this patch test
ymax=2.;
xmax=2.;

    node = [                % list of xy "node" coordinates
        0, 0                % outer square
        2, 0
        2, 2
        0, 2 
        ] ;
    
    edge = [                % list of "edges" between nodes
        1, 2                % outer square 
        2, 3
        3, 4
        4, 1 
          ] ;

%------------------------------------------- call mesh-gen.
   [vert,etri, tria,tnum] = refine2(node,edge,[],[],1.0) ; % do not touch
  
  the_coor=vert; 
% just perturb the location of the middle node
pert=0.2;
the_coor(9,:)=the_coor(9,:)+[random('norm',0.,pert) random('norm',0.,pert)]

vert=the_coor

    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;

connect=tria;
 
ne_t=length(connect(:,1))
% simple grid plot (note should work for all linear-like mesh) .... very
% naive, lots of duplicates.... should be a loop on mesh edges.
 scatter(the_coor(:,1),the_coor(:,2)); hold on;
 for e=1:ne_t
     line(the_coor(connect(e,:),1),the_coor(connect(e,:),2)); hold on;
 end

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

% impose displacement  everywhere except at the mid point
%block bottom y_dof 

 
% impose displacement .... block bottom y_dof and
kti=find(the_coor(:,2)==0.);
Imp_displacement=[ ];
for i=1:length(kti)
    Imp_displacement =[Imp_displacement; kti(i) 2 0.];
end

% block x_dof for x==0.
kti=find(the_coor(:,1)==0.);
for i=1:length(kti)
    Imp_displacement =[Imp_displacement; kti(i) 1 0.];
end


% top vertical load, we set the load by segment 
ktl=find(the_coor(:,2)==ymax);
Boundary_loads = [ ]; 
 Boundary_loads=[ Boundary_loads;  2 3 1 8 1];
 Boundary_loads=[ Boundary_loads;  2 8 1 4 1];
 
 Boundary_loads=[ Boundary_loads;  1 2 1 7 1];
 Boundary_loads=[ Boundary_loads;  1 7 1 3 1];
 
%  
% %[ dof_dim  node1 value1 node2 value2]
% for i=1:length(ktl)-1
%     % aie, needs to find edge to be general !
%     
%     Boundary_loads=[ Boundary_loads; 2 ktl(i) 1 ktl(i+1) 1];
% 
% end
% 
% % right side horizontal load, we set the load by segment 
% ktl=find(the_coor(:,1)==xmax);
% %[ dof_dim  node1 value1 node2 value2]
% for i=1:length(ktl)-1
%     Boundary_loads=[ Boundary_loads; 1 ktl(i) 1 ktl(i+1) 1];
% end
% 

% no Initial stress field 
mySig_o= zeros(mesh_fem.Nelts,3); 

% create elasticity Block
obj_elas=Elasticity_Block(Config,mesh_fem,propObject,Imp_displacement,...
    Boundary_loads,mySig_o);

[K,dof_aux]=BuildStiffness(obj_elas); 

[F]=BuildBoundaryLoad(obj_elas); 

Usol=full(Solve(obj_elas)); 

[StressG,StrainG,AvgCoor]=Stress_And_Strain(obj_elas,Usol,'Gauss');
 
[Stress,Strain,AvgCoor]=Stress_And_Strain(obj_elas,Usol);

Stress

 

%figure(2)
%voronoi(objN.coor(:,1),objN.coor(:,2))

