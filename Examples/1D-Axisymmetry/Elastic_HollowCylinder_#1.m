% TEST #1: Quasi-Static ELASTICITY
% Isotopric Hollow CYLINDER
% Inner AND Outer pressure
% Lame solution
% 

%% Geometry and mesh
R_b=1.; % outer radius
R_i=0.2;% inner radius
coor=[R_i:0.02:R_b]';
% connectivity table
Ien=[];
for i=1:length(coor)-1;
    Ien(i,:)=[i i+1]; 
end;
% Mesh object
objN=FEnode(coor); 
mat_zones=ones(length(Ien(:,1))); %single material
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones); % linear FE

% MATERIAL PROPERTIES
% all stiffness in MPa, 
k=5e3;
g=4.1e3; 

propObject=Properties_Elastic_Isotropic(k,g,1.);

% Elasticity problem  1D plane-strain axisymmetry 
% ProblemType='Elasticity';
Config='Axis';
dof_h=DOF_handle(mesh_fem,1,'Matrix');  % the dof_handle

Imp_displacement=[   ]; % no imposed displacement 

% inner and outer pressure  in MPa
p_inner=0.2; 
p_outer=1;
%  BOUNDARY LOADS descriptor 
Boundary_loads=[1 1  p_inner ;dof_h.nrow 1   -p_outer ];

% outer   force  (0d "hands"integration)
F_outer=2*pi*R_b*Boundary_loads(2,3); %'by hand' intg for axis
% inner   force (0d "hands"integration)
F_inner=2*pi*R_i*Boundary_loads(1,3);

% Load Vector 
F=zeros(dof_h.nrow,1);
F(dof_h.nrow)=F_outer;
F(1)=F_inner;

% no Initial stress field 
mySig_o=-zeros(length(coor)-1,3); 

% create elasticity Block
obj_elas=Elasticity_Block('Axis',mesh_fem,propObject,Imp_displacement,...
    Boundary_loads,mySig_o);

[K,dof_aux]=BuildStiffness(obj_elas);   %  Assemble stiffness matrix...

% SOLUTION
U=K\F;

%  GET stress and strain  (ct over elememt)
[Stress ,Strain ,AvgCoor]=Stress_And_Strain(obj_elas,U);

% Analytical solution (see Theory of Elasticity Timoshenko 1935 pg 57)
sr=(1./(AvgCoor.^2))*((R_i^2)*(R_b^2)*(p_outer-p_inner))...
    /(R_b^2-R_i^2)+(p_inner*R_i^2-p_outer*R_b^2)/(R_b^2-R_i^2);
st=-(1./(AvgCoor.^2))*((R_i^2)*(R_b^2)*(p_outer-p_inner))/...
    (R_b^2-R_i^2)+(p_inner*R_i^2-p_outer*R_b^2)/(R_b^2-R_i^2);
% Comparisons Numerical / Analytical 
% for radial stress
plot(AvgCoor,sr,'r'); hold on;
plot(AvgCoor,Stress(:,1),'.')
% for hoop stress
plot(AvgCoor,st,'r'); hold on;
plot(AvgCoor,Stress(:,2),'.')
title('Radial and orthoradial stress profile');

