% TEST #2 : Quasi-Static Elasticity (plane-strain, axisymmetry)
% 3 co-axial cylinders with an initial stress state in the middle cylinder
% only
% 
%% non-unform mesh, 3 zones
coor=[[1.:0.025:1.5]'; [1.51:0.1:5]'; [5.5:0.5:16]' ];

objN=FEnode(coor); 
Ien=[];
for i=1:length(coor)-1;
    Ien(i,:)=[i i+1];
    ID_array(i)=[i]';
    ID_array(i+1)=[i+1]';
end;

% vector of index for elements in the 3 different material zones
% zone 1   
R_12=1.5;
R_23=11;
kt1=find(coor(Ien(:,1),1)<R_12);
kt2=find(coor(Ien(:,1),1)>=R_12 & coor(Ien(:,1),1)<R_23 );
kt3=find(coor(Ien(:,1),1)>=R_23) ;

% List of material zones and properties
mat_zones=ones(length(Ien(:,1)),1);
mat_zones(kt1)=1;
mat_zones(kt2)=2;
mat_zones(kt3)=3;

propObject(1)=Properties_Elastic_Isotropic(20e3,15e3,1);
propObject(2)=Properties_Elastic_Isotropic(5.5e3,4.1e3,1);
propObject(3)=Properties_Elastic_Isotropic(200e3,175e3,1);

% mesh_fem object with material zone, mesh etc.
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones);

% stress free 
Imp_displacement=[];
Boundary_loads=[];
 
%  Unit & isotropic initial stress state in middle cylinder
mySig=-ones(length(coor)-1,3);   %% sigma_o equal -1
mySig(kt1,:)=0*mySig(kt1,:);
mySig(kt3,:)=0*mySig(kt3,:);

obj=Elasticity_Block('Axis',mesh_fem,propObject,...
    Imp_displacement,Boundary_loads,mySig);

U=Solve(obj);

% Average stress and strain within elts.
[Stress,Strain,AvgCoor]=Stress_And_Strain(obj,U);
%plot(AvgCoor,Stress);

load CoaxCylindersIniStressCenter.dat % analytical solution 
plot(CoaxCylindersIniStressCenter(:,1),...
    CoaxCylindersIniStressCenter(:,2),'-')
hold on
plot(coor,U,'.r' );
title('Radial displacement profile');
