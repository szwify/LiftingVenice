% TEST #3 for the Laplace equation
% Hollow cylinder  with given inner field variable
% (Dirichlet BC)

%%  

R_b=1;
coor=[R_b:0.5:25]';
Ien=[];
for i=1:length(coor)-1;
    Ien(i,:)=[i i+1]; 
end;
% mesh
objN=FEnode(coor); 
mat_zones=ones(length(Ien(:,1)));
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones);

% material properties
propObject=Properties_Laplacian_Isotropic(1,1);

% problem set-up
ProblemType='Laplacian';
Config='Axis';
dof_h=DOF_handle(mesh_fem,1,'Matrix');
obj_mat=PDE_part_matrix(ProblemType,Config,mesh_fem,propObject,dof_h);

% stiffness matrix
tic; 
K=Assembly(obj_mat,3);  
toc;

% Applied Temp on first nodes
Tappl=1;
Faux=-K(2:dof_h.nrow-1,1)*Tappl;
Kaux=K(2:dof_h.nrow-1,2:dof_h.nrow-1);

Ur=Kaux\Faux;
U=[Tappl;Ur;0]; % zero BC at r_outer

%%% analytical solution
T_solu=(log(coor(length(coor)))-log(coor))/log(coor(length(coor)));

plot(coor,T_solu); hold on;
plot(coor,U,'.r');
title('Field variable profile');
