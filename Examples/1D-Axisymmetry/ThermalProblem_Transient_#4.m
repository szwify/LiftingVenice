% TEST
% Wellbore in an infinite medium with given constant inner Temperature 
% Time-dependent thermal diffusion problem

%% simple 1d mesh
R_b=1;
h=0.25;
coor=[R_b:h:25]';
Ien=[];
for i=1:length(coor)-1;
    Ien(i,:)=[i i+1]; 
end;
 
% Mesh fem object
objN=FEnode(coor); 
mat_zones=ones(length(Ien(:,1)));
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones);

% properties
rhoC=10;
cond=1;
propObject=Properties_Laplacian_Isotropic(cond,rhoC);

ProblemType='Laplacian';
Config='Axis';

dof_h=DOF_handle(mesh_fem,1,'Matrix');
obj_cond=PDE_part_matrix(ProblemType,Config,mesh_fem,propObject,dof_h);
% stiffness matrix
tic; 
K=Assembly(obj_cond,3);  
toc;
% mass matrix
obj_mass=PDE_part_matrix('Mass',Config,mesh_fem,propObject,dof_h);
tic; 
M=Assembly(obj_mass,3);  
toc;

% Applied Temp on first nodes
Tappl=1;
Faux=-K(2:dof_h.nrow-1,1)*Tappl;
Kaux=K(2:dof_h.nrow-1,2:dof_h.nrow-1);

%%% steady state solution.....
Ur=Kaux\Faux;
U=[Tappl;Ur;0];

%%%% Implicit algorithm for the transient evolution
dt=1./((cond/rhoC)/h^2); % the basic CFL dt

T_o=zeros(length(dof_h.ID_array),1);
T_o(1)=Tappl;

%%%%% 
T_k=T_o;
tk=0;
allTs=zeros(length(dof_h.ID_array),3000);

for i=1:3000;
    
    tk=tk+dt;
    Faux=-dt*K*T_k;    
    KK=(M+dt*K);
    DeltaT=KK(2:dof_h.nrow-1,2:dof_h.nrow-1)\Faux(2:dof_h.nrow-1);
    T_k=T_k+[0;DeltaT;0];
    % store solution
    allTs(:,i)=T_k;
%    if (i/10 == 1 );
%        plot(coor,T_k);hold on;
%   end;
    
end
% Steady-state analytical solution
T_solu=(log(coor(length(coor)))-log(coor))/log(coor(length(coor)));

plot(coor,T_solu); hold on;

% the  few time selected for comparisons of T profile with the analytical
% solution 
kt=[2;50;200;600;2000;3000];
the_ts=kt*dt;
 
 
%  LOAD Mathematica Solution for transient
load HollowCylinderTransientThermal.dat;

plot(coor,HollowCylinderTransientThermal); hold on;
plot(coor(1:2:length(coor)),allTs(1:2:length(coor),kt),'.')
title('Temperature profile at different time');

