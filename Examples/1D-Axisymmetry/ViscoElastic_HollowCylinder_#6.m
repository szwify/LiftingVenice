% TEST #6: 
% Isotopric VISCOELASTIC KELVIN-VOIGT  HOLLOW CYLINDER  
%  
% ct presssure applied at t=0 on outer cylinder
% 
% There is 2 concentric cylinders,
% only the inner cylinder is viscoelastic
% we compare with the infinite solution (which has moduli k*kv/(k+kv) )
%%
R_b=5;
coor=[[1.:0.2:R_b]';  ];

objN=FEnode(coor);
Ien=[];
for i=1:length(coor)-1;
    Ien(i,:)=[i i+1];
    ID_array(i)=[i]';
    ID_array(i+1)=[i+1]';
end;

% vector of index for elements in the 3 different material zones
% zone 1
 

% List of material zones and properties
mat_zones=ones(length(Ien(:,1)),1);
R_1=3;
kt1=find(coor(Ien(:,1),1)<=R_1);
kt2=find(coor(Ien(:,1),1)>R_1);

mat_zones(kt2)=2;

k=5.5e3;
g=4.1e3;

propObject(1)=Properties_Elastic_Isotropic(k,g,1);
propObject(2)=Properties_Elastic_Isotropic(20e3,12e3,1);

k_v=2e3;
g_v=1.3e3;
vis_eta=1e3;

CV=Properties_Viscoelastic_Isotropic(propObject(1),vis_eta,k_v,g_v);

% mesh_fem object with material zone, mesh etc.
myPk=1;
mesh_fem=Mesh_with_a_FEM(myPk,objN,Ien,mat_zones);
dof_h=DOF_handle(mesh_fem,1,'Matrix');  % the dof_handle
Imp_displacement=[];

 % ct applied pressure
 p_outer=1; % 1 Mpa

Boundary_loads=[ dof_h.nrow 1  -p_outer ];

% outer   force  (0d "hands"integration)
F_outer=2*pi*R_b*Boundary_loads(1,3); %'by hand' intg for axis
 
% no Initial stress field 
mySig=-zeros(length(coor)-1,3); 


% SOLUTION @ t=0+
obj_elas=Elasticity_Block('Axis',mesh_fem,propObject,Imp_displacement,...
    Boundary_loads,mySig);
[K_1,dof_h]=BuildStiffness(obj_elas);

% Load Vector :: ct
F=zeros(dof_h.nrow,1);
F(dof_h.nrow)=F_outer;

[U_o]=K_1\F;

% GET stress and strain  (ct over elememt)
[Stress_o ,Strain_o, AvgCoor]=Stress_And_Strain(obj_elas,U_o);

% SOLUTION @ t=inf
%%%%% solution @ t=inf
propObject_inf=propObject;
propObject_inf(1)=Properties_Elastic_Isotropic(k*k_v/(k+k_v),g*g_v/(g+g_v),1);
obj_elas_inf=Elasticity_Block('Axis',mesh_fem,propObject_inf,Imp_displacement,...
    Boundary_loads,mySig);
[K_inf,dof_h]=BuildStiffness(obj_elas_inf);

[U_inf]=K_inf\F;
% GET stress and strain  (ct over elememt)
[Stress_inf ,Strain_inf ,AvgCoor]=Stress_And_Strain(obj_elas_inf,U_inf);


% for radial stress
plot(AvgCoor,Stress_o(:,1),'-k');       hold on;
plot(AvgCoor,Stress_inf(:,1),'-r');     hold on;
% for hoop stress
plot(AvgCoor,Stress_o(:,2),'-k');       hold on;
plot(AvgCoor,Stress_inf(:,2),'-r') ;    hold on;

title('Radial and orthoradial stress profile');


%%     
%      LOOP FOR VISCOELASTICITY
% 
dt=0.1;             % the scheme is second order accurate, so we can have large time-steps.

imax=18;

D_Sig_star=mySig*0;
Stress_i=Stress_o;
Strain_i=Strain_o;
Strain_visc_i=0.*Strain_o(kt1,:);


for i=1:imax,
   
    % can test both integration (no ageing and ageing)
%    [D_Strain,D_Stress,D_Strain_visc]=ViscoelasticStep_NonAgeing(dt,obj_elas,1,CV,Strain_visc_i,Stress_i,Strain_visc_i(:,1)*0,D_Sig_star);
    [D_Strain,D_Stress,D_Strain_visc]=ViscoelasticStep_Ageing(dt,obj_elas,1,CV,CV,Strain_visc_i,Stress_i,Strain_visc_i(:,1)*0,D_Sig_star);    

    Strain_visc_i=Strain_visc_i+D_Strain_visc ;
    Stress_i=Stress_i + D_Stress;
    Strain_i=Strain_i +D_Strain;
    
% Comparisons Numerical / Analytical
% for radial stress
    plot(AvgCoor,Stress_i(:,1),'.b'); hold on;

    plot(AvgCoor,Stress_i(:,2),'.b'); hold on;

end

% the algorithm stress should tend toward the infinite one 


