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
connect=tria;

% add a FEM On it  
mesh = FEmesh(the_coor,connect);

% MATERIAL PROPERTIES
% all stiffness in MPa, 

k=4.2e3;
g=3.1e3; 

L_elas=Elastic_Isotropic_Stiffness(k,g,'PlaneStrain');


% Elasticity problem  1D plane-strain axisymmetry 
  

% impose z displacement bottom  (dof 2) 
klt_z = find(mesh.XY(:,2)==0.);
Imp_displacement=[klt_z  2*ones(length(klt_z),1)  0.*ones(length(klt_z),1) ];

% impose r displacement left  (dof 1)
klt_r = find(round(mesh.XY(:,1),4)==0.);
Imp_displacement=[Imp_displacement;
    klt_r  1*ones(length(klt_r),1)  0.*ones(length(klt_r),1) ];
 
hold on;
plot(mesh.XY(klt_z,1),mesh.XY(klt_z,2),'or')
hold on;
plot(mesh.XY(klt_r,1),mesh.XY(klt_r,2),'or')

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

ktl_rad=find(round(sqrt(mesh.XY(:,1).^2.+mesh.XY(:,2).^2),4)==Radius )
Boundary_loads =[];

%    switch from fr,0=  to----- fx fy    - to be recheck
for seg=1:length(ktl_rad)-1
    
    theta_1 = atan2(mesh.XY(ktl_rad(seg),2),mesh.XY(ktl_rad(seg),1));
    theta_2 = atan2(mesh.XY(ktl_rad(seg+1),2),mesh.XY(ktl_rad(seg+1),1));
    
    Boundary_loads = [Boundary_loads ;...
                   1 ktl_rad(seg) cos(theta_1) ktl_rad(seg+1) cos(theta_2); ...
                   2 ktl_rad(seg) sin(theta_1) ktl_rad(seg+1) sin(theta_2); 
                  ]; 
end
%
% s Initial stress field 
mySig_o= zeros(mesh.Nelts,3);  % 4 components in 2D axi-symmetry, 3 in plane-strain
  mySig_o(:,1)=0.;
  mySig_o(:,2)=0.;

 %% assemble stiffness
 proplist={L_elas};
 [K,ID_array]=AssembleMatrix(mesh,'2D','Elasticity',proplist,3);
 
   
%%  
[Fload]=AssembleVector(mesh,'2D','BoundaryLoads',Boundary_loads,ID_array,3);

Fbody=Fload;
%% solution of the system

[eq_free,fix_nonZero,eq_fix]=PrepareDirichletBC(Imp_displacement,ID_array);

if (isempty(fix_nonZero))
    Ur=K(eq_free,eq_free)\(Fbody(eq_free)+Fload(eq_free));
else
    eq_fix_nonZero=[];
    for imp=1:length(fix_nonZero)
        eq_fix_nonZero=[eq_fix_nonZero ; ID_array(Imp_displacement(fix_nonZero(imp),1),Imp_displacement(fix_nonZero(imp),2)) ];
    end
    %  disp(eq_fix_nonZero);
    %  disp(size(obj.Imp_displacement(fix_nonZero,3)));
    F_disp=-K(eq_free,eq_fix_nonZero)*obj.Imp_displacement(fix_nonZero,3);
    Ur=K(eq_free,eq_free)\(Fbody(eq_free)+Fload(eq_free)+F_disp);
end

% glue back solution for all nodes
Usol(eq_free)=Ur;
if (isempty(eq_fix)==0)
    Usol(eq_fix)=0.;
end
if (isempty(fix_nonZero)==0)
    Usol(eq_fix_nonZero)=Imp_displacement(fix_nonZero,3);
end

%% 
% radial displacement at z=0 and at at r=0 -> should be the same....
%

figure(2)
plot(mesh.XY(klt_z,1), Usol(ID_array(klt_z,1)),'ok' )
hold on
plot(mesh.XY(klt_r,2), Usol(ID_array(klt_r,2)),'*k' )
hold on
% analytical solution for the case of wellbore pressure...
myr=linspace(Radius,Lend,60);
ur =1.*Radius./myr/(2*g);
plot(myr,ur,'-r');

%%
 [Stress,Strain,AvgCoor]=Compute_Stress_And_Strain(mesh,'2D',proplist,3,Usol,ID_array,mySig_o)
% 
% [Stress,Strain,AvgCoor]=Compute_Stress_And_Strain(mesh,'2D',proplist,3,Usol,ID_array,mySig_o,'Gauss')

  

