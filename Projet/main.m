L = 120;
Htot = 60;
H1 = 25;
H2 = 40;
H3 = 45;
H4 = 60;

Hs = [0,5;
    -H1,5;
    -H2,0.5; 
    -H3,0.5;
    -H4,5;
    ];

[vert,etri, tria,tnum] = stratimesher(120, Hs);
the_coor=vert; 
connect=tria;

figure(1);
plotmesh(the_coor,connect,[.2 .2 .2],'w')

mesh = FEmesh(the_coor,connect);


% MATERIAL PROPERTIES
% all stiffness in MPa, 

k=4.2e3;
g=3.1e3; 

L_elas=Elastic_Isotropic_Stiffness(k,g,'Axis');

% Elasticity problem  1D plane-strain axisymmetry 
% ProblemType='Elasticity';
% Config='Axis';

% dof_h=DOF_handle(mesh_fem,2,'Matrix');  % the dof_handle

% impose z displacement bottom  (dof 2) 
% impose r displacement left  (dof 1)

klt_down = find(mesh.XY(:,2)==-H4);
Imp_displacement=[klt_down  ones(length(klt_down),1)  ones(length(klt_down),1) ];

klt_left = find(mesh.XY(:,1)==0);
Imp_displacement=[Imp_displacement;
    klt_left  zeros(length(klt_left),1)  ones(length(klt_left),1) ];

klt_right = find(mesh.XY(:,1)==L);
Imp_displacement=[Imp_displacement;
    klt_right  zeros(length(klt_right),1)  ones(length(klt_right),1) ];

 
hold on;
plot(mesh.XY(klt_down,1),mesh.XY(klt_down,2),'or')
hold on;
plot(mesh.XY(klt_left,1),mesh.XY(klt_left,2),'or')
hold on;
plot(mesh.XY(klt_right,1),mesh.XY(klt_right,2),'or')





%%%% create boundary loads matrix ....format is dof node1 intensity1 node2
%%%% intensity 2   (where node 1 node 2 define a segment)
%%%% 
deltap=1;
ktl_rad=find(round(sqrt(mesh.XY(:,1).^2.+mesh.XY(:,2).^2),4)==Radius )
Boundary_loads =[];
hold on;
plot(mesh.XY(ktl_rad,1),mesh.XY(ktl_rad,2),'og')
hold on;
%    switch from fr,0=  to----- fx fy    - to be recheck
for seg=1:length(ktl_rad)-1
    
    theta_1 = atan2(mesh.XY(ktl_rad(seg),2),mesh.XY(ktl_rad(seg),1));
    theta_2 = atan2(mesh.XY(ktl_rad(seg+1),2),mesh.XY(ktl_rad(seg+1),1));
    
    Boundary_loads =[Boundary_loads ;...
                   1 ktl_rad(seg) deltap*cos(theta_1) ktl_rad(seg+1) deltap*cos(theta_2); ...
                   2 ktl_rad(seg) deltap*sin(theta_1) ktl_rad(seg+1) deltap*sin(theta_2); 
                  ]; 
end

%Boundary_loads =[];
%
% no Initial stress field 
mySig_o= zeros(mesh.Nelts,4);  % 4 components in 2D axi-symmetry (srr, szz, srz, stt
% mySig_o(:,1)=1; 
% mySig_o(:,2)=1;
% mySig_o(:,4)=1;

 %% assemble stiffness
 proplist={L_elas};
 [K,ID_array]=AssembleMatrix(mesh,'Axis','Elasticity',proplist,3);

[Fload]=AssembleVectorBoundaryTerm(mesh,'Axis','BoundaryLoads',Boundary_loads,ID_array,2);
%%
Fbody=Fload*0.;
% solution of the system

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


nu =  (1-2*g/(3*k))/(2+2*g/(3*k));
%
AnalyticSolUr= @(r) (r*(1-2*nu)/(2.*g*(1.+nu)));

% radial displacement at z=0 and at at r=0 -> should be the same....
% 
figure(2)
plot(mesh.XY(klt_y,1), Usol(ID_array(klt_y,1)),'ok' )
hold on
plot(mesh.XY(klt_r,2), Usol(ID_array(klt_r,2)),'*k' )
hold on
plot(mesh.XY(klt_r,2),AnalyticSolUr(mesh.XY(klt_r,2)),'-r');

%% Reshape UsolF
udisp =reshape(Usol,[length(ID_array(:,1)) 2 ]);
udisp(:,1)=Usol(ID_array(:,1));
udisp(:,2)=Usol(ID_array(:,2));


% plot deformed mesh
figure(3)
plotmesh(the_coor,connect,[.2 .2 .2],'w')
hold on;
plotmesh(the_coor+udisp*1e3,connect,[.8 .2 .2],'none')

%% Stress
[StressG,StrainG,AvgCoor]=Compute_Stress_And_Strain(mesh,'Axis',proplist,3,Usol,ID_array,mySig_o)
 
%[Stress,Strain,AvgCoor]=Stress_And_Strain(obj_elas,Usol);

% stress should be uniform & equal to applied pressure (srr,szz,stt) and
% srz=0

 
