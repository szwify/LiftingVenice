% implementation of TRI3 element for axi-symmetry elasticity
%    
%     SPHERE
% 
% ---- there is a bug somewhere in boundary load integration 
%  - works ok with inital stress state.

Radius = 1.;
nrad = 40;
thetas=linspace(0.,pi/2.,nrad); 
node = [ 0., 0. ;
   Radius*cos(thetas)' , Radius*sin(thetas)' 
];
edge = [];
for e=1:nrad,
edge = [ edge ; e e+1 ];
end
edge = [ edge ; nrad+1 1 ];
 
%------------------------------------------- call mesh-gen.
   [vert,etri, tria,tnum] = refine2(node,edge,[],[],0.1) ; % do not touch
  
  the_coor=vert; 
  connect=tria;
 
  plotmesh(the_coor,connect,[.2 .2 .2],'w')
  % boundary 
   patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    
ne_t=length(connect(:,1)) 

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

klt_z = find(mesh.XY(:,2)==0.);
Imp_displacement=[klt_z  2*ones(length(klt_z),1)  0.*ones(length(klt_z),1) ];

klt_r = find(round(mesh.XY(:,1),4) ==0.);
Imp_displacement=[Imp_displacement;
    klt_r  1*ones(length(klt_r),1)  0.*ones(length(klt_r),1) ];
 
hold on;
plot(mesh.XY(klt_z,1),mesh.XY(klt_z,2),'or')
hold on;
plot(mesh.XY(klt_r,1),mesh.XY(klt_r,2),'or')

%%%% create boundary loads matrix ....format is dof node1 intensity1 node2
%%%% intensity 2   (where node 1 node 2 define a segment)
%%%% 
deltap=1;
ktl_rad=find(round(sqrt(mesh.XY(:,1).^2.+mesh.XY(:,2).^2),4)==Radius );
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

Boundary_loads =[];
%
% no Initial stress field 
mySig_o= zeros(mesh.Nelts,4);  % 4 components in 2D axi-symmetry (srr, szz, srz, stt
 mySig_o(:,1)=1; 
 mySig_o(:,2)=1;
 mySig_o(:,4)=1;

 %% assemble stiffness
 
 proplist={L_elas};
 
 [K,ID_array]=AssembleMatrix(mesh,'Axis','Elasticity',proplist,3);
  
[Fbody]=AssembleVectorVolumeTerm(mesh,'Axis','InitialStress',mySig_o,ID_array,3);
 Fload=Fbody*0.;
 
 
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

%
AnalyticSolUr= @(r) (r*(1-2*nu)/(2.*g*(1.+nu)));
%dissplacement at z=0 and at at r=0 -> should be the same....
% 
figure(2)
plot(mesh.XY(klt_z,1), Usol(ID_array(klt_z,1)),'ok' )
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
