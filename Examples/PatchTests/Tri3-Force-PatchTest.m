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
 
mesh=FEmesh(the_coor,connect);

% MATERIAL PROPERTIES
% all stiffness in MPa, 

k=4.2e3;
g=3.1e3; 
 L_elas=Elastic_Isotropic_Stiffness(k,g,'PlaneStrain');
 
 
% Elasticity problem  1D plane-strain axisymmetry 
 Config='2D';

 
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
  
% no Initial stress field 
mySig_o= zeros(mesh.Nelts,3); 

%  
proplist={L_elas};
 
 [K,ID_array]=AssembleMatrix(mesh,'2D','Elasticity',proplist,3);
 
[Fload]=AssembleVectorBoundaryTerm(mesh,'2D','BoundaryLoads',Boundary_loads,ID_array,3);

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

 
[Stress,Strain,AvgCoor]=Compute_Stress_And_Strain(mesh,'2D',proplist,3,Usol,ID_array,mySig_o,'Gauss')

[Stress,Strain,AvgCoor]=Compute_Stress_And_Strain(mesh,'2D',proplist,3,Usol,ID_array,mySig_o)

Stress



udisp =reshape(Usol,[length(ID_array(:,1)) 2 ]);
udisp(:,1)=Usol(ID_array(:,1));
udisp(:,2)=Usol(ID_array(:,2));


% plot deformed mesh
figure(3)
plotmesh(the_coor,connect,[.2 .2 .2],'w')
hold on;
plotmesh(the_coor+udisp*1e3,connect,[.8 .2 .2],'none')


