%
%
%  Laplacian problem on the Hollow sphere
%  2D axisymmetric modeling
%   Dirichlet Boundary conditions at r=ra and r=rb 

Ra = 0.4;   Ta = 0.3;
Rb= 1.;    Tb = 1.3 ;

nrad = 30;
thetas=linspace(0.,pi/2.,nrad); 
node = [    Ra*cos(thetas)' , Ra*sin(thetas)' ;
        Rb*cos(sort(thetas,'descend'))' , Rb*sin(sort(thetas,'descend'))' 
];
edge = [];
for e=1:2*nrad-1,
edge = [ edge ; e e+1 ];
end
edge = [ edge ; 2*nrad  1 ];
 
%------------------------------------------- call mesh-gen.
[vert,etri, tria,tnum] = refine2(node,edge,[],[],0.05) ; % do not touch
  
  the_coor=vert; 
  connect=tria;
 
  figure(1);
  plotmesh(the_coor,connect,[.2 .2 .2],'w')
 ne_t=length(connect(:,1)) 
 
 mesh = FEmesh(the_coor,connect);
 %
 pK=1;
 
% material properties
propList={1};

% Assemble problem matrix 
  
[K,ID_array]=AssembleMatrix(mesh,'Axis','Laplacian',propList,pK+1);
 
 % Dirichlet boundary conditions.
ktl_ra=find(round(sqrt(mesh.XY(:,1).^2.+mesh.XY(:,2).^2),4)==Ra );
ktl_rb=find(round(sqrt(mesh.XY(:,1).^2.+mesh.XY(:,2).^2),4)==Rb );


eq_fix=[ID_array(ktl_ra) ;ID_array(ktl_rb) ];
 
Timp=[ones(length(ktl_ra),1)*Ta ;ones(length(ktl_rb),1)*Tb];

eq_free=setdiff(ID_array(:),eq_fix(:));


F=-K(eq_free,eq_fix)*Timp; 
solT = K(eq_free,eq_free)\F;



T_sol = zeros(length(ID_array),1 ) ;
T_sol(eq_fix)=Timp;
T_sol(eq_free)=solT;


% Plotting solution and compare with analytical solution

AnalyticSol= @(r) (r.*Ra*Ta - Ra*Rb*Ta - r.*Rb*Tb + Ra*Rb*Tb)./(r.*Ra - r.*Rb);
% 

klt_z = find(mesh.XY(:,2)==0.);
klt_r = find(round(mesh.XY(:,1),4) ==0.);


figure(2)
plot(mesh.XY(klt_z,1), T_sol(ID_array(klt_z,1)),'ok' )
hold on
plot(mesh.XY(klt_r,2), T_sol(ID_array(klt_r,1)),'*k' )
hold on
plot(sort(mesh.XY(klt_r,2)), AnalyticSol(sort(mesh.XY(klt_r,2))),'r.-' )
hold on 
title('Field variable profile');
legend('Numer. along x=0','Numer. along y=0','Analytical','Location','south')
legend('boxoff');

% compute 2 rel error along the x-0 and y=0 axis

re_err_x = abs(T_sol( ID_array(klt_z,1))-AnalyticSol( (mesh.XY(klt_z,1))))./AnalyticSol( (mesh.XY(klt_z,1)))
mean(re_err_x)

re_err_y = abs(T_sol( ID_array(klt_r,1))-AnalyticSol( (mesh.XY(klt_r,2))))./AnalyticSol( (mesh.XY(klt_r,2)))

mean(re_err_y)

figure(3)
trisurf(mesh.conn,mesh.XY(:,1),mesh.XY(:,2),T_sol)

