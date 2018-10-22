
% 1D diffusion equation - with a constant initial condition 
% zero dirichlet boundary at the bottom, zero flux at the top

ymax=1.;
xmax=1.;

    node = [                % list of xy "node" coordinates
        0, 0                % outer square
        xmax, 0
        xmax, ymax
        0, ymax 
        ] ;
    
    edge = [                % list of "edges" between nodes
        1, 2                % outer square 
        2, 3
        3, 4
        4, 1 
          ] ;

%------------------------------------------- call mesh-gen.
   [vert,etri, tria,tnum] = refine2(node,edge,[],[],0.1) ; % do not touch
  
  the_coor=vert; 
  

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
  
 ne_t=length(connect(:,1));
  
% simple grid plot (note should work for all linear-like mesh) .... very
% naive, lots of duplicates.... should be a loop on mesh edges.
 scatter(the_coor(:,1),the_coor(:,2)); hold on;
 for e=1:ne_t
     line(the_coor(connect(e,:),1),the_coor(connect(e,:),2)); hold on;
 end

mesh=FEmesh(the_coor,connect);
%%

propList={1.};

[L,ID_array]=AssembleMatrix(mesh,'2D','Laplacian',propList,3);

[M,ID_array2]=AssembleMatrix(mesh,'2D','Mass',propList,3);

%%


% put bottom dof to zero
kl=find(mesh.XY(:,2)==0);
eq_fix=ID_array(kl);
eq_free= setdiff(ID_array(:),eq_fix(:))

np_unknowns=length(eq_free);

L_f=L(eq_free,eq_free);
M_f=M(eq_free,eq_free);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  simple drainage with uniform initial pp %%%%%%%%%%%%%%%
yns=unique(mesh.XY(2:end,2)); % 

% implicit scheme .....
pn =ones(np_unknowns,1);  % unit initial pore pressure
time_step=0.01;%*hy^2.; % Pe h^2/6
n_step=100;
pns=[ ];
pns(1,:)=ones(length(eq_free)+length(eq_fix),1);
tn=0.;
taus=[tn];
res =[];
res(1,:)=ones(length(yns),1);
theta=0.5;

for j=1:n_step
    tn=tn+time_step;
    dp=(M_f+theta*time_step*L_f)\(-time_step*L_f*pn) ;
    pn=pn+dp;
    pns(j+1,eq_free)=pn;
    pns(j+1,eq_fix)=0.;
    
    taus(j+1)=tn;
    res(j+1,:)=pp_terzaghi(yns,tn/4.);   % be careful the domain is of size 2 so dimensionless time = time/ 2^2
    
end

kl=find((mesh.XY(:,1)==0)) ; % && ( reg_mesh.FEnode.coor(:,2)~=0)));
kl_id=ID_array(kl);

figure(2)
plot(mesh.XY(kl,2),pns(1:10:n_step ,kl_id),'.k'); hold on;
plot(yns,res(1:10:n_step,:),'-r')

figure(3)
plot(taus,pns(:,kl_id(2))) ; hold on;
plot(taus,res(:,end),'-r');
 
 
