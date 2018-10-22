
% 1D diffusion equation - with a constant initial condition 
% zero dirichlet boundary at the bottom, zero flux at the top

ne_x=1; %  
ne_y=40; %  

ne_t=ne_x*ne_y;
ymax=1.;
xmax=1.;

hx=xmax/ne_x;
hy =ymax/ne_y;

xs=linspace(0.,xmax,ne_x+1)
ys=linspace(0.,ymax,ne_y+1)

the_coor=zeros((ne_x+1)*(ne_y+1),2);
k=1;
for i=1:length(ys)
    
    the_coor(k:(k+length(xs)-1),1)=xs';
    the_coor(k:(k+length(xs)-1),2)=ys(i);
    
    k=k+length(xs);
    
end
% create the corresponding connectivity table
connect=zeros(ne_t,4);
e=1;
for ey=1:ne_y
    for ex=1:ne_x
       connect(e,1:4)=[ex+(ey-1)*(ne_x+1), ex+1+(ey-1)*(ne_x+1), (ex+1)+(ey)*(ne_x+1), (ex)+(ey)*(ne_x+1)];   
       e=e+1;
    end
end

% simple grid plot (note should work for all linear-like mesh) .... very
% naive, lots of duplicates.... should be a loop on mesh edges.
 scatter(the_coor(:,1),the_coor(:,2)); hold on;
 for e=1:ne_t
     line(the_coor(connect(e,:),1),the_coor(connect(e,:),2)); hold on;
 end

%%

 
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
time_step=0.1*(1/ne_y)^2.; % Pe h^2/6
n_step=5000;
pns=[ ];
pns(1,:)=pn; tn=0.;
taus=[tn];
res =[];
res(1,:)=ones(length(yns),1);
theta=0.5;

for j=1:n_step
    tn=tn+time_step;
    dp=(M_f+theta*time_step*L_f)\(-time_step*L_f*pn) ;
    pn=pn+dp;
    pns(j+1,:)=pn;
    taus(j+1)=tn;
    res(j+1,:)=pp_terzaghi(yns,tn/4.);   % be careful the domain is of size 2 so dimensionless time = time/ 2^2
end

%%

figure(2)
plot(mesh.XY(3:end,2),pns(1:10:n_step,:),'.-k'); hold on;
plot(yns,res(1:20:n_step,:),'-r')
xlabel(' x');
ylabel('Pressure');

figure(3)

plot(taus,pns(:,2*30)) ; hold on;
plot(taus,res(:,30),'-r');
xlabel(' time');
ylabel('Pressure at 30th nodes');

%%

rel0=abs(res(:,30)-pns(:,2*30))./res(:,30)
abs0=abs(res(:,30)-pns(:,2*30))

figure(4);
plot(taus,rel0); hold on;
xlabel(' time');
ylabel('relative error at 30th nodes');

%
% analytical solution 
% tau=[0.0001 0.001 0.01 0.05 0.1   1.];
% y_fine=linspace(0.,1.,100);
% [res]=pp_terzaghi(y_fine,tau,100);
% 
% figure(2)
% for i=1:length(tau)
%     plot(y_fine,res(:,i)); hold on;
% end
% 

