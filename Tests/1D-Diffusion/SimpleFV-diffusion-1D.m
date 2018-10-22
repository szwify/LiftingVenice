
% 1D diffusion equation - with a constant initial condition 
% simple FV 
% zero dirichlet boundary at the bottom, zero flux at the top
% domain is between 0 and 1 

ne =200;
h = 1./ne;

xend =linspace(0,1.,ne);
np_unknowns=ne-1; % 

%%% create FV matrix L 
L=sparse(np_unknowns,np_unknowns);
for e=1:np_unknowns
    
       L(e,e)=-2.; % diagonal term 
    
    switch e
        case 1    % first cell has no cell to its left
           L(e,e+1)=1;
        case  np_unknowns  % the last cell has no outgoing flux (v_i+1/2=0)
            L(e,e-1)=1;
            L(e,e)=-1;
        otherwise       % all the other cells have left and right flux
             L(e,e+1)=1;
             L(e,e-1)=1;
    end
    
end

L=L/h^2.;
Id =speye(np_unknowns,np_unknowns);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  simple drainage with uniform initial pp %%%%%%%%%%%%%%%
xs=xend(1:ne-1)+h/2;   % computing the center of each cell 
% 
pn =ones(np_unknowns,1);  % unit initial pore pressure

time_step=60*h^2./2.; % 

tMax=1.;  %  maximum time up to which we seek the solution
n_nax=1000;    % max nunber of steps to compute
pns=[ ];       % 
pns(1,:)=pn;   % pns will be a matrix with the solution of each time steps (1 step = 1 row)
tn=0.;
taus=[tn];     % storing in a vector all the time where the solution is computed.
res =[];       % a matrix for the corresponding analytical solution results 
res(1,:)=pn;   

theta=1;    % theta - time integration scheme choice (theta \in [0,1] )

j=0;
while tn<tMax && j<n_nax
    j=j+1;
    tn=tn+time_step;    % increase the time by delta t
    dp=(Id-theta*time_step*L)\(time_step*L*pn) ;  % I compute the increment of pressure dp
    pn=pn+dp;      % adding to the previous pressure
    pns(j+1,:)=pn;   % store the results in a 'big' matrix
    taus(j+1)=tn;    % store the corresponding time
    res(j+1,:)=pp_terzaghi(xs,tn/4.,900)';   % compute the corresponding analytical
    %solution at each center of cells for the given time tn
    %  be careful the domain is of size 2 so the dimensionless time = time/ 2^2
end

%%
figure(2)
title('Pressure profile at different time '); hold on;
plot(xs,pns(1:20:n_nax,:),'.-k'); hold on;
plot(xs,res(1:20:n_nax,:),'-r'); hold on;
xlabel(' x');
ylabel('Pressure');


%%
figure(3)
cell_id = 20;

title(strcat('Pressure versus time at  cell #',string(cell_id))); hold on;
plot(taus,pns(:,cell_id)) ; hold on;
plot(taus,res(:,cell_id),'-r');
xlabel(' time');
ylabel(strcat('Pressure at cell #', string(cell_id)));


 % corresponding relative error between the numerical results and the analytical solution
rel0=abs(res(:,10)-pns(:,10))./res(:,10)  
abs0=abs(res(:,10)-pns(:,10))
 
figure(4)
plot(taus,rel0); hold on;
xlabel(' time');
ylabel(strcat('relative error on Pressure at cell #', string(cell_id)));


