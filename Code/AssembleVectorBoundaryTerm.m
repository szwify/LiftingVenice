function [F]=AssembleVectorBoundaryTerm(meshObj,Geometry,Operator,ValueList,ID_array,IntegrationOrder)

% ValueList
% either
% [dof node1 intensity1 node2   intensity2 ] -- we assume therefore that the Load varies as the element order
% we assume that we have adjacent nodes (node1 node2)... (such that we get it define an element edge)..

nnodes=meshObj.n_nodes;

switch Operator
    
    case 'BoundaryLoads'  % elastic boundary loads.
        ndof = 2; % number of dof per nodes
   
        tot_dof = length(ID_array(:));
        n_row=tot_dof;
        if (n_row~=nnodes*ndof)
            disp('error in id_array given');
            abort();
        end
        
        [nr , ncv]=size(ValueList);
        if (ncv ~=5)
            disp('Error in boundary value input');
            abort()
        end
        
end

disp('number of segments with applied Neumann BC: ' );
disp(nr);

F=sparse(n_row,1);

for i=1:nr
    
    % first node of segment
    [nrow1,ncol1]=find(meshObj.conn(:,:)==ValueList(i,2)) ; % get element with that node
    
    %            disp(nrow1);
    % note the test below was working only for Qua4 type mesh
    if length(nrow1)>2
        %                error(' this is not a boundary node !');
    end
    xyn1=meshObj.XY(ValueList(i,2),:);
    
    % second node of segment
    [nrow2,ncol2]=find(meshObj.conn(:,:)==ValueList(i,4)) ; % get element with that node
    
    % note the test below was working only for Qua4 type mesh
    if length(nrow2)>2
        %                error(' this is not a boundary node !');
    end
    
    xyn2=meshObj.XY(ValueList(i,4),:);
    e = intersect(nrow1,nrow2);
    
    switch meshObj.dim
        
        case 1% body mesh is 1D -> for now the integration of boundary nodes is not implemented and should be done manually
            % Not implemented
            error(' F load integration not implemented for 1D problems');
        case 2
            % this is ok for plane-strain,
            
            switch Geometry
                case '2D'
                    
                    mycoor= [0 ; norm(xyn2-xyn1)] ; % we create a dummy seg2D element
                    local_elt=Element_Seg2(Geometry,mycoor,[1 2]);
                    
                case 'Axis' %    axisymmetry..... with a hack for the 2pi r factor
                    rc_aux = (xyn2+xyn1)/2.;
                    rc = rc_aux(1);
                    % disp(rc);
                    xdiff = (xyn2-xyn1);s = xdiff / norm(xdiff);
                    theta = acos(s(1)/norm(s));
                    if (s(2)<0)
                        theta=-theta;
                    end
                    
                    mycoor=[-norm(xyn2-xyn1)/2.; norm(xyn2-xyn1)/2.] ; % we create a dummy seg2D element that needs to get
                    local_elt=Element_Seg2(Geometry,mycoor,[1 2],rc,theta);
                    
            end
            
            % 3D not implemented yet....
            
    end
    
    switch Operator
        
        case 'BoundaryLoads'
            tt=Operator_Load_Boundary(local_elt,ValueList(i,[3, 5]));
            
    end
    
    % element matrix integration
    elt_F=@(xi) K_Elt_Int(xi,tt); % create function handle  
    
    [Kel]=GaussIntegration(elt_F,1,IntegrationOrder);
    
    [nr nc]=size(Kel);
    
    %      disp([nr nc]);
    %      disp([length(neq_row) length(neq_col)]);
    neq_row =ID_array(ValueList(i,[2;4]),ValueList(i,1));
    
    F(neq_row,1)=F(neq_row,1)+Kel;
    
end

end