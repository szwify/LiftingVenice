function [K,ID_array_u,ID_array_p]=AssembleCouplingMatrix(meshObj_u,meshObj_p,Geometry,PropertiesList,IntegrationOrder)
% assembly of H-M or T-M coupling matrices
% Inputs
%
% meshObj_u :: a mesh object with nodes coor and connectivity etc. (created
% with FEmesh) on which the displacement   
% Geometry :: string either '2D' (for plane strain), 'Axis' for
% axisymmetric problem
% Operator :: String depicting the corresponding PDE operator 
% PropertiesList: list of either stiffness matrix, or relevant properties -
% of length equal to the number of different material in the mesh
% IntegrationOrder :: gauss integration order
%
% Outputs:
% 
% K :: the global stiffness matrix (sparse)
% ID_array :: the id array (not in the case of coupling, we output only the
% elastic ID_array)
%
%
% Restrictions : 2D problem supported only 


nnodes_u=meshObj_u.n_nodes;
nnodes_p=meshObj_p.n_nodes;
if (nnodes_u == nnodes_p)
    
    disp('Same interpolation order for both coupling field - violating LBB'); % 
end
if (meshObj_u.Nelts ~= meshObj_p.Nelts)
    disp('Error the 2 given mesh do not have the same number of elements');
    return
end


ndof_u = 2; % number of dof per nodes for mech
ndof_p = 1; % number of dof per nodes for mech

ID_array_u=reshape([1:ndof_u*nnodes_u]',ndof_u,nnodes_u)';
ID_array_p=reshape([1:ndof_p*nnodes_p]',ndof_p,nnodes_p)';

n_row=length(ID_array_u(:));
n_col=length(ID_array_p(:));
        

K=sparse(n_row,n_col); 

for e=1:meshObj_u.Nelts
    
    kt=find(meshObj_u.conn(e,:)~=-1);
    local_ien_u= meshObj_u.conn(e,kt);
    neq_row=[];
    for i=1:length(local_ien_u)
        neq_row=[neq_row ID_array_u(local_ien_u(i),:)];
    end
       
    mycoor=meshObj_u.XY(local_ien_u,:);
    switch meshObj_u.EltType{e}
        case 'Seg2'           
            local_elt_u=Element_Seg2(Geometry,mycoor,local_ien_u); %,e
        case 'Qua4'           
            local_elt_u=Element_Qua4(Geometry,mycoor,local_ien_u);
        case 'Tri3'         
            local_elt_u=Element_Tri3(Geometry,mycoor,local_ien_u);   
            % add here case for Quadratic element....
    end

    kt=find(meshObj_p.conn(e,:)~=-1);
    local_ien_p= meshObj_p.conn(e,kt);
    neq_col=[];
    for i=1:length(local_ien_p)
        neq_col=[neq_col ID_array_p(local_ien_p(i),:)];
    end
    
    mycoor=meshObj_p.XY(local_ien_p,:);
    switch meshObj_p.EltType{e}
        case 'Seg2'           
            local_elt_p=Element_Seg2(Geometry,mycoor,local_ien_p); %,e
        case 'Qua4'           
            local_elt_p=Element_Qua4(Geometry,mycoor,local_ien_p);
        case 'Tri3'         
            local_elt_p=Element_Tri3(Geometry,mycoor,local_ien_p);            
    end
    
     
    Prop_local=PropertiesList{meshObj_u.matID(e)};
    
    
    tt=Operator_Coupling_Matrices(local_elt_u,local_elt_p,Prop_local);
    
    % element matrix integration
    elt_F=@(xi) K_Elt_Int(xi,tt); % create function handle
    Kel= GaussIntegration(elt_F,meshObj_u.dim,IntegrationOrder);
    
    [nr nc]=size(Kel);
    
%      disp([nr nc]);
%      disp([length(neq_row) length(neq_col)]);

    if ( (nr~=length(neq_row)) || (nc~=length(neq_col))  )
        disp('error : incompatible ');
        return
    end
    
    K(neq_row,neq_col)=K(neq_row,neq_col)+Kel;
    
end

end