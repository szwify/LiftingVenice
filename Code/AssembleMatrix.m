function [K,ID_array]=AssembleMatrix(meshObj,Geometry,Operator,PropertiesList,IntegrationOrder)
% Inputs
%
% meshObj :: a mesh object with nodes coor and connectivity etc. (created
% with FEmesh
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


switch Operator
    case 'Laplacian'
        nnodes=meshObj.n_nodes;

        ndof = 1;
        max_eqn_numb=ndof*nnodes;
        ID_array=reshape([1:max_eqn_numb]',ndof,nnodes)';
        tot_dof = length(ID_array);
        n_row=tot_dof;
        n_col=tot_dof;
        
    case 'Mass'
        nnodes=meshObj.n_nodes;

        ndof = 1;
        max_eqn_numb=ndof*nnodes;
        ID_array=reshape([1:max_eqn_numb]',ndof,nnodes)';
        tot_dof = length(ID_array);
        n_row=tot_dof;
        n_col=tot_dof;
        
    case 'Elasticity'
        nnodes=meshObj.n_nodes;

        ndof = 2; % number of dof per nodes
        max_eqn_numb=ndof*nnodes;
        ID_array=reshape([1:max_eqn_numb]',ndof,nnodes)';
        tot_dof = length(ID_array(:));
        n_row=tot_dof;
        n_col=tot_dof;
        
end

K=sparse(n_row,n_col); 

for e=1:meshObj.Nelts
    
    kt=find(meshObj.conn(e,:)~=-1);
    local_ien= meshObj.conn(e,kt);
    
    neq_row=[];
    for i=1:length(local_ien)
        neq_row=[neq_row ID_array(local_ien(i),:)];
    end
    
    
    mycoor=meshObj.XY(local_ien,:);
    
    switch meshObj.EltType{e}
        case 'Seg2'
            
            local_elt=Element_Seg2(Geometry,mycoor,local_ien); %,e
            
        case 'Qua4'
            
            local_elt=Element_Qua4(Geometry,mycoor,local_ien);
            
        case 'Tri3'
%             
            local_elt=Element_Tri3(Geometry,mycoor,local_ien);
            
   end
    
    Prop_local=PropertiesList{meshObj.matID(e)};
    
    switch Operator
        case 'Laplacian'
            tt=Operator_Laplacian(local_elt,Prop_local);
            neq_col=neq_row;
        case 'Mass'
            tt=Operator_Mass(local_elt,Prop_local);
            neq_col=neq_row;
        case 'Elasticity'
            tt=Operator_Elasticity(local_elt,Prop_local);
            neq_col=neq_row;

    end
    
    % element matrix integration
    elt_F=@(xi) K_Elt_Int(xi,tt); % create function handle
    Kel= GaussIntegration(elt_F,meshObj.dim,IntegrationOrder);
    
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