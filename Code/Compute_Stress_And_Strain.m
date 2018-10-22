%-------------------------------------------------------------------
%> @brief Compute element stress and strain from nodal Displacement solution
%> @param obj: an Elasticity_Block object
%> @param  U: solution nodal displacement vector
%> @retval Stress : 6 by n_elt matrix Stresses in each element
%> @retval Strain : 6 by n_elt matrix Strain in each element
%-------------------------------------------------------------------
function [Stress,Strain,AvgCoor]=Compute_Stress_And_Strain(mesh,Geometry,PropertiesList,IntegrationOrder,U,ID_Array,Initial_Stresses,varargin)
% Compute stress strain on solution
if (isempty(U))
    disp('Solve the problem first ...');
    return
end

[nr, ~]=size(Initial_Stresses);
if (nr==1)
    Sig0=Initial_Stresses ;
else
    if (nr~=mesh.Nelts) 
        disp('size of Initial_stresses incompatible with number of elements');
        abort();
    end
    
end

%%% loop over the elts.....
Stress=zeros(mesh.Nelts,6);
Strain=zeros(mesh.Nelts,6);
AvgCoor=zeros(mesh.Nelts,mesh.dim);

for e=1:mesh.Nelts
    
    kt=find(mesh.conn(e,:)~=-1);
    local_ien= mesh.conn(e,kt);
    
    mycoor=mesh.XY(local_ien,:);

    switch mesh.EltType{e}
        
        case 'Seg'
             local_elt=Element_Seg2(Geometry,mycoor,local_ien); %,e
            
        case 'Qua'
            
             local_elt=Element_Qua4(Geometry,mycoor,local_ien);
            
        case 'Tri'
            
            local_elt=Element_Tri3(Geometry,mycoor,local_ien);
            
    end
    
   Prop_local=PropertiesList{mesh.matID(e)};

    obj_aux=Operator_Elasticity(local_elt,Prop_local);
    
    % get the correct eq number
    neq_local=ID_Array(local_ien,:);
        
    if (issparse(U))
        Ue=full(U(neq_local(:)))';  
    else
        Ue=U(neq_local(:))';       
    end
     
    % initial stress field
    if (nr~=1) % cases where it is different in each element. 
     Sig0=Initial_Stresses(e,:);
    end
    
    if (~isempty(varargin))
        if strcmp(varargin{1},'Gauss')
            % if there is a
            [Sg,St]=StressStrainAtGaussPoints(obj_aux,Ue,Sig0,IntegrationOrder);
            Stress(e,1:length(Sg(:)))=Sg(:) ;
            Strain(e,1:length(St(:)))=St(:) ;
        end

    else
        
        [Sg,St,coorAux]=AverageStressStrainElt(obj_aux,Ue,Sig0,IntegrationOrder);
        Stress(e,1:length(Sg))=Sg ;
        Strain(e,1:length(St))=St ;
        
        AvgCoor(e,1:mesh.dim)=coorAux;
        
    end
    
end

end
