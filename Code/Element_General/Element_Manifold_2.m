%------------------------------------------------------------------
%> @brief A   class for the description of 1D element
%>
%> This class provides basic properties of a 1D element
%------------------------------------------------------------------
classdef Element_Manifold_2 < Element_General
   
    %  2D element class derived from Element_General
    
    properties
        %> type of the element 
        type ; % 'PlaneStrain', 'Axis' etc.
        %> nodal coordinates of the element
        xae; 
    end
    
    
    methods
        
        %----------------------------------------------------------------
        % @brief Class Constructor
        %
        % @param type : element type
        % @param  xae : element nodal coordinates
        % @param the_ien : elt connectivity table  -> unused
        % @param ncon_nodes : number of connected nodes with adjacent elt
        % -> unused
        % @param varargin : if exist element id number
        %----------------------------------------------------------------
        function obj=Element_Manifold_2(type,xae,the_ien,varargin) %ncon_nodes,
            
            obj = obj@Element_General(the_ien,2,varargin); % 2D elt ncon_nodes,
            
            if (size(xae,2)~=2)
                disp('error not a valid 2D elemt');
                %   clear obj;
                return
            end
            
            obj.type=type;
            obj.xae=xae;
            
        end
        
    end
    
    
end