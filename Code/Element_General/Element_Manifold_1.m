%------------------------------------------------------------------
%> @brief A   class for the description of 1D element
%>
%> This class provides basic properties of a 1D element
%------------------------------------------------------------------
classdef Element_Manifold_1 < Element_General
   
    %  1D element class derived from Element_General
    
    properties
        %> type of the element 
        type ; % Plane_strain, axisym etc.
        %> nodal coordinates of the element
        xae; 
    end 
       
    methods
        
        %----------------------------------------------------------------
        % @brief Class Constructor
        %
        % @param type : element type
        % @param  xae : element nodal coordinates
        % @param the_ien : elt connectivity table
        % @param ncon_nodes : number of connected nodes with adjacent elt
        % @param varargin : if exist element id number
        %----------------------------------------------------------------
        function obj=Element_Manifold_1(type,xae,the_ien,varargin)
            
            obj = obj@Element_General(the_ien,1,varargin);
            
            if (size(xae,2)~=1)
                disp('error not a valid 1D elemt');
                %   clear obj;
                return
            end
            
            obj.type=type;
            obj.xae=xae;
            
        end
        
    end
    
    
end