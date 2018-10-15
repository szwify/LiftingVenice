%------------------------------------------------------------------
%> @brief A high-level class for description of a single element
%>
%> This class provides basic properties of a mesh
%------------------------------------------------------------------
classdef Element_General
   
    % this is a very high level class, description of a single element (or
    % cell)
    
    properties
        %> Number of nodes the element should connect
       % ncon_nodes;
        %> Local connectivity table of the elt
        ien ;    % local connectivity table
        %> Element id number
        elt_id=0 ; %%% an integer (e.g. : element number)
        %> spatial dimension of the element
        dim;     % spatial dimension of the element
    end
    
    methods
        %------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param the_ien : connectivity table of the element
      %  %> @param ncon_nodes   : number of connected nodes (= number of nodes for linear FE)
        %> @param dim : spatial dimension
        %> @param varargin : elt_id if present
        %------------------------------------------------------------------        
        function obj=Element_General(the_ien,dim,varargin) %ncon_nodes,
            
%             if (length(the_ien)~=ncon_nodes)
%                 disp('Error between connectivity length and number of connected nodes');
%                 return;
%             end
            obj.ien=the_ien;
         %   obj.ncon_nodes=ncon_nodes;
            obj.dim=dim;
            if (nargin>0)
                obj.elt_id=varargin{1}{1};
            end
            
            
        end
        
        
    end
    
    
end