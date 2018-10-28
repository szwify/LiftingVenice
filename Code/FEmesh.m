%> @brief Class for a Finite Element mesh
%
%>  A class for the finite element mesh with material zones (by element)
%>  Author: Brice Lecampion
%>  First date : October 5 2008
%> modified Sept 2018

classdef  FEmesh
    
    properties
        %> The connectivity table of the mesh
        conn ; 
        %> the nodes coordinates
        XY; 
        
        %> a vector containting the material zone id (integer) for each element
        matID; 
        
    end
    
    properties (Dependent)
        %> Total number of element
        Nelts;
        %> Number of material zones in the mesh
        Nmats;
        %> A vector of size Nelts containing element type segment, triangle, square, etc.....
        EltType;
         %> spatial dimension, scalar 
		dim ; % 
		%> total number of nodes
        n_nodes; 
    end
        
    methods
        %---------------------------------------------------------------
        %> @brief Class Constructor
        %> @param obj_fenode FEnode : object (grid)
        %> @param ien : connectivity table
        %> @param varargin : optional input : material zones
        %> @return  an instance of the FEmesh class
        %---------------------------------------------------------------
        function obj=FEmesh(XY,ien,varargin)
            
            optargin = size(varargin,2);
            if (optargin==0)
                matID = ones(size(ien,1),1);
            else
                matID = varargin{1}{1};
            end
            
            if (size(matID,1)~=size(ien,1))
                error('ien and mat_zones not consistent: mesh object is not created' );
                %              disp(size(mat_zones,1));disp(size(ien,1));
            end
            
            obj.conn=ien;
           
            obj.XY=XY;
            
            obj.matID=matID;
            
        end
       
        %---------------------------------------------------------------
        %> @brief Function  getting number of elt in the object
        %> @param obj : FEmesh object
        %> @retval Nelts : Number of elts in the mesh
        %---------------------------------------------------------------
        function Nelts=get.Nelts(obj)
            % number of elements
            Nelts=size(obj.conn,1);
        end
         %---------------------------------------------------------------
        % dependent properties
		%> @brief get the spatial dimension of the mesh
		%> @param obj : a FEmesh object
		%> @retval dim: spatial dimension of the mesh
         %---------------------------------------------------------------
        function dim=get.dim(obj)
              dim=size(obj.XY,2);
        end
        
        %---------------------------------------------------------------
		%> @brief Get the total number of nodes
		%> @param obj : a FEnode object
		%> @retval : n_tot tot number of nodes
        %---------------------------------------------------------------
        function n_nodes=get.n_nodes(obj)
         n_nodes=size(obj.XY,1);
        end

        %---------------------------------------------------------------
        %> @brief Function  getting number of different material in the Mesh object
        %> @param obj : FEmesh object
        %> @retval Nmats : Number of materials in the mesh
        %---------------------------------------------------------------
        function Nmats=get.Nmats(obj)
            % a loop on the element needed
            aux=sort(obj.matID);
            z1=aux(1);Nmats=1;
            for i=2:obj.Nelts
                if(aux(i)~=z1)
                    z1=aux(i);
                    Nmats=Nmats+1;
                end
            end
        end
        
        %---------------------------------------------------------------
        %> @brief Function getting the type of each element in the object
        %> @param obj : FEmesh object
        %> @retval GeomType : array of string containing the type of each
        %> elt
        %---------------------------------------------------------------
        function EltType=get.EltType(obj)
            %%% get the geometry type of each element.....
            % GeomType={};
            for e=1:obj.Nelts
                aux=find(obj.conn(e,:)~=-1); %%% -1 denotes no node there in ien
                switch length(aux)
                    case 2
                        EltType{e}='Seg2';
                    case 3
                        EltType{e}='Tri3';
                    case 4
                        if (obj.dim==2)
                            EltType{e}='Qua4';
                        else
                            if (obj.dim==3)
                                EltType{e}='Tetra';
                            end
                        end
                end
            end
            
        end
        
        %---------------------------------------------------------------
        %> @brief : Plot the mesh via Delaunay (seg in 1d, triangle in 2D)
        %> @param obj : FEmesh object
        % work only for a FE with similar element type
        %---------------------------------------------------------------
        function []=plot_mesh(obj)
            % plot the tessalation via Delaunay :: Triangle
            % test the dimension
            switch obj.dim
                case 1
                    plot(obj.XY.coor,0.*(obj.XY.coor),'.-');
                case 2

                    if (length(obj.conn(1,:))==3)
                        triplot(obj.conn,obj.XY(:,1),obj.XY(:,2));
                    else % case Quadrilateral
                        scatter(obj.XY(:,1),obj.XY(:,2)); hold on;
                        for e=1:ne_t
                            line(obj.XY(obj.conn(e,:),1),obj.XY(obj.conn(e,:),2)); hold on;
                        end

                    end
                    
                case 3
                    %                    T=delaunayn(obj.FEnode.coor);
                    tetramesh(obj.conn,obj.XY);
            end
            
        end
        
        %---------------------------------------------------------------
        %> @brief : Compute the middle point coordinates of the element in the mesh
        %> @param obj : FEmesh object
        %---------------------------------------------------------------
        function [MidCoor]=Element_MidPoint(obj)
            
            %%% loop over the element ...
            for e=1:obj.Nelts
                MidCoor(e,:)=mean(obj.XY(obj.conn(e,:)));
            end
            
        end
        
    end
    
end