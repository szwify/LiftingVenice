%----------------------------------------------------------------------
%> @brief Class for the description of the linear 2 nodes 1D element Seg2
%>
%----------------------------------------------------------------------
classdef Element_Seg2 < Element_Manifold_1
        
    % shall we add a properties in there for the 2d axisymmetric case... ?
    % or we derive another class ?
    
    properties
          % 2D axis cases property dependence (for boundary load
          % integration)
    theta = 0.;
    rc = 0.;
    end
    
    methods
        
        %----------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param type : problem type (axis, planestrain ...)
        %> @param  xae : element nodal coordinates
        %> @param the_ien : elt connectivity table
    %    %> @param ncon_nodes : number of connected nodes with adjacent elt
        %- not really needed
        %> @param varargin : if exist element id number
        %----------------------------------------------------------------
        
        function obj=Element_Seg2(type,xae,the_ien,varargin) %ncon_nodes,
            
            obj=obj@Element_Manifold_1(type,xae,the_ien,1,varargin); %ncon_nodes,
            
            if (size(xae,1)~=2)
                disp('error not a valid Seg2 elemt');
                return
            end
            
            if (length(varargin)>0)
                obj.rc = varargin{1};
                obj.theta = varargin{2};
            end
            
        end
        
        %----------------------------------------------------------------
        %> @brief Method to estimate the strain matrix of the element at a
        %> gauss pt
        %>
        %> only 1d axisymm implemented
        %
        %> @param xil : gauss pt coordinate (on the reference [-1,1] elt)
        %> @param obj: instance of the Element_Seg2 class
        %> @retval B : the strain matrix (size depends on problem
        %>  configuration)
        %> @retval j : element jacobian
        %----------------------------------------------------------------
        function [B,j]=Bamat(xil,obj)
            switch obj.type
                
                case 'Axis'  %% One dimensional One dofs axi-symmetry elasticity e_r e_t (e_z=0)
                    
                    DNaDxi=[-0.5  0.5];
                    DxDxi=DNaDxi*(obj.xae);
                    j=DxDxi;
                    
                    DNaDx=(1./j)*DNaDxi;
                    
                    raux=[ 0.5*(1-xil) 0.5*(1+xil)]*(obj.xae);
                    
                    B=[ DNaDx(1) 0 DNaDx(2) 0 ;...
                        0.5*(1-xil)/raux 0 0.5*(1+xil)/raux 0 ; ...
                        0 0 0 0 ];
                    
                otherwise    % 1D stuff... to be checked
                
                    DNaDxi=[-0.5  0.5];
                    DxDxi=DNaDxi*(obj.xae);
                    j=DxDxi;
                    
                    DNaDx=(1./j)*DNaDxi;
                    
                    B=[ DNaDx(1)  DNaDx(2) ;];
                    
            end
            
        end
        
        %----------------------------------------------------------------
        %> @brief Method to compute the element jacobian
        %>
        %> only 1d axisymm implemented
        %
        %> @param obj: instance of the Element_Seg2 class
        %> @retval j : element jacobian
        %----------------------------------------------------------------
        function [j]=Jacobian(xil,obj)
            
            DNaDxi=[-0.5  0.5];
            DxDxi=DNaDxi*(obj.xae);
            j=DxDxi;
            
        end
        
        %----------------------------------------------------------------
        %> @brief Method to estimate global   coordinate from local
        %> coordinates (ref element) and global nodes coordinates of the elt
        %
        %
        %> @param xil : gauss pt coordinate (on the reference [-1,1] elt)
        %> @param obj: instance of the Element_Seg2 class
        %> @retval Xaux : corresponding global coordinates of xil
        %----------------------------------------------------------------
        function [Xaux]=Mapx(xil,obj)
            
            Na_xi=[0.5*(1-xil)  0.5*(1+xil)];
            Xaux=Na_xi*(obj.xae);
            
        end
        
        
        %----------------------------------------------------------------
        %> @brief Method to estimate elt shape function
        %
        %
        %> @param xil : gauss pt coordinate (on the reference [-1,1] elt)
        %> @param obj: instance of the Element_Seg2 class
        %> @retval Na : vector containing value of the different nodal shape
        %> function at xil
        %----------------------------------------------------------------
        function [Na]=Na(xil,obj)
            
            Na=[0.5*(1-xil)  0.5*(1+xil)];
            
        end
        
    end
    
end