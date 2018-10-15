%----------------------------------------------------------------------
%> @brief Class for the description of the Quad 4 nodes 2D element Qua4
%>
%    [-1 1]*[-1 1]
%----------------------------------------------------------------------
classdef Element_Qua4 < Element_Manifold_2
    
    methods
        
        %----------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param type : element type
        %> @param  xae : element nodal coordinates
        %> @param the_ien : elt connectivity table
        %> @param ncon_nodes : number of connected nodes with adjacent elt
        %> @param varargin : if exist element id number
        %----------------------------------------------------------------
        
        function obj=Element_Qua4(type,xae,the_ien,varargin) %,ncon_nodes
            
            obj=obj@Element_Manifold_2(type,xae,the_ien,varargin); %,ncon_nodes
            
            if (size(xae,1)~=4)
                disp('error not a valid Qua4 elemt');
                return
            end
            
        end
        
        %----------------------------------------------------------------
        %> @brief Method to estimate the strain matrix of the element at a
        %> gauss pt
        %>
        %> only 1d axisymm implemented
        %
        %> @param xil : gauss pt coordinate (on the reference [-1,1]*[-1,1] elt)
        %> @param obj: instance of the Element_Qua4 class
        %> @retval B : the strain matrix (size depends on problem
        %>  configuration)
        %> @retval j : element jacobian
        %----------------------------------------------------------------
        function [B,j]=Bamat(xil,obj)
            
            DNaDxi=0.25*[ -1.*(1-xil(2)) 1.*(1-xil(2)) (1+xil(2))  -1.*(1+xil(2))   ;... %w.r. to xi
                -1.*(1-xil(1)) -1*(1+xil(1)) (1+xil(1)) (1-xil(1)) ];   % w.r. to eta
            
            DxDxi=DNaDxi*(obj.xae);
            j=det(DxDxi);
            
            DxiDx = inv(DxDxi); % could be optimized  - (2,2) mat.
            DNaDx =  DxiDx*(DNaDxi);
            
            switch obj.type
                
                case '2D'  %'PlaneStrain'
                    
                    % 3 by 8 mat
                    B = [DNaDx(1,1)   0 DNaDx(1,2) 0 DNaDx(1,3) 0 DNaDx(1,4) 0 ;...
                        0   DNaDx(2,1) 0 DNaDx(2,2) 0 DNaDx(2,3) 0 DNaDx(2,4) ;...
                        DNaDx(2,1) DNaDx(1,1) DNaDx(2,2) DNaDx(1,2) DNaDx(2,3) DNaDx(1,3) DNaDx(2,4) DNaDx(1,4)
                        ];
                    
                case 'Axis'
                    % err, ezz, 2 erz, ett
                    xaux=Mapx(xil,obj);
                    [Nar]=Na(xil,obj);
                    
                    % 4 by 8 mat
                    B = [DNaDx(1,1)   0 DNaDx(1,2) 0 DNaDx(1,3) 0 DNaDx(1,4) 0 ;...
                        0   DNaDx(2,1) 0 DNaDx(2,2) 0 DNaDx(2,3) 0 DNaDx(2,4) ;...
                        DNaDx(2,1) DNaDx(1,1) DNaDx(2,2) DNaDx(1,2) DNaDx(2,3) DNaDx(1,3) DNaDx(2,4) DNaDx(1,4);...
                        Nar(1)/xaux(1) 0. Nar(2)/xaux(1) 0. Nar(1)/xaux(3) 0. Nar(4)/xaux(1) 0. ...
                        ];
                     
            end
            
        end
        
        %----------------------------------------------------------------
        %> @brief Method to estimate the gradient of shape fction matrix  
        %>
        %
        %> @param xil : gauss pt coordinate (on the reference [-1,1]*[-1,1] elt)
        %> @param obj : instance of the Element_Qua4 class
        %> @retval GN : the gradient of shape fction matrix (size depends on problem
        %>  configuration)
        %> @retval j : element jacobian
        %----------------------------------------------------------------
        function [DNaDx,j]=GradN(xil,obj)
            
            DNaDxi=0.25*[ -1.*(1-xil(2)) 1.*(1-xil(2)) (1+xil(2))  -1.*(1+xil(2))   ;... %w.r. to xi
                -1.*(1-xil(1)) -1*(1+xil(1)) (1+xil(1)) (1-xil(1)) ];   % w.r. to eta
            
            DxDxi=DNaDxi*(obj.xae);
            j=det(DxDxi);
            
            %DxiDx = inv(DxDxi); % could be optimized  - (2,2) mat.
            DNaDx =DxDxi\DNaDxi;  %DxiDx*(DNaDxi);
            
        end
        
        %----------------------------------------------------------------
        %> @brief Method to compute the element jacobian
        %>
        %>
        %
        %> @param obj: instance of the Element_Qua4 class
        %> @retval j : element jacobian
        %----------------------------------------------------------------
        function [j]=Jacobian(xil,obj)
            
            DNaDxi=0.25*[ -1.*(1-xil(2)) 1.*(1-xil(2)) (1+xil(2))  -1.*(1+xil(2))   ;...
                -1.*(1-xil(1)) -1*(1+xil(1)) (1+xil(1)) (1-xil(1)) ];
            
            DxDxi=DNaDxi*(obj.xae);
            j=det(DxDxi);
            
        end
        
        %----------------------------------------------------------------
        %> @brief Method to estimate global   coordinate from local
        %> coordinates (ref element) and global nodes coordinates of the elt
        %
        %>
        %
        %> @param xil : gauss pt coordinate (on the reference [-1,1]*[-1,1] elt)
        %> @param obj: instance of the Element_Qua4 class
        %> @retval Xaux : corresponding global coordinates of xil
        %----------------------------------------------------------------
        function [Xaux]=Mapx(xil,obj)
            
            Na_xi=0.25*[(1-xil(1))*(1-xil(2)) (1+xil(1))*(1-xil(2)) (1+xil(1))*(1+xil(2)) (1-xil(1))*(1+xil(2))];
            Xaux=[Na_xi*(obj.xae(:,1)) Na_xi*(obj.xae(:,2))]; % [x y]
            
        end
        
        
        %----------------------------------------------------------------
        %> @brief Method to estimate elt shape function
        %
        %>
        %
        %> @param xil : gauss pt coordinate (on the reference [-1,1]*[-1,1] elt)
        %> @param obj: instance of the Element_Qua4 class
        %> @retval Na : vector containing value of the different nodal shape
        %> function at xil
        %----------------------------------------------------------------
        function [Na]=Na(xil,obj)
            
            Na=0.25*[(1-xil(1))*(1-xil(2)) (1+xil(1))*(1-xil(2)) (1+xil(1))*(1+xil(2)) (1-xil(1))*(1+xil(2))];
            
        end
        
        %----------------------------------------------------------------
        %> @brief Method to estimate elt volume / area
        %
        %>
        %
        %> @param obj: instance of the Element_Qua4 class
        %> @retval volume: scalar with volume
        %----------------------------------------------------------------
        function [vol]=Volume(obj)
            [j]=Jacobian(obj);
            vol = 4.*j;
            
        end
        
        
    end
    
end