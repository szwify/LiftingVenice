%----------------------------------------------------------------------
%> @brief Class for the description of the Tri 3 nodes 2D element Tri3
%>  degenerated from the QU4 element
%>    [-1 1]*[-1 1]
%> the barycenter is at [0. -1./3] on this reference elt
%----------------------------------------------------------------------
classdef Element_Tri3 < Element_Manifold_2
    
    methods
        
        %----------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param type : element type
        %> @param  xae : element nodal coordinates
        %> @param the_ien : elt connectivity table
     %   %> @param ncon_nodes : number of connected nodes with adjacent elt
        %> @param varargin : if exist element id number
        %----------------------------------------------------------------
        
        function obj=Element_Tri3(type,xae,the_ien,varargin) %ncon_nodes,
            
            obj=obj@Element_Manifold_2(type,xae,the_ien,varargin); %ncon_nodes,
            
            if (size(xae,1)~=3)
                disp('error not a valid Tri3 elemt');
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
            
            % 2-3 mat
            DNaDxi=0.25*[ -1.*(1-xil(2)) 1.* (1-xil(2)) 0.   ;... %w.r. to xi
                -1.*(1-xil(1)) -1*(1+xil(1))  2. ];   % w.r. to eta
            
            DxDxi=DNaDxi*(obj.xae);
            j=det(DxDxi);
            
            DxiDx = inv(DxDxi); % could be optimized  - (2,2) mat.
            DNaDx = DxiDx*(DNaDxi);
            
            switch obj.type
                
                case '2D'
                    % exx, eyy, 2 exy
                    
                    % 3 by 6 mat
                    B = [DNaDx(1,1)  0 DNaDx(1,2) 0 DNaDx(1,3) 0  ;...
                        0   DNaDx(2,1) 0 DNaDx(2,2) 0 DNaDx(2,3) ;...
                        DNaDx(2,1) DNaDx(1,1) DNaDx(2,2) DNaDx(1,2) DNaDx(2,3) DNaDx(1,3)
                        ];
                    
                    % easy to add the axisymmetry here ! ....
                case 'Axis'
                 % err, ezz, 2 erz, ett

                    xaux=Mapx(xil,obj);
                    [Nax]=Na(xil,obj);
                    
                    % 4 by 6 mat
                    B = [DNaDx(1,1)  0 DNaDx(1,2) 0 DNaDx(1,3) 0  ;...
                        0   DNaDx(2,1) 0 DNaDx(2,2) 0 DNaDx(2,3) ;...
                        DNaDx(2,1) DNaDx(1,1) DNaDx(2,2) DNaDx(1,2) DNaDx(2,3) DNaDx(1,3)
                        Nax(1)/xaux(1)  0.  Nax(2)/xaux(1) 0 Nax(3)/xaux(1) 0. ];                    
 
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
            
            % 2-3 mat
            DNaDxi=0.25*[ -1.*(1-xil(2)) 1.* (1-xil(2)) 0.   ;... %w.r. to xi
                -1.*(1-xil(1)) -1*(1+xil(1))  2. ];   % w.r. to eta
            
            DxDxi=DNaDxi*(obj.xae);
            j=det(DxDxi);
            
            DxiDx = inv(DxDxi); % could be optimized  - (2,2) mat.
            DNaDx = DxiDx*(DNaDxi);
                     
        end
        %----------------------------------------------------------------
        %> @brief Method to compute the element jacobian
        %>
        %>
        %
        %> @param obj: instance of the Element_Tri3 class (degenerate Qua4)
        %> @retval j : element jacobian
        %----------------------------------------------------------------
        function [j]=Jacobian(xil,obj)
            
            DNaDxi=0.25*[ -1.*(1-xil(2)) 1.* (1-xil(2)) 0.   ;... %w.r. to xi
                -1.*(1-xil(1)) -1*(1+xil(1))  2. ];   % w.r. to eta
            
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
        %> @param obj: instance of the Element_Tri3 class (degenerate Qua4)
        %> @retval Xaux : corresponding global coordinates of xil
        %----------------------------------------------------------------
        function [Xaux]=Mapx(xil,obj)
            
            Na_xi=0.25*[(1-xil(1))*(1-xil(2)) (1+xil(1))*(1-xil(2)) 2*(1+xil(2))];
            Xaux=[Na_xi*(obj.xae(:,1)) Na_xi*(obj.xae(:,2))]; % [x y]
            
        end
        
        %----------------------------------------------------------------
        %> @brief Method to estimate elt shape function
        %
        %>
        %
        %> @param xil : gauss pt coordinate (on the reference [-1,1]*[-1,1] elt)
        %> @param obj: instance of the Element_Tri3 class (degenerate Qua4)
        %> @retval Na : vector containing value of the different nodal shape
        %> function at xil
        %----------------------------------------------------------------
        function [Na]=Na(xil,obj)
            
            Na=0.25*[(1-xil(1))*(1-xil(2)) (1+xil(1))*(1-xil(2)) 2*(1+xil(2))];
            
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
            [j]=Jacobian([0. 0. ],obj);
            vol = 4.*abs(j);  %hummm this is weird.... check todo 
            
        end

    end
    
end