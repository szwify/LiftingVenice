%------------------------------------------------------------------
%> @brief Laplacian operator class at the element level
%
%> This class provides the FE implementation of the Laplacian Operator at
%> the element level.
%>
%>  Bi-linear form
%>
%>   \int_Elt \nabla u \times Cond \times \nabla v    d Elt
%>
%>  Element matrix
%>
%>  k_el = \int_elt Grad N'* Cond * Grad N d Elt
%------------------------------------------------------------------
classdef Operator_Laplacian 
    
    properties
        %> Geometry
        Geometry; 
        %> EltObject
        EltObject;
        %> conductivit proper
        Cond;
    end
    
    methods
        %------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param EltObject : Elemenet object
        %> @param Cond conductivity properties
       
        %------------------------------------------------------------------
        function obj=Operator_Laplacian(EltObject,Cond)
            
            obj.EltObject=EltObject;
            obj.Geometry = EltObject.type;
            obj.Cond=Cond;
            
        end
        
        %------------------------------------------------------------------
        %> @brief Compute the integrand of the  Element ' stiffness' matrix
        %>entry
        %
        %> B'*Cond*B @ xil.
        %>  For Laplacian problem,
        %
        %> @param xil : gauss points coordinates (on the reference unit
        %>element)
        %> @param obj : instance of the class
        %> @retval K : element stiffness matrix integrand part
        %------------------------------------------------------------------
        function [K]=K_Elt_Int(xil,obj)
            % Integrand part of the element stiffness matrix
            % Input
            % obj :: object laplacian operator for the elt
            % xil :: gaus quadrature point  (local coordinates)
            % OUTPUT
            %  K : element operator matrix
            %
            
         %   [B,j]=Bamat(xil,obj.EltObject);
            [DNaDx,j]=GradN(xil,obj.EltObject);
            
            Cond=obj.Cond;
            
            switch obj.Geometry
                
                case 'Axis'
                    
                    xaux=Mapx(xil,obj.EltObject);
                    
                    switch obj.EltObject.dim % 
                        case 1
                            
                           % Bcond=B(1,[ 1  3]);
                            K=(2*pi*xaux)*j*(DNaDx'*Cond*DNaDx);
                        case 2 % 2D axis
                            
%                             nn=length(B(1,:));
%                             
%                             Bcond=B(1,1:2:nn);
%                             
%                             Bcond=[Bcond; B(2,2:2:nn)];
%                             
%                             K=(2*pi*xaux(1))*j*(Bcond'*Cond*Bcond);
                            K=(2*pi*xaux(1))*j*(DNaDx'*Cond*DNaDx);
                    end
                    
                case '2D' %  
                    
%                     nn=length(B(1,:));
%                     
%                     Bcond=B(1,1:2:nn);
%                     
%                     Bcond=[Bcond; B(2,2:2:nn)];
%                     
%                     K=j*(Bcond'*Cond*Bcond);
                    
                    K=j*DNaDx'*Cond*DNaDx;
                    
            end
            
        end
        
    end
    
end

