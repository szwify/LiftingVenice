%------------------------------------------------------------------
%> @brief High level class at the element level for a coupling operator
%> 
%> The Bi-Linear form is :
%>
%> \int_Elt \eps(u) \times c \times p d Elt
%>
%> Elt Matrix entry is:
%>
%> k_el= \int_Elt B'* c * Na d Elt
%>it requires 2 elts object (case of different interpolation)....
%------------------------------------------------------------------
classdef  Operator_Coupling_Matrices
    
    % class
    
    properties
        %> string Laplacian, Elasticity, Mass etc.
        ProblemType ;
        %> Element obj.
        EltObject_a ;   %  element obj.   field a
        EltObject_b ;   %  element obj.   field b (same geometry but not necessary the same interpolation)
        %> Properties object (properties of the corresponding elt)
        coupling_parameter;   %
        
    end
    
    methods
        
        %------------------------------------------------------------------
        %> @brief Class Constructor
        %>
        %> @param EltObject element object (instance of the class ...)
        %> @param ProblemType
        %> @param alpha :: coupling parameter  SCALAR only - no anisotropy.
        %------------------------------------------------------------------
        function obj=Operator_Coupling_Matrices(EltObjecta,EltObjectb,alpha)
            
            obj.EltObject_a=EltObjecta;
            obj.EltObject_b=EltObjectb;
            
            if (size(alpha)~= [1 1])
                disp('error anisotropy not yet in coupling operator');
                return;
                
            else
                
                obj.coupling_parameter=alpha;
            end
            
        end
        
        
        %------------------------------------------------------------------
        %> @brief Compute the integrand of the  Element ' stiffness' matrix
        %>entry   - 1st order coupling matrix
        %
        %> B'*b*Na @ xil.
        %>   for first order matrix....
        %
        %> @param xil : gauss points coordinates (on the reference unit
        %>element)
        %> @param obj : instance of the class
        %> @retval K : element stiffness matrix integrand part
        %------------------------------------------------------------------
        function [K]=K_Elt_Int(xil,obj)
            % Integrand part of the element stiffness matrix
            % Input
            % obj :: object  operator for the elt
            % xil :: gaus quadrature point  (local coordinates)
            % OUTPUT
            %  K : element operator matrix
            %
            
            [B,j]=Bamat(xil,obj.EltObject_a);
            
            [Naxi]=Na(xil,obj.EltObject_b);
            
            c_const=obj.coupling_parameter;
            
            switch obj.EltObject_a.type
                
                case 'Axis'
                    
                    xaux=Mapx(xil,obj.EltObject_b);
                    
                    switch obj.EltObject_a.dim % todo switch should be on dimension here not dof ! -
                        case 1
                            Baux=B(1,:)+B(2,:);
                            K=(2*pi*xaux)*c_const*j*(Baux'*Naxi);
                            
                        case 2 % 2D axis
                            
                            Baux=B(1,:)+B(2,:)+B(4,:);
                            K=(2*pi*xaux(1)*c_const)*j*(Baux'*Naxi);
                            
                    end
                    
                case '2D' %
                    
                    Baux=B(1,:)+B(2,:);
                    K=j*c_const*(Baux'*Naxi);
                    
            end
            
        end
        
    end
    
end


