%------------------------------------------------------------------
%> @brief Mass operator class at the element level
%
%> This class provides the FE implementation of the Mass Operator at
%> the element level.
%>
%>  Bi-linear form
%>
%>   \int_Elt u \times rho \times v d Elt
%>
%>  Element matrix
%>
%>  k_el = \int_Elt  N'* rho N   d Elt
%------------------------------------------------------------------
classdef Operator_Mass 
    
    properties
        %> Geometry
        Geometry; 
        %> EltObject
        EltObject;
        %> density like property (scalar)
        rho;

    end
    
    methods
        %------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param EltObject : Elemenet object
        %> @param PropObject   : a properties object containing 'density'
        %> properties
        %> @param ndof : number of dof per nodes
        %------------------------------------------------------------------
        function obj=Operator_Mass(EltObject,rho)

            obj.EltObject=EltObject;
            obj.Geometry = EltObject.type;
            obj.rho=rho;
            
        end
        
        %------------------------------------------------------------------
        %> @brief Compute the integrand of the  Element ' mass' matrix
        %>        entry
        %
        %> N'*Cond*N gauss pts @ xil.
        %>
        %
        %> @param xil : gauss points coordinates (on the reference unit
        %>element)
        %> @param obj : instance of the class
        %> @retval K : element mass matrix integrand part
        %------------------------------------------------------------------
        function [K]=K_Elt_Int(xil,obj)
            
            [j]=Jacobian(xil,obj.EltObject);
             
            [N]=Na(xil,obj.EltObject);
            
            switch obj.EltObject.type
                
                case 'Axis'
                    
                    switch obj.EltObject.dim %
                        % % %
                        case 1
                            % %
                            xaux=Mapx(xil,obj.EltObject);
                            
                            K=(2*pi*xaux)*j*(N'*obj.rho*N);
                        
                        case 2
                            
                            xaux=Mapx(xil,obj.EltObject);
                            
                            K=(2*pi*xaux(1))*j*(N'*obj.rho*N);
                            
                    end
                    
                case '2D' %
                    
                    K = j*(N'*obj.rho*N);
                    
            end
            
        end
        
    end
    
end

