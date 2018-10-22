%---------------------------------------------------------------------
%> @brief Provide nodal force terms associated with body source / sink for
%> Laplacian
%
%> It overload  the K_Elt_Int method to provide elt nodal force  vector
%>
%>
%---------------------------------------------------------------------
classdef Operator_Load_FluxSource < Operator_Laplacian
     %
     %  Class deriving from the Operator_Laplacian class
     %  provides way of computing  the nodals force terms associated
     %  with a body source / sink terms
     %  by overloading the K_Elt_Int method to provide elt nodal forces
     %  vector
     
    properties
       %> Body sink / source term (Laplacian problem)
        Qsource ; % body sink/source 
    end
    
    methods

        %------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param EltObject : Element object
        %> @param PropObject   : a properties object containing 'cond'
        %> properties
        %> @param ndof : number of dof per nodes
        %> @param Qsource : body source / sink value
        %------------------------------------------------------------------
        function obj=Operator_Load_FluxSource(EltObject,Cond,Qsource)
            
            obj=obj@Operator_Laplacian(EltObject,Cond);
            obj.Qsource=Qsource;
        end
        
        %------------------------------------------------------------------
        %> @brief Compute the integrand of the  element nodal forces
        %> contribution for a body sink/source
        %>
        %>
        %> @param xil : gauss points coordinates (on the reference unit
        %>element)
        %> @param obj : instance of the class
        %> @retval K : vector containing the corresponding nodal forces
        %> contribution for this element and the
        %>  given source
        %------------------------------------------------------------------
        function [K]=K_Elt_Int(xil,obj)
            
            [j]=Jacobian(xil,obj.EltObject);
            N=Na(xil,obj.EltObject);
            [nr nc]=size(N');
            q_flux=obj.Qsource;
            
            if (nc~=length(q_flux))
                disp('fluxes field does not coincides with options and dim');
                return;
            end
    
            switch   obj.EltObject.type
                
                case 'Axis'
                    
                    xaux=Mapx(xil,obj.EltObject);
                    
                     switch obj.EltObject.dim  
                        case 1 % 1d axis
                            K=(2*pi*xaux)*j*N'*q_flux;
                        case 2  % 2d axis
                            K=(2*pi*xaux(1))*N'*q_flux;
                    end
                    
                case '2D'
                    
                    K = j*N'*q_flux;  
                    
            end
        end
        
        
    end
    
    
end
