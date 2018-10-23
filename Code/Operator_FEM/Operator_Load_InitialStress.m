%---------------------------------------------------------------------
%> @brief Provide nodal force terms associated with a pre-stress for
%> the elasticity eq.
%
%>
%>
%> The Bi-Linear form is :
%>
%> \int_Elt \sigma_o \times \eps(v) d Elt
%>
%> Elt Matrix entry is :
%>
%> f_el= \int_Elt B'* sigma_o d Elt
%---------------------------------------------------------------------
classdef Operator_Load_InitialStress 
    
    properties
        %> A given pre-stress field on the element (ct on the element)
        Stress_field;
      
        EltObject;
        % todo - make modification to have different value at each Gauss
        % Points...
    end
    
    
    methods
        
        %------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param EltObject : Element object
        %> @param Sig : element pre-stress field
        %------------------------------------------------------------------
        function obj=Operator_Load_InitialStress(EltObject,Sig)
            
            obj.EltObject=EltObject;

            obj.Stress_field=Sig;
            
        end
        
        %------------------------------------------------------------------
        %> @brief Compute the integrand of the  element nodal forces
        %> contribution for a pre-stress field
        %>
        %>
        %> @param xil : gauss points coordinates (on the reference unit
        %>element)
        %> @param obj : instance of the class
        %> @retval K : vector containing the corresponding nodal forces
        %> contribution for this element and the
        %>  given pre-stress field
        %------------------------------------------------------------------
        function [K]=K_Elt_Int(xil,obj)
            
            [B,j]=Bamat(xil,obj.EltObject);
            [nr nc]=size(B');
            Sig=obj.Stress_field; % ct stress field per element ...
            
            if (nc~=length(Sig))
                disp([nc length(Sig)]);
                error('error in initial stress field component');       
            end
            
            switch   obj.EltObject.type
                
                case 'Axis'
                    
                    xaux=Mapx(xil,obj.EltObject);
                    
                    switch obj.EltObject.dim
                        case 1  % unidimensional case
                            Faux=(2*pi*xaux)*j*B'*Sig;
                            K=Faux([1;3]);
                            
                        case 2
                            
                             Faux=(2*pi*xaux(1))*j*B'*Sig;
                             K=Faux;
                              
                    end
                    
                case '2D'
                    
                    Faux= j*B'*Sig;
                    K=Faux; % column vector of length 2*nnodes
                   
            end
            
        end
        
    end 
    
end
