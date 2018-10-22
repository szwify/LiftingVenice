%---------------------------------------------------------------------
%> @brief Provide nodal force terms associated with boundary load for
%> the elasticity eq.
%
%> It overload  the K_Elt_Int method of the Operator_Elasticity class to provide elt nodal force  vector
%>
%>
%> The Bi-Linear form is :
%>
%> \int_\partial_Elt  f \times \eps(v) d Elt
%>
%> Elt Matrix entry is :
%>
%> f_el= \int_Elt   f_i \times Na   d Elt
%
%> note that for a 2D problem with Qua4 element, the EltObject require here
%> is Seg2 for the EltProObjsect
%---------------------------------------------------------------------
classdef Operator_Load_Boundary  
    %
    
    properties
        
        %> a given force (Fx or Fy at the  nodes of an edge element....), vector
        %> of size nnodes
        Force_i;
        EltObject;
        
    end
    
    methods
        
        %------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param EltObject : Element object (note Seg2 if elasticity with
        %> Qua4)
        %> @param force : force applied at nodes (matrix with Fdof1 Fdof2
        %> for each nodes where a force is applied
        %------------------------------------------------------------------
        function obj=Operator_Load_Boundary(EltObject,force)
            
            obj.EltObject=EltObject;
            if (size(force(:,1))~=1)
                error(' Force vector is not a vector' );
            end
            
            obj.Force_i=force;
            
        end
        
        %------------------------------------------------------------------
        %> @brief Compute the integrand of the  element nodal forces
        %> contribution for point force
        %>
        %>
        %> @param xil : gauss points coordinates (on the reference unit
        %>element)
        %> @param obj : instance of the class
        %> @retval K : vector containing the corresponding nodal forces
        %> contribution for the segment element and the
        %------------------------------------------------------------------
        function [K]=K_Elt_Int(xil,obj)
            
            [j]=Jacobian(xil,obj.EltObject);
            [Na_xi]=Na(xil,obj.EltObject);
            nc = length(Na_xi);
            
            f=obj.Force_i; % vector 
            
            if (nc~=length(f(1,:)))
                disp([nc length(f(:,1))]);
                error('error in force component #of nodes');
            end
            
            
            switch   obj.EltObject.type
                
                case 'Axis'
                    % todo implementation of boundary force for the axisymmetric case.....
                    %                    TO BE CHECKED - Not working
                    %                    currently.....
                   
                    xaux=Mapx(xil,obj.EltObject);

                   % disp('----'); % there is a problem for axisymmm.... HERE.... 
                   % disp(xaux); 
                   % disp(obj.EltObject.dim);
                   % switch obj.EltObject.dim
                   %                          case 1 % unidimensional case % none integration by hand...  2 pi R * F
                   % 
                   %     case 2
%                      
%                      disp('val theta'); disp(obj.EltObject.theta); 
%                      disp(',,,');
                     Naf = Na_xi*f'; 
                     Faux= 2.*pi*(obj.EltObject.rc+j*xil/cos(obj.EltObject.theta))*j*Naf*Na_xi';  % 2 pi r(s) = 2 pi (rc+ he \xi cos \theta ) 
                     K=Faux;
                        
                   % end
                    
                case '2D'
                    
                    Naf = Na_xi*f';
                    Faux= j*Naf*Na_xi'; %
                    K=Faux; % column vector of length  nnodes
                    
            end
            
        end
        
    end
    
end
