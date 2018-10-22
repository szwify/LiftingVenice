%------------------------------------------------------------------
%> @brief FE implementation of the Elasticity Operator
%
%> This class provides the FE implementation of the Elasticity Operator at
%> the element level.
%>
%> The Bi-Linear form is :
%>
%> \int_Elt \eps(u) \times L \times \eps(v) d Elt
%>
%> Elt Matrix entry is :
%>
%> k_el= \int_Elt B'* L * B d Elt
%------------------------------------------------------------------
classdef Operator_Elasticity  
    
    properties
        %> Geometry
        Geometry; 
        %> EltObject
        EltObject;
        %> Stiffness tensor
        L;
    end
    
    methods
        %------------------------------------------------------------------
        %> @brief Class Constructor
        %> @param EltObject : Elemenet object
        %> @param PropObject   : a properties object containing elastic
        %> properties
        %> @param ndof : number of dof per nodes
        %------------------------------------------------------------------
        function obj=Operator_Elasticity(EltObject,L)
            
            obj.EltObject=EltObject;
            obj.Geometry = EltObject.type;
            obj.L=L;
            
        end
        %------------------------------------------------------------------
        %>  @brief Compute the integrand of the  Element stiffness matrix
        %>entry
        %
        %> B'.L.B @ gauss pts xil.
        %> For elastic problem, only plane-strain axisymmetric problem are
        %> implemented for now
        %>
        %> @param xil : gauss points coordinates (on the reference unit
        %>element)
        %> @param obj : instance of the class
        %> @retval K : element stiffness matrix integrand part
        %------------------------------------------------------------------
        function [K]=K_Elt_Int(xil,obj)
            
            [B,j]=Bamat(xil,obj.EltObject);
            L=obj.L;
            
            switch obj.EltObject.type
                
                case 'Axis'
                    
                    xaux=Mapx(xil,obj.EltObject);
                    
                    switch obj.EltObject.dim
                        case 1
                            
                            K_aux=(2*pi*xaux)*j*(B'*L(1:3,1:3)*B);
                            K=K_aux([1; 3],[1 ;3]);
                            
                        case 2 % 2d axis.
                            %Lp=obj.PropObject.Lp;
                            Lax=[L(1,1) L(1,2) L(1,4) L(1,3) ;...
                                L(2,1) L(2,2) L(2,4)    L(2,3) ; ...
                                L(4,1)  L(4,2)  L(4,4)/2.  L(4,3)  ;...
                                L(3,1)  L(3,2)  L(3,4)      L(3,3);];
                            
                            %                             Lp=obj.PropObject.Lp;
                            %                            % Lax=[    ]; % to tune up
                            K=(2*pi*xaux(1))*j*(B'*Lax(1:4,1:4)*B);
                            
                    end
                    
                case '2D'
                    %
                    %                     Lp = L(1:3,1:3) ;% valid for isotropy only.
                    %                     Lp(3,:)=0.;L(:,3)=0.;Lp(3,3)=L(4,4)/2.;
                     
                    K =j*(B'*obj.L*B);
                    
            end
            
        end
        
        
        %------------------------------------------------------------------
        %> @brief Compute stress and strain from displacement @ point
        %> inside the elemet
        %
        %> This fction compute at a point inside the element (local
        %> coordinates on the reference element) the stress and strain for
        %> a given displacement and initial stress field
        %> @param xil: pt inside the elt (local coordinates on reference elt)
        %> @param obj : object elasticity operator for the elt
        %> @param  Ue: displacement vector of the elt nodes
        %> @param Sig_o : initial stress field on the element
        %> @retval Sg :  stress @ xil in vector form
        %> @retval  St :  strain @ xil in vector form
        %------------------------------------------------------------------
        function [Sg,St]=StressStrainInsideElt(xil,obj,Ue,Sig_o)
            
            [B,j]=Bamat(xil,obj.EltObject);
            L=obj.L;
            
            switch obj.EltObject.type
                
                case 'Axis' % axysym
                    
                    switch obj.EltObject.dim
                        
                        case 1   %1D axisym
                            
                            [nr nc] = size(Ue);
                            if nr ==1;
                                Ue=Ue';
                            end
                            
                            aux=L(1:3,1:3)*B;
                            St=(B(1:3,[1; 3])*Ue)'; %%% strain
                            
                            [nr nc]=size(Sig_o(1:3));
                            
                            Sg =(aux(1:3,[1 ;3])*Ue)'+Sig_o(1:3);
                            
                        case 2    % 2D axi-sym here
                            Lax=[L(1,1) L(1,2) L(1,4) L(1,3) ;...
                                L(2,1) L(2,2) L(2,4)    L(2,3) ; ...
                                L(4,1)  L(4,2)  L(4,4)/2.  L(4,3)  ;...
                                L(3,1)  L(3,2)  L(3,4)      L(3,3);];
                            [nr nc] = size(Ue);
                            if nr ==1;
                                Ue=Ue';
                            end
                            
                            %  aux=Lax*B;
                            St=(B*Ue); %%% strain
                            
                            [nr nc]=size(Sig_o(1:4));
                            
                            Sg =(Lax*St)'+Sig_o(1:4);
                            
                    end
                    
                case '2D'
                    
                    Lp=obj.L; % plane-strain elastic stiff matrix
                    
                    [nr nc] = size(Ue);
                    if nr ==1;
                        Ue=Ue';
                    end
                    
                    aux=Lp*B;
                    St=(B*Ue); %%% strain
                    
                    [nr nc]=size(Sig_o(1:3));
                    
                    Sg =(Lp*St)'+Sig_o(1:3);
                    
            end
            
            
        end
        
        %------------------------------------------------------------------
        %> @brief Compute stress and strain at gauss pts in a given elt
        %>from solution vector
        %
        %> This method computes the Stress and strain at all gauss pts in
        %> an element from the displacement and initial stress field
        %>
        %> @param obj: elasticity operator object for the elt
        %> @param Ue:  Displacement at element nodes
        %> @param Sig_o : Initial stress (ct over the elt)
        %> @param IntegrationOrder : order of intg -> furnish gauss points
        %> @retval Sg : Stress at all gauss pts  matrix of size nb of gauss pts
        %> \times nstress components (3 for axis)
        %> @retval St : Strain at all gauss pts matrix of size nb of gauss pts
        %> \times nstrain components (3 for axis)
        %------------------------------------------------------------------
        function [Sg,St]=StressStrainAtGaussPoints(obj,Ue,Sig_o,IntegrationOrder)
            
            
            dim=obj.EltObject.dim;
            
            % create the right object
            
            switch dim
                case 1
                    
                    %one dimensional integration
                    [xeta,Wl]=GaussQuadratureCoefs(IntegrationOrder);
                    Sg=zeros(length(Wl),3);
                    St=zeros(length(Wl),3);
                    for i=1:length(Wl),
                        [Sg(i,:) St(i,:)]=StressStrainInsideElt(xeta(i),obj,Ue,Sig_o);
                    end
                    
                case 2
                    %2 dimensional integration
                    [xeta,Wl]=GaussQuadratureCoefs(IntegrationOrder);
                    
                    switch obj.EltObject.type
                        case '2D'
                            Sg=zeros(IntegrationOrder*IntegrationOrder,3);
                            St=zeros(IntegrationOrder*IntegrationOrder,3);
                        case 'Axis'
                            
                            Sg=zeros(IntegrationOrder*IntegrationOrder,4);
                            St=zeros(IntegrationOrder*IntegrationOrder,4);
                    end
                    k=0;
                    for i=1:IntegrationOrder
                        for j=1:IntegrationOrder
                            k=k+1;
                            [Sg(k,:) St(k,:)]=StressStrainInsideElt([xeta(i) xeta(j)],obj,Ue,Sig_o);
                        end
                    end
                    
                    
            end
            
        end
        
        %------------------------------------------------------------------
        %> @brief This method compute the average of GaussPts stress and
        %> strain within a given elt
        %
        %> This method compute the average of GaussPts stress and
        %> strain within a given elt from the displacement and initial
        %> stress field
        %> @param obj: elasticity operator object for the elt
        %> @param Ue:  Displacement at element nodes
        %> @param Sig_o : Initial stress (ct over the elt)
        %> @param IntegrationOrder : order of intg -> furnish gauss points
        %> @retval Sg :Average Stress in the elt  matrix of size nb of gauss pts
        %> \times nstress components (3 for axis)
        %> @retval St :  Average strain in the elt matrix of size nb of gauss pts
        %> \times nstrain components (3 for axis)
        %> @retval avgCoor : average global coordinates (i.e. center) of the elt
        %------------------------------------------------------------------
        function [Sg,St,avgCoor]=AverageStressStrainElt(obj,Ue,Sig_o,IntegrationOrder)
            
            [SgG,StG]=StressStrainAtGaussPoints(obj,Ue,Sig_o,IntegrationOrder);
            Sg=mean(SgG);
            St=mean(StG);
            
            avgCoor=mean(obj.EltObject.xae);
            
        end
    end
    
end

