classdef CoreShell_1DL < handle
    %Class describing a 1-D linear finite element
    
    properties
        % this is all of the stuff that defines the element
        num 
        nodes
        const
        nnpe
        dx
        force
        fhandle
        yp
        eleM
        eleC
        eleK
        ypOld
        qOld
        qdOld
        qNew
        qdNew
    end
    
    methods
        % Constructor... send in nodes and material properties
        function obj = CoreShell_1DL(numIn, nodeInfo, constIn, fhandle)
            obj.num = numIn;
            obj.nodes = nodeInfo;
            obj.const = constIn;
            obj.nnpe = 2;
            obj.dx = abs(obj.nodes(1).loc - obj.nodes(2).loc );
            obj.yp = 0;
            
            r1 = obj.nodes(1).loc;
            r2 = obj.nodes(2).loc;
            
            rcut = constIn(1);
            rhocp_c = constIn(2);
            rhocp_s = constIn(3);
            k_c = constIn(4);
            k_s = constIn(5);
            
            
            if obj.nodes(2).loc <= rcut
                rhocp = rhocp_c;
                k = k_c;
                disp('found a core element')
            elseif obj.nodes(2).loc > rcut
                rhocp = rhocp_s;
                k = k_s;
                disp('found a shell element')
            else
                moop
            end
                 
            %1/r^2 d/dr [ r^2 k dT/dr ]
            KGd11 = -(r1^2+r1*r2+r2^2)/(3*r1-3*r2);
            KGd12 = (r1^2+r1*r2+r2^2)/(3*r1-3*r2);
            KGd21 = KGd12;
            KGd22 = KGd11; 
            
            obj.eleK = k*[KGd11, KGd12 ; ...
                          KGd21, KGd22];
            
            C11 = -((r1-r2)*(6*r1^2 + 3*r1*r2 + r2^2))/30.;
            C12 = -((r1-r2)*(3*r1^2 + 4*r1*r2 + 3*r2^2))/60.;
            C21 = C12;
            C22 = -((r1-r2)*(r1^2 + 3*r1*r2 + 6*r2^2))/30.;
                                 
            obj.eleC = rhocp * [ C11 C12 ; ...
                                 C21 C22 ];
            obj.eleM = 0.0*obj.dx/6 * [ 0 0 ; ...
                                        0 0 ];
                                          
            obj.fhandle = fhandle;
            obj.force = [0;0];
 
        end
        
        function [] = view(obj)
            % displays the relevant information about the element
            
            fprintf('--------------------------------------------------------\n')
            fprintf('I am a 1-D linear element with the following properties:\n')
            fprintf('--------------------------------------------------------\n')
            fprintf('Element #: %i\n', obj.num)
            fprintf('Nodes: %i, %i\n', obj.nodes(1).num, obj.nodes(2).num)
            fprintf('M^(e): [ %f, %f ]\n', obj.eleK(1,1), obj.eleK(1,2) )
            fprintf('       [ %f, %f ]\n', obj.eleK(2,1), obj.eleK(2,2) )
            fprintf('K^(e): [ %f, %f ]\n', obj.eleK(1,1), obj.eleK(1,2) )
            fprintf('       [ %f, %f ]\n', obj.eleK(2,1), obj.eleK(2,2) )
            fprintf('C^(e): [ %f, %f ]\n', obj.eleC(1,1), obj.eleC(1,2) )
            fprintf('       [ %f, %f ]\n', obj.eleC(2,1), obj.eleC(2,2) )
            fprintf('y'': %f\n', obj.yp)
            fprintf('Force: [ %e ]\n', obj.force(1) )
            fprintf('       [ %e ]\n', obj.force(2) )
            
        end
        
        function [] = computeForce(obj,t)
            % evaluate integrals to get forcing function contributions to
            % each node
            % note the need for the 'ArrayValued' flag in order to handle a
            % call to integral() in the force term definition

            x1 = obj.nodes(1).loc;
            x2 = obj.nodes(2).loc;
            ff1 = @(x) (1-(x2-x)./obj.dx).*x*x*obj.fhandle(x,t);
            ff2 = @(x) (x2-x)./obj.dx.*x*x*obj.fhandle(x,t);
            obj.force = [ integral(ff1,x1,x2,'ArrayValued', true);
                          integral(ff2,x1,x2,'ArrayValued', true) ];
            
        end
        
        function [] = computeDir(obj)
            % define spatial derivative within element
            obj.yp = (obj.nodes(2).y - obj.nodes(1).y)/obj.dx;
            
        end
        
        function [] = computeFlux(obj,gs4,Kn)
            % this funtion computes the heat flux for models with a time
            % derivative appearing in their definition
        
            obj.qdNew = ( -Kn/3*(obj.ypOld+gs4.w1*(obj.yp-obj.ypOld))...
                     +(gs4.lam6w1+gs4.lam5w2*gs4.dt-gs4.w1*gs4.dt-1)*obj.qdOld ...
                     - obj.qOld)/(gs4.lam6w1+gs4.lam5w2*gs4.dt);
                 
            obj.qNew = obj.qOld + gs4.dt*obj.qdOld ...
                + gs4.lam5*gs4.dt*(obj.qdNew-obj.qdOld); 
            
            obj.qdOld = obj.qdNew;
            obj.qOld = obj.qNew;
            obj.ypOld = obj.yp;
            
        end
            
    end
    
end

