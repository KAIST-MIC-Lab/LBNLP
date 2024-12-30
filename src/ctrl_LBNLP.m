classdef ctrl_LBNLP
    properties
        Nc = 1     % control horizon
        ctrl_Ts

        A0
        A 
        B 
        C
        xp 
        d 
        Q 
        y_ref 

        isActive  = 0;

        w = 1e-0; %1e-4, 1e-2 ,1e-0
        w_ref
        up % preserve for next step

        k1_nom
        k2_nom

        Vdc
        Vmax
        Imax

        Binv
        xn0

        sc1
        sc2

        lam = [0 ; 0];
    end

    methods
        function obj = ctrl_LBNLP(Ts)
            obj.ctrl_Ts = Ts;

            Vdc = 98;
            obj.Vmax = Vdc/sqrt(3);
            obj.Imax = 50;
            
            obj.up = [0;0];

            obj.sc1 = 40/obj.Vmax^2 * 1;
            obj.sc2 = 1000/obj.Imax^2;
        end

        function obj = preControl(obj, nom_param)

            Ld_nom = nom_param.Ld_nom;
            Lq_nom = nom_param.Lq_nom;
            LAMpm_nom = nom_param.LAMpm_nom;
            Rs_nom = nom_param.Rs_nom;
            obj.k1_nom = nom_param.k1_nom;
            obj.k2_nom = nom_param.k2_nom;

            Ts = obj.ctrl_Ts;
            wr = obj.w_ref;

            obj.d = [0 ; obj.k1_nom];
            obj.Q = [0 , obj.k2_nom/2 ; obj.k2_nom/2 , 0];

            obj.A0 = Ts*[-Rs_nom/Ld_nom , Lq_nom/Ld_nom*wr ; -Ld_nom/Lq_nom*wr , -Rs_nom/Lq_nom];
            obj.A = eye(2) + obj.A0;
            obj.B = Ts*diag([1/Ld_nom,1/Lq_nom]);
            obj.Binv = diag([Ld_nom/Ts , Lq_nom/Ts]);
            obj.C = Ts*[0 ; -LAMpm_nom/Lq_nom*wr];
        
            obj.xn0 = obj.A*obj.xp + obj.C;
        end

        function obj = controlIteration(obj)

            xn = obj.xn0 + obj.B*obj.up;

            ceq = obj.k1_nom*xn(2) + obj.k2_nom*xn(1)*xn(2) - obj.y_ref;
            uss = obj.Binv*(-obj.A0*xn - obj.C);
            cin = -(obj.Vmax^2 - uss'*uss)*obj.sc1;
    
            if obj.isActive  > 0
                ca = [ceq ; cin];
                Ca = [
                        obj.B'*(obj.d + 2*obj.Q*xn), ...
                        2*obj.B'*(-obj.A0')*(obj.Binv)'*uss*obj.sc1
                    ];
                grad_L = 2*obj.B'*(xn) + (obj.lam(1)*Ca(:,1) + obj.lam(2)*Ca(:,2));
                W = 2*obj.B'*(eye(2) + obj.lam(1)*obj.Q + ...
                    obj.lam(2)*((-obj.A0')*(obj.Binv)'*obj.Binv*(-obj.A0)*obj.sc1))*obj.B;
            else
                ca = ceq;
                Ca = obj.B'*(obj.d + 2*obj.Q*xn);
                grad_L = 2*obj.B'*(xn) + (obj.lam(1)*Ca(:,1));
                W = 2*obj.B'*(eye(2) + obj.lam(1)*(obj.Q))*obj.B;
            end
    
            temp = eig(W);
            % temp2 = min(abs(temp(1)),abs(temp(2)));
            temp3 = eig(inv(W));            
    
            if min(temp) >= 0
                beta = min(temp);
            elseif max(temp) <= 0
                beta = (max(temp));
            else
                beta = 0;
            end
            % beta = 0
        
            % beta
            % beta = - norm(W) - 100;
            % beta = -beta
            
            alp = 0.5;

            del_x0 = -alp*(W'*grad_L + (Ca)*ca);
            del_lam0 = beta*ca - 1*Ca'*grad_L/2;
            
            del_g0 = (W*del_x0 + Ca*del_lam0);
            del_c0 = Ca'*del_x0;

            l = -1*(del_g0'*grad_L + del_c0'*ca)/(del_g0'*del_g0 + del_c0'*del_c0);
            
            K = [W , Ca ; Ca', - 2/(temp3(1) + temp3(2))*eye(length(ca))];
            % eig(K)
            diff = l*[del_x0 ; del_lam0];

            if obj.isActive  > 0
                sol = [obj.up ; obj.lam] + diff;
                obj.up = sol(1:2);
                obj.lam = sol(3:4);
            else
                sol = [obj.up ; obj.lam(1)] + diff;
                obj.up = sol(1:2);
                obj.lam(1) = sol(3);
            end
    
            %lam(2:end) = max(lam(2:end),0);
            obj.lam(2) = max(obj.lam(2),0);
    
            if obj.isActive  > 0 && obj.lam(2) == 0
                obj.isActive  = 0;
            elseif obj.isActive  == 0 && cin > 0
                obj.isActive  = 1;
            end
        end

        function obj = postControl(obj)

        end 

        function [obj, comp_time, u] = getControl(obj, current, w_ref, y_ref, nom_param)
            obj.xp = current;
            obj.w_ref = w_ref;
            obj.y_ref = y_ref;
            
            obj = obj.preControl(nom_param);

            tic
            for iter = 1:1:obj.Nc
                obj = obj.controlIteration();
            end
            comp_time = toc;

            u = obj.up;

            obj.postControl();
        end


    end
end