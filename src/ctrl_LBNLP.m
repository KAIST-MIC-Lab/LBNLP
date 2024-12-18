classdef ctrl_LBNLP
    properties
        ctrl_Ts

        A 
        B 
        C
        xp 
        d 
        Q 
        y_ref 
        
        w = 1e-0; %1e-4, 1e-2 ,1e-0
        w_ref
        u_prev % preserve for next step

        Vdc
        Vmax
        Imax

        Binv
        xn0
    end

    methods
        function obj = ctrl_LBNLP(Ts)
            obj.ctrl_Ts = Ts;

        end

        function obj = preControl(obj, nom_param)

            Ld_nom = nom_param.Ld_nom;
            Lq_nom = nom_param.Lq_nom;
            LAMpm_nom = nom_param.LAMpm_nom;
            Rs_nom = nom_param.Rs_nom;
            k1_nom = nom_param.k1_nom;
            k2_nom = nom_param.k2_nom;

            Ts = obj.ctrl_Ts;
            wr = obj.w_ref;

            obj.d = [0 ; k1_nom];
            obj.Q = [0 , k2_nom/2 ; k2_nom/2 , 0];

            A0 = Ts*[-Rs_nom/Ld_nom , Lq_nom/Ld_nom*wr ; -Ld_nom/Lq_nom*wr , -Rs_nom/Lq_nom];
            obj.A = eye(2) + A0;
            obj.B = Ts*diag([1/Ld_nom,1/Lq_nom]);
            obj.Binv = diag([Ld_nom/Ts , Lq_nom/Ts]);
            obj.C = Ts*[0 ; -LAMpm_nom/Lq_nom*wr];
        
            obj.xn0 = obj.A*obj.xp + obj.C;
        end

        function obj = controlIteration(obj)

            xn0 = obj.xn0;


            xn = xn0 + B*up;
            %ceq = (d'*xn + xn'*Q*xn - y_ref);
            ceq = k1_nom*xn(2) + k2_nom*xn(1)*xn(2) - y_ref;
            %uss = B\((eye(2) - A)*xn - C);
            uss = Binv*(-A0*xn - C);
            %cin = [(Vmax^2 - uss'*uss)*sc1 ; (Imax^2 - xn'*xn)*sc2];
            cin = -(Vmax^2 - uss'*uss)*sc1;
    
            if as > 0
                ca = [ceq ; cin];
                %Ca = [B'*(d + 2*Q*xn) , -2*B'*(eye(2) - A')*(Binv)'*uss*sc1];
                Ca = [B'*(d + 2*Q*xn) , 2*B'*(-A0')*(Binv)'*uss*sc1];
                grad_L = 2*B'*(xn) + (lam(1)*Ca(:,1) + lam(2)*Ca(:,2));
                %W = 2*B'*(eye(2) - lam(1)*Q - lam(2)*(-(eye(2) - A')*(Binv)'*Binv*(eye(2) - A)*sc1))*B;
                W = 2*B'*(eye(2) + lam(1)*Q + lam(2)*((-A0')*(Binv)'*Binv*(-A0)*sc1))*B;
            else
                ca = ceq;
                Ca = B'*(d + 2*Q*xn);
                grad_L = 2*B'*(xn) + (lam(1)*Ca(:,1));
                W = 2*B'*(eye(2) + lam(1)*(Q))*B;
            end
    
            temp = eig(W);
            temp2 = min(abs(temp(1)),abs(temp(2)));
            temp3 = eig(inv(W));
            %temp2 = norm(W);
    
            if min(temp) >= 0
                beta = min(temp);
            elseif max(temp) <= 0
                beta = (max(temp));
            else
                beta = 0;
            end
            %beta = 0
            del_x0 = -0.5*(W'*grad_L + (Ca)*ca);
            del_lam0 = beta*ca - 1*Ca'*grad_L/2;
            
            del_g0 = (W*del_x0 + Ca*del_lam0);
            del_c0 = Ca'*del_x0;
            l = -1*(del_g0'*grad_L + del_c0'*ca)/(del_g0'*del_g0 + del_c0'*del_c0);
            
            K = [W , Ca ; Ca', - 2/(temp3(1) + temp3(2))*eye(length(ca))];
            eig(K)
            diff = l*[del_x0 ; del_lam0];
            if as > 0
                sol = [up ; lam] + diff;
                up = sol(1:2);
                lam = sol(3:4);
            else
                sol = [up ; lam(1)] + diff;
                up = sol(1:2);
                lam(1) = sol(3);
            end
    
            %lam(2:end) = max(lam(2:end),0);
            lam(2) = max(lam(2),0);
    
            if as > 0 && lam(2) == 0
                as = 0;
            elseif as == 0 && cin > 0
                as = 1;
            end
        end

        function obj = getControl(obj, current, w_ref, y_ref, nom_param)
            obj.xp = current;
            obj.w_ref = w_ref;
            obj.y_ref = y_ref;
            
            obj = obj.preControl(nom_param);

            for iter = 1:1
                obj = obj.controlIteration();
            end
        end

    end
end