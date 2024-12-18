classdef ctrl_SQP
% THe code is from Problem1_SQP_ver1.m


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
    end

    methods
        function obj = ctrl_SQP(Ts, w) 
            obj.ctrl_Ts = Ts;
            obj.w = w;

            obj.u_prev = [0; 0];

            obj.Vdc = 98;
            obj.Vmax = obj.Vdc/sqrt(3);
            obj.Imax = 50;
        end

        function obj = preControl(obj, nom_param)

            Ld_nom = nom_param.Ld_nom;
            Lq_nom = nom_param.Lq_nom;
            LAMpm_nom = nom_param.LAMpm_nom;
            Rs_nom = nom_param.Rs_nom;
            k1_nom = nom_param.k1_nom;
            k2_nom = nom_param.k2_nom;

            obj.d = [0 ; k1_nom];
            obj.Q = [0 , k2_nom/2 ; k2_nom/2 , 0];

            obj.A = eye(2) + obj.ctrl_Ts*[-Rs_nom/Ld_nom , Lq_nom/Ld_nom*obj.w_ref ; -Ld_nom/Lq_nom*obj.w_ref , -Rs_nom/Lq_nom];
            obj.B = obj.ctrl_Ts*diag([1/Ld_nom,1/Lq_nom]);
            obj.C = obj.ctrl_Ts*[0 ; -LAMpm_nom/Lq_nom*obj.w_ref];
        end

    
        function [obj, comp_time, u] = getControl(obj, current, w_ref, y_ref, nom_param)
            obj.xp = current;
            obj.w_ref = w_ref;
            obj.y_ref = y_ref;

            obj = obj.preControl(nom_param);

            up = obj.u_prev;
        
            options = optimset('Algorithm','sqp','Display','off','GradObj','on','GradConstr','on','MaxFunEvals',25,'MaxIter',10,'Hessian','bfgs'); %25 / 10
            
            tStart = tic;
            u = fmincon(@(x) objfun(obj, x), ...
                up,[],[],[],[],[],[], ...
                @(x) nonlcon(obj, x), ...
                options());
            comp_time = toc(tStart);
            
            obj.u_prev = u;

        end

    end
    

end

%% LOCAL FUNCTIONS
function [f, g] = objfun(obj, u)        
    A = obj.A;
    B = obj.B;
    C = obj.C;
    xp = obj.xp;
    d = obj.d;
    Q = obj.Q;
    y_ref = obj.y_ref;
    w = obj.w;

    xn = A*xp + B*u + C;
    f = (d'*xn + xn'*Q*xn - y_ref)^2 + w*(xn'*xn);
    if nargout > 1 % gradient required
        g = 2*B'*((d'*xn + xn'*Q*xn - y_ref)*(d + 2*Q*xn) + w*xn);
    end
end

function [c, ceq, DC, DCeq] = nonlcon(obj, u)
    A = obj.A;
    B = obj.B;
    C = obj.C;
    xp = obj.xp;
    % d = obj.d;
    % Q = obj.Q;
    % y_ref = obj.y_ref;
    % w = obj.w;

    Imax = obj.Imax;
    Vmax = obj.Vmax;
    
    
    xn = A*xp + B*u + C;
    u_ss = B\(eye(2) - A)*(xn) - B\C;
    c = [xn'*xn - Imax^2 ; u_ss'*u_ss - Vmax^2];
    ceq = [];
    
    if nargout > 2
        DC = [2*B'*xn 2*(B'*(eye(2) - A')*(inv(B))')*u_ss];
        DCeq = [];
    end
end