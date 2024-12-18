classdef env_SM
    properties
        Ts 

        PP
        Ld
        Lq
        LAMpm
        Rs
        k1
        k2
        k3
        Vdc
        Vmax
        Imax

        i % current
    end

    methods
        function obj = env_SM(Ts)
            obj.Ts = Ts;

            % PMSM parameter
            obj.PP = 8;
            obj.Ld = 0.45*10^-3;
            obj.Lq = 0.66*10^-3;
            obj.LAMpm = 0.0563;
            obj.Rs = 0.025;

            obj.k1 = 1.5*obj.PP*obj.LAMpm;
            obj.k2 = 1.5*obj.PP*(obj.Ld - obj.Lq);
            obj.k3 = obj.k1/(2*obj.k2);

            obj.Vdc = 98;
            obj.Vmax = obj.Vdc/sqrt(3);
            obj.Imax = 50;

            obj.i = [0; 0];
        end

        function nom_param = getNomParam(obj)
            nom_param.Ld_nom = obj.Ld;
            nom_param.Lq_nom = obj.Lq;
            nom_param.LAMpm_nom = obj.LAMpm;
            nom_param.Rs_nom = obj.Rs;
            
            nom_param.k1_nom = 1.5*obj.PP*nom_param.LAMpm_nom;
            nom_param.k2_nom = 1.5*obj.PP*(nom_param.Ld_nom - nom_param.Lq_nom);    
        end

        function obs_current = getObsCurrent(obj)
            obs_current = [
                obj.i(1) + 0*2*(rand - 0.5)*obj.Imax*0.015 ; 
                obj.i(2) + 0*2*(rand - 0.5)*obj.Imax*0.015
            ];
        end

        function real_info = getRealObs(obj)
            real_info.i = obj.i;
            real_info.Te = (obj.k1+obj.k2*obj.i(1))*obj.i(2);
        end

        function obj = step(obj, u_opt, wr)
            % Step Forward           
            id = obj.i(1) + obj.Ts/obj.Ld*(-obj.Rs*obj.i(1) ...
                + wr*obj.Lq*obj.i(2) + u_opt(1));
            iq = obj.i(2) + obj.Ts/obj.Lq*(-obj.Rs*obj.i(2) ...
                - wr*(obj.Ld*obj.i(1) + obj.LAMpm) + u_opt(2));

            obj.i = [id; iq];

        end

    end

end