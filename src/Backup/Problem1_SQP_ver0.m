global A B C xp d Q y_ref Vmax Imax w

% PMSM parameter
PP = 8;
Ld = 0.45*10^-3;
Lq = 0.66*10^-3;
LAMpm = 0.0563;
Rs = 0.025;

k1 = 1.5*PP*LAMpm;
k2 = 1.5*PP*(Ld - Lq);
k3 = k1/(2*k2);

Vdc = 98;
Vmax = Vdc/sqrt(3);
Imax = 50;
fs = 10*10^3;
Ts = 1/fs;

% Simulation setup
t_end = 0.02; % 0.02
time = (Ts:Ts:t_end)';
Ns = t_end*fs;

% PMSM variables
id = zeros(Ns,1);
iq = zeros(Ns,1);
u_opt_list = zeros(Ns,2);
Te = zeros(Ns,1);

% Operating conditions
%Te_ref = [60*ones(Ns/3,1); 120*ones(Ns/3,1) ; 60*ones(Ns/3,1)];
Te_ref = 34*ones(Ns,1);
Te_ref_LPF = zeros(Ns,1);
wc_LPF = 2*pi*200;
a_LPF = 1/(1 + wc_LPF*Ts);
Te_ref_lim = zeros(Ns,1);

wr = 1000*pi/30*PP + 2*5*min((2000*time - 0*4000*max(time - 3,0) + 0*2000*max(time - 4.8,0)),8000)*pi/30*PP;

% Control variables
w = 0.05*0.1;
computation_time = zeros(Ns,1);

%Te_ref_LPF(1) = 20; id(1) = -3; iq(1) = 28; u_opt_list(1,1) = -18; u_opt_list(1,2) = 47;
Te_ref_LPF(1) = 0; id(1) = 0; iq(1) = 0; u_opt_list(1,1) = 0; u_opt_list(1,2) = wr(1) * LAMpm;

for i = 2:Ns
    
    xp = [id(i - 1) + 0*2*(rand - 0.5)*Imax*0.015 ; iq(i - 1) + 0*2*(rand - 0.5)*Imax*0.015];
    up = u_opt_list(i - 1,:)';

%     Ld_nom = 1.2*Ld;
%     Lq_nom = 1.1*Lq;
%     LAMpm_nom = 0.9*LAMpm;
%     Rs_nom = Rs;    
    
    Ld_nom = Ld;
    Lq_nom = Lq;
    LAMpm_nom = LAMpm;
    Rs_nom = Rs;
    
    k1_nom = 1.5*PP*LAMpm_nom;
    k2_nom = 1.5*PP*(Ld_nom - Lq_nom);    
    
    d = [0 ; k1_nom];
    Q = [0 , k2_nom/2 ; k2_nom/2 , 0];
    
    Te_ref_LPF(i) = a_LPF*Te_ref_LPF(i - 1) + (1 - a_LPF)*Te_ref(i);
    y_ref = Te_ref_LPF(i);
    
    A = eye(2) + Ts*[-Rs_nom/Ld_nom , Lq_nom/Ld_nom*wr(i) ; -Ld_nom/Lq_nom*wr(i) , -Rs_nom/Lq_nom];
    B = Ts*diag([1/Ld_nom,1/Lq_nom]);
    C = Ts*[0 ; -LAMpm_nom/Lq_nom*wr(i)];

    tStart = tic;
    u_opt = myproblem(up);
    computation_time(i) = toc(tStart);
    u_opt_list(i,:) = u_opt;
    
    % PMSM modeling
    id(i) = id(i - 1) + Ts/Ld*(-Rs*id(i - 1) + wr(i)*Lq*iq(i - 1) + u_opt(1));
    iq(i) = iq(i - 1) + Ts/Lq*(-Rs*iq(i - 1) - wr(i)*(Ld*id(i - 1) + LAMpm) + u_opt(2));
    Te(i) = (k1 + k2*id(i))*iq(i);
end

% MTPA & MTPV trajactory
iq_MTPA = 0:1:Imax;
id_MTPA = - k3 - sqrt(k3^2 + iq_MTPA.^2);

MTPV_pos = zeros(2,Imax);
MTPV_neg = zeros(2,Imax);
for j = 1:Imax*10
    wr_MTPV = 3000*j/Imax;
    LAM_d = (-Lq*LAMpm + sqrt((Lq*LAMpm)^2 + 8*((Ld - Lq)*Vmax/wr_MTPV)^2))/(4*(Ld - Lq));
    MTPV_pos(:,j) = [(LAM_d - LAMpm)/Ld ; sqrt((Vmax/wr_MTPV)^2 - LAM_d^2)/Lq];
    MTPV_neg(:,j) = [(LAM_d - LAMpm)/Ld ; -sqrt((Vmax/wr_MTPV)^2 - LAM_d^2)/Lq];
end
MTPV = [MTPV_pos [-LAMpm/Ld ; 0] fliplr(MTPV_neg)];

figure(1)
plot(time,id,'b')
hold on
plot(time,iq,'r')
%plot(time,wr,'k')
plot(time,Te_ref_lim,'g--')
plot(time,Te,'c')
plot(time,Te_ref_LPF,'k--')
hold off

axis = [-Imax Imax -Imax Imax];

figure(2)
plot(id_MTPA,iq_MTPA,'k--')
hold on
plot(MTPV(1,:),MTPV(2,:),'k--')
plot(id(Ns - 1),iq(Ns - 1),'ro')
plot(id(Ns),iq(Ns),'go')
fimplicit(@(x,y) ((Rs*x - wr(Ns)*Lq*y).^2 + (Rs*y + wr(Ns)*(Ld*x + LAMpm)).^2 - Vmax^2), axis,'Color','g','Linestyle',':','LineWidth',1.5)
fimplicit(@(x,y) ((k1 + k2*x).*y - Te_ref_LPF(end)), axis,'Color','k','Linestyle','--','LineWidth',1.5)
fimplicit(@(x,y) (x.^2 + y.^2 - Imax^2), axis,'Color','b','Linestyle','-.','LineWidth',1.5)
plot(id,iq,'c')
hold off
xlim([-Imax 0.5*Imax])
ylim([-Imax Imax])

axis = [-Vmax Vmax -Vmax Vmax];

figure(3)
scatter(u_opt_list(:,1),u_opt_list(:,2),'g')
hold on
fimplicit(@(vd,vq) ((vd).^2 + (vq).^2 - Vmax^2), axis,'Color','g','Linestyle',':','LineWidth',1.5)
fimplicit(@(vd,vq) (a1*(vd) + a2*(vq) + a0), axis,'Color','r','Linestyle',':','LineWidth',1.5)
%fimplicit(@(vd,vq) (g1(1)*(vd) + g1(2)*(vq) + g0), axis,'Color','b','Linestyle',':','LineWidth',1.5)
hold off
grid on
% xlim([-2000 2000])
% ylim([-2000 2000])

figure(4)
plot(time,u_opt_list(:,1),'b')
hold on
plot(time,u_opt_list(:,2),'r')
hold off

figure(5)
plot(time,computation_time,'b')
ylim([0 0.01])

figure(6)
plot(time,wr)

Te_MPC1 = Te;
id_MPC1 = id;
iq_MPC1 = iq;
vd_MPC1 = u_opt_list(:,1);
vq_MPC1 = u_opt_list(:,2);
time_MPC1 = computation_time;

% save Te_MPC1.mat Te_MPC1
% save id_MPC1.mat id_MPC1
% save iq_MPC1.mat iq_MPC1
% save vd_MPC1.mat vd_MPC1
% save vq_MPC1.mat vq_MPC1
% save time_MPC1.mat time_MPC1

mean(computation_time) * 10^3

function u = myproblem(up)

options = optimset('Algorithm','sqp','Display','off','GradObj','on','GradConstr','on','MaxFunEvals',25,'MaxIter',10,'Hessian','bfgs'); %25 / 10
u = fmincon(@objfun,up,[],[],[],[],[],[],@nonlcon,options);

    function [f, g] = objfun(u)        
        global A B C xp d Q y_ref w
        xn = A*xp + B*u + C;
        f = (d'*xn + xn'*Q*xn - y_ref)^2 + w*(xn'*xn);
        if nargout > 1 % gradient required
            g = 2*B'*((d'*xn + xn'*Q*xn - y_ref)*(d + 2*Q*xn) + w*xn);
        end
    end

    function [c, ceq, DC, DCeq] = nonlcon(u)
        global A B C xp Vmax Imax
        xn = A*xp + B*u + C;
        u_ss = B\(eye(2) - A)*(xn) - B\C;
        c = [xn'*xn - Imax^2 ; u_ss'*u_ss - Vmax^2];
        ceq = [];
        
        if nargout > 2
            DC = [2*B'*xn 2*(B'*(eye(2) - A')*(inv(B))')*u_ss];
            DCeq = [];
        end
    end
end