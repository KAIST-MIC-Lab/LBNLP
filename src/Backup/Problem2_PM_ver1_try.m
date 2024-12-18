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
t_end = 0.002; % 0.02
time = (Ts:Ts:t_end)';
Ns = t_end*fs;

% PMSM variables
id = zeros(Ns,1);
iq = zeros(Ns,1);
u_opt_list = zeros(Ns,2);
Te = zeros(Ns,1);

% Operating conditions
%Te_ref = [60*ones(Ns/3,1); 120*ones(Ns/3,1) ; 60*ones(Ns/3,1)];
%Te_ref = 34*ones(Ns,1);
Te_ref = 30*ones(Ns,1);
%Te_ref = 28.9*ones(Ns,1);
Te_ref_LPF = zeros(Ns,1);
wc_LPF = 2*pi*500;
a_LPF = 1/(1 + wc_LPF*Ts);
Te_ref_lim = zeros(Ns,1);

wr = 1000*pi/30*PP + 2*1.5*0.1*5*min((2000*time - 0*4000*max(time - 3,0) + 0*2000*max(time - 4.8,0)),8000)*pi/30*PP;
%wr = min(wr,1150);
wr = 1089.1*ones(Ns,1);

% Control variables
computation_time = zeros(Ns,1);

%Te_ref_LPF(1) = 20; id(1) = -3; iq(1) = 28; u_opt_list(1,1) = -18; u_opt_list(1,2) = 47;
Te_ref_LPF(1) = 0; id(1) = 0; iq(1) = 0; u_opt_list(1,1) = 0; u_opt_list(1,2) = wr(1) * LAMpm;

% PM design parameters
lam = [0 ; 0];
u_opt_list = zeros(Ns,2);
lam_list = zeros(Ns,2);
V_list = zeros(Ns,3);
Log_list = zeros(Ns,3);

%sc1 = 1000/Vmax^2;
%sc2 = 1000/Imax^2;
sc1 = 40/Vmax^2 * 1;
sc2 = 1000/Imax^2;

grad_L = [0 ; 1];
ca = [0; 0];
as = 0;

for i = 2:Ns*0.8
    %i
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

    A0 = Ts*[-Rs_nom/Ld_nom , Lq_nom/Ld_nom*wr(i) ; -Ld_nom/Lq_nom*wr(i) , -Rs_nom/Lq_nom];
    A = eye(2) + A0;
    B = Ts*diag([1/Ld_nom,1/Lq_nom]);
    Binv = diag([Ld_nom/Ts , Lq_nom/Ts]);
    C = Ts*[0 ; -LAMpm_nom/Lq_nom*wr(i)];

    xn0 = A*xp + C;

    tStart = tic;
    %StartTime = cputime;
    for iter = 1:2
        xn = xn0 + B*up;
        %ceq = (d'*xn + xn'*Q*xn - y_ref);
        ceq = k1_nom*xn(2) + k2_nom*xn(1)*xn(2) - y_ref;
        %uss = B\((eye(2) - A)*xn - C);
        uss = Binv*(-A0*xn - C);
        %cin = [(Vmax^2 - uss'*uss)*sc1 ; (Imax^2 - xn'*xn)*sc2];
        cin = (Vmax^2 - uss'*uss)*sc1;

        if as > 0
            ca = [ceq ; cin];
            %Ca = [B'*(d + 2*Q*xn) , -2*B'*(eye(2) - A')*(Binv)'*uss*sc1];
            Ca = [B'*(d + 2*Q*xn) , -2*B'*(-A0')*(Binv)'*uss*sc1];
            grad_L = 2*B'*(xn) - (lam(1)*Ca(:,1) + lam(2)*Ca(:,2));
            %W = 2*B'*(eye(2) - lam(1)*Q - lam(2)*(-(eye(2) - A')*(Binv)'*Binv*(eye(2) - A)*sc1))*B;
            W = 2*B'*(eye(2) - lam(1)*Q - lam(2)*(-(-A0')*(Binv)'*Binv*(-A0)*sc1))*B;
        else
            ca = ceq;
            Ca = B'*(d + 2*Q*xn);
            grad_L = 2*B'*(xn) - (lam(1)*Ca(:,1));
            W = 2*B'*(eye(2) - lam(1)*(Q))*B;
        end

        del_x0 = -0.5*(W'*grad_L + (Ca)*ca);
        del_lam0 = -pinv(Ca)*W*Ca*ca + Ca'*grad_L/2;
        del_g0 = (W*del_x0 - Ca*del_lam0);
        del_c0 = Ca'*del_x0;
        l = -(del_g0'*grad_L + del_c0'*ca)/(del_g0'*del_g0 + del_c0'*del_c0);

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
        elseif as == 0 && cin < 0
            as = 1;
        end
    end
    computation_time(i) = toc(tStart);
    %computation_time(i) = cputime - StartTime;

    V = 0.5*(grad_L'*grad_L) + 0.5*(ca'*ca);

    u_opt = up;
    u_opt_list(i,:) = u_opt;
    lam_list(i,:) = lam;
    V_list(i,:) = [V ; 0.5*grad_L'*grad_L ; 0.5*ca'*ca];
    %Log_list(i,:) = [lmax ; del_g0'*grad_L + del_c0'*ca ; del_g0'*del_g0 + del_c0'*del_c0];
    Log_list(i,:) = [cin ; l ; Imax];

    % PMSM modeling
    id(i) = id(i - 1) + Ts/Ld*(-Rs*id(i - 1) + wr(i)*Lq*iq(i - 1) + u_opt(1));
    iq(i) = iq(i - 1) + Ts/Lq*(-Rs*iq(i - 1) - wr(i)*(Ld*id(i - 1) + LAMpm) + u_opt(2));
    Te(i) = (k1 + k2*id(i))*iq(i);
end

color_blue = [0 0.4470 0.7410];
color_yellow = [0.9290 0.6940 0.1250];
color_sky = [0.3010 0.7450 0.9330];

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

figure(21)
%plot(MTPV(1,:),MTPV(2,:),'k--')
%plot(id(i - 1),iq(i - 1),'ro')
%plot(id(i),iq(i),'go')
plot(x1,Te_ref_LPF(i - 1)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
hold on
fimplicit(@(x,y) ((Rs*x - wr(i)*Lq*y).^2 + (Rs*y + wr(i)*(Ld*x + LAMpm)).^2 - Vmax^2), axis,'Color',color_yellow,'Linestyle','-','LineWidth',1.5)
plot(id_MTPA,iq_MTPA,'k--')%,'LineWidth',1.5)
%fimplicit(@(x,y) (x.^2 + y.^2 - Imax^2), axis,'Color','b','Linestyle','-.','LineWidth',1.5)
plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
x1 = -50:30;
for ii = 1:1:i - 2
%fimplicit(@(x,y) ((k1 + k2*x).*y - Te_ref_LPF(ii)), axis,'Color',[0 0 0 0.9],'Linestyle','-','LineWidth',2.5)
plot(x1,Te_ref_LPF(ii)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
end
plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
p1 = [-5 15];                         % First Point
p2 = [-10 35];                         % Second Point
dp = p2-p1;                         % Difference
quiver(p1(1),p1(2),dp(1),dp(2),0,'Color',color_blue)
t1 = text(-30,25,'$t = 0,1,2,\cdots$','FontSize',16);
set(t1,'Interpreter','latex')
hold off
xlim([-Imax 0.5*Imax])
ylim([-0.2*Imax Imax])
%title('Vehicle operation')
xl1 = xlabel('$x_1$');
yl1 = ylabel('$x_2$');
%xticks(0:100:500)
%title('Loss of Optimality (%)')
set(xl1,'Interpreter','latex')
set(yl1,'Interpreter','latex')
set(gcf,'color','w')
l1 = legend('$e_{t+N} = 0$','$h(x,u) = 0$','MTPA line','$x_t$');
set(l1,'Interpreter','latex','location','southwest','Orientation','vertical');
grid on
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;

figure(22)
%plot(MTPV(1,:),MTPV(2,:),'k--')
%plot(id(i - 1),iq(i - 1),'ro')
%plot(id(i),iq(i),'go')
plot(x1,Te_ref_LPF(i - 1)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
hold on
fimplicit(@(x,y) ((Rs*x - wr(i)*Lq*y).^2 + (Rs*y + wr(i)*(Ld*x + LAMpm)).^2 - Vmax^2), axis,'Color',color_yellow,'Linestyle','-','LineWidth',1.5)
plot(id_MTPA,iq_MTPA,'k--')
%fimplicit(@(x,y) (x.^2 + y.^2 - Imax^2), axis,'Color','b','Linestyle','-.','LineWidth',1.5)
plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
x1 = -50:30;
for ii = 1:1:i - 2
%fimplicit(@(x,y) ((k1 + k2*x).*y - Te_ref_LPF(ii)), axis,'Color',[0 0 0 0.9],'Linestyle','-','LineWidth',2.5)
plot(x1,Te_ref_LPF(ii)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
end
plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
p1 = [5 15];                         % First Point
p2 = [-15 35];                         % Second Point
dp = p2-p1;                         % Difference
quiver(p1(1),p1(2),dp(1),dp(2),0,'Color',color_blue)
t1 = text(0,25,'$t = 0,1,2,\cdots$','FontSize',16);
set(t1,'Interpreter','latex')
hold off
xlim([-Imax 0.5*Imax])
ylim([-0.2*Imax Imax])
%title('Vehicle operation')
xl1 = xlabel('$x_1$');
yl1 = ylabel('$x_2$');
%xticks(0:100:500)
%title('Loss of Optimality (%)')
set(xl1,'Interpreter','latex')
set(yl1,'Interpreter','latex')
set(gcf,'color','w')
l1 = legend('$e_{t+N} = 0$','$h(x,u) = 0$','MTPA line','$x_t$');
set(l1,'Interpreter','latex','location','southwest','Orientation','vertical');
grid on
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;


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
plot(time,computation_time*10^3,'b')
ylim([0 1])

figure(6)
plot(time,wr)

figure(7)
plot(time,lam_list)
hold off

figure(8)
plot(time,V_list)

figure(9)
plot(time,Log_list(:,1))

figure(10)
plot(time,Log_list(:,2:3))

mean(computation_time) * 10^3