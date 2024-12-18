color_yellow = [0.9290, 0.6940, 0.1250];
color_blue = [0, 0.4470, 0.7410];

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
Te_ref_LPF = zeros(Ns,1);
wc_LPF = 2*pi*500;
a_LPF = 1/(1 + wc_LPF*Ts);
Te_ref_lim = zeros(Ns,1);

wr = 1000*pi/30*PP + 2*1.5*0.1*5*min((2000*time - 0*4000*max(time - 3,0) + 0*2000*max(time - 4.8,0)),8000)*pi/30*PP;
%wr = min(wr,950);
wr = 1089.1*ones(Ns,1);

% Control variables
w = 1e-0; %1e-4, 1e-2 ,1e-0
computation_time = zeros(Ns,1);

%Te_ref_LPF(1) = 20; id(1) = -3; iq(1) = 28; u_opt_list(1,1) = -18; u_opt_list(1,2) = 47;
Te_ref_LPF(1) = 0; id(1) = 0; iq(1) = 0; u_opt_list(1,1) = 0; u_opt_list(1,2) = wr(1) * LAMpm;

for i = 2:Ns*0.8
    
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

% id_SQP3_MTPA = id;
% iq_SQP3_MTPA = iq;
% save id_SQP3_MTPA.mat id_SQP3_MTPA
% save iq_SQP3_MTPA.mat iq_SQP3_MTPA

% MTPA & MTPV trajactory
iq_MTPA = 0:1:2*Imax;
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
x1 = -50:30;
plot(x1,Te_ref_LPF(i - 1)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
hold on
fimplicit(@(x,y) ((Rs*x - wr(i)*Lq*y).^2 + (Rs*y + wr(i)*(Ld*x + LAMpm)).^2 - Vmax^2), axis,'Color',color_yellow,'Linestyle','-','LineWidth',1.5)
plot(id_MTPA,iq_MTPA,'k--')%,'LineWidth',1.5)
%fimplicit(@(x,y) (x.^2 + y.^2 - Imax^2), axis,'Color','b','Linestyle','-.','LineWidth',1.5)
plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
for ii = 1:1:i - 2
    %fimplicit(@(x,y) ((k1 + k2*x).*y - Te_ref_LPF(ii)), axis,'Color',[0 0 0 0.9],'Linestyle','-','LineWidth',2.5)
    plot(x1,Te_ref_LPF(ii)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
end
plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
p1 = [-5 15];                         % First Point
p2 = [-10 35];                         % Second Point
dp = p2-p1;                         % Difference
%quiver(p1(1),p1(2),dp(1),dp(2),0,'Color',color_blue)
%t1 = text(-30,25,'$t = 0,1,2,\cdots$','FontSize',16);
%set(t1,'Interpreter','latex')
hold off
xlim([-Imax 0.5*Imax])
ylim([-0.2*Imax Imax])
tit1 = title('SQP $(w = 10^{-4})$');
xl1 = xlabel('$x_1$');
yl1 = ylabel('$x_2$');
%xticks(0:100:500)
set(xl1,'Interpreter','latex')
set(yl1,'Interpreter','latex')
set(tit1,'Interpreter','latex')
set(gcf,'color','w')
l1 = legend('$e_{t+N} = 0$','$h(x,u) = 0$','MTPA line','$x_t$');
set(l1,'Interpreter','latex','location','northeast','Orientation','vertical');
grid on
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;

figure(22)
%plot(MTPV(1,:),MTPV(2,:),'k--')
%plot(id(i - 1),iq(i - 1),'ro')
%plot(id(i),iq(i),'go')
x1 = -50:30;
plot(x1,Te_ref_LPF(i - 1)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
hold on
fimplicit(@(x,y) ((Rs*x - wr(i)*Lq*y).^2 + (Rs*y + wr(i)*(Ld*x + LAMpm)).^2 - Vmax^2), axis,'Color',color_yellow,'Linestyle','-','LineWidth',1.5)
plot(id_MTPA,iq_MTPA,'k--')
%fimplicit(@(x,y) (x.^2 + y.^2 - Imax^2), axis,'Color','b','Linestyle','-.','LineWidth',1.5)
plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
for ii = 1:1:i - 2
    %fimplicit(@(x,y) ((k1 + k2*x).*y - Te_ref_LPF(ii)), axis,'Color',[0 0 0 0.9],'Linestyle','-','LineWidth',2.5)
    plot(x1,Te_ref_LPF(ii)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
end
plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
p1 = [5 15];                         % First Point
p2 = [-15 35];                         % Second Point
dp = p2-p1;                         % Difference
%quiver(p1(1),p1(2),dp(1),dp(2),0,'Color',color_blue)
%t1 = text(0,25,'$t = 0,1,2,\cdots$','FontSize',16);
%set(t1,'Interpreter','latex')
hold off
xlim([-Imax 0.5*Imax])
ylim([-0.2*Imax Imax])
tit1 = title('SQP $(w = 10^{-4})$');
xl1 = xlabel('$x_1$');
yl1 = ylabel('$x_2$');
%xticks(0:100:500)
set(xl1,'Interpreter','latex')
set(yl1,'Interpreter','latex')
set(tit1,'Interpreter','latex')
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
plot(time,computation_time,'b')
ylim([0 0.01])

figure(6)
plot(time,wr)


%%

Optimal_MTPA = flipud([
    [-6.6,42.6184]
    [-6.35,41.7427]
[-5.77000000000000,39.7510527336279]
[-5.04000000000000,37.1467649564978]
[-4.52000000000000,35.1909439783787]
[-3.91000000000000,32.6005146183656]
[-3.13000000000000,29.1762308035465]
[-2.23000000000000,24.6347138151198]
[-1.29000000000000,18.6034732990339]
[-0.420000000000002,10.5987256026389]
[0,0]]);

Optimal_FW = flipud([
[-28.13,39.5208]
[-27.9500000000000,39.3350565871342]
[-27.7200000000000,39.0897233546424]
[-27.4100000000000,38.7677515549208]
[-27.0200000000000,38.3413413530990]
[-26.5100000000000,37.7789302985536]
[-25.8300000000000,37.0380971432740]
[-24.9900000000000,36.0526980648474]
[-23.9100000000000,34.7462619405194]
[-22.5300000000000,33.0101667506180]
[-20.8300000000000,30.6913677713154]
[-18.7400000000000,27.5884169565877]
[-16.3400000000000,23.4126577356443]
[-13.7200000000000,17.7829315331949]
[-11.2300000000000,10.1885506233000]
[0,0]]);

time_seq = [1,2,3,4,5,6,7,8,10,13,16];

load('id_SQP1_MTPA.mat')
load('iq_SQP1_MTPA.mat')
load('id_SQP2_MTPA.mat')
load('iq_SQP2_MTPA.mat')
load('id_SQP3_MTPA.mat')
load('iq_SQP3_MTPA.mat')

axis = [-Imax Imax -Imax Imax];

figure(101)
%plot(MTPV(1,:),MTPV(2,:),'k--')
%plot(id(i - 1),iq(i - 1),'ro')
%plot(id(i),iq(i),'go')
x1 = -50:0.01:30;
plot(x1,Te_ref_LPF(i)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
plot(Optimal_MTPA(:,1),Optimal_MTPA(:,2),'ks','MarkerSize',15,'MarkerFaceColor','w','LineWidth',1)
hold on
plot(id_SQP1_MTPA(time_seq),iq_SQP1_MTPA(time_seq),'-o','Color',color_green,'MarkerFaceColor',color_green)
plot(id_SQP2_MTPA(time_seq),iq_SQP2_MTPA(time_seq),'-o','Color',color_blue,'MarkerFaceColor',color_blue)
plot(id_SQP3_MTPA(time_seq),iq_SQP3_MTPA(time_seq),'-o','Color',color_red,'MarkerFaceColor',color_red)
%fimplicit(@(x,y) ((Rs*x - wr(i)*Lq*y).^2 + (Rs*y + wr(i)*(Ld*x + LAMpm)).^2 - Vmax^2), axis,'Color',color_yellow,'Linestyle','-','LineWidth',1.5)
plot(x1,Te_ref_LPF(ii)./(k1 + k2.*x1),'Color','k','Linestyle','-','LineWidth',1.5)
% 전압원 그리기
theta = linspace(0, 2*pi, 1000);
xy_fill = [Rs , -wr(i)*Lq ; wr(i)*Ld , Rs]\(Vmax*[cos(theta) ; sin(theta)] + [0 ; -wr(i)*LAMpm]);
x_fill = xy_fill(1,:);
y_fill = xy_fill(2,:);
% fill 함수를 사용하여 타원 내부 색칠
fill(x_fill, y_fill,color_yellow,'FaceAlpha',0.2,'EdgeColor',color_yellow); % 'r'은 빨간색으로 채우기를 의미
plot(Optimal_MTPA(:,1),Optimal_MTPA(:,2),'ks','MarkerSize',15,'MarkerFaceColor','w','LineWidth',1)
plot(id_MTPA,iq_MTPA,'Color',color_purple,'LineWidth',2)
%fimplicit(@(x,y) (x.^2 + y.^2 - Imax^2), axis,'Color','b','Linestyle','-.','LineWidth',1.5)
%plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
% for ii = time_seq
%     %fimplicit(@(x,y) ((k1 + k2*x).*y - Te_ref_LPF(ii)), axis,'Color',[0 0 0 0.9],'Linestyle','-','LineWidth',2.5)
%     plot(x1,Te_ref_LPF(ii)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
% end
plot(id_SQP1_MTPA(time_seq),iq_SQP1_MTPA(time_seq),'-o','Color',color_green,'MarkerFaceColor',color_green)
plot(id_SQP2_MTPA(time_seq),iq_SQP2_MTPA(time_seq),'-o','Color',color_blue,'MarkerFaceColor',color_blue)
plot(id_SQP3_MTPA(time_seq),iq_SQP3_MTPA(time_seq),'-o','Color',color_red,'MarkerFaceColor',color_red)

t1 = text(-24.5,50,'$h({x},u) \ge 0$','FontSize',16,'Color',color_yellow);
set(t1,'Interpreter','latex')
%quiver(-11,47,10,0,0,'Color','k')
t2 = text(-8,50,'MTPC','FontSize',16,'Color',color_purple);
set(t2,'Interpreter','latex')
t3 = text(3.5,50,'$m_m = 30$ Nm','FontSize',16,'Color','k');
set(t3,'Interpreter','latex')
hold off
xlim([-0.5*Imax 0.3*Imax])
ylim([-0.2*Imax 1.2*Imax])
%tit1 = title('SQP $(w = 10^{-4})$');
xl1 = xlabel('$x_1(=i_s^d)$');
yl1 = ylabel('$x_2(=i_s^q)$');
%xticks(0:100:500)
set(xl1,'Interpreter','latex')
set(yl1,'Interpreter','latex')
%set(tit1,'Interpreter','latex')
set(gcf,'color','w')
l1 = legend('Optimal (MTPC)','SQP ($w=10^{-4}$)','SQP ($w=10^{-2}$)','SQP ($w=10^0$)');%,'$e[k+1] = 0$');
set(l1,'Interpreter','latex','location','southwest','Orientation','vertical','FontSize',14);
grid on
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;

load('id_SQP1_FW.mat')
load('iq_SQP1_FW.mat')
load('id_SQP2_FW.mat')
load('iq_SQP2_FW.mat')
load('id_SQP3_FW.mat')
load('iq_SQP3_FW.mat')

figure(102)
%plot(MTPV(1,:),MTPV(2,:),'k--')
%plot(id(i - 1),iq(i - 1),'ro')
%plot(id(i),iq(i),'go')
x1 = -50:0.01:30;
%plot(x1,Te_ref_LPF(i)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
plot(Optimal_FW(time_seq,1),Optimal_FW(time_seq,2),'ks','MarkerSize',15,'MarkerFaceColor','w','LineWidth',1)
hold on
plot(id_SQP1_FW(time_seq),iq_SQP1_FW(time_seq),'-o','Color',color_green,'MarkerFaceColor',color_green)
plot(id_SQP2_FW(time_seq),iq_SQP2_FW(time_seq),'-o','Color',color_blue,'MarkerFaceColor',color_blue)
plot(id_SQP3_FW(time_seq),iq_SQP3_FW(time_seq),'-o','Color',color_red,'MarkerFaceColor',color_red)
%fimplicit(@(x,y) ((Rs*x - wr(i)*Lq*y).^2 + (Rs*y + wr(i)*(Ld*x + LAMpm)).^2 - Vmax^2), axis,'Color',color_yellow,'Linestyle','-','LineWidth',1.5)
plot(x1,Te_ref_LPF(ii)./(k1 + k2.*x1),'Color','k','Linestyle','-','LineWidth',1.5)
% 전압원 그리기
theta = linspace(0, 2*pi, 1000);
xy_fill = [Rs , -wr(i)*Lq ; wr(i)*Ld , Rs]\(Vmax*[cos(theta) ; sin(theta)] + [0 ; -wr(i)*LAMpm]);
x_fill = xy_fill(1,:);
y_fill = xy_fill(2,:);
% fill 함수를 사용하여 타원 내부 색칠
fill(x_fill, y_fill,color_yellow,'FaceAlpha',0.2,'EdgeColor',color_yellow); % 'r'은 빨간색으로 채우기를 의미
plot(Optimal_FW(time_seq,1),Optimal_FW(time_seq,2),'ks','MarkerSize',15,'MarkerFaceColor','w','LineWidth',1)
plot(id_MTPA,iq_MTPA,'Color',color_purple,'LineWidth',2)
%fimplicit(@(x,y) (x.^2 + y.^2 - Imax^2), axis,'Color','b','Linestyle','-.','LineWidth',1.5)
%plot(id(1:1:i),iq(1:1:i),'-o','Color',color_blue)
% for ii = time_seq
%     %fimplicit(@(x,y) ((k1 + k2*x).*y - Te_ref_LPF(ii)), axis,'Color',[0 0 0 0.9],'Linestyle','-','LineWidth',2.5)
%     plot(x1,Te_ref_LPF(ii)./(k1 + k2.*x1),'Color',[0 0 0 0.2],'Linestyle','-','LineWidth',1.5)
% end
plot(id_SQP1_FW(time_seq),iq_SQP1_FW(time_seq),'-o','Color',color_green,'MarkerFaceColor',color_green)
plot(id_SQP2_FW(time_seq),iq_SQP2_FW(time_seq),'-o','Color',color_blue,'MarkerFaceColor',color_blue)
plot(id_SQP3_FW(time_seq),iq_SQP3_FW(time_seq),'-o','Color',color_red,'MarkerFaceColor',color_red)

t1 = text(-49.5,44.5,'$h({x},u) \ge 0$','FontSize',16,'Color',color_yellow);
set(t1,'Interpreter','latex')
%quiver(-11,47,10,0,0,'Color','k')
t2 = text(-16,46,'MTPC','FontSize',16,'Color',color_purple);
set(t2,'Interpreter','latex')
t3 = text(-49.5,33,'$m_m = 30$ Nm','FontSize',16,'Color','k');
set(t3,'Interpreter','latex')
hold off
xlim([-Imax 0.0*Imax])
ylim([-0.2*Imax Imax])
%tit1 = title('SQP $(w = 10^{-4})$');
xl1 = xlabel('$x_1(=i_s^d)$');
yl1 = ylabel('$x_2(=i_s^q)$');
%xticks(0:100:500)
set(xl1,'Interpreter','latex')
set(yl1,'Interpreter','latex')
%set(tit1,'Interpreter','latex')
set(gcf,'color','w')
l1 = legend('Optimal (FW)','SQP ($w=10^{-4}$)','SQP ($w=10^{-2}$)','SQP ($w=10^0$)');%,'$e[k+1] = 0$');
set(l1,'Interpreter','latex','location','southwest','Orientation','vertical','FontSize',14);
grid on
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;

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
