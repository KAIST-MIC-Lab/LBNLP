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

% Ns = 98;
% time = (Ts:Ts:Ts*(Ns))';

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
%wr = min(wr,950);
wr = 1089.1*ones(Ns,1);

% Control variables
m = 100;
computation_time = zeros(Ns,m);

%Te_ref_LPF(1) = 20; id(1) = -3; iq(1) = 28; u_opt_list(1,1) = -18; u_opt_list(1,2) = 47;
Te_ref_LPF(1) = 0; id(1) = 0; iq(1) = 0; u_opt_list(1,1) = 0; u_opt_list(1,2) = wr(1) * LAMpm;

% ALM design parameters

%mut = 1; muvt = Vmax^2*100; muc = Imax^2*0.1;
%mut = 0.1 * 10 * 0.001; muvt = mu; muc = mu;
mut = 1e-0; muvt = mut; muc = mut;

%alpha = 1e-2;

for n = 1:m

    lambdat = 0; lambdavt = 0; lambdac = 0;
    %mut = 0.1; muvt = Vmax^2; muc = Imax^2 * 10;

    u_opt_list = zeros(Ns,2);
    lam_list = zeros(Ns,2);

    for i = 2:Ns*0.8
        xp = [id(i - 1) + 0*2*(rand - 0.5)*Imax*0.015 ; iq(i - 1) + 0*2*(rand - 0.5)*Imax*0.015];
        up = u_opt_list(i - 1,:)';

        %     Ld_nom = 1.2*Ld;
        %     Lq_nom = 1.1*Lq;
        %     LAMpm_nom = 0.9*LAMpm;
        %     Rs_nom = Rs;

        %Ld_nom = Ld;
        %Lq_nom = Lq;
        %LAMpm_nom = LAMpm;
        %Rs_nom = Rs;

        %k1_nom = 1.5*PP*LAMpm_nom;
        %k2_nom = 1.5*PP*(Ld_nom - Lq_nom);

        %d = [0 ; k1_nom];
        %Q = [0 , k2_nom/2 ; k2_nom/2 , 0];

        Te_ref_LPF(i) = a_LPF*Te_ref_LPF(i - 1) + (1 - a_LPF)*Te_ref(i);
        y_ref = Te_ref_LPF(i);

        idq = xp;
        LAMdq = [Ld * idq(1) + LAMpm ; Lq * idq(2)];
        p = PP;
        vdqRef = up;

        Linv = [Ld , 0 ; 0 , Lq]\eye(2);
        J = [0 , 1 ; -1 , 0];
        A11 = eye(2) + Ts * wr(i) * J;
        A12 = -Ts * Rs * eye(2);
        A21 = Ts * wr(i) * Linv * J;
        A22 = eye(2) - Ts * Rs * Linv;
        B1 = Ts * eye(2);
        B2 = Ts * Linv;

        D2 = 1.5 * Rs * (B2' * B2);
        D1 = 1.5 * Rs * 2 * (A21*LAMdq + A22*idq)' * B2;
        E2 = -1.5 * p * B1' * J * B2;
        E1 = -1.5 * p * (A21*LAMdq + A22*idq)' * J' * B1 - 1.5 * p * (A11*LAMdq + A12*idq)' * J * B2;
        E0 = y_ref - 1.5 * p * (A11*LAMdq + A12*idq)' * J * (A21*LAMdq + A22*idq);
        Ftemp1 = -wr(i) * J * B1 + Rs * B2;
        Ftemp2 = -wr(i) * J * (A11*LAMdq + A12*idq) + Rs * (A21*LAMdq + A22*idq);
        F2 = -Ftemp1' * Ftemp1;
        F1 = -2 * Ftemp2' * Ftemp1;
        F0 = Vmax^2 - Ftemp2' * Ftemp2;
        G2 = -(B2' * B2);
        G1 = -2 * (A21*LAMdq + A22*idq)' * B2;
        G0 = Imax^2 - (A21*LAMdq + A22*idq)' * (A21*LAMdq + A22*idq);
        %H2 = -eye(2);
        %H1 = 0;
        %H0 = Vmax^2;
        ct = vdqRef'*E2*vdqRef + E1*vdqRef + E0;
        cvt = vdqRef'*F2*vdqRef + F1*vdqRef + F0;
        cc = vdqRef'*G2*vdqRef + G1*vdqRef + G0;

        tStart = tic;
        for iter = 1:1
            %temp1 = (E2 + E2') * vdqRef + E1';
            %temp2 =  (F2 + F2') * vdqRef + F1';
            %flagMTPV = temp1' * temp2 / (norm(temp1) * norm(temp2));
            %cond = ((((flagMTPV > 1 - alpha) || (flagMTPV < -1 + alpha)) && (lambdavt > 0)) || (lambdac > 0)) && ...
            %    ((TqRef > 0) && (TqRef > Te) || (TqRef <= 0) && (TqRef < Te));
            cond = 0;
            gradL = [10 ; 10]; cnt = 0;
            while ((norm(gradL) > 0.0001) && (cnt < 10))
                gradct = (E2 + E2') * vdqRef + E1';
                condcvt = cvt - lambdavt * muvt;
                condcc = cc - lambdac * muc;
                if cond
                    gradL = - lambdat * gradct;
                    grad2L = - lambdat * (E2 + E2');
                else
                    gradL = (D2 + D2') * vdqRef + D1' - lambdat * gradct + 1/mut * ct * gradct;
                    grad2L = (D2 + D2') - lambdat * (E2 + E2') + 1/mut * (ct * (E2 + E2') + gradct * gradct');
                end
                % grad2L
                % eig(grad2L)
                % grad2L(1,1) = grad2L(1,1) + 5;
                % grad2L
                % eig(grad2L)
                if condcvt <= 0
                    gradcvt = (F2 + F2') * vdqRef + F1';
                    gradL = gradL - lambdavt * gradcvt + 1/muvt * cvt * gradcvt;
                    grad2L = grad2L - lambdavt * (F2 + F2') + 1/muvt * (cvt * (F2 + F2') + gradcvt * gradcvt');
                end

                % if condcc <= 0
                %     gradcc = (G2 + G2') * vdqRef + G1';
                %     gradL = gradL - lambdac * gradcc + 1/muc * cc * gradcc;
                %     grad2L = grad2L - lambdac * (G2 + G2') + 1/muc * (cc * (G2 + G2') + gradcc * gradcc');
                % end
                vdqRef = vdqRef - grad2L \ gradL;

                ct = vdqRef'*E2*vdqRef + E1*vdqRef + E0;
                cvt = vdqRef'*F2*vdqRef + F1*vdqRef + F0;
                cc = vdqRef'*G2*vdqRef + G1*vdqRef + G0;
                cnt = cnt + 1;
            end
            %cnt
            if ~cond
                lambdat = lambdat - 1/mut * ct;
            end
            lambdavt = max(lambdavt - 1/muvt * cvt,0);
            lambdac = max(lambdac - 1/muc * cc,0);


            %diff = [-l/2*(grad_L + (W\Ca)*ca) ; -l*(ca)];
            %stepsize = 1;
            %sol0 = [up ; lam];
            % for k = 1:10
            %     sol = sol0 + stepsize*diff*Ts;
            %     up = sol(1:2);
            %     lam = sol(3:4);
            %
            %     xn = A*xp + B*up + C;
            %     ceq = (d'*xn + xn'*Q*xn - y_ref);
            %     uss = B\((eye(2) - A)*xn - C);
            %     %cin = Vmax^2 - uss'*uss;
            %     cin = (Vmax^2 - uss'*uss)*sc;
            %
            %     ca = [ceq ; min(cin,0)];
            %     %Ca = [B'*(d + 2*Q*xn) , -2*B'*(eye(2) - A')*(inv(B))'*uss];
            %     Ca = [B'*(d + 2*Q*xn) , -2*B'*(eye(2) - A')*(inv(B))'*uss*sc];
            %
            %     grad_L = 2*B'*(xn) - (lam(1)*Ca(:,1) + lam(2)*Ca(:,2));
            %
            %     V_temp = 0.5*(grad_L'*grad_L) + 0.5*(ca'*ca);
            %     if (V_temp <= V)
            %         break;
            %     end
            %     stepsize = stepsize*0.5;
            % end
            % k
        end

        % if norm(up) > Vmax
        %     u_opt = up * Vmax/norm(up);
        % else
        %     u_opt = up;
        % end
        u_opt = vdqRef;
        computation_time(i,n) = toc(tStart);
        u_opt_list(i,:) = u_opt;
        %lam_list(i,:) = lam;

        % PMSM modeling
        id(i) = id(i - 1) + Ts/Ld*(-Rs*id(i - 1) + wr(i)*Lq*iq(i - 1) + u_opt(1));
        iq(i) = iq(i - 1) + Ts/Lq*(-Rs*iq(i - 1) - wr(i)*(Ld*id(i - 1) + LAMpm) + u_opt(2));
        Te(i) = (k1 + k2*id(i))*iq(i);
    end
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

% idq_MTPA_m2 = [id , iq];
% save idq_MTPA_m2.mat idq_MTPA_m2
% idq_MTPA_m0 = [id , iq];
% save idq_MTPA_m0.mat idq_MTPA_m0
% idq_MTPA_p2 = [id , iq];
% save idq_MTPA_p2.mat idq_MTPA_p2

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

x = 1:length(computation_time);
CM_avg = mean(computation_time) * 10^6;
CM_std = std(computation_time) * 10^6;
CM_up = CM_avg + 1 * CM_std;
CM_lo = CM_avg - 1 * CM_std;
x_CM = [x fliplr(x)];
inBetween_CM = [CM_up fliplr(CM_lo)];

sc = 1;
figure(5)
%plot(1:sc:length(CM_avg),CM_avg(1:sc:end),'Color',color_blue,'LineStyle','-.','LineWidth',2)%,'Marker','square','MarkerIndices',1:sc:length(SOC10),'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',4)
plot(1:sc:length(CM_avg),CM_avg(1:sc:end),'Color',color_yellow,'LineWidth',2)
hold on
%fill(x_CM,inBetween_CM,'w','EdgeColor','none','FaceColor',color_blue,'FaceAlpha',0.4)
fill(x_CM,inBetween_CM,'w','EdgeColor','none','FaceColor',color_yellow,'FaceAlpha',0.4)
hold off
xlim([0 length(CM_avg)])
ylim([0 40])
%title('Vehicle operation')
xl1 = xlabel('$t$');
yl1 = ylabel('Computation time ($\mu$s)');
%xticks(0:100:500)
title('ALM')
set(xl1,'Interpreter','latex')
set(yl1,'Interpreter','latex')
set(gcf,'color','w')
%l1 = legend('CS','PS');
%set(l1,'Orientation','horizontal','Location','northeast')
grid on
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;

figure(6)
plot(time,wr)

% figure(7)
% plot(time,lam_list(:,1),'b')
% hold on
% plot(time,lam_list(:,2),'r')
% hold off


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

%mean(computation_time) * 10^3

computation_time_array_temp = computation_time(2:16,:);
computation_time_array = reshape(computation_time_array_temp,numel(computation_time_array_temp),1);
mean(computation_time_array*10^3)
std(computation_time_array*10^3)