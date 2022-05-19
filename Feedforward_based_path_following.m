%% ME593 Feedward Control Project %%
% Bob Hu, Toby Ouyang, May 2022
% Mechanical Engineering Departmen,University of Washington - Seattle
%
% Abstract: 
% Implemente inversion-based feedforward output tracking method to a
% path-following problem from published paper. Compare and evalue results.
%
% Reference: 
% A. P. Aguiar, J. P. Hespanha and P. V. Kokotovic, "Path-following for 
% nonminimum phase systems removes performance limitations," in IEEE 
% Transactions on Automatic Control, vol. 50, no. 2, pp. 234-239, 
% Feb. 2005, doi: 10.1109/TAC.2004.841924.

clc
close all 
clear all
format short

%% state-space representation
% state:[y1 y2 y1d y2d z1 z2 z1d z2d]'
A = [ 0 0 1 0 0 0 0 0;
      0 0 0 1 0 0 0 0;
      0 0 15 0 0 0 15 0;
      0 0 0 10 0 0 0 10;
      0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 1;
      -15 0 0 0 15 0 -15 0;
      0 -1 0 0 0 10 0 -10];
B = [0 0; 0 0; 1 0; 0 1; 0 0; 0 0; 0 0; 0 0];
C1  = [1 0 0 0 0 0 0 0]; 
C2 = [0 1 0 0 0 0 0 0];
C = [C1;C2];
D = [0];
sys = ss(A,B,C,D); 
% C*A*B %non zero relative degree = [2;2]

%% generate reference signal (circle with radius of 1 centered at origin)
radius = 1;
T = 5;
Tmax = 50;
dt  = 0.01;
time = 0:dt:Tmax-dt;
t = time;
T_pre = 10;
T_pos = 10;

for ix=1:length(time)
    if ix < T_pre/dt
    yd1(ix) = (1/(T_pre)^2)*(ix*dt)^2;
    yd2(ix) = 0;
    elseif ix < (Tmax-dt-T_pos)/dt
    wt = 2*pi*time(ix-T_pre/dt+1)/T;
    yd1(ix) = radius*cos(wt);
    yd2(ix) = radius*sin(wt);
    yd1final = yd1(ix);
    yd2final = yd2(ix);
    else
    yd1(ix) = yd1final;
    yd2(ix) = yd2final;
    end
end

wf = 100*2*pi/T; % cutoff frequency (10x original) 
sys_filter = ss(-wf, wf, -wf, wf);  % low-pass filter
vd1 = lsim(sys_filter, yd1, time);
vd1 = vd1';
ad1 = lsim(sys_filter, vd1, time);
ad1 = ad1';
vd2 = lsim(sys_filter, yd2, time);
vd2 = vd2';
ad2 = lsim(sys_filter, vd2, time);
ad2 = ad2';

figure
plot(yd1, yd2,'LineStyle',':','LineWidth',1.5,'Color','g');
set(gca,'XTick',(-1:0.2:1))
set(gca,'YTick',(-1:0.2:1))
title('Desired Trajectory/Reference')
axis equal

figure
subplot(311), plot(t,ad1,'-',t,ad2,'--');
xlabel('time (s)')
ylabel('acc')
grid ; set(gca,'FontSize',10)
subplot(312), plot(t,vd1,'-',t,vd2,'--');
xlabel('time (s)')
ylabel('vel')
grid; set(gca,'FontSize',10)
subplot(313), plot(t,yd1,'-',t,yd2,'--');
xlabel('time (s)')
ylabel('pos')
grid; set(gca,'FontSize',10)
% return

%% Case 1: cheap optimal control with LQR
Q = eye(10)
R = 0.001;
[Kx,S,E] = lqi(sys,Q,R);
Acl=A-B*Kx(:,1:8);
[V_Acl,Pn_Acl] = eig(Acl); % closed-loop system is stable as all eig are negative
sys_cl = ss(Acl,B,C,D);
figure
plot(yd1, yd2,'LineStyle',':','LineWidth',1.5,'Color','g')
hold on

DC=dcgain(sys_cl)
 rr1 =yd1/DC(1,1);
 rr2 =yd2/DC(2,2);
[ycheap,tn_cl] = lsim(sys_cl,[rr1;rr2]',t,[0 0 0 0 0 0 0 0 ]');
ycheap = ycheap';
plot(ycheap(1,:),ycheap(2,:),'LineStyle','-','Color','r','LineWidth',1.2)
set(gca,'XTick',(-1:0.2:1))
set(gca,'YTick',(-1:0.2:1))
axis equal
legend('y_d','y_{actual}')
title('Desired Trajectory VS Actual Trajectory using fb ')
%return
%% Case 2: inverse of the original linear MIMO nonminimum phase system
A = Acl;
sys = ss(A,B,C,D); 
% Check relative degree
check_CB_CAB = [C1*B ;C1*A*B;C2*B ;C2*A*B ];
Tt = [C1; C1*A; C2; C2*A];
Tb = (null(Tt))';
T = [Tt; Tb];
rank_of_T = rank(T);
invT = inv(T);
TinvL = invT(:,1:4);
TinvR = invT(:,5:8);
Ay = [C1*A*A;C2*A*A]; % relative degree is 2 2 
By = [C1*A*B;C2*A*B]; % relative degree is 2 2 
YY = [yd1; vd1; yd2; vd2; ad1; ad2];

Ainv = Tb*(A - B*(inv(By))*Ay)*TinvR;
Binv = [(Tb*(A - B*(inv(By))*Ay)*TinvL) Tb*B*(inv(By))];
Cinv = -(inv(By))*Ay*TinvR;
Dinv = [-((inv(By))*Ay*TinvL)  inv(By)];
sys_inv = ss(Ainv,Binv,Cinv,Dinv);
%return

%% split the internal dynamics into stable and unstable subsystem
[Vn_inv,Pn_inv] = eig(Ainv)
Tsplit = [Vn_inv(:,2) Vn_inv(:,4) Vn_inv(:,1) Vn_inv(:,3)]
invTsplit = inv(Tsplit); 
check_identity = Tsplit * invTsplit;
%return
% stable and unstable system
A_split = invTsplit*Ainv*Tsplit;
B_split = invTsplit*Binv;

A_s = A_split(1:2,1:2);  B_s = B_split(1:2,:);
A_u = A_split(3:4,3:4);  B_u = B_split(3:4,:);
Stable_part = A_s;
Unstable_part = A_u;
%change
YYd = YY;
YYd_ini = [0;0;0;0;0;0]; 
YYd_fin = [yd1(length(time)); 
           vd1(length(time)); 
           yd2(length(time)); 
           vd2(length(time));
           ad1(length(time));
           ad2(length(time))];
eta_s_initial = -inv(A_s)*B_s*YYd_ini;
eta_u_final   = -inv(A_u)*B_u*YYd_fin;
%return

%% Solve the stable internal dynamics 
sys_s =ss(A_s,B_s,eye(2,2),[0]);
[eta_s, tn, x_s] = lsim(sys_s, YYd, time, eta_s_initial);

figure
plot(time, eta_s)
xlabel('time(s)')
ylabel('\eta_s')
title('Stable Internal States')
Time_for_IC_to_be_zero = 4/abs(min(eig(A_s)));
%return

%% Solve the unstable internal dynamics
sys_u_backward =ss(-A_u,-B_u,eye(2,2), [0]);
YYdbackward = fliplr(YYd); 
[eta_u_backward,tm,x_u] = lsim(sys_u_backward,YYdbackward,time,eta_u_final); 

% flip the unstable solution to the correct time direction
eta_u = fliplr(eta_u_backward');

figure
plot(time, eta_u)
xlabel('time (s)')
ylabel('\eta_u')
title('Unstable Internal States')
time_for_error_to_become_zero = 4/abs(min(eig(A_u)));

% find the original internal dynamics eta_ref 
eta_d = Tsplit*[eta_s'; eta_u];

figure
plot(time, eta_d)
xlabel('time (s)')
ylabel('\eta_d')
title('Internal States')

xi_d = [yd1;vd1;yd2;vd2];
xd = invT*[xi_d; eta_d];
uff = Cinv*eta_d + Dinv*YYd;
figure
plot(time,uff)
xlabel('time (s)')
ylabel('u_{ff}')
title('Feedforward Input')
pole(sys_inv) % check the poles of the inverse system
tzero(sys)  % check the zeros of the original closed-loop system
%return

%% apply inverse input to the dynamic system
[yff] = lsim(sys_cl,uff,time);
yff=yff';
figure
plot(time, yd1)
hold on
plot(time, yff(1,:),'--')
xlabel('time (s)')
ylabel('y')
title('Output Tracking with Feedforward Only')
legend('y_{desired}','y_{actual}')
hold off

% use LQR method to find Kx
% Try changing R and Q
R_ff = 1; Q_ff = C'*C; [Kx_ff,S_ff,E_ff] = lqr(sys,Q_ff,R_ff) ;
Kx_ff
Acl_ff = A -B*Kx_ff; 
%Ut = uff +Kx_ff*xd; % feedback + feedforward 
sysCLff = ss(Acl_ff,B,C,D);
xo = [0;0;0;0;0;0;0;0];

[yffcl,time,xcl] = lsim(sys_cl,uff,t,xo);
yffcl = yffcl';

% the actual input to the system 
Ufb = -Kx_ff*(xcl'-xd);
Utotal = uff +Ufb; 

% plot the outputs and the inputs 
figure
plot(t,yffcl(1,:),t,yffcl(2,:),t,yd1,t,yd2)
legend('y_{cl1}','y_{cl2}','y_{d1}','y_{d2}')
xlabel('time')
ylabel('output, ff+fb')
set(gca,'FontSize',10)
figure
plot(yd1, yd2,'LineStyle',':','LineWidth',1.5,'Color','g')
hold on
plot(yffcl(1,:),yffcl(2,:),'LineStyle','-','Color','r','LineWidth',1.2)
set(gca,'XTick',(-1:0.2:1))
set(gca,'YTick',(-1:0.2:1))
axis equal
legend('y_d','y_{actual}')
title('Desired Trajectory VS Actual Trajectory')
%return

figure 
curve1 = animatedline('Color','g','LineStyle',':','LineWidth',1.5)
curve2 = animatedline('Color','r','LineWidth',1.3)
axis equal
legend('y_d','y_{actual}')
title('Desired Trajectory VS Actual Trajectory fb+ff')
for i =  1:length(yd1)
    addpoints(curve1,yd1(i),yd2(i))
    addpoints(curve2,yffcl(1,i),yffcl(2,i))
    drawnow
    grid on

end

figure
plot(t,Ufb(1,:),'m:',t,Ufb(2,:),'r:',t,uff(1,:),'g',t,uff(2,:),'b','linewidth',1.5)
title('Input')
legend('U_{fb1}','U_{fb2}','U_{ff1}','U_{ff2}','location','best')
xlabel('Time(sec)')
ylabel('Amplitude')
title('Amount of U_{ff} and U_{fb}')