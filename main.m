%% Written by: Diego S. Dantonio

%% INITIALZING WORKSPACE

close all; clc; clear all;


addpath('./Geometry-Toolbox/');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALZING SYSTEM PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% moment of Inertia in Kg(m^2)
    params.mass = 0.3 ;
% moment of Inertia in Kg(m^2)
    params.J = diag([0.557, 0.557, 1.05]*10e-2);
% acceleration via gravity contant
    params.g = 9.81 ;
% interial fram axis
    params.e1 = [1;0;0] ;
    params.e2 = [0;1;0] ;
    params.e3 = [0;0;1] ;
% distance of center of mass from fram center in m
    params.d = 0.315;
% fixed constant in m
    params.c = 8.004*10e-4;
%defining parameters different for each trajectories
params.x_desired =  nan;
params.gen_traj = 1;        %change here
params.vel = nan;
params.acc = nan;
params.b1d = nan;
params.w_desired =  [0;0;0];
params.k1 = diag([15, 15 ,15]);
params.k2 = diag([5, 5, 10]);
params.kR = 200;
params.kOm = 1;

params.l = 2.0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTIALIZING - INTIAL PARAMETERS x,v,R,w %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intial position
    x_quad_A = [-0.08; 0.8; 1.];
    x_quad_B = [0.08; 0.8; 1.];
%     xQ0 = [1;3;2];
% Intial velocity
    v_quad_0 = zeros(3,1);
% Initial orientation
R0 = RPYtoRot_ZXY(0*pi/180, 0*pi/180, 180*pi/180);
%R0 = eye(3);
% Intial angular velocity
    w0= zeros(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenating the entire initial condition into a single vector
    qA = [x_quad_A; v_quad_0; reshape(R0,9,1); w0];
    qB = [x_quad_B; v_quad_0; reshape(R0,9,1); w0];
    x0 = [qA ; qB];
    x0 = [qA; qB];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% SIMULATION
odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% odeopts = [] ;
tend = 25;
tspan =  linspace(0,tend,tend/0.05);
[t, x] = ode15s(@odefun_quadDynamics, [0 25], x0, odeopts, params) ;

% Computing Various Quantities
disp('Evaluating...') ;
index = round(linspace(1, length(t), round(1*length(t))));
% ind = 0:length(t);
for i = index
   [~,xdA_,fA_,MA_,xdB_,fB_,MB_, xC_, xdC_] =  odefun_quadDynamics(t(i),x(i,:)',params);
   xdA(i,:) = xdA_';
   posA_err_fx(i) = norm(x(i,1:3)-xdA(i,1:3));
   velA_err_fx(i) = norm(x(i,4:6)-xdA(i,4:6));
   fA(i,1)= fA_;
   MA(i,:)= MA_';
   
   xdB(i,:) = xdB_';
   posB_err_fx(i) = norm(x(i,1:3)-xdB(i,1:3));
   velB_err_fx(i) = norm(x(i,4:6)-xdB(i,4:6));
   fB(i,1)= fB_;
   MB(i,:)= MB_';
   
   xC(i,:) = xC_';
   xdC(i,:) = xdC_';
end

%%% Plotting graphs
plott(t,x,xdA,posA_err_fx,velA_err_fx, xdB,posB_err_fx,velB_err_fx, xC, xdC,params);


            
%%   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% INTIALIZING - INTIAL PARAMETERS x,v,R,w %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Intial position
%     x_quad_0 = [0;0;0];
% %     xQ0 = [1;3;2];
% % Intial velocity
%     v_quad_0 = zeros(3,1);
% % Initial orientation
% %     R0 = RPYtoRot_ZXY(0*pi/180, 10*pi/180, 20*pi/180);
%     R0 = RPYtoRot_ZXY(0*pi/180, 0*pi/180, 0*pi/180);
% % Intial angular velocity
%     w0= zeros(3,1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Concatenating the entire initial condition into a single vector
%     x0 = [x_quad_0; v_quad_0; reshape(R0,9,1); w0];
    
    
    
% %%%%%%%% SIMULATION
% odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% % odeopts = [] ;
% [t, x] = ode15s(@odefun_quadDynamics, [0 20], x0, odeopts, quad_params) ;
% 
% % Computing Various Quantities
% disp('Evaluating...') ;
% index = round(linspace(1, length(t), round(1*length(t))));
% % ind = 0:length(t);
% for i = index
%    [~,xd_,f_,M_] =  odefun_quadDynamics(t(i),x(i,:)',quad_params);
%    xd(i,:) = xd_';
%    pos_err_fx(i) = norm(x(i,1:3)-xd(i,1:3));
%    vel_err_fx(i) = norm(x(i,4:6)-xd(i,4:6));
%    f(i,1)= f_;
%    M(i,:)= M_';
% end
% 
% %%% Plotting graphs
% plott(t,x,xd,pos_err_fx,vel_err_fx);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Function definitions %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Catenary robot simulation.
function[dx, xdA, fA,MA, xdB, fB,MB, xC, xdC] = odefun_quadDynamics(t,x,params)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fetching desired states %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xA = x(1:18);
    xB = x(19:end);

    %[traj_desired] = return_traj(t);
    dtj = desired_trajectory(t,params);
    
    x_bar = norm(xA(1:3) - xB(1:3));
    
    fun = @(c)  (- params.l/2 + c * sinh(x_bar/(2*c))); x0 = 0.15;
    c = abs(fzero(fun,x0));
    
    
    xC = xA(1:3) - Rot(dtj.yaw,'z')*[ - x_bar/2;...
                                     0;...
                                     c*(cosh(x_bar/(2*c))-1)];
    xdC = dtj.xCd;
    
    traj_desiredA = [dtj.xAd; dtj.dxAd; dtj.ddxAd; dtj.vBA1];
    traj_desiredB = [dtj.xBd; dtj.dxBd; dtj.ddxBd; dtj.vBB1];
    
    
    [dxA, xdA, fA,MA] = quadrotor(xA, traj_desiredA,params);
    [dxB, xdB, fB,MB] = quadrotor(xB, traj_desiredB,params);
    
    dx = [dxA;
          dxB];
      
end

%% Local simulation
function [dx, xd, f,M] = quadrotor(x, traj_desired,params)

    mass = params.mass;
    J = params.J;
    g = params.g;
    e1 = params.e1;
    e2 = params.e2;
    e3 = params.e3;

    %x_des = traj_desired.pos;
    x_des = traj_desired(1:3);
    v_des = traj_desired(4:6);
    a_des = traj_desired(7:9);
    b1d = traj_desired(10:12);
    w_desired = params.w_desired;
    w_desired_dot = zeros(3,1);


    x_curr = x(1:3);
    v_curr = x(4:6);
    R = reshape(x(7:15),3,3);
    w = x(16:18);
    b3 = R(:,3);
    dx = [];
    %%%%%%%%%%%
    % CONTROL %
    %%%%%%%%%%%
    % Position Control
    err_x = x_curr - x_des;
    err_v = v_curr - v_des;
    % constants for force and moment
    k1 = params.k1;
    k2 = params.k2;
    A = (-k1*err_x - k2*err_v + mass*a_des + mass*g*e3);
    normA = norm(A);
    b3_desired = A/normA;

    f = vec_dot(A,b3);

    %%%%%%%%%%%%%%%%%%%%
    % Attitude Control %
    %%%%%%%%%%%%%%%%%%%%
    %b1d definded above;
    b2d = vec_cross(b3_desired,b1d);
    norm_b2d = norm(b2d);
    b2d = b2d/norm_b2d;
    projection_b1d = -vec_cross(b3_desired,b2d);
    Rd = [vec_cross(b2d,b3_desired) b2d b3_desired];    
    %angluar velocity and rotation matrix constants
    kR = params.kR;
    kOm = params.kOm;
    %calculating error in angular velocity and attitude
%         psi_R = ;
    err_R = 1/2 * vee_map(Rd'*R - R'*Rd) ;
    err_Om = w - R'*Rd*w_desired ;
    M = -kR*err_R - kOm*err_Om + vec_cross(w, J*w)...
        - J*(hat_map(w)*R'*Rd*w_desired - R'*Rd*w_desired_dot) ;

    % Equations of Motion
    % -------------------
    xQ_dot = v_curr;
    vQ_dot = -g*e3 + (f/mass)*R*e3;


    R_dot = R*hat_map(w) ;
    Omega_dot = J\( -vec_cross(w, J*w) + M ) ;
    % Computing xd
    xd = [x_des; v_des; reshape(Rd, 9,1); w_desired];

    % Computing dx
    %-------------
    dx = [xQ_dot;
          vQ_dot;
          reshape(R_dot, 9,1);
          Omega_dot;];
end

%% plot function
function plott(t,x,xdA,posA_err_fx,velA_err_fx, xdB,posB_err_fx,velB_err_fx, xC, xdC,params)
    disp('Plotting graphs...');
    index = round(linspace(1, length(t), round(1*length(t))));
    figure(1);
    subplot(3,1,1);
        plot(t(index),x(index,1),'Color',[0.3010 0.7450 0.9330],'LineWidth',2); hold on;
%         plot(t(index),xdA(index,1),':r','LineWidth',2); 
        plot(t(index),x(index,19),'Color',[0.9290 0.6940 0.1250],'LineWidth',2); 
%         plot(t(index),xdB(index,1),':r','LineWidth',2); 
        plot(t(index),xC(index,1),'b','LineWidth',2); 
        plot(t(index),xdC(index,1),':r','LineWidth',2); 
              
        
        hold off;
        %axis equal;
        grid on;
        legend('$x_A$','$x_B$','$x_C$','$x_C$ desired','Interpreter','latex');%axis equal;
        ylabel('$x(m)$','Interpreter','latex');set(gca,'xticklabel',[])
        
    subplot(3,1,2);
        plot(t(index),x(index,2),'Color',[0.3010 0.7450 0.9330],'LineWidth',2); hold on;
%         plot(t(index),xdA(index,2),':r','LineWidth',2); 
        plot(t(index),x(index,20),'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
%         plot(t(index),xdB(index,2),':r','LineWidth',2); 
        plot(t(index),xC(index,2),'b','LineWidth',2); 
        plot(t(index),xdC(index,2),':r','LineWidth',2);hold off;
        %axis equal;
        grid on;
        legend('$x_A$','$x_B$','$x_C$','$x_C$ desired','Interpreter','latex');
        ylabel('$y(m)$','Interpreter','latex');set(gca,'xticklabel',[])
        
    subplot(3,1,3);
        plot(t(index),x(index,3),'Color',[0.3010 0.7450 0.9330],'LineWidth',2); hold on;
        plot(t(index),x(index,21),'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
        plot(t(index),xC(index,3),'b','LineWidth',2); 
        plot(t(index),xdC(index,3),':r','LineWidth',2);hold off;
        %axis equal;
        
        grid on;
        legend('$x_A$','$x_B$','$x_C$','$x_C$ desired','Interpreter','latex');
        xlabel('$t(s)$','Interpreter','latex');ylabel('$z(m)$','Interpreter','latex');
        
    figure(2)
        plot3(x(index,1),x(index,2),x(index,3),'-b','LineWidth',2); hold on;
        plot3(xdA(index,1),xdA(index,2),xdA(index,3),'-.r','LineWidth',3.5);
        
        plot3(x(index,19),x(index,20),x(index,21),'-k','LineWidth',2);
        plot3(xdB(index,1),xdB(index,2),xdB(index,3),':r','LineWidth',3.5); 
        
        plot3(xC(index,1),xC(index,2),xC(index,3),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
        plot3(xdC(index,1),xdC(index,2),xdC(index,3),'or','LineWidth',5);
        plot3(x(end,1),x(end,2),x(end,3),'ok','LineWidth',4);
        plot3(x(end,19),x(end,20),x(end,21),'ok','LineWidth',4);
        
        plot3(x(1,1),x(1,2),x(1,3),'or','LineWidth',2);
        plot3(x(1,19),x(1,20),x(1,21),'or','LineWidth',2);
        
        XA = x(1,1:3)';
        XB = x(1,19:21)';
        XC = xC(1,1:3)';
        vCA = XA - XC;
        vCB = XB - XC;
        vD = - vec_cross(vCA,vCB);
        vBA1 = (vD/norm(vD));
        vBA3 = params.e3;
        vBA2 = vec_cross(vBA1,vBA3);
        vBA2 = vBA2/norm(vBA2);
        
        vBA1 = XA + vBA1/8;
        vBA2 = XA + vBA2/8;
        vBA3 = XA + vBA3/8;
        
        
        vBB1 = (vD/norm(vD));
        
        vBB3 = params.e3;
        vBB2 = vec_cross(vBB1,vBB3);
        vBB2 = vBB2/norm(vBB2);
        
        vBB1 = XB + vBB1/8;
        vBB2 = XB + vBB2/8;
        vBB3 = XB + vBB3/8;
        
        
        arrow3(XA',vBA1','r'); 
        arrow3(XA',vBA2','g');
        arrow3(XA',vBA3','b');
        
        arrow3(XB',vBB1','r'); 
        arrow3(XB',vBB2','g');
        arrow3(XB',vBB3','b');

        XA = x(end,1:3)';
        XB = x(end,19:21)';
        XC = xC(end,1:3)';
        vCA = XA - XC;
        vCB = XB - XC;
        vD = -vec_cross(vCA,vCB);
        vBA1 = (vD/norm(vD));
        vBA3 = params.e3;
        vBA2 = vec_cross(vBA1,vBA3);
        vBA2 = vBA2/norm(vBA2);
        
        vBA1 = XA + vBA1/8;
        vBA2 = XA + vBA2/8;
        vBA3 = XA + vBA3/8;
        
        
        vBB1 = (vD/norm(vD));
        
        vBB3 = params.e3;
        vBB2 = vec_cross(vBB1,vBB3);
        vBB2 = vBB2/norm(vBB2);
        
        vBB1 = XB + vBB1/8;
        vBB2 = XB + vBB2/8;
        vBB3 = XB + vBB3/8;
        
        
        arrow3(XA',vBA1','r'); 
        arrow3(XA',vBA2','g');
        arrow3(XA',vBA3','b');
        
        arrow3(XB',vBB1','r'); 
        arrow3(XB',vBB2','g');
        arrow3(XB',vBB3','b');
        

        
        xcat = catenarycurve(t,x,xC,1,params);        
        plot3(xcat(:,1),xcat(:,2),xcat(:,3),'-k','LineWidth',1);  
        j = length(x);
        xcat = catenarycurve(t,x,xC,j,params);        
        plot3(xcat(:,1),xcat(:,2),xcat(:,3),'-k','LineWidth',1);                             

        axis equal;
        grid on; 
        legend('$x_A$','$x_A$ desired','$x_B$','$x_B$ desired','$x_C$','$x_C$ desired','Interpreter','latex');
        xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');
    view(2300,30)
%     %%%Plotting error functions
% %     figure;
% %     subplot(2,1,1);
% %         plot(t(index),posA_err_fx(index),'-b','LineWidth',2);
% %         grid on; title('error in position');
% %         legend('{e}_x');
% %         xlabel('Time');ylabel('{x}-{x}_d');
% %     subplot(2,1,2);
% %         plot(t(index),velA_err_fx(index),'-b','LineWidth',2);
% %         grid on; title('error in velocity');legend('{e}_v');
% %         xlabel('Time');ylabel('{v}-{v}_d');
%     
% %     %%Plotting force and Moments
% %     figure;
% %     subplot(2,2,1);
% %     plot(t(index),f(index,1),'-b','LineWidth',2); %hold on;
% % %     plot(t(index),xd(index,1),':r','LineWidth',2); hold off;
% %     grid on; title('f');%legend('x','x_d');%axis equal;
% %     xlabel('time');ylabel('x [m]');
% %     subplot(2,2,2);
% %     plot(t(index),M(index,1),'-b','LineWidth',2); %hold on;
% % %     plot(t(index),xd(index,2),':m','LineWidth',2); hold off;
% %     grid on; title('M_x');%legend('y','y_d');%axis equal;
% %     xlabel('time');ylabel('M in x direction');
% %     subplot(2,2,3);
% %     plot(t(index),M(index,2),'-b','LineWidth',2); %hold on;
% % %     plot(t(index),xd(index,3),':m','LineWidth',2); hold off;
% %     grid on; title('z');%legend('z','z_d');%axis equal;
% %     xlabel('time');ylabel('M in y direction');
% %     subplot(2,2,4);
% %     plot(t(index),M(index,3),'-b','LineWidth',2); %hold on;
% % %     plot(t(index),xd(index,3),':m','LineWidth',2); hold off;
% %     grid on; title('M_z');%legend('z','z_d');%axis equal;
% %     xlabel('time');ylabel('M in z direction');


end


%% Numerical solution 
function xcat = catenarycurve(t,x,xC,j,params)
        xbar = norm(x(j,1:3)-x(j,19:21));
        fun = @(c)  (- params.l/2 + c * sinh(xbar/(2*c))); x0 = 0.1;
        c = abs(fzero(fun,x0));
        
        dtj = desired_trajectory(t(j),params);

        xc = linspace(-xbar/2,xbar/2,100);
        xcat=[];
        for i=1:1:100
            xcat(i,:) = xC(j,1:3)' + Rot(dtj.yaw,'z')*[ xc(i);...
                                         0;...
                                         c*(cosh(xc(i)/(c))-1)];
        end
end
