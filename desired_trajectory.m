%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computing our desired trajectory%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Written by: Diego S. Dantonio

function dtj = desired_trajectory(t,params)
%     dtj.xCd =  [0.1*sin(t/5);... %xCd
%                          0.002*t;... %yCd
%                         0.1*cos(t/5) + 0.5]; %zCd    
%     
%     dtj.dxCd = [(0.1*cos(t/5))/5;... %dxCd
%                          0.002;... %dyCd    
%                          -(0.1*sin(t/5))/5]; %dzCd
%     
%     dtj.ddxCd = [-(0.1*sin(t/5))/(25);... %ddxCd
%                           0;... %ddzCd    
%                           -(0.1*cos(t/5))/(25)]; %ddzCd 

%     dtj.xCd =  [-0.2,0.2,t/200 + 0.5]'; %zCd    
%     
%     dtj.dxCd = [0,0,1/200]'; %dzCd
%     
%     dtj.ddxCd = [0,0,0]'; %ddzCd 

    t = t; 

    dtj.xCd =  [0,0,0.4]'; %zCd    
    dtj.dxCd = [0,0,0]'; %dzCd
    dtj.ddxCd = [0,0,0]'; %ddzCd 

    
    
    
     dtj.yaw = t/10;
%    dtj.yaw = 0;
    dtj.Omega = [0, 0, 1/10]';
    dtj.dOmega = [0, 0, 0]';
    
    dtj.R = Rot(dtj.yaw,'z');
    dtj.dR = dtj.R * hat(dtj.Omega);
    dtj.ddR = dtj.dR * hat(dtj.dOmega);
                  
    
    %% xbar desired
    dtj.x_bar_d = 0.7 + 0.3*cos(t/1);
    dtj.dx_bar_d = - 0.3*sin(t/1)/1;
    dtj.ddx_bar_d = - 0.3*cos(t/1)/(1*1);
 
    fun = @(c)  (- params.l/2 + c * sinh(dtj.x_bar_d/(2*c))); x0 = 0.15;
    dtj.c = fzero(fun,x0);    

    aux = dtj.x_bar_d/(2*dtj.c);
    
    fun = @(dc) (dc * sinh(aux)...
        + (((dtj.dx_bar_d)/2 - (dtj.x_bar_d * dc)/(2*dtj.c))...
        * cosh(aux)));  x0 = 0;
    dtj.dc = fzero(fun,x0);    
    

fun = @(ddc) ((( dtj.dx_bar_d)/(2) ...
    - (dtj.x_bar_d * dtj.dc)/(2 * dtj.c))^2 * sinh(aux) ...
    ...
    + 2*dtj.dc * (( dtj.dx_bar_d)/(2 * dtj.c) ...
    - (dtj.x_bar_d * dtj.dc)/(2 * dtj.c^2)) * cosh(aux) ... 
    ...
    + ((dtj.x_bar_d * ( dtj.dc)^2)/dtj.c^2 ...
    - ( dtj.dc * dtj.dx_bar_d)/dtj.c ...
    - (dtj.x_bar_d * ddc)/(2 * dtj.c) ...
    + ( dtj.ddx_bar_d)/(2)) * cosh(aux) ... 
    ...
    + ddc * sinh(aux)); x0 = 0;

    x0 = 0;

    dtj.ddc = fzero(fun,x0);  
 
    
    %% xA, xB and zAB desired
    % xA and xB desired
    % xA = xC_d - x_bar/2
    % xB = xC_d + x_bar/2


    % yA and yB desired  
    
    % zAB desired 
    % zABd = zCd + c*(cosh(x_bar_d/(2*c)) - 1);
    dtj.zABd = dtj.c*(cosh(dtj.x_bar_d/(2*dtj.c)) - 1);
    
    dtj.dzABd =  dtj.dc*(cosh(aux) - 1)...
        + (dtj.dx_bar_d/(2) - (dtj.x_bar_d*dtj.dc)/(2*dtj.c))*sinh(aux);
    
    

    dtj.ddzABd = ((( dtj.dx_bar_d)/(2) ...
                    - (dtj.x_bar_d * dtj.dc)/(2 * dtj.c))^2 * cosh(aux) ...
                    ...
                    + 2*dtj.dc * (( dtj.dx_bar_d)/(2 * dtj.c) ...
                    - (dtj.x_bar_d * dtj.dc)/(2 * dtj.c^2)) * sinh(aux) ... 
                    ...
                    + ((dtj.x_bar_d * ( dtj.dc)^2)/dtj.c^2 ...
                    - ( dtj.dc * dtj.dx_bar_d)/dtj.c ...
                    - (dtj.x_bar_d * dtj.ddc)/(2 * dtj.c) ...
                    + ( dtj.ddx_bar_d)/(2)) * sinh(aux) ... 
                    ...
                    + dtj.ddc * cosh(aux));
    
    xAdv = [- dtj.x_bar_d/2;...
            0;...
            dtj.zABd];
            
    xBdv = [dtj.x_bar_d/2;...
            0;...
            dtj.zABd];
        
    dxAdv = [- dtj.dx_bar_d/2;...
            0;...
            dtj.dzABd];
            
    dxBdv = [dtj.dx_bar_d/2;...
            0;...
            dtj.dzABd]; 
        
    ddxAdv = [- dtj.ddx_bar_d/2;...
            0;...
            dtj.ddzABd];
            
    ddxBdv = [dtj.ddx_bar_d/2;...
            0;...
            dtj.ddzABd]; 
        
    %% Desired positions to vectors 
    
    dtj.xAd = dtj.xCd + dtj.R * xAdv;                                                              
    dtj.xBd = dtj.xCd + dtj.R * xBdv;

    dtj.dxAd = dtj.dxCd + dtj.dR * xAdv + dtj.R * dxAdv;  
    dtj.dxBd = dtj.dxCd + dtj.dR * xBdv + dtj.R * dxBdv;

    dtj.ddxAd = dtj.ddxCd + dtj.ddR*xAdv + 2*dtj.dR*dxAdv + dtj.R*ddxAdv;  
    dtj.ddxBd = dtj.ddxCd + dtj.ddR*xBdv + 2*dtj.dR*dxBdv + dtj.R*ddxBdv;

    %% Desired Rotation
    vCA = dtj.xAd - dtj.xCd;
    vCB = dtj.xBd - dtj.xCd;
    vD = vec_cross(vCA,vCB);
    vBA1 = (vD/norm(vD));
    vBA3 = params.e3;
    vBA2 = vec_cross(vBA1,vBA3);
%    vBA2 = vBA2/norm(vBA2);

    dtj.vBA1 = vBA1;
    dtj.vBA2 = vBA2;
    dtj.vBA3 = vBA3;
    
    vBB1 = (vD/norm(vD));

    vBB3 = params.e3;
    vBB2 = vec_cross(vBB1,vBB3);
%    vBB2 = vBB2/norm(vBB2);

    dtj.vBB1 = vBB1;
    dtj.vBB2 = vBB2;
    dtj.vBB3 = vBB3;    
    

end