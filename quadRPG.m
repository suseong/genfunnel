classdef quadRPG < DrakeSystem
% Defines the dynamics for the powered plane.
  
  properties
    m = 1.0; % kg
    g = -9.81; % m/s^2
  end
  
  methods
    function obj = quadRPG()
      obj = obj@DrakeSystem(9,0,4,9,0,1);
      obj = setStateFrame(obj,CoordinateFrame('quadState',9,'x',{'x','y','z','dx','dy','dz','phi','theta','psi'}));
      obj = setInputFrame(obj,CoordinateFrame('quadInput',4,'u',{'T','p','q','r'}));
      obj = setOutputFrame(obj,getStateFrame(obj));  % allow full state feedback
    end
    
    function [xdot, df, d2f, d3f] = dynamics(obj,t,x,u)
%     function xdot = dynamics(obj,t,x,u)
        
        phi = x(7); theta = x(8); psi = x(9);
        Rr = [1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
        Rp = [cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
        Ry = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];
        R = Ry*Rp*Rr;   
        invQ = [ 1 (sin(theta)*sin(phi))/cos(theta) (cos(phi)*sin(theta))/cos(theta);
                 0                 cos(phi)                -sin(phi);
                 0          sin(phi)/cos(theta)          cos(phi)/cos(theta)];
        
        ddx = obj.g*[0;0;1] + u(1)*R*[0;0;1]/obj.m;
        dw = invQ*[u(2);u(3);u(4)];
        
        xdot = [x(4);x(5);x(6);ddx(1);ddx(2);ddx(3);dw(1);dw(2);dw(3)];
        
        if (nargout>1)
            [df,d2f,d3f]= dynamicsGradients(obj,t,x,u,nargout-1);
        end
    end
    
    function y = output(obj,t,x,u)
      y = x;
    end
    
    function x = getInitialState(obj)
      x = [0 0 0 0 0 0 0 0 0]';
    end
        
  end
  
end
