clear;

t_data=0 :0.001 :20;
u1=14*sin(7*t_data);
u = @(t) interp1(t_data, u1, t, 'linear', 'extrap');
[t,x]=ode45(@(t,x)odefun(t,x,u),[0 :0.001 :20],[0;0]);

t=t_data;

y2 = @(t) interp1(t_data, x(:,2), t, 'linear', 'extrap');
y1 = @(t) interp1(t_data, x(:,1), t, 'linear', 'extrap');

sys1=tf([0 1],[1 1]);

[phi11,~]=lsim(sys1,x(:,1),t);
[phi22,~]=lsim(sys1,x(:,2),t);
[phi33,~]=lsim(sys1,u1,t);

phi1 = @(t) interp1(t_data, phi11, t, 'linear', 'extrap');
phi2 = @(t) interp1(t_data, phi22, t, 'linear', 'extrap');
phi3 = @(t) interp1(t_data, phi33, t, 'linear', 'extrap');

gamma=1*[0.4,0.4,0.4,0.4,0.4,0.4];   %0.6*[2,0.2,0.4,0.5,0.24,0.3]; %0.6*[1,0.2,0.4,0.5,0.3,0.2];
sigma=0.02;
[t,Z]=ode45(@(t,Z)odefun3(t,Z,y1,y2,phi1,phi2,phi3,u,gamma,sigma),[0 :0.001 :20],[-1.5;0.6;-0.9;-1.2;1;10]);

theta1=Z(:,1)-1;
theta2=Z(:,2);
theta3=Z(:,3);
theta4=Z(:,4);
theta5=Z(:,5);
theta6=Z(:,6)-1;

a11=-2.15*ones(length(t),1);
a12=0.25*ones(length(t),1);
a21=-0.75*ones(length(t),1);
a22=-2*ones(length(t),1);
b1=0*ones(length(t),1);
b2=1.5*ones(length(t),1);

figure();
plot(t,theta1,t,theta2,t,theta3,t,theta4,t,theta5,t,theta6);
hold on;
plot(t,a11,'Color','b','LineStyle','--');
hold on;
plot(t,a12,'Color','r','LineStyle','--');
hold on;
plot(t,a21,'Color','y','LineStyle','--');
hold on;
plot(t,a22,'Color','magenta','LineStyle','--');
hold on;
plot(t,b1,'Color','green','LineStyle','--');
hold on;
plot(t,b2,'Color','cyan','LineStyle','--');
title('Εκτιμήσεις παραμέτρων');
legend({'$\hat{a_{11}}$','$\hat{a_{12}}$','$\hat{a_{21}}$','$\hat{a_{22}}$','$\hat{b_{1}}$','$\hat{b_{2}}$'},'Interpreter','latex');
grid on;

x_hat1=phi11.*theta1+phi22.*theta2+phi33.*theta5;
x_hat2=phi11.*theta3+phi22.*theta4+phi33.*theta6;

 figure();
 plot(t,x(:,1),t,x_hat1);
 legend({'x1(t)','$\hat{x1}(t)$'},'Interpreter','latex');
 grid on;
 
 error=x(:,1)-x_hat1;
 
 figure();
 plot(t,error);
 title('e1=x1(t)-$\hat{x1}$(t)','Interpreter','latex');
 grid on;
 
 figure();
 plot(t,x(:,2),t,x_hat2);
 legend({'x2(t)','$\hat{x2}(t)$'},'Interpreter','latex');
 grid on;
 
 error2=x(:,2)-x_hat2;
 
 figure();
 plot(t,error);
 title('e2=x2(t)-$\hat{x2}$(t)','Interpreter','latex');
 grid on;
 
%2,10,100,0.01

function dx=odefun(t,x,u)
    % 
    x1=quantize_value(x(1),0.001);
    x2=quantize_value(x(2),0.001);
    dx=[-2.15*x1+0.25*x2;-0.75*x1-2*x2+1.5*u(t)];

    function x1_quant = quantize_value(x1, Delta)
    % quantize_value: κβαντίζει την τιμή x1 με βήμα Delta
    % x1: τιμή προς κβάντιση
    % Delta: βήμα κβαντισμού
    % x1_quant: κβαντισμένη τιμή
    
        x1_quant = Delta * round(x1 / Delta);
    end
    disp(t);
end
function dZ=odefun3(t,Z,y1,y2,phi1,phi2,phi3,u,gamma,sigma)

    theta=Z(1:6);

    x1_hat=[theta(1),theta(2),theta(5)]*[phi1(t);phi2(t);phi3(t)];
    x2_hat=[theta(3),theta(4),theta(6)]*[phi1(t);phi2(t);phi3(t)];
    ey1=y1(t)-x1_hat;
    ey2=y2(t)-x2_hat;
   
    dtheta(1)=gamma(1)*(ey1)*phi1(t)-gamma(1)*sigma*theta(1);
    if (theta(1)>-2 && theta(1)<0) || (theta(1)==0 && dtheta(1)<=0) || (theta(1)==-2 && dtheta(1)>=0)
        i=1;
    else
        dtheta(1)=0;
    end
    dtheta(2)=gamma(2)*(ey1)*phi2(t)-gamma(2)*sigma*theta(2);
    dtheta(5)=gamma(5)*(ey1)*phi3(t)-gamma(5)*sigma*theta(5);
    dtheta(3)=gamma(3)*(ey2)*phi1(t)-gamma(3)*sigma*theta(3);
    dtheta(4)=gamma(4)*(ey2)*phi2(t)-gamma(4)*sigma*theta(4);
    dtheta(6)=gamma(6)*(ey2)*phi3(t)-gamma(6)*sigma*theta(6);
    if (theta(6)>2) || (theta(6)==2 && dtheta(6)>=0)
        i=1;
    else
        dtheta(6)=0;
    end
    
    disp(t);
    dZ=[dtheta'];

end





function dZ=odefun2(t,Z,y1,y2,phi1,phi2,phi3,gamma,dtheta)
    theta = Z(1:6);
 


    ey1=y1(t)-theta(1)*phi1(t)-theta(2)*phi2(t)-theta(5)*phi3(t);
    ey2=y2(t)-theta(3)*phi1(t)-theta(4)*phi2(t)-theta(6)*phi3(t);

    [d1,d2]=dis(ey1,ey2,2);

    dtheta(1)=gamma(1)*(ey1+d1)*phi1(t);
    if (theta(1)<0 && theta(1)>-2) || (theta(1)==0 && dtheta(1)<=0) || (theta(1)==-2 && dtheta(1)>=0)
        i=1;
    else
        dtheta(1)=0;
    end
    dtheta(2)=gamma(2)*(ey1+d1)*phi2(t);
    dtheta(5)=gamma(5)*(ey1+d1)*phi3(t);
    dtheta(3)=gamma(3)*(ey2+d2)*phi1(t);
    dtheta(4)=gamma(4)*(ey2+d2)*phi2(t);
    dtheta(6)=gamma(6)*(ey2+d2)*phi3(t);
    if (theta(6)>1) || (theta(6)==1 && dtheta(6)>=0)
        i=1;
    else
        dtheta(6)=0;
    end
    
    disp(t);
    dZ=[dtheta'];

    function [d1,d2]=dis(ey1,ey2,omega)
        if(ey1<-omega)
            d1=omega;
        elseif(ey1>omega)
            d1=-omega;
        else
            d1=-ey1;
        end

        if(ey2<-omega)
            d2=omega;
        elseif(ey2>omega)
            d2=-omega;
        else
            d2=-ey2;
        end

    end

end