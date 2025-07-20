clear;

t_data=0 :0.001 :30;
u1=5*sin(4*t_data);
u = @(t) interp1(t_data, u1, t, 'linear', 'extrap');
[t,x]=ode45(@(t,x)odefun(t,x,u),[0 :0.001 :30],[0;0]);

y2 = @(t) interp1(t_data, x(:,2), t, 'linear', 'extrap');
y1 = @(t) interp1(t_data, x(:,1), t, 'linear', 'extrap');

gamma=[4.5,10,4,4.8,10,0.2];

[t,Z]=ode45(@(t,Z)odefun2(t,Z,u,y1,y2,gamma,[0,0,0,0,0,0]),[0 :0.001 :30],[0;0;-2;0.6;-0.84;0;1;2]);

theta1=Z(:,3);
theta2=Z(:,4);
theta3=Z(:,5);
theta4=Z(:,6);
theta5=Z(:,7);
theta6=Z(:,8);

a11=-2.15*ones(length(t),1);
a12=0.25*ones(length(t),1);
a21=-0.75*ones(length(t),1);
a22=-2*ones(length(t),1);
b1=0*ones(length(t),1);
b2=1.5*ones(length(t),1);

x_hat1=Z(:,1);
x_hat2=Z(:,2);

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

function dx=odefun(t,x,u)

    dx=[-2.15*x(1)+0.25*x(2);-0.75*x(1)-2*x(2)+1.5*u(t)];

end

function dZ=odefun2(t,Z,u,y1,y2,gamma,dtheta)
    x = Z(1:2);
    theta = Z(3:8);
    
    dx=[theta(1)*x(1)+theta(2)*x(2)+theta(5)*u(t);theta(3)*x(1)+theta(4)*x(2)+theta(6)*u(t)];
    dtheta(1)=gamma(1)*(y1(t)-x(1))*x(1);
    if (theta(1)<-1 && theta(1)>-3) || (theta(1)==-1 && dtheta(1)<=0) || (theta(1)==-3 && dtheta(1)>=0)
        i=1;
    else
        dtheta(1)=0;
    end
    dtheta(2)=gamma(2)*(y1(t)-x(1))*x(2);
    dtheta(3)=gamma(3)*(y2(t)-x(2))*x(1);
    dtheta(4)=gamma(4)*(y2(t)-x(2))*x(2);
    dtheta(5)=gamma(5)*(y1(t)-x(1))*u(t);
    dtheta(6)=gamma(6)*(y2(t)-x(2))*u(t);
    if (theta(6)>1) || (theta(6)==1 && dtheta(6)>=0)
        i=1;
    else
        dtheta(6)=0;
    end
    
    disp(t);
    dZ=[dx;dtheta'];
end

function dZ=odefun3(t,Z,u,y1,y2,gamma,dtheta)
    x = Z(1:2);
    theta = Z(3:8);
    dx=[theta(1)*y1(t)+theta(2)*y2(t)+theta(5)*u(t)+4*(y1(t)-x(1));theta(3)*y1(t)+theta(4)*y2(t)+theta(6)*u(t)+4*(y2(t)-x(2))];

    dtheta(1)=gamma(1)*(y1(t)-x(1))*y1(t);
    if (theta(1)<-1 && theta(1)>-3) || (theta(1)==-1 && dtheta(1)<=0) || (theta(1)==-3 && dtheta(1)>=0)
        i=1;
    else
        dtheta(1)=0;
    end
    dtheta(2)=gamma(2)*(y1(t)-x(1))*y2(t);
    dtheta(3)=gamma(3)*(y2(t)-x(2))*y1(t);
    dtheta(4)=gamma(4)*(y2(t)-x(2))*y2(t);
    dtheta(5)=gamma(5)*(y1(t)-x(1))*u(t);
    dtheta(6)=gamma(6)*(y2(t)-x(2))*u(t);
    if (theta(6)>1) || (theta(6)==1 && dtheta(6)>=0)
        i=1;
    else
        dtheta(6)=0;
    end

    dZ=[dx;dtheta'];
    disp(t);
end

function [omega1,omega2]=disturbances(y1,y2,omega)

    omega1_max=omega*sqrt(2)/2;
    omega2_max=omega*sqrt(2)/2;

    omega1_list=linspace(-omega1_max,omega1_max,8);
    omega2_list=linspace(-omega2_max,omega2_max,8);

    if (y1>omega1_list(end))
       omega1=omega1_list(end);
    elseif (y1<omega1_list(1))
        omega1=omega1_list(1);
    else
        for i=1:length(omega1_list)
            if y1>omega1_list(i)
                omega1=omega1_list(i);
            end
        end
    end

    if (y2>omega2_list(end))
       omega2=omega2_list(end);
    elseif (y1<omega2_list(1))
        omega2=omega2_list(1);
    else
        for i=1:length(omega2_list)
            if y2>omega2_list(i)
                omega2=omega2_list(i);
            end
        end
    end
end