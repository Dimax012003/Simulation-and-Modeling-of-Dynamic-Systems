clear;

phi_inf=0.01;
phi_0=100;
lamda=10;
rho=100;
k1=10;
k2=20;

t1=[0:0.01:20];
r1=-(pi/1000)*t1.^2+(pi/50)*t1;
rd= @(t) interp1(t1, r1, t, 'linear', 'extrap');

[t,x]=ode15s(@(t,x)odefun(t,x,phi_inf,phi_0,lamda,rho,k1,k2,rd),[0:0.01:20],[0;0]);

phi=(phi_0-phi_inf)*exp(-lamda*t)+phi_inf;
z1=(x(:,1)-r1')./phi;
T1=log((1+z1)./(1-z1));
a=-k1*T1;
z2=(x(:,2)-a)./rho;
T2=log((1+z2)./(1-z2));
u=-k2*T2;

t_data=t;

x2 = @(t) interp1(t_data, x(:,2), t, 'linear', 'extrap');
x1 = @(t) interp1(t_data, x(:,1), t, 'linear', 'extrap');
u  = @(t) interp1(t_data, u, t, 'linear', 'extrap');

gamma=[180000,3090000,150000,9000];

[t,z]=ode15s(@(t,Z)odefun2(t,Z,x1,x2,u,gamma),[0:0.01:20],[0;0;0.01;0.01;0.25;0.01]);

x_hat=z(:,1);
a1_hat=z(:,3);
a2_hat=z(:,4);
a3_hat=z(:,5);
b_hat=z(:,6);

error=x(:,1)-x_hat;

a1=1.315*ones(length(t),1);
a2=0.725*ones(length(t),1);
a3=0.225*ones(length(t),1);
b=1.175*ones(length(t),1);

figure();
plot(t,a1_hat,t,a1);
title('$hat{a}_1$','Interpreter','latex');
legend({'$\hat{a}_1(t)$','a_1'},'Interpreter','latex');
grid on;

figure();
plot(t,a2_hat,t,a2);
title('$\hat{a}_2$','Interpreter','latex');
legend({'$\hat{a}_2(t)$','a_2'},'Interpreter','latex');
grid on;

figure();
plot(t,a3_hat,t,a3);
title('$\hat{a}_3$','Interpreter','latex');
legend({'$\hat{a}_3(t)$','a_3'},'Interpreter','latex');
grid on;

figure();
plot(t,b_hat,t,b);
title('$\hat{b}$','Interpreter','latex');
legend({'$\hat{b}(t)$','b'},'Interpreter','latex');
grid on;

figure();
plot(t,error);
title('$r(t)-\hat{r}(t)$','Interpreter','latex');
grid on;

figure();
plot(t,x(:,1),t,x_hat);
legend({'r(t)','$\hat{r}(t)$'},'Interpreter','latex');
grid on;

function dx=odefun(t,x,phi_inf,phi_0,lamda,rho,k1,k2,rd)

    a1=1.315;
    a2=0.725;
    a3=0.225;
    b1=1.175;

    phi=(phi_0-phi_inf)*exp(-lamda*t)+phi_inf;
    z1=(x(1)-rd(t))/phi;

    T1=log((1+z1)/(1-z1));
    a=-k1*T1;
   
    z2=(x(2)-a)/rho;
    T2=log((1+z2)/(1-z2));

    u=-k2*T2;
    dx=[x(2);-a1*x(2)-a2*sin(x(1))+a3*x(2)^2*sin(2*x(1))+b1*u+0.15*sin(0.5*t)];
    disp(t);
end

function dZ=odefun2(t,Z,x1,x2,u,gamma)
    x=Z(1:2);
    theta=Z(3:6);
    dx=[x2(t)+1000*(x1(t)-x(1));-theta(1)*x2(t)-theta(2)*sin(x1(t))+theta(3)*x2(t)^2*sin(2*x1(t))+theta(4)*u(t)+1000*(x2(t)-x(2))];
    dtheta(1)=-gamma(1)*(x2(t)-x(2))*x2(t);
    dtheta(2)=-gamma(2)*sin(x1(t))*(x2(t)-x(2));
    dtheta(3)=gamma(3)*(x2(t)^2)*sin(2*x1(t))*(x2(t)-x(2));
    dtheta(4)=gamma(4)*u(t)*(x2(t)-x(2));
    disp(t);
    dZ=[dx;dtheta'];
end