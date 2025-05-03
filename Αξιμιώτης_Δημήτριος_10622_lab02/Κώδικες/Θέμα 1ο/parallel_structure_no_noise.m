[t,x]=ode45(@odefun,[0 :0.0001 :50],[0;0]);

t_data=t;
u=2.5*sin(t_data);

y_dot= @(t) interp1(t_data, x(:,2), t, 'linear', 'extrap');
y = @(t) interp1(t_data, x(:,1), t, 'linear', 'extrap');
u = @(t) interp1(t_data, u, t, 'linear', 'extrap');

%[0.1,0.8,1]
%[0.1,3,4]
%[10,10,10]

gamma=[0.1,0.86,1.27];

[t,Z]=ode15s(@(t,Z)odefun2(t,Z,u,y,y_dot,gamma),[0 :0.0001 :50],[0;0;0.2;0.2;0.2]);

theta1=Z(:,3);
theta2=Z(:,4);
theta3=Z(:,5);

m_hat=1./theta3;
k_hat=m_hat.*theta2;
b_hat=m_hat.*theta1;

x_hat1=Z(:,1);
x_hat2=Z(:,2);

m=1.315*ones(length(t_data),1);
k=0.725*ones(length(t_data),1);
b=0.225*ones(length(t_data),1);

figure();
plot(t,m_hat,t,k_hat,t,b_hat);
hold on;
plot(t,m,'Color','b','LineStyle','--');
hold on;
plot(t,k,'Color','r','LineStyle','--');
hold on;
plot(t,b,'Color','y','LineStyle','--');
title('Εκτιμήσεις παραμέτρων');
legend({'$\hat{m}$','$\hat{k}$','$\hat{b}$','m','k','b'},'Interpreter','latex');
grid on;

figure();
plot(t,x(:,1),t,x_hat1);
legend({'x(t)','$\hat{x}(t)$'},'Interpreter','latex');
grid on;

error=x(:,2)-x_hat2;

figure();
plot(t,error);
title('e=x(t)-$\hat{x}$(t)','Interpreter','latex');
grid on;

function dx=odefun(t,x)

    % ορισμός παραμέτρων 
    b=0.225;
    k=0.725;
    m=1.315;
    % υλοποίηση διαφορικής εξίσωσης
    dx=[x(2);-(b/m)*x(2)-(k/m)*x(1)+(1/m)*2.5*sin(t)];

end

function dZ=odefun2(t,Z,u,y,y_dot,gamma)
    x = Z(1:2);
    theta = Z(3:5);
    dx=[x(2);-theta(1)*x(2)-theta(2)*x(1)+theta(3)*2.5*sin(t)];
    dtheta(1)=-gamma(1)*((y_dot(t)-x(2))*x(2));
    dtheta(2)=-gamma(2)*((y_dot(t)-x(2))*x(1));
    dtheta(3)=gamma(3)*((y_dot(t)-x(2))*u(t));
    dZ=[dx;dtheta'];
    disp(t);
end