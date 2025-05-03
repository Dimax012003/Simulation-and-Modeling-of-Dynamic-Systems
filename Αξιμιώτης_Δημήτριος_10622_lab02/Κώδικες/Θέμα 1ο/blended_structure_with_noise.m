clear;

n0=0:0.1:1;

% 0.25
% 10
% 25
% 100

[t,error,e_m,e_k,e_b]=noise_data(0.25);
figure();
plot(t,e_m,t,e_b,t,e_k);
legend({'$\tilde{m}$','$\tilde{b}$','$\tilde{k}$'},'Interpreter','latex');
title('Σφάλματα εκτιμήσεων παραμέτρων');
grid on;

figure();
plot(t,error);
title('Σφάλμα εκτίμησης e(t)=x(t)-$\hat{x}$(t)','Interpreter','latex');
grid on;

 E_m=zeros(200001,length(n0));
 E_k=zeros(200001,length(n0));
 E_b=zeros(200001,length(n0));
 
 
 for i=1:length(n0)
     [t,error,e_m,e_k,e_b]=noise_data(n0(i));
     E_m(:,i)=e_m;
     E_k(:,i)=e_k;
     E_b(:,i)=e_b;
 
 end
 
 figure(1);
 for i = 1:length(n0)
     plot(t, E_m(:,i), 'DisplayName', num2str(n0(i)));
     hold on;
 end
 title('$\tilde{m}$','Interpreter','latex');
 grid on;
 legend show;
 
 figure(2);
 for i = 1:length(n0)
     plot(t, E_k(:,i), 'DisplayName', num2str(n0(i)));
     hold on;
 end
 title('$\tilde{k}$','Interpreter','latex');
 grid on;
 legend show;
 
 figure(3);
 for i = 1:length(n0)
     plot(t, E_b(:,i), 'DisplayName', num2str(n0(i)));
     hold on;
 end
 title('$\tilde{b}$','Interpreter','latex');
 grid on;
 legend show;


function [t,error,e_m,e_k,e_b]=noise_data(n0)
    [t,x]=ode45(@odefun,[0 :0.0001 :20],[0;0]);

    t_data=t;
    u=2.5*sin(t_data);

    y_dot= @(t) interp1(t_data, x(:,2), t, 'linear', 'extrap');
    y = @(t) interp1(t_data, x(:,1), t, 'linear', 'extrap');
    u = @(t) interp1(t_data, u, t, 'linear', 'extrap');

    gamma=[0.1,0.4,0.7];


    [t,Z]=ode45(@(t,Z)odefun2(t,Z,u,y,y_dot,gamma,n0),[0 :0.0001 :20],[0;0;0.2;0.3;0.2]);

    theta1=Z(:,3);
    theta2=Z(:,4);
    theta3=Z(:,5);
    x_hat1=Z(:,1);
    error=x(:,1)-x_hat1;

    m_hat=1./theta3;
    k_hat=m_hat.*theta2;
    b_hat=m_hat.*theta1;

    m=1.315*ones(length(t_data),1);
    k=0.725*ones(length(t_data),1);
    b=0.225*ones(length(t_data),1);

    e_m=abs(m_hat-m)/1.315;
    e_k=abs(k_hat-k)/0.725;
    e_b=abs(b_hat-b)/0.225;

end

function dx=odefun(t,x)
        % ορισμός παραμέτρων 
        b=0.225;
        k=0.725;
        m=1.315;
        % υλοποίηση διαφορικής εξίσωσης
        dx=[x(2);-(b/m)*x(2)-(k/m)*x(1)+(1/m)*2.5*sin(t)];
    
end
function dZ=odefun2(t,Z,u,y,y_dot,gamma,n0)
        x = Z(1:2);
        theta = Z(3:5);
        n=n0*sin(2*pi*20*t);
        dx=[y_dot(t)+2*(y(t)+n-x(1));-theta(1)*y_dot(t)-theta(2)*(y(t)+n)+theta(3)*2.5*sin(t)+2*(y_dot(t)-x(2))];
        dtheta(1)=-gamma(1)*(y_dot(t)-x(2))*y_dot(t);
        dtheta(2)=-gamma(2)*(y_dot(t)-x(2))*(y(t)+n);
        dtheta(3)=gamma(3)*(y_dot(t)-x(2))*u(t);
        dZ=[dx;dtheta'];
        disp(t);
end