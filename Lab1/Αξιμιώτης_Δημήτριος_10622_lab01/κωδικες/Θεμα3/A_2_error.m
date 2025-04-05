
A=[0.01,0.1,0.2,0.5,2.5,5,10,50,100];
mse=[];
for i=1:length(A)   
    mse=[mse create_error_graph_8(A(i))];
    disp(i);
end

figure();
plot(A,mse,LineWidth=0.8);
xlabel('Α');
ylabel('MSE');
grid on;


function [mse,e]=create_error_graph_8(A)
[t,x,u]=getsignals_A(A);

sys1=tf([-1 0],[1 3 2]);
sys2=tf([0 -1],[1 3 2]);
sys3=tf([0 1],[1 3 2]);

Ts=0.1;
Dt=0.0001;
step=Ts/Dt;

k=[1:step:length(t)]';
t_new=[0:Ts:20]';

[x1,~]=lsim(sys1,x(k,1),t_new);
[x2,~]=lsim(sys2,x(k,1),t_new);
[x3,~]=lsim(sys3,u(k),t_new);

X=[x1,x2,x3];
y=x(k,1);

[theta,~]=least_squares(X,y);

L_hat = 9.81 / (theta(2)+2);
m_hat = 1 / (theta(3)*L_hat^2);
c_hat = (3+theta(1))*(m_hat*L_hat^2);

[t,q_hat]=ode45(@(t,q_hat)odefun(t,q_hat,m_hat,L_hat,c_hat,A),[0 :Dt :20],[0;0]);

disp([L_hat,m_hat,c_hat]);

e=x(:,1)-q_hat(:,1);
mse=mean(e.^(2));

figure();
plot(t,e);
title(["Error plot ",'A=',num2str(A)]);
xlabel('t');
ylabel('e');
legend(['L=',num2str(L_hat),' ','m=',num2str(m_hat),' ','c=',num2str(c_hat),' ']);
grid on;

function dq=odefun(t,q,m_hat,L_hat,c_hat,A0)

    % ορισμός παραμέτρων 
    g=9.81;
    omega=2;
    % υλοποίηση διαφορικής εξίσωσης
    dq=[q(2);(-c_hat/(m_hat*L_hat^2))*q(2)-(g/L_hat)*q(1)+(1/(m_hat*L_hat^2))*A0*cos(omega*t)];

end
end