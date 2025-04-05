function [t,x,u]=getsignals_A(A)
% χρονικές εξισώσεις
[t,x]=ode45(@(t,x)odefun(t,x,A),[0 :0.0001 :20],[0;0]);
u=A*cos(2*t);

function dx=odefun(t,x,A0)

    % ορισμός παραμέτρων 
    m=0.75;
    L=1.25;
    c=0.15;
    g=9.81;
    omega=2;

    % υλοποίηση διαφορικής εξίσωσης
    dx=[x(2);(-c/(m*L^2))*x(2)-(g/L)*x(1)+(1/(m*L^2))*A0*cos(omega*t)];

end

end