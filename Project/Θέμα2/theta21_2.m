clear;
N = 40;
t_data = 0:0.001:N;

%% 1. Προσομοίωση συστήματος για διάφορες εισόδους
u1 = 1*sin(20*t_data)';
u2 = 2*sin(5*t_data)';
u3 = 1*cos(10*t_data)' + 5*sin(0.4*t_data)';
u4 = 2*ones(length(t_data), 1)+2*cos(10*t_data).^2';
u5 = 0.01*ones(length(t_data), 1) + cos(0.5*t_data)';
u6 = 5*ones(length(t_data), 1) + 2*cos(0.05*t_data)';

% Συνάρτηση για προσομοίωση συστήματος
odefun = @(t, x, u) -x^3 + tanh(x) + (1/(1+x^2)) + u(t);

% Δημιουργία δεδομένων
U = [u1, u2, u3, u4, u5, u6];
X = zeros(length(t_data), 6);

for i = 1:6
    u_fun = @(t) interp1(t_data, U(:,i), t, 'linear', 'extrap');
    [~, x] = ode45(@(t,x) odefun(t,x,u_fun), t_data, 0);
    X(:,i) = x;
end

%% 2. Cross-Validation με φιλτραρισμένο μοντέλο (όρος 1/(s+1)x με συντελεστή 1)
k_folds = 6;
theta = zeros(6, k_folds);  % 5 παράμετροι + σφάλμα
Em=0;
for fold = 1:k_folds
    theta0 = zeros(6, 1);  % 5 παράμετροι (χωρίς τον όρο 1/(s+1)x)
   
    x0 = [theta0];
    
    theta_temp = zeros(5, 1);
    gamma=10*[1,1,1,1,1];
    for i = 1:6
        if i ~= fold
            y = @(t) interp1(t_data, X(:,i), t, 'linear', 'extrap');
            u = @(t) interp1(t_data, U(:,i), t, 'linear', 'extrap');
            [~, x] = ode15s(@(t,x) odefun3(t, x, u,y, gamma), t_data, x0);
            theta_temp = theta_temp + 0.2 * mean(x(:, 2:6))';
        end
    end
        
    u = @(t) interp1(t_data, U(:,fold), t, 'linear', 'extrap');
    [~, Y] = ode15s(@(t,x) odefun2(t, x, u,theta_temp), t_data, 0);
    total_error = mean((Y - X(:,fold)).^2);
    disp(['Fold-',num2str(fold),' mse:',num2str(total_error)]);
    Em=Em+(1/6)*total_error;
    theta(:, fold) = [theta_temp; total_error];

    figure();
    plot(t_data,Y,t_data,X(:,fold));
    title(['Εκπαίδευση fold:',num2str(fold)]);
    legend({'Εκτιμώμενη έξοδος','Πραγματική έξοδος'});
    grid on;

end

disp(['Μέσο τετραγωνικό σφάλμα κατά την εκπαίδευση:',num2str(Em)]);

%% 3. Επιλογή καλύτερου μοντέλου (ελάχιστο σφάλμα)
[~, idx] = min(theta(end, :));
theta_final = theta(1:5, idx);  % Οι 5 παράμετροι (χωρίς το σφάλμα)

%% 4. Επικύρωση με νέα δεδομένα
u1=ones(length(t_data),1);
u2=4*(cos(0.2*t_data).*sin(4*t_data))'+2*ones(length(t_data),1);
u3=ones(length(t_data),1)+1.5*sin(0.5*t_data)'+0.4*cos(12*t_data)';
u4=2*sin(20*t_data).^(3)'+0.5*ones(length(t_data),1);

U1=[u1,u2,u3,u4];
for i=1:4
u_val = U1(:,i);
[~, x_val] = ode45(@(t,x) odefun(t, x, @(t) interp1(t_data, u_val, t)), t_data, 0);

% Πρόβλεψη με τελικό μοντέλο


[~,y_pred] = ode45(@(t,x) odefun2(t, x, @(t) interp1(t_data, u_val, t),theta_final), t_data, 0);
val_error = mean((y_pred - x_val).^2);
fprintf('Σφάλμα Επικύρωσης: %.4f\n', val_error);

%% 5. Γραφήματα
figure();
plot(t_data, x_val, 'b', t_data, y_pred, 'r');
legend({'Πραγματικό x(t)', 'Προβλεπόμενο x(t)'});
title('Αποτελέσματα Μοντέλου με Όρο 1/(s+1)x');
grid on;

end


function dx=odefun2(t,x,u,theta)

    phi1 = exp(-(x - u(t))^2 / 1);
    phi2 = exp(-(x - u(t))^2 / 10);
    phi3 = exp(-(x - u(t))^2 / 2);
    phi4 = u(t);
    phi5 = cos(x);
    
    dx=theta(1)*phi1+theta(2)*phi2+theta(3)*phi3+theta(4)*phi4+theta(5)*phi5;


end

function dZ=odefun3(t,Z,u,y,gamma)
    x = Z(1);
    theta = Z(2:6);

    phi1 = exp(-(y(t) - u(t)).^2 / 1);
    phi2 = exp(-(y(t) - u(t)).^2 / 10);
    phi3 = exp(-(y(t)- u(t)).^2 / 2);
    phi4 = u(t);
    phi5 = cos(y(t));

    dx=theta(1)*phi1+theta(2)*phi2+theta(3)*phi3+theta(4)*phi4+theta(5)*phi5+6.5*(y(t)-x);
    dtheta(1)=gamma(1)*(y(t)-x)*phi1;
    dtheta(2)=gamma(2)*(y(t)-x)*phi2;
    dtheta(3)=gamma(3)*(y(t)-x)*phi3;
    dtheta(4)=gamma(4)*(y(t)-x)*phi4;
    dtheta(5)=gamma(5)*(y(t)-x)*phi5;

    dZ=[dx;dtheta'];

end