% Define the system matrix
A = [1.01 0.01 0; 0.01 1.01 0.01; 0 0.01 1.01];
B = [1 0 0 ;0 1 0; 0 0 1];
% C = [1 1];
% D = 0;

% Define the LQR weight matrix
Q = [0.001 0 0; 0 0.001 0; 0 0 0.001];
R = [1 0 0; 0 1 0; 0 0 1];

% Calculate the LQR gain matrix K
[K, S, e] = dlqr(A, B, Q, R);

disp('is A stable?:');
disp(isStable(A));

% Print results
disp('The LQR gain matrix K:');
disp(K);
disp('Solution of the Riccati equation P*:');
disp(S);

% max iteration step
stepnumber = 40+1;
results = zeros(10,stepnumber);

%% Intial point P0=0
P0 = [0 0 0;0 0 0;0 0 0];

%VI algorithm
Pi_VI = P0;
for i = 1:stepnumber
    error_VI = norm(Pi_VI - S)/norm(S);
    results(1, i) = error_VI;
    Piplus1_VI = A'* Pi_VI * A - A'*Pi_VI*B*inv(R+B'*Pi_VI*B)*B'*Pi_VI*A + Q;
    Pi_VI = Piplus1_VI;
end
%PI algorithm
Pi_PI = P0;
Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;

for i = 1:stepnumber
    error_PI = norm(Pi_PI - S)/norm(S);
    results(2, i) = error_PI;
    % Calculate the closed-loop system matrix A - BK
    A_cl = A - B*Ki_PI;
    % Calculation Q_eff = K'RK + Q
    Q_eff = Ki_PI'*R*Ki_PI + Q;
    % Solve lyapunov eq and update Ki
    Pi_PI = dlyap(A_cl',Q_eff);
    Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;
end
%% Intial point P0=2P*

P0 = S + S;

%VI algorithm
Pi_VI = P0;
for i = 1:stepnumber
    error_VI = norm(Pi_VI - S)/norm(S);
    results(3, i) = error_VI;
    Piplus1_VI = A'* Pi_VI * A - A'*Pi_VI*B*inv(R+B'*Pi_VI*B)*B'*Pi_VI*A + Q;
    Pi_VI = Piplus1_VI;
end

%PI algorithm
Pi_PI = P0;
Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;

for i = 1:stepnumber
    error_PI = norm(Pi_PI - S)/norm(S);
    results(4, i) = error_PI;
    A_cl = A - B*Ki_PI;
    Q_eff = Ki_PI'*R*Ki_PI + Q;

    Pi_PI = dlyap(A_cl',Q_eff);
    Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;
end

%% Intial point P0=0.5P*

P0 = 0.5*S;
%The range 0.3–0.4 is the boundary between convergence and non-convergence, while 0.5–0.6 is the boundary for whether the first step converges.
%PI algorithm
Pi_PI = P0;
Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;

for i = 1:stepnumber
    error_PI = norm(Pi_PI - S)/norm(S);
    results(5, i) = error_PI;
    A_cl = A - B*Ki_PI;
    Q_eff = Ki_PI'*R*Ki_PI + Q;

    Pi_PI = dlyap(A_cl',Q_eff);
    Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;
end

%% Intial point P0=0.7P*

P0 = 0.7*S;
%PI algorithm
Pi_PI = P0;
Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;

for i = 1:stepnumber
    error_PI = norm(Pi_PI - S)/norm(S);
    results(6, i) = error_PI;
    A_cl = A - B*Ki_PI;
    Q_eff = Ki_PI'*R*Ki_PI + Q;

    Pi_PI = dlyap(A_cl',Q_eff);
    Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;
end

%% settings for plot
figure;

x = 0:(size(results, 2) - 1);  
plot(x,results(1, :),'--','Color','red','LineWidth',1.5)
hold on
plot(x,results(3, :),'-','Color','red','LineWidth',1)
hold on
plot(x,results(2, :),'--','Color','blue','LineWidth',1.5)
hold on
plot(x,results(4, :),'-','Color','blue','LineWidth',1)
hold on
plot(x,results(5, :),':','Color','magenta','LineWidth',2)
hold on
plot(x,results(6, :),'-.','Color','black','LineWidth',2)
hold on

ylim([0,1.5])
ylabel('$\frac{\|P_i-P^*\|}{\|P^*\|}$','interpreter','latex','FontSize',12)
xlabel('$\mathrm{iteration}~i$','interpreter','latex','FontSize',12)
legend('VI with $P_{0}$ = 0','VI with $P_{0} = 2P^*$','PI with $P_{0}$ = 0','PI with $P_{0} = 2P^*$','PI with $P_{0} = 0.5P^* \in \mathcal{B}_{\delta_0}(P^*)$','PI with $P_{0} = 0.7P^*\in \mathcal{B}_{\delta_1}(P^*)$','interpreter','latex','FontSize',9)

% save as .png
saveas(gcf, 'Comparision_Convergence.png');

function stability = isStable(A)
    % Calculating eigenvalues
    eig_values = eig(A);
    
    % stability
    if all(abs(eig_values) < 1)
        stability = true; % System stable
    else
        stability = false; % System not stable
    end
end



