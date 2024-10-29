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
stepnumber = 60+1;
results = zeros(10,stepnumber);
%% Exact VI/PI with Intial point P0=P*+0.01E

P0 = S+0.01*eye(3);

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
    A_cl = A - B*Ki_PI;
    Q_eff = Ki_PI'*R*Ki_PI + Q;
    Pi_PI = dlyap(A_cl',Q_eff);
    Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;
end


%% Inxact VI/PI with Intial point P0=P*+0.01E and errors = 0.01^i 

error_A0 = 0.01*eye(3); 
error_B0 = 0.01*eye(3);
%VI algorithm
Ai = A + error_A0 ;
Bi= B + error_B0 ;
error_Ai = error_A0;
error_Bi = error_B0;
Pi_VI = P0;
for i = 1:stepnumber
    error_VI = norm(Pi_VI - S)/norm(S);
    results(3, i) = error_VI;
    Piplus1_VI = Ai'* Pi_VI * Ai - Ai'*Pi_VI*Bi*inv(R+Bi'*Pi_VI*Bi)*Bi'*Pi_VI*Ai + Q;
    Pi_VI = Piplus1_VI;
    error_Ai = 0.01*error_Ai; 
    error_Bi = 0.01*error_Bi; 
    Ai = A + error_Ai ;
    Bi = B + error_Bi ;
end

%PI algorithm

Ai = A + error_A0 ;
Bi= B + error_B0 ;
error_Ai = error_A0;
error_Bi = error_B0;
Pi_PI = P0;
Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;

for i = 1:stepnumber
    error_PI = norm(Pi_PI - S)/norm(S);
    results(4, i) = error_PI;
    A_cl = Ai - Bi*Ki_PI;
    Q_eff = Ki_PI'*R*Ki_PI + Q;

    Pi_PI = dlyap(A_cl',Q_eff);

    error_Ai = 0.01*error_Ai; 
    error_Bi = 0.01*error_Bi; 
    Ai = A + error_Ai ;
    Bi = B + error_Bi ;
    Ki_PI = inv(R+Bi'*Pi_PI*Bi)*Bi'*Pi_PI*Ai;
end

%% Inxact VI/PI with Intial point P0=P*+0.01E and errors = 0.01^i 

error_A0 = 0.01*eye(3); 
error_B0 = 0.01*eye(3);
%VI algorithm
 
Ai = A + error_A0;
Bi= B + error_B0;
error_Ai = error_A0;
error_Bi = error_B0;
Pi_VI = P0;
for i = 1:stepnumber
    error_VI = norm(Pi_VI - S)/norm(S);
    results(5, i) = error_VI;
    Piplus1_VI = Ai'* Pi_VI * Ai - Ai'*Pi_VI*Bi*inv(R+Bi'*Pi_VI*Bi)*Bi'*Pi_VI*Ai + Q;
    Pi_VI = Piplus1_VI;
    error_Ai = 0.9*error_Ai; 
    error_Bi = 0.9*error_Bi; 
    Ai = A + error_Ai;
    Bi = B + error_Bi;
end

%PI algorithm

Ai = A + error_A0;
Bi= B + error_B0;
error_Ai = error_A0;
error_Bi = error_B0;
Pi_PI = P0;
Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;

for i = 1:stepnumber
    error_PI = norm(Pi_PI - S)/norm(S);
    results(6, i) = error_PI;
    A_cl = Ai - Bi*Ki_PI;
    Q_eff = Ki_PI'*R*Ki_PI + Q;

    Pi_PI = dlyap(A_cl',Q_eff);

    error_Ai = 0.9*error_Ai; 
    error_Bi = 0.9*error_Bi; 
    Ai = A + error_Ai;
    Bi = B + error_Bi;
    Ki_PI = inv(R+Bi'*Pi_PI*Bi)*Bi'*Pi_PI*Ai;
end

%% Inxact VI/PI with Intial point P0=P*+0.01E and errors = 0.01 * 0.9^(i-1) 

error_A0 = 0.01*eye(3); 
error_B0 = 0.01*eye(3);
%VI algorithm
 
Ai = A + error_A0;
Bi= B + error_B0;
error_Ai = error_A0;
error_Bi = error_B0;
Pi_VI = P0;
for i = 1:stepnumber
    error_VI = norm(Pi_VI - S)/norm(S);
    results(5, i) = error_VI;
    Piplus1_VI = Ai'* Pi_VI * Ai - Ai'*Pi_VI*Bi*inv(R+Bi'*Pi_VI*Bi)*Bi'*Pi_VI*Ai + Q;
    Pi_VI = Piplus1_VI;
    error_Ai = 0.9*error_Ai; 
    error_Bi = 0.9*error_Bi; 
    Ai = A + error_Ai;
    Bi = B + error_Bi;
end

%PI algorithm

Ai = A + error_A0;
Bi= B + error_B0;
error_Ai = error_A0;
error_Bi = error_B0;
Pi_PI = P0;
Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;

for i = 1:stepnumber
    error_PI = norm(Pi_PI - S)/norm(S);
    results(6, i) = error_PI;
    A_cl = Ai - Bi*Ki_PI;
    Q_eff = Ki_PI'*R*Ki_PI + Q;

    Pi_PI = dlyap(A_cl',Q_eff);

    error_Ai = 0.9*error_Ai; 
    error_Bi = 0.9*error_Bi; 
    Ai = A + error_Ai;
    Bi = B + error_Bi;
    Ki_PI = inv(R+Bi'*Pi_PI*Bi)*Bi'*Pi_PI*Ai;
end

%% Inxact VI/PI with Intial point P0=P*+0.01E and errors = 0.01 * 0.9^(i-1) +0.001I

error_A0 = 0.01*eye(3); 
error_B0 = 0.01*eye(3);
%VI algorithm
 
Ai = A + error_A0 + 0.001*eye(3);
Bi= B + error_B0 + 0.001*eye(3);
error_Ai = error_A0;
error_Bi = error_B0;
Pi_VI = P0;
for i = 1:stepnumber
    error_VI = norm(Pi_VI - S)/norm(S);
    results(7, i) = error_VI;
    Piplus1_VI = Ai'* Pi_VI * Ai - Ai'*Pi_VI*Bi*inv(R+Bi'*Pi_VI*Bi)*Bi'*Pi_VI*Ai + Q;
    Pi_VI = Piplus1_VI;
    error_Ai = 0.6*error_Ai; 
    error_Bi = 0.6*error_Bi; 
    Ai = A + error_Ai + 0.001*eye(3);
    Bi = B + error_Bi + 0.001*eye(3);
end

%PI algorithm

Ai = A + error_A0;
Bi= B + error_B0;
error_Ai = error_A0;
error_Bi = error_B0;
Pi_PI = P0;
Ki_PI = inv(R+B'*Pi_PI*B)*B'*Pi_PI*A;

for i = 1:stepnumber
    error_PI = norm(Pi_PI - S)/norm(S);
    results(8, i) = error_PI;
    A_cl = Ai - Bi*Ki_PI;
    Q_eff = Ki_PI'*R*Ki_PI + Q;

    Pi_PI = dlyap(A_cl',Q_eff);

    error_Ai = 0.6*error_Ai; 
    error_Bi = 0.6*error_Bi; 
    Ai = A + error_Ai + 0.001*eye(3);
    Bi = B + error_Bi + 0.001*eye(3);
    Ki_PI = inv(R+Bi'*Pi_PI*Bi)*Bi'*Pi_PI*Ai;
end
%% settings for plot
figure;
x = 0:(size(results, 2) - 1);  

plot(x,results(1, :),'-','Color','red','LineWidth',1)
hold on
plot(x,results(2, :),'-','Color','blue','LineWidth',1)
hold on
%plot(x,results(3, :),'LineStyle',':','Color','red','LineWidth',1.5)
%hold on
%plot(x,results(4, :),'LineStyle',':','Color','blue','LineWidth',1.5)
%hold on
plot(x,results(5, :),'--','Color','red','LineWidth',2)
hold on
plot(x,results(6, :),'--','Color','blue','LineWidth',2)
hold on
plot(x,results(7, :),':','Color','red','LineWidth',1.5)
hold on
plot(x,results(8, :),':','Color','blue','LineWidth',1.5)
hold on
ylim([0 .27])
ylabel('$\frac{\|\hat{P}_i-P^*\|}{\|P^*\|}$','interpreter','latex','FontSize',12)
xlabel('$\mathrm{iteration}~i$','interpreter','latex','FontSize',12)
%legend('exact VI with $A,B$','exact PI with $A,B$',...
%    'inexact VI with $\hat{A}_i^1 = A+0.01 \times 0.9^iI,\tilde{B} = B+0.01 \times 0.9^iI$','inexact PI with $\tilde{A} = A+0.01 \times 0.9iI,\tilde{B} = B+0.01 \times 0.9^iI$',...
 %   'inexact VI with $\hat{A}_i^2 = A+0.01 \times 0.9^iI +0.001I,\tilde{B} = B+0.01 \times 0.9^iI +0.001I$','inexact PI with $\tilde{A} = A+0.01 \times 0.9^iI +0.001I,\tilde{B} = B+0.01 \times 0.9^iI +0.001I$',...
 %   'interpreter','latex','FontSize',9)
legend('exact VI with $A,B$','exact PI with $A,B$',...
   'inexact VI with $\hat{A}_i^{m},\hat{B}_i^m$','inexact PI with $\hat{A}_i^{m},\hat{B}_i^m$',...
   'inexact VI with $\hat{A}_i^{n},\hat{B}_i^n$','inexact PI with $\hat{A}_i^{n},\hat{B}_i^n$',...
   'interpreter','latex','FontSize',9)
% save as .png
saveas(gcf, 'Comparision_Robustness.png');

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

