%%%% Applying Lagrangian Relaxation in Integer Programming

%%%% Subgradient Method
numiter = 10;
Aineq = [1, 1, 0, 0;...
         0, 0, 1, 1];
bineq = [1; 1];
lb = zeros(4, 1);
ub = ones(4, 1);
ctype = repmat("U", 1, 2);
vartype = repmat("I", 1, 4);

% Constant alphak
lambdak = 0; % initial value
alphak = 1.0; % initial value
lambdavals = zeros(numiter, 1);
for iter = 1:numiter
    c = [ 16 - 8*lambdak; 10 - 2*lambdak; -lambdak; 4 - 4*lambdak];
    xk = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
    b_Ax = 10 - [8, 2, 2, 4]*xk;
    lambdak = max(0, lambdak - alphak*b_Ax);
    lambdavals(iter) = lambdak;
end
plot(lambdavals, 'ro-');
hold on;

% alphak = alphak/2
lambdak = 0;
alphak = 1.0;
lambdavals = zeros(numiter,1);
for iter=1:numiter
    c = [16 - 8*lambdak; 10 - 2*lambdak; -lambdak; 4 - 4*lambdak];
    xk = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
    b_Ax = 10 - [8, 2, 2, 4]*xk;
    lambdak = max(0, lambdak - alphak*b_Ax);
    lambdavals(iter) = lambdak;
    alphak = alphak/2.0;
end
plot(lambdavals, 'b*-');
hold on;

% alphak = alphak/3
lambdak = 0;
alphak = 1.0;
lambdavals = zeros(numiter,1);
for iter=1:numiter
    c = [16 - 8*lambdak; 10 - 2*lambdak; -lambdak; 4 - 4*lambdak];
    xk = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
    b_Ax = 10 - [8, 2, 2, 4]*xk;
    lambdak = max(0, lambdak - alphak*b_Ax);
    lambdavals(iter) = lambdak;
    alphak = alphak/3.0;
end

plot(lambdavals, 'g+-'); xlabel('Iterations (k)'); ylabel('lambda^{(k)}'); 
legend('alfa_k = 1', 'alfa_k = alfa_{k-1}/2', 'alfa_k = alfa_{k-1}/3', ...
    'Location','northoutside', 'Orientation', 'horizontal');
grid on;

%%%% Changing step length

lambdak = 0;
alphak = 1.0;
Zstar = 0; % Initial value
betak = 2; % Initial value
Zstarvals = zeros(numiter,1);
lambdavals = zeros(numiter,1);
ZLkprev = realmax;
lb = zeros(4, 1); ub = ones(4, 1);
ctype = repmat("U", 1, 2); vartype = repmat("I", 1, 4);
for iter=1:numiter
    c = [16 - 8*lambdak; 10 - 2*lambdak; -lambdak; 4 - 4*lambdak];
    [xk, ZLk] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
    ZLk = ZLk + 10*lambdak;
    if (ZLkprev < ZLk)
        rediter = rediter + 1;
    else
        rediter = 1;
    end
    ZLkprev = ZLk;
    b_Ax = 10 - [8, 2, 2, 4]*xk;
    lambdak = max(0, lambdak - alphak*b_Ax);
    lambdavals(iter) = lambdak;
    if (b_Ax >= 0)
        ZL = [16, 10, 0, 4]*xk;
        if (ZL > Zstar)
            Zstar = ZL;
        end
    end
    if (rediter == 3) % Objective does not decrease for 3 consecutive iterations
        betak = betak/2.0;
        rediter = 1;
    end
    alphak = (betak*(ZLk - Zstar))/(b_Ax^2);
    Zstarvals(iter) = Zstar;
end
subplot(2,1,1);
plot(Zstarvals, 'ro-'); xlabel('Iterations (k)'); ylabel('Z*'); 
grid on;
subplot(2,1,2);
plot(lambdavals, 'ro-'); xlabel('Iterations (k)'); ylabel('lambda^{(k)}'); 
grid on;



%%%%%% LP Relaxation

lambdak = 0;
alphak = 1.0;
zlkvals = zeros(numiter,1);
lb = zeros(4, 1); ub = ones(4, 1);
ctype = repmat("U", 1, 2); vartype = repmat("I", 1, 4);
for iter=1:numiter
    c = [16 - 8*lambdak; 10 - 2*lambdak; -lambdak; 4 - 4*lambdak];
    [xk, ZLk] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
    b_Ax = 10 - [8, 2, 1, 4]*xk;
    lambdak = max(0, lambdak - alphak*b_Ax);
    zlkvals(iter) = ZLk + 10*lambdak;
    alphak = alphak/2.0;
end
plot(zlkvals, 'b*-'); xlabel('Iterations (k)'); ylabel('Z_L(lambda^{(k)})'); 
grid on;

%LP relaxation of the original problem
Aineqlp = [8, 2, 2, 4; ...
         1, 1, 0, 0; ...
         0, 0, 1, 1];
bineqlp = [10; 1; 1]; clp = [16; 10; 0; 4];

ctype = repmat("U", 1, 3); vartype = repmat("C", 1, 4);
[x, fval] = glpk(clp, Aineqlp, bineqlp, lb, ub, ctype, vartype, -1);
fprintf("Optimal objective function value of LP relaxation: %f", fval);


%%%%% solve the knapsack problem with an IP solver here
numiter = 15;
lambdak = 0; gammak = 0; alphak = 1.0;
zlkvals = zeros(numiter,1);
Aineq = [8, 2, 1, 4]; bineq = 10;
lb = zeros(4, 1); ub = ones(4, 1);
ctype = "U"; vartype = repmat("I", 1, 4);
for iter=1:numiter
    c = [16 - lambdak; 10 - lambdak; -gammak; 4 - gammak];
    [xk, ZLk] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
    b_Ax = 1 - [1, 1, 0, 0]*xk;
    lambdak = max(0, lambdak - alphak*b_Ax);
    b_Ax = 1 - [0, 0, 1, 1]*xk;
    gammak = max(0, gammak - alphak*b_Ax);
    zlkvals(iter) = ZLk + lambdak + gammak;
end
plot(zlkvals, 'b*-'); xlabel('Iterations (k)'); ylabel('Z_L(lambda^{(k)})'); 
grid on;
fprintf('The last Lagrangian solution = '); disp(xk');

