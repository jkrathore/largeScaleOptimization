%%%%% Column Generation
A = eye(3);
n = length(A);
b = [32; 28; 16];
c = ones(n, 1);

lb = zeros(n, 1);
ub = [];
ctype = repmat("L", 1, 3);
vartype = repmat("C", 1, n);
[x, fval, ~, extra] = glpk(c, A, b, lb, ub, ctype, vartype, 1);
fprintf('Optimal objective function value : %f\n', fval);
fprintf('Optimal solution : ');
disp(x');
fprintf('Dual optimal solution: ');
extra.lambda'

ck = 1;
ak = [4; 2; 1];
cBBinv = extra.lambda';
ckbar = ck - cBBinv*ak

A = [A, ak];
n = length(A);
b = [32; 28; 16];
c = ones(n, 1);
lb = zeros(n,1);
ub = [];
ctype = repmat("L", 1, 3);
vartype = repmat("C", 1, n);
[x, fval, ~, extra] = glpk(c, A, b, lb, ub, ctype, vartype, 1);
fprintf('Optimal objective function value: %f\n', fval);
fprintf('Optimal solution: '); disp(x');

A = eye(3); b = [32; 28; 16];
D = [2, 5, 7]; e = 26;
n = length(A);
c = ones(n, 1);

lb = zeros(n, 1); ub = [];
ctype = repmat("L", 1, 3); vartype = repmat("C", 1, n);
lbks = zeros(3, 1); ubks = [];
ctypeks = "U"; vartypeks = repmat("I", 1, 3);
    
barcs = -100; tol = -0.001;
while (barcs+1 < tol)
    [xopt, fval, ~, extra] = glpk(c, A, b, lb, ub, ctype, vartype, 1);
    fprintf('Optimal objective function value: %f\n', fval);
    d = -extra.lambda';
    [ak, barcs] = glpk(d, D, e, lbks, ubks, ctypeks, vartypeks, 1);
    fprintf('New column''s reduced cost: %f\n', barcs+1);
    if (barcs+1 < tol)
        A = [A, ak];
        n = length(A); c = ones(n, 1);
        vartype = repmat("C", 1, n);
        lb = zeros(n, 1);
    end
end
fprintf('\nOptimal solution: '); disp(xopt(xopt>0)');
fprintf('Selected combinations: \n'); disp(A(:, xopt>0));



%%%%% Column Generation Applied to Dantzig-Wolfe Decomposition

A = [-1, 2]; c = [1; -3]; D = [1, 1]; d = 5;
xp1 = [0; 0]; xp2 = [0; 5];
Asub = [A*xp1, A*xp2; 1, 1]; 
csub = [c'*xp1; c'*xp2];
b = [6; 1]; 
lb = zeros(2, 1); ub = [];
ctype = repmat("U", 1, 2); vartype = repmat("C", 1, 2);
lbps = zeros(2, 1); ubps = [];
ctypeps = "U"; vartypeps = repmat("C", 1, 2);

barcs = -100; tol = -0.001;
while (barcs < tol)
    [xopt, fval, ~, extra] = glpk(csub, Asub, b, lb, ub, ctype, vartype, 1);
    pi = extra.lambda(1); alpha = extra.lambda(2); 
    fprintf('Objective function value: %f\n', fval);
    % Pricing subproblem
    [xp, barcs] = glpk(c - A'*pi, D, d, lbps, ubps, ctypeps, vartypeps);
    barcs = barcs - alpha;
    if (barcs < tol)
        addcol = [A*xp; 1];
        Asub = [Asub, addcol];
        n = length(Asub); csub = [csub; c'*xp];
        vartype = repmat("C", 1, n);
        lb = zeros(n, 1);
    end
end

# Solving the original problem
A = [-1, 2; 1, 1]; c = [1; -3]; b = [6; 5]; 
lb = zeros(2, 1); ub = [];
ctype = repmat("U", 1, 2); vartype = repmat("C", 1, 2);
[xopt, fval, ~, extra] = glpk(c, A, b, lb, ub, ctype, vartype, 1);
fprintf("\nOptimal objective function value of the original problem: %f", fval);

















