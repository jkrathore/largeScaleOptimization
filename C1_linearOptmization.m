%%%%%%% Linear Programming

c = [2; -1; 1; 0; 0; 0]; % Objective function coefficients
Aeq = [3,  1,  1,  1,  0,  0; ...
       1, -1,  2,  0,  1,  0; ...
       1,  1, -1,  0,  0,  1]; % Coefficient matrix for equalities
beq = [6; 1; 2]; % RHS for equalities
lb = zeros(6, 1); ub = [];
ctype = repmat("S", 1, 3); vartype = repmat("C", 1, 6);
[x, fval, errnum] = glpk(c, Aeq, beq, lb, ub, ctype, vartype, -1);
fprintf('Solution Status = %d\n', errnum);
fprintf('Optimal Objective Function Value = %f\n', fval);
fprintf('Optimal Solution = '); disp(x')

%%%%%%% Simplex Method 

bv = [1,2,4];
nbv = [3,5];

B = Aeq(:, bv);
N = Aeq(:, nbv);
cB = c(bv);
cN = c(nbv);
fprintf('Primal feasibility: ');
disp(inv(B)*beq);
fprintf('Dual fesibility: ');
disp(cN' - cB'*inv(B)*N);
                                                                                  
%%%%%% Duality

% First the "primal" problem
c = [5; 11; 8];
Aineq = [2, 3, 3;...
         3, 7, 4;...
         1, 2, 2];
bineq = [7; 11; 5];
lb = zeros(3,1);
ub = [];
ctype = repmat("L", 1, 3);
vartype = repmat("C", 1, 3);
[x, fval] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, 1);
fprintf('Optimal Objective Function Value of The Primal problem = %f\n', fval);
fprintf('Optimal Solution of the Primal Problem = ');
disp(x')

% Next the "dual" problem
ctype = repmat("U", 1, 3);
[y, fval] = glpk(bineq, Aineq, c, lb, ub, ctype, vartype, -1);
fprintf('Optimal Objective Function Value of Dual Problem = %f\n', fval);
fprintf('Optimal Solution of the Dual Problem = ');
disp(y')

%%%%% Mixed Integer Programming

c = [0;1];
Aineq = [-2, 2;...
          2, 2];
bineq = [1; 7];
lb = zeros(2, 1);
ub = [];
ctype = repmat("U", 1, 2);
vartype = repmat("C", 1, 2);


%first LP solution
[x, fval_lp] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Optimal objective function value of LP = %f\n', fval_lp);
fprintf('Optional Solution of LP = ');
disp(x')


vartype = repmat("I", 1, 2);
[x, fval_ip] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Optimal objective function value of IP = %f\n', fval_ip);
fprintf('Optimal Solution of IP = ');
disp(x')

fprintf('Relative Gap = %.0f%%', 100*(fval_lp/fval_ip));

%%%%%%%%%%%%% branch and bound method 

% Root problem
c = [4; -2; 7; -1];
Aineq = [ 1, 0, 5, 0;...
          1, 1, -1, 0;...
          6, -5, 0, 0;...
         -1,  2, 0, -2];
bineq = [10; 1; 0; 3];
lb = zeros(4, 1);
ub = [];
ctype = repmat("U", 1, 4);
vartype = repmat("C", 1, 4);
[x, fval, errnum] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Solution Status = %d\n', errnum);
fprintf('Optimal Objective Function Value of Root Problem = %f\n', fval);
fprintf('Optimal Solution of the Root Problem = ');
disp(x');


% splitting the exiting problem into two subproblems
% x1 <= 1 and x1 >= 2


c = [ 4; -2; 7; -1];
% subproblem 1
Aineq = [ 1, 0, 5, 0;...
          1, 1, -1, 0;...
          6, -5, 0, 0;...
         -1, 0,  2, -2;...
          1, 0,  0, 0];
bineq = [10; 1; 0; 3; 1];
lb = zeros(4,1);
ub = [];
ctype = repmat("U", 1, 5); 
vartype = repmat("C", 1, 4);
[x, fval, errnum] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Solution Status = %d\n', errnum);
fprintf('Optimal Objective Function Value of subproblem 1 = %f\n', fval);
fprintf('Optimal Solution of Subproblem 1 = ');
disp(x')

% subproblem 2
Aineq = [1, 0, 5, 0;...
         1, 1, -1, 0;...
         6, -5, 0, 0;...
        -1, 0, 2, -2;...
        -1, 0, 0, 0];
bineq = [10; 1; 0; 3; -2];
[x, fval, errnum] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Solution Status = %d\n', errnum);
fprintf('Optimal Objective Function Value of subproblem 1 = %f\n', fval);
fprintf('Optimal Solution of Subproblem 1 = ');
disp(x')

% subproblem 2 is infeasible and hence it is fathomed, subproblem 1 is further branched
% subproblem 3 : original problem plus additional constraints x1 <= 1 and x2 <= 1
% subproblem 4 : original problem plus additional constraints x1 <= 1 and x2 <= 2

c = [4; -2; 7; -1];
% Subproblem 3
Aineq = [ 1,  0,  5,  0; ...
          1,  1, -1,  0; ...
          6, -5,  0,  0; ...
         -1,  0,  2, -2; ...
          1,  0,  0,  0; ...
          0,  1,  0,  0]; 
bineq = [10; 1; 0; 3; 1; 1];

lb = zeros(4, 1); ub = [];
ctype = repmat("U", 1, 6);
vartype = repmat("C", 1, 4);
[x, fval, errnum] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Solution Status = %d\n', errnum);
fprintf('Optimal Objective Function Value of Subproblem 3 = %f\n', fval);
fprintf('Optimal Solution of Subproblem 3 = ');
disp(x')

% Subproblem 4
Aineq = [ 1,  0,  5,  0; ...
          1,  1, -1,  0; ...
          6, -5,  0,  0; ...
         -1,  0,  2, -2; ...
          1,  0,  0,  0; ...
          0, -1,  0,  0]; 
bineq = [10; 1; 0; 3; 1; -2];

[x, fval, errnum] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Solution Status = %d\n', errnum);
fprintf('Optimal Objective Function Value of Subproblem 4 = %f\n', fval);
fprintf('Optimal Solution of Subproblem 4 = ');
disp(x')

%Subproblem 5: Original problem plus additional constraints  x1≤1 ,  x2≤1  and  x1≤0  ( x1=0 )
%Subproblem 6: Original problem plus additional constraints  x1≤1 ,  x2≤1  and  x1≥1  ( x1=1 )

c = [4; -2; 7; -1];
% Subproblem 5
Aineq = [ 1,  0,  5,  0; ...
          1,  1, -1,  0; ...
          6, -5,  0,  0; ...
         -1,  0,  2, -2; ...
          1,  0,  0,  0; ...
          0,  1,  0,  0; ...
          1,  0,  0,  0]; 
bineq = [10; 1; 0; 3; 1; 1; 0];

lb = zeros(4, 1); ub = [];
ctype = repmat("U", 1, 7); vartype = repmat("C", 1, 4);
[x, fval, errnum] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Solution Status = %d\n', errnum);
fprintf('Optimal Objective Function Value of Subproblem 5 = %f\n', fval);
fprintf('Optimal Solution of Subproblem 5 = '); disp(x')

% Subproblem 6
Aineq = [ 1,  0,  5,  0; ...
          1,  1, -1,  0; ...
          6, -5,  0,  0; ...
         -1,  0,  2, -2; ...
          1,  0,  0,  0; ...
          0,  1,  0,  0; ...
         -1,  0,  0,  0]; 
bineq = [10; 1; 0; 3; 1; 1; -1];

[x, fval, errnum] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Solution Status = %s\n', errnum);
fprintf('Optimal Objective Function Value of Subproblem 6 = %f\n', fval);
fprintf('Optimal Solution of Subproblem 6 = '); disp(x')



%Let's double-check our result and conclude this lecture.

% Root problem
c = [4; -2; 7; -1];
Aineq = [ 1,  0,  5,  0; ...
          1,  1, -1,  0; ...
          6, -5,  0,  0; ...
         -1,  0,  2, -2]; 
bineq = [10; 1; 0; 3];

lb = zeros(4, 1); ub = [];
ctype = repmat("U", 1, 4); vartype = "IIIC";
[x, fval_ip, errnum] = glpk(c, Aineq, bineq, lb, ub, ctype, vartype, -1);
fprintf('Optimal Objective Function Value of IP = %f\n', fval_ip);
fprintf('Optimal Solution of IP = '); disp(x')





















