# gomory cutting plan , lattice tree cuts, crooked crossed cuts, splits cuts, cove inequilities

# Cutting Plane Methods

c = [3; 4];
A = [2, 5; 2, -2];
b = [15; 5];

lb = zeros(2, 1); ub = [];
ctype = repmat("U", 1, 2); vartype = repmat("C", 1, 2);
[xopt, objval] = glpk(c, A, b, lb, ub, ctype, vartype, -1);
fprintf('\nOptimal objective function value: %f', objval); 
fprintf('\nOptimal solution: '); disp(xopt');


A = [A; 1, 0];
b = [b; 3];
ctype = repmat("U", 1, 3);
[xopt, objval] = glpk(c, A, b, lb, ub, ctype, vartype, -1);
fprintf('\nOptimal objective function value: %f', objval); 
fprintf('\nOptimal solution: '); xopt'


A = [A; 1, 1];
b = [b; 4];
ctype = repmat("U", 1, 4);
[xopt, objval] = glpk(c, A, b, lb, ub, ctype, vartype, -1);
fprintf('\nOptimal objective function value: %f', objval); 
fprintf('\nOptimal solution: '); disp(xopt')

A = [A; 1, 2];
b = [b; 6];
ctype = repmat("U", 1, 5);
[xopt, objval] = glpk(c, A, b, lb, ub, ctype, vartype, -1);
fprintf('\nOptimal objective function value: %f', objval); 
fprintf('\nOptimal solution: '); disp(xopt')