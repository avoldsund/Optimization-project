function SQP

rho = 1; E = 1; 
M = 10; mu = 10^-4;
A_top = 1; A_bottom = 10^-4; 




% 1 Tower

n = 13;
m = 58;
ns = 4;     % number of fixed nodes (z = 0)

% initial values:
A = 10^-3*ones(m,1);
q = ones(m,1);
f_supp = ones(3*ns,1);
lambda = ones(3*n,1);

% tau and eta in algorithm
t = 0.5;
eta = 0.25;


x = [0, 1]; y = [0, 1]; z = [0, 1, 2];
% Coordinates v:
v = zeros(3,n);    
v = generate_v(v,x,y,z);
v(:,end) = [0.5; 0.5; 2.5];

% Length of bars [first node, second node, lenght]
l = generate_l(v, m, n);

% Tau, vector of the bar's directional cosines:
tau = generate_tau(m,v,l);

% B-matrix:
B = generate_B(n, m, tau, l);
 
% Force on each node
f_ext = spalloc(3 * (n - ns), 1, 0);
% only non-zero ext load: top node [0.5, 0.5, 2.5] which is node n
f_ext(end-2 : end) = [0; 0; -1];

I_supp = [speye(3 * ns); spalloc(3 * (n - ns), 3 * ns, 0)];
I_ext = [spalloc(3 * ns, 3 * (n - ns), 0); speye(3 * (n - ns))];




% D-matrix as function of A
D = @(A) spdiags(l(:,3)./(E*A),0,m,m);

% Last part of function f
A_sum = @(A) sum(log((A-A_bottom).*(A_top-A)));

% Function evaluated at q, A:
f = @(A,q) 0.5 * q'*D(A)*q - mu*log(M - rho*l(:,3)'*A) - mu*A_sum(A);


% Calculation of gradient of f, dim 2m + 3ns:
    function gradient_f = grad_f(A,q)
        grad_f_A = -0.5*D(A)*(q.^2./A.^2) - (mu/(M-sum(rho*l(:,3).*A)))*rho*l(:,3)-mu*(2*A+A_bottom+A_top)./((A-A_bottom).*(A_top-A));
        grad_f_q = D(A)*q;
        grad_f_fsupp = zeros(3*ns,1);
        gradient_f = [grad_f_A; grad_f_q; grad_f_fsupp];
    end

% Constraint
c = @(q, f_supp) B*q - I_supp*f_supp - I_ext*f_ext;

% Jacobian matrix of the constraints (called A(x) in algorithm):
grad_c = [spalloc(m, 3*n, 0); B'; -I_supp'];

% Calculation of Hessian-components
    function hessian = hess_L(A, q)
        hess_L_AA = 1/4 *spdiags(D(A)*(q.^2./A.^3),0,m,m) - (mu/(M-sum(rho*l(:,3).*A))^2)*rho^2*l(:,3)*l(:,3)' + mu*(spdiags((1./(A-A_bottom).^2)+(1./(A_top-A).^2),0,m,m));
        hess_L_qA = -spdiags(D(A)*(q/A.^2),0,m,m);
        hess_L_fsuppA = zeros(m,3*ns);
        
        hess_L_Aq = hess_L_qA';
        hess_L_qq = D(A);
        hess_L_fsuppq = zeros(m,3*ns);
        
        hess_L_Afsupp = hess_L_fsuppA';
        hess_L_qfsupp = hess_L_fsuppq';
        hess_L_fsuppfsupp = zeros(3*ns,3*ns);
        
        hessian = [hess_L_AA hess_L_qA hess_L_fsuppA;
            hess_L_Aq hess_L_qq hess_L_fsuppq;
            hess_L_Afsupp hess_L_qfsupp hess_L_fsuppfsupp];
    end

% -----------------------------------------------
% ALGORITHM


% (18.9) newton: matrix, s: right-hand side, sol: solution
sol = solve_newton(A, q, f_supp);

% step length and lagrange multiplier from (18.9):
p = sol(1:2*m+3*ns);
lambda_hat = sol(2*m+3*ns+1:end);

p_lambda = lambda_hat - lambda;


% Penalty Parameter (mu in book)
pen = 0.5;
% Parameter in (0,1)
par = 0.5;

% Sigma, determines if Hessian of L is pos. def.
if all(eigs(hess_L(A,q))) > 0
    sigma = 1;
else
    sigma = 0;
end

while pen*1.01 < (grad_f(A,q)'*p + sigma/2*p'*hess_L(A,q)*p)/((1-par)*norm(c(q, f_supp),1))
    pen = pen*1.2;
end


merit_D = @(A,q,f_supp,k) f(A,q) + k*abs(c(q,f_supp));
dir_D = @(A,q,f_supp,k,p) grad_f(A,q,f_supp)'*p - k*abs(c(q,f_supp));
alpha = 1;







end

%%
function v = generate_v(v,x,y,z)
it = 1;
for k = 1:length(z)
    for j = 1:length(y)
        for i = 1:length(x)
            v(:,it) = [x(i); y(j); z(k)];
            it = it + 1;
        end
    end
end
end


%%
function l = generate_l(v, m, n)
l = zeros(m,3); 
it = 1;
for i = 1:n
    for j = i:n
        if i ~= j && norm(v(:,i)-v(:,j)) <= sqrt(3)
            l(it,:) = [i, j, norm(v(:,i)-v(:,j))];
            it = it + 1;
        end
    end
end
end

%%
function tau = generate_tau(m,v,l)
tau = zeros(3,m);
for j = 1:m
    tau(:,j) = (v(:,l(j,1))-v(:,l(j,2)))/l(j,3);
end
end

%%
function B = generate_B(n, m, tau, l)

B = spalloc(3*n, m, 0);
for i = 1:n
    for j = 1:m
        for k = 1:3
            if l(j,1) == i
                B(3*(i-1)+k, j) = -tau(k,j);     
            elseif l(j,2) == i
                B(3*(i-1)+k, j) = tau(k,j);
            end
        end
    end
end
end

%%    
function sol = solve_newton(A,q,f_supp)
newton = [hess_L(A,q) -grad_c; grad_c' zeros(3*n)];
s = [-grad_f(A,q); -c(q,fsupp)];
sol = newton(A,q) \ s(A,q,f_supp);
end

