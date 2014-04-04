function [A, q, f_supp] = SQP()

rho = 1; E = 1; 
M = 10; mu = 10^-4;
A_top = 1; A_bottom = 10^-4; 

problem = input('1 = Tower, 2 = Cantilever, 3 = Bridge\n');

switch problem
    case 1

    n = 13;
    m = 58;
    ns = 4;     % number of fixed nodes (z = 0)
    x = [0, 1]; y = [0, 1]; z = [0, 1, 2];    
    
    
    % Coordinates v:
    v = zeros(3,n);    
    v = generate_v(v,x,y,z);
    v(:,end) = [0.5; 0.5; 2.5];
    
    % initial values:
    A = 0.12*ones(m,1);
    q = 0.1*ones(m,1);
    f_supp = 10000*ones(3*ns,1);
    lambda = 10*ones(3*n,1);


    [l, B, I_supp, I_ext] = generate_truss(m, n, ns, v);


    % Force on each node
    f_ext = zeros(3*(n-ns),1); %spalloc(3 * (n - ns), 1, 0);
    % only non-zero ext load: top node [0.5, 0.5, 2.5] which is node n
    f_ext(end-2 : end) = [0; 0; -1];

end

%% LOCAL FUNCTIONS

        % D-matrix as function of A
        D = @(A) spdiags(l(:,3)./(E*A),0,m,m);
        
        % Last part of function f
        function tot = A_sum(A)
            tot = 0;
            for j = 1:m
                tot = tot + log((A(j)-A_bottom)*(A_top-A(j)));
            end
        end

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
    grad_c = [zeros(m, 3*n); B'; -I_supp']';

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

        function sol = solve_newton(A,q,f_supp)
            newton = [hess_L(A,q), -grad_c'; grad_c, zeros(3*n)];
            s = [-grad_f(A,q); -c(q,f_supp)];
            sol = newton \ s;
        end

    
        function res = isValid(alpha,eta)
            
            if merit(A + alpha*p(1:m),q + alpha*p(m+1:2*m),f_supp + alpha*p(2*m+1:2*m+3*ns)) > merit(A,q,f_supp) + eta*alpha*merit_D(A,q,f_supp)
                res = 0;
            else res = 1;
            end
        end
    

        % Function for norm of KKT conditions
        function tol = convergence(A, q, f_supp)

            % KKT conditions must be close to zero for convergence
            tol1 = hess_L(A,q)*p + grad_f(A,q) - grad_c'*lambda;
            tol2 = grad_c*p + c(q, f_supp);
            tol = norm(tol1) + norm(tol2);
        end

    
        function out = merit_D(A,q,f_supp)
%             disp('grad_f')
%             grad_f(A,q)'*p
            out = grad_f(A,q)'*p - pen*norm(c(q,f_supp),1);
        end

        function out = merit(A,q,f_supp)
            out = f(A,q) + pen*norm(c(q,f_supp),1);
        end
    
        function [alpha_new alpha_store] = interpolation(i,alpha_store)
            switch i
                case 1
                    alpha_new = quadInterpolation(alpha_store(1));
                    alpha_store(2) = alpha_new;
                otherwise 
%                     alpha_new = cubicInterpolation(alpha_store(i-1),alpha_store(i));
                    alpha_new = alpha_store(i)*0.01;
                    alpha_store(i+1) = alpha_new;
            end
        end



        function alpha_new = quadInterpolation(alpha)
            alpha_new = - merit_D(A,q,f_supp)*alpha^2/(2*(merit(A + alpha*p(1:m),q + alpha*p(m+1:2*m),f_supp + alpha*p(2*m+1:2*m+3*ns))...
                - merit(A,q,f_supp) - merit_D(A,q,f_supp)*alpha));
        end

        function alpha_new = cubicInterpolation(alpha_0,alpha_1)
            coeff = 1/(alpha_0^2*alpha_1^2*(alpha_1-alpha_0))*[alpha_0^2 -alpha_1^2; -alpha_0^3 alpha_1^3]*...
                [(merit(A + alpha_1*p(1:m),q + alpha_1*p(m+1:2*m),f_supp + alpha_1*p(2*m+1:2*m+3*ns)) - merit(A,q,f_supp) ...
                - merit_D(A,q,f_supp)*alpha_1); (merit(A + alpha_0*p(1:m),q + alpha_0*p(m+1:2*m),f_supp + alpha_0*p(2*m+1:2*m+3*ns)) ...
                - merit(A,q,f_supp) - merit_D(A,q,f_supp)*alpha_0)];
            alpha_new = (-coeff(2) + sqrt(coeff(2)^2 - 3*coeff(1)*merit_D(A,q,f_supp)))/(3*coeff(1));
        end


    %%

    % ALGORITHM

    % tau and eta in algorithm
%     t = 0.5;
    eta = 0.25;
    
    p = zeros(2*m+3*ns,1);
    iter = 0;

    while convergence(A, q, f_supp) > 5 %% iter < 1000
        convergence(A, q, f_supp);
        iter = iter + 1
        % Compute p by solving (18.9)
        % newton: matrix, s: right-hand side, sol: solution
        sol = solve_newton(A, q, f_supp);
        
        p = sol(1:2*m+3*ns);
        lambda_hat = sol(2*m+3*ns+1:end);

        p_lambda = lambda_hat - lambda;

        % Choose penalty parameter to satisfy (18.36) with sigma evaluated:
        pen = 0.5;
        % Parameter in (0,1)
        par = 0.5;

        % Sigma, determines if Hessian of L is pos. def.
        [~,posDef] = chol(hess_L(A,q));
        if posDef == 0
            sigma = 1;
        else
            sigma = 0;
        end

        % Penalty parameter must satisfy (18.36):
        while pen*1.1 < (grad_f(A,q)'*p + sigma/2*p'*hess_L(A,q)*p)/((1-par)*norm(c(q, f_supp),1))
            pen = pen*1.1;
        end

        alpha = 1;
        
        i = 1;
            
%         alpha_store(1) = alpha;
        alpha_new = alpha;
        
        while ~isValid(alpha_new,eta) || (M-sum(rho*l(:,3).*(A+alpha_new*p(1:m)))) < 0 || ~isempty(find(((A+alpha_new*p(1:m))-A_bottom).*(A_top-(A+alpha_new*p(1:m))) < 0))
%             [alpha_new alpha_store] = interpolation(i,alpha_store); 
            alpha_new = alpha_new*0.5;
            if i > 1000
                error('initial conditions outside feasible set')
            end
            i = i + 1;
        end
        
        
        alpha = alpha_new


            A = A + alpha*p(1:m);
            q = q + alpha*p(m+1:2*m);
            f_supp = f_supp + alpha*p(2*m+1:2*m+3*ns);
            
            lambda = lambda + alpha*p_lambda;
            
            
    end
        
    plot_structure(A,l,v);

end
%%
    
    


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
function [l, B, I_supp, I_ext] = generate_truss(m, n, ns, v)
    % Length of bars [first node, second node, lenght]
    l = generate_l(v, m, n);

    % Tau, vector of the bar's directional cosines:
    tau = generate_tau(m,v,l);

    % B-matrix:
    B = generate_B(n, m, tau, l);

    I_supp = [eye(3 * ns); zeros(3*(n-ns), 3*ns)];%spalloc(3 * (n - ns), 3 * ns, 0)];
    I_ext = [zeros(3*ns, 3*(n-ns)); eye(3*(n-ns))];%[spalloc(3 * ns, 3 * (n - ns), 0); speye(3 * (n - ns))];
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

function plot_structure(A, l, v)

m = length(A);
[~, n] = size(v);

figure
hold on

for i = 1:n
    plot3(v(1,i), v(2,i), v(3,i), 'bo');
end


thickness = A/pi;
for j = 1:m
    points = 
    plot3([
    
end



end



plot3(pts(:,1), pts(:,2), pts(:,3))


