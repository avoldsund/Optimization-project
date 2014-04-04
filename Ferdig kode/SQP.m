function [A, q, f_supp, energy] = SQP()


global rho E M mu A_top A_bottom
global m n ns v l B I_supp I_ext f_ext lambda forc

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
        dist = sqrt(3);
        
        % initial values:
        A = 0.12*ones(m,1);
        q = 0.1*ones(m,1);
        f_supp = 10000*ones(3*ns,1);
        lambda = 10*ones(3*n,1);
        
        
        [l, B] = generate_truss(m, n, v, dist);
        
        I_supp = [eye(3 * ns); zeros(3*(n-ns), 3*ns)];%spalloc(3 * (n - ns), 3 * ns, 0)];
        I_ext = [zeros(3*ns, 3*(n-ns)); eye(3*(n-ns))];%[spalloc(3 * ns, 3 * (n - ns), 0); speye(3 * (n - ns))];
        
        % Force on each node
        f_ext = zeros(3*(n-ns),1); %spalloc(3 * (n - ns), 1, 0);
        % only non-zero ext load: top node [0.5, 0.5, 2.5] which is node n
        f_ext(end-2 : end) = [0; 0; -1];
        
        forc = I_supp*f_supp + I_ext*f_ext;
        
        
    case 2
            
        n = 13;
        m = 58;
        ns = 4;     % number of fixed nodes (z = 0)
        x = [0, 1]; y = [0, 1]; z = [0, 1, 2];
        
        
        % Coordinates v:
        v = zeros(3,n);
        v = generate_v(v,x,y,z);
        v(:,end) = [0.5; 0.5; 2.5];
        dist = sqrt(3);
        
        % initial values:
        A = 0.12*ones(m,1);
        q = 0.1*ones(m,1);
        f_supp = 10000*ones(3*ns,1);
        lambda = 10*ones(3*n,1);
        
        
        [l, B] = generate_truss(m, n, v, dist);
        
        I_supp = [eye(3 * ns); zeros(3*(n-ns), 3*ns)];%spalloc(3 * (n - ns), 3 * ns, 0)];
        I_ext = [zeros(3*ns, 3*(n-ns)); eye(3*(n-ns))];%[spalloc(3 * ns, 3 * (n - ns), 0); speye(3 * (n - ns))];
        
        % Force on each node
        f_ext = zeros(3*(n-ns),1); %spalloc(3 * (n - ns), 1, 0);
        % only non-zero ext load: top node [0.5, 0.5, 2.5] which is node n
        f_ext(end-2 : end) = [0; 0; 1];
        
        forc = I_supp*f_supp + I_ext*f_ext;
        
        
    case 3
        
        n = 20;
        ns = 4;
        m = 94;
        x = [0, 1, 2, 3, 4]; y = [0, 1]; z = [0, 1];
        
        v = zeros(3,n);
        v = generate_v(v,x,y,z);
        dist = sqrt(3);
        
        % initial values:
        A = 0.12*ones(m,1);
        q = 9*ones(m,1);
        f_supp = 100000*ones(3*ns,1);
        lambda = 100*ones(3*n,1);
        
        %fixed nodes: node number 1, 5, 6 and 10
        % fixed = [0,4,0,4;0,0,1,1;0,0,0,0];
        
        [l, B] = generate_truss(m, n, v, dist);
        
        
        I_supp = blkdiag(eye(3), [zeros(9,6); eye(6)], [zeros(9,3); eye(3); zeros(30,3)]);
        I_ext = blkdiag([zeros(3,9); eye(9)], [zeros(6,9); eye(9)], [zeros(3,30); eye(30)]);
        
        f_ext = zeros(3*(n-ns),1);
        for i = 1:n
            if v(3,i) == 0
                f_ext((i-1)*3+3) = -1;
            end
        end
        
        forc = I_supp*f_supp + I_ext*f_ext; 
end
   

[A, q, f_supp] = SQP_lineSearch(A, q, f_supp);
plot_structure(A)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function [A, q, f_supp] = SQP_lineSearch(A,q,f_supp)
  
    global m ns pen par rho l lambda p M A_top A_bottom
    eta = 0.25;
    iter = 0;
    alpha = 1;
    p = zeros(2*m+3*ns,1);

    while convergence(A, q, f_supp, lambda, p) > 1 && iter < 2000 && alpha > 10^-16%% iter < 1000
        
        iter = iter + 1;
        
        % Compute p by solving (18.9)
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
            
%       alpha_store(1) = alpha;
        alpha_new = alpha;
        
        while ~isValid(alpha_new,eta, A, q, f_supp) || (M-sum(rho*l(:,3).*(A+alpha_new*p(1:m)))) < 0 || ~isempty(find(((A+alpha_new*p(1:m))-A_bottom).*(A_top-(A+alpha_new*p(1:m))) < 0))
%           [alpha_new alpha_store] = interpolation(i,alpha_store, A, q, f_supp); 
            alpha_new = alpha_new*0.5;
            if i > 1000
                error('initial conditions outside feasible set')
            end
            i = i + 1;
        end
        alpha = alpha_new;

        A = A + alpha*p(1:m);
        q = q + alpha*p(m+1:2*m);
        f_supp = f_supp + alpha*p(2*m+1:2*m+3*ns);

        lambda = lambda + alpha*p_lambda;
            
            
    end
    
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = D(A)
    global l E m
    d = spdiags(l(:,3)./(E*A),0,m,m);
end

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

function [l, B] = generate_truss(m, n, v, dist)
    % Length of bars [first node, second node, lenght]
    l = generate_l(v, m, n, dist);
    % Tau, vector of the bar's directional cosines:
    tau = generate_tau(m,v,l);
    B = generate_B(n, m, tau, l);
end

function l = generate_l(v, m, n, dist)
l = zeros(m,3); 
it = 1;
for i = 1:n
    for j = i:n
        if i ~= j && norm(v(:,i)-v(:,j)) <= dist
            l(it,:) = [i, j, norm(v(:,i)-v(:,j))];
            it = it + 1;
        end
    end
end
end

function tau = generate_tau(m,v,l)
    tau = zeros(3,m);
    for j = 1:m
        tau(:,j) = (v(:,l(j,1))-v(:,l(j,2)))/l(j,3);
    end
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function tot = A_sum(A)
    global A_bottom A_top m
    tot = 0;
    for j = 1:m
        tot = tot + log((A(j)-A_bottom)*(A_top-A(j)));
    end
end

function func = f(A,q)
    global mu M rho l 
    func = 0.5 * q'*D(A)*q - mu*log(M - rho*l(:,3)'*A) - mu*A_sum(A);
end

function gradient_f = grad_f(A,q)
    global mu M rho l A_bottom A_top ns
    grad_f_A = -0.5*D(A)*(q.^2./A.^2) - (mu/(M-sum(rho*l(:,3).*A)))*rho*l(:,3)-mu*(2*A+A_bottom+A_top)./((A-A_bottom).*(A_top-A));
    grad_f_q = D(A)*q;
    grad_f_fsupp = zeros(3*ns,1);
    gradient_f = [grad_f_A; grad_f_q; grad_f_fsupp];
end

function constraints = c(q, f_supp)
    global B I_supp I_ext f_ext
    constraints = B*q - I_supp*f_supp - I_ext*f_ext;
end

function gradient = grad_c
    global B I_supp m n
    gradient = [zeros(m, 3*n); B'; -I_supp']';
end

function hessian = hess_L(A, q)
    global m mu M rho l A_bottom A_top ns 
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
    global n
    newton = [hess_L(A,q), -grad_c'; grad_c, zeros(3*n)];
    s = [-grad_f(A,q); -c(q,f_supp)];
    sol = newton \ s;
end

function res = isValid(alpha, eta, A, q, f_supp)
    global m ns p
    if merit(A + alpha*p(1:m),q + alpha*p(m+1:2*m),f_supp + alpha*p(2*m+1:2*m+3*ns)) > merit(A,q,f_supp) + eta*alpha*merit_D(A,q,f_supp)
        res = 0;
    else res = 1;
    end
end

function tol = convergence(A, q, f_supp, lambda, p)
    % KKT conditions must be close to zero for convergence
    tol1 = hess_L(A,q)*p + grad_f(A,q) - grad_c'*lambda;
    tol2 = grad_c*p + c(q, f_supp);
    tol = norm(tol1) + norm(tol2);
end

function out = merit_D(A,q,f_supp)
    global pen p
    out = grad_f(A,q)'*p - pen*norm(c(q,f_supp),1);
end

function out = merit(A,q,f_supp)
    global pen
    out = f(A,q) + pen*norm(c(q,f_supp),1);
end

function [alpha_new alpha_store] = interpolation(i,alpha_store, A, q, f_supp)
    switch i
        case 1
            alpha_new = quadInterpolation(alpha_store(1), p, A, q, f_supp);
            alpha_store(2) = alpha_new;
        otherwise 
            alpha_new = cubicInterpolation(alpha_store(i-1),alpha_store(i), p, A, q, f_supp);
            alpha_store(i+1) = alpha_new;
    end
end

function alpha_new = quadInterpolation(alpha, p, A, q, f_supp)
    global m ns 
    alpha_new = - merit_D(A,q,f_supp)*alpha^2/(2*(merit(A + alpha*p(1:m),q + alpha*p(m+1:2*m),f_supp + alpha*p(2*m+1:2*m+3*ns))...
        - merit(A,q,f_supp) - merit_D(A,q,f_supp)*alpha));
end

function alpha_new = cubicInterpolation(alpha_0,alpha_1, p, A, q, f_supp)
    global m ns 
    coeff = 1/(alpha_0^2*alpha_1^2*(alpha_1-alpha_0))*[alpha_0^2 -alpha_1^2; -alpha_0^3 alpha_1^3]*...
        [(merit(A + alpha_1*p(1:m),q + alpha_1*p(m+1:2*m),f_supp + alpha_1*p(2*m+1:2*m+3*ns)) - merit(A,q,f_supp) ...
        - merit_D(A,q,f_supp)*alpha_1); (merit(A + alpha_0*p(1:m),q + alpha_0*p(m+1:2*m),f_supp + alpha_0*p(2*m+1:2*m+3*ns)) ...
        - merit(A,q,f_supp) - merit_D(A,q,f_supp)*alpha_0)];
    alpha_new = (-coeff(2) + sqrt(coeff(2)^2 - 3*coeff(1)*merit_D(A,q,f_supp)))/(3*coeff(1));
end

function plot_structure(A)
    global l v m n forc

figure(1)
hold on

thickness = sqrt(A)/pi;
for j = 1:m
    pts = [v(:,l(j,1))'; v(:,l(j,2))'];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth',30*thickness(j));
end
hold on

% for i = 1:n
%     x = [v(1,i); v(1,i) + forc((i-1)*3 + 1)];
%     y = [v(2,i); v(2,i) + forc((i-1)*3 + 2)];
%     z = [v(3,i); v(3,i) + forc((i-1)*3 + 3)];
%     plot3(x,y,z, 'r');
% end


end



















