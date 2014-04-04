  function [A, q, f_supp] = SQP_lineSearch(A,q,f_supp)
  
    global m ns pen par rho l 
    eta = 0.25;
    iter = 0;
    alpha = 1;

    while convergence(A, q, f_supp, lambda) > 1 && iter < 2000 && alpha > 10^-16%% iter < 1000
        
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
        
        while ~isValid(alpha_new,eta, p, f_supp) || (M-sum(rho*l(:,3).*(A+alpha_new*p(1:m)))) < 0 || ~isempty(find(((A+alpha_new*p(1:m))-A_bottom).*(A_top-(A+alpha_new*p(1:m))) < 0))
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