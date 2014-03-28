% Project Tandoori


rho = 1; E = 1; 
M = 10; mu = 10^-4;
A_top = 1; A_bottom = 10^-4; 


problem = input('1 = Tower, 2 = Cantilever, 3 = Bridge');

switch problem
    case 1
        n = 13;
        m = 58;
        ns = 4;     % number of fixed nodes (z = 0)

        % initial values:
        A = 0.8*ones(m,1);
        q = ones(m,1);
        f_supp = ones(3*ns,1);

        v = zeros(3,n);    % coordinates of node
        l = zeros(m,3);     % length of bars with index of nodes in v


        % Coordinates of v:
        x = [0, 1]; y = [0, 1]; z = [0, 1, 2];
        it = 1;
        for k = 1:length(z)
            for j = 1:length(y)
                for i = 1:length(x)
                    v(:,it) = [x(i); y(j); z(k)];
                    it = it + 1;
                end
            end
        end
        v(:,end) = [0.5; 0.5; 2.5];



        % Force on each node
        force = zeros(3*n,1);
        f_ext = zeros(3*(n-ns),1);
        % only non-zero ext load: top node [0.5, 0.5, 2.5] which is node n
        f_ext(end-2:end) = [0; 0; -1];



        % Length of bars [first node, second node, lenght]
        it = 1;
        for i = 1:n
            for j = i:n
                if i ~= j && norm(v(:,i)-v(:,j)) <= sqrt(3)
                    l(it,:) = [i, j, norm(v(:,i)-v(:,j))];
                    it = it + 1;
                end
            end
        end


        % Tau, vector of the bar's directional cosines:
        tau = zeros(3,m);
        for j = 1:m
            tau(:,j) = (v(:,l(j,1))-v(:,l(j,2)))/l(j,3);
        end



        % B-matrix:
        B = zeros(3*n, m);
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



        % D-matrix
        
        D = @(A) spdiags(l(j,3)./(E*A),0,m,m);

        A_sum = sum((A-A_bottom).*(A_top-A));

        % Function evaluated at q, A:
        f = @(A,q) 0.5 * q'*D(A)*q - mu*log(M - rho*l(:,3)'*A) - mu*A_sum;

        A_inv = 1./A;


        % Calculation of gradient of f, dim 2m + 3ns:
        grad_f_A = @(A,q) -0.5*D(A)*(q.^2./A.^2) - (mu/(M-sum(rho*l(:,3).*A)))*rho*l(:,3)-mu*(2*A+A_bottom+A_top)./((A-A_bottom).*(A_top-A));
        grad_f_q = @(A,q) D(A)*q;
        grad_f_fsupp = zeros(3*ns,1);

        grad_f = @(A,q) [grad_f_A(A,q); grad_f_q(A,q); grad_f_fsupp];
end