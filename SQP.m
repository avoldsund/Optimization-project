

rho = 1; E = 1; 
M = 10; mu = 10^-4;
A_top = 1; A_bottom = 10^-4; 




% 1 Tower

n = 13;
m = 58;
ns = 4;     % number of fixed nodes (z = 0)

% initial values:
A = ones(m,1);
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
f_ext = spalloc(3 * (n - ns), 1, 0);
% only non-zero ext load: top node [0.5, 0.5, 2.5] which is node n
f_ext(end-2 : end) = [0; 0; -1];

I_supp = [speye(3 * ns); spalloc(3 * (n - ns), 3 * ns, 0)];
I_ext = [spalloc(3 * ns, 3 * (n - ns), 0); speye(3 * (n - ns))];

force = I_supp * f_supp + I_ext * f_ext;



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
D = zeros(m,m);
for j = 1:m
    D(j,j) = l(j,3)/(E*A(j));
end



A_sum = 0;
for j = 1:m
    A_sum = A_sum + (A(j)-A_bottom)*(A_top-A(j));
end 

% Function evaluated at q, A:
f = 0.5 * q'*D*q - mu*log(M - rho*l(:,3)'*A) - mu*A_sum;

A_inv = zeros(m,1);
for j = 1:m
    A_inv(j) = 1/(A(j));
end


% Calculation of gradient of f, dim 2m + 3ns:
grad_f_A = -0.5*D*q.*q.*A_inv.*A_inv;
grad_f_q = D*q;
grad_f_fsupp = zeros(3*ns,1);

grad_f = [grad_f_A; grad_f_q; grad_f_fsupp];



