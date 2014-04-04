function plot_structure(A, l, v)

m = length(A);
[~, n] = size(v);

figure
hold on

for i = 1:n
    plot(v(1,i), v(2,i), v(3,i));
end
