% Project Tandoori


rho = 1;
E = 1;
A_h = 1;
A_l = 10^-4;
M = 10;
my = 10^-4;


problem = input('1 = Tower, 2 = Cantilever, 3 = Bridge');

switch problem
    case 1
        n = 13;
        m = 58;
        