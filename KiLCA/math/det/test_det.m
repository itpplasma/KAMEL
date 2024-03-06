clear all

N = 11;

pi    = 3.141592653589793238462643383279502884197
euler = 0.5772156649015328606065120900824024310422

for k = 1:N,
    for l = 1:N,
        %A(k,l) = 1e8*(sqrt(k*l)+i*l^2/k^2);
        A(k,l) = k^(l+1)-i*l^k;
    end
end

A

D=det(A);

[L,U,P] = lu(A)

D
