function [real, imag] = NewtonIteration(r_0,i_0,fr,fi,dfr,dfi,iters,tol)
%NewtonIteration Newton's Iteration method
%   Newtons iteration method for solving roots

real = r_0;
imag = i_0;
for i = 1:iters
    real(i+1) = real(i) - fr(real(i),imag(i))/dfr(real(i));
    imag(i+1) = imag(i) - fi(real(i),imag(i))/dfi(imag(i));
    if abs(fr(real(i+1),imag(i+1))) < tol && abs(fi(real(i+1),imag(i+1))) < tol
        pass = pass + 1;
        if pass >= 5
            % The requirement is the solution converges within tolerance
            % and stays converged for 5 iterations
            fprintf('The initial guess %.2f converges to %.2f in %i iterations',x(1),x(end),i-5)
            break
        end
    else
        pass = 0;
    end
end
