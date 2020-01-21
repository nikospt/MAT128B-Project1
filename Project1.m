% MAT 128B: Project 1
% UC Davis Winter 2020
% Nikos Trembois and

%% Fractals
fr = @(a,b) a^2 - b^2;
fi = @(a,b) 2*a*b;
dfr = @(a) 2*a;
dfi = @(b) 2*b;
real = linspace(-1,1,100);
imag = linspace(-1,1,100);
for i = 1:length(real)
    for j = 1:length(imag)
        [Re{i,j}, Im{i,j}] = NewtonIteration(real(i),real(j),fr,fi,dfr,dfi,100,1e-10);
    end
end

x = Re{30,70}; y = Im{30,70};
plot(x,y)

%% 
phi = @(z) z^2;
fp1 = 0;
fp2 = 0;
a = linspace(-1,1,100);
b = linspace(-1,1,100);
M = ones(100,100);

for r = 1:100
    for i = 1:100
        pass = 0;
        clear z;
        z = a(r) + 1i*b(i);
        for j = 1:100
            z(j+1) = phi(z(j));
%             if z(j+1) < 1e-6
%                 pass = pass + 1;
%                 if pass >= 5
%                     fprintf('The initial guess %.2f converges to %.2f in %i iterations\n',z(1),z(end),i-5)
%                     break;
%                 end
%             else
%                 pass = 0;
%             end
            if abs(z(j+1)) > 2
                M(r,i) = 2;
                break;
            end
        end
    end
end

colormap([1 0 0; 1 1 1]);
image( [-1 1], [-1 1], M)
axis xy
axis('equal')