% MAT 128B: Project 1
% UC Davis Winter 2020
% Nikos Trembois, Caitlin Brown, and

%% Fractals
phi = @(z) z^2;
fp1 = 0;
fp2 = 0;
a = linspace(-1,1,500);
b = linspace(-1,1,500);
M = ones(length(a),length(b));

for r = 1:length(a)
    for i = 1:length(b)
        clear z;
        z = a(r) + 1i*b(i);
        for j = 1:100
            z(j+1) = phi(z(j));
            if abs(z(j+1)) > 2
                M(r,i) = 2;
                break;
            end
        end
    end
end

% Here is a change

colormap([1 0 0; 1 1 1]);
image( [-1 1], [-1 1], M)
axis xy
axis('equal')

%% Fractals
phi = @(z) z^2 - 1.25;
fp1 = (1+sqrt(6))/2;
fp2 = (1-sqrt(6))/2;
a = linspace(-1.8,1.8,500);
b = linspace(-.7,.7,500);
M = ones(length(a),length(b));

for r = 1:length(a)
    for i = 1:length(b)
        clear z;
        z = a(r) + 1i*b(i);
        for j = 1:100
            z(j+1) = phi(z(j));
            if abs(z(j+1)) > 2
                M(r,i) = 2;
                break;
            end
        end
    end
end

colormap([1 0 0; 1 1 1]);
image( [-1.8 1.8], [-.7 .7], M)
axis xy
axis('equal')
