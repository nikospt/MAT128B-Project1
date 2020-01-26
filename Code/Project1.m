% MAT 128B: Project 1
% UC Davis Winter 2020
% Nikos Trembois, Caitlin Brown, and

%% Part 1: Fractals
phi = @(z) z^2;
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

colormap([1 0 0; 1 1 1]);
image( [-1 1], [-1 1], M)
axis xy
axis('equal')

%% Part 2:  Fractals
clear M
phi = @(z,c) z^2 - c;
a = linspace(-1,1,100);
b = linspace(-1,1,100);
c = [0.36 + 0.1i, -.123 - .745i];
for k = 1:length(c)
    M{k} = ones(length(a),length(b));
    for r = 1:length(a)
        for i = 1:length(b)
            clear z;
            z = a(r) + 1i*b(i);
            for j = 1:100
                z(j+1) = phi(z(j),c(k));
                if abs(z(j+1)) > 2
                    M{k}(r,i) = 2;
                    break;
                end
            end
        end
    end
end

for i = 1:length(c)
    figure(); hold on
    colormap([1 0 0; 1 1 1]);
    image( [-1 1], [-1 1], M{i})
    axis xy
    axis([ -1 1 -1 1])
    hold off
end

%% Part 3: Julia Sets
