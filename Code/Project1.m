% MAT 128B: Project 1
% UC Davis Winter 2020
% Nikos Trembois, Caitlin Brown, and Shuai Zhi

%% Part 1: Fractals
phi = @(z) z^2;
a = linspace(-1,1,100);
b = linspace(-1,1,100);
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

figure(); hold on
title('Filled Julia Set of $\phi = z^2$','Interpreter','Latex')
xlabel('Real')
ylabel('Imaginary')
colormap([1 0 0; 1 1 1]);
image( [-1 1], [-1 1], M)
axis xy
axis('equal')
axis([ -1 1 -1 1])
hold off
saveas(gcf,'../Figures/UnitDisk.png')

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
    axis('equal')
    axis([ -1 1 -1 1])
    hold off
end

%% Part 3: Julia Sets

%% Part 8: Mandelbrot Set
phi = @(z,c) z^2 + c;
a = linspace(-1,1,500);
b = linspace(-1,1,500);
M = ones(length(a),length(b));

for r = 1:length(a)
    for i = 1:length(b)
        clear z;
        z = 0;
        c = a(r) + b(i)*1i;
        for j = 1:100
            z(j+1) = phi(z(j),c);
            if abs(z(j+1)) > 2
                M(r,i) = j;
                break;
            end
        end
    end
end

figure(); hold on
title('Mandelbrot Set')
xlabel('Real')
ylabel('Imaginary')
colormap(jet(100))
image( [-1 1], [-1 1], M)
colorbar
axis xy
axis('equal')
axis([ -1 1 -1 1])
%% Zoom in 
% zoom in on a fractal by changing limits
% Weird it seems to change significantly when using higher resolution in
% this section from what is shown in the previous plot
phi = @(z,c) z^2 + c;
a = linspace(-.8,-.4,500);
b = linspace(.2,.5,500);
M = ones(length(a),length(b));

for r = 1:length(a)
    for i = 1:length(b)
        clear z;
        z = 0;
        c = a(r) + b(i)*1i;
        for j = 1:100
            z(j+1) = phi(z(j),c);
            if abs(z(j+1)) > 2
                M(r,i) = j;
                break;
            end
        end
    end
end

figure(); hold on
title('Mandelbrot Set')
xlabel('Real')
ylabel('Imaginary')
colormap(jet(100))
image( [-.7 -.4], [.2 .45], M)
colorbar
axis xy
axis('equal')
axis([-.7 -.4 .2 .45])