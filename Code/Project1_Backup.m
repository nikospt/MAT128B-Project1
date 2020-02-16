% MAT 128B: Project 1
% UC Davis Winter 2020
% Nikos Trembois, Caitlin Brown, and Shuai Zhi

global c xRange yRange pts bsave
c = [0.36 + 0.1i, -.123 - .745i,-.749, -1.25];
xRange = [1, 1, 1, 1];
yRange = [1, 1, 1, 1];
pts = 100;
bsave = 0;

%% Part 1: Fractals
clearvars -except c xRange yRange pts bsave
phi = @(z,c) z^2;
M = FilledJuliaSet(phi,1,1,100,0);

figure(); hold on
title('Filled Julia Set of $\phi = z^2$','Fontsize',16,'Interpreter','Latex')
xlabel('\Re','Fontsize',18)
ylabel('\Im','Fontsize',18)
colormap([1 0 0; 1 1 1]);
image( [-1 1], [-1 1], M')
axis xy
axis equal
axis([ -1 1 -1 1])
hold off
if bsave == 1
    saveas(gcf,'../Figures/UnitDisk.png')
end

%% Part 2:  Fractals
clearvars -except c xRange yRange pts bsave
phi = @(z,c) z^2 + c;
xrange = 2; yrange = 2;
M = cell(length(c),1);
for i = 1:length(c)
    M{i} = FilledJuliaSet(phi,xrange,yrange,100,c(i));
end

for i = 1:length(c)
    figure(); hold on
    if (real(c(i) > 0 ))
        stitle = strcat('Filled Julia Set of $z^2 + ',num2str(c(i)),'$');
    else
        stitle = strcat('Filled Julia Set of $z^2 ',num2str(c(i)),'$');
    end
    title(stitle,'Interpreter','Latex','FontSize',24)
    colormap([1 0 0; 1 1 1]);
    image( [-xrange xrange], [-yrange yrange], M{i}')
    axis xy
    axis equal
    ax = gca;
    plot(ax.XLim,[0,0],'LineStyle','--','Color',[.5,.5,.5])
    plot([0,0],ax.YLim,'LineStyle','--','Color',[.5,.5,.5])
    xlabel('\Re','Fontsize',18)
    ylabel('\Im','Fontsize',18)
    hold off
    if bsave == 1
        ssave = strcat('../Figures/FilledJulia',num2str(i),'.png');
        saveas(gcf,ssave)
    end
end

%% Part 3: Julia Sets
clearvars -except c xRange yRange pts bsave
psi = @(z,c) sqrt(z - c);
x = zeros(100,4);
y = zeros(100,4);
for k = 1:length(c)
    clear z;
    z = c(k);
    for j = 1:1000
        x(j,k) = real(z(j));
        y(j,k) = imag(z(j));
        if randi([0 1],1,1) == 1
            z(j+1) = psi(z(j),c(k));
        else
            z(j+1) = -psi(z(j),c(k));
        end
    end
end

for i = 1:length(c)
    figure(); hold on
    if (real(c(i) > 0 ))
        stitle = strcat('Julia Set of $z^2 + ',num2str(c(i)),'$');
    else
        stitle = strcat('Julia Set of $z^2 ',num2str(c(i)),'$');
    end
    title(stitle,'Interpreter','Latex','FontSize',24) 
    scatter(x(:,i),y(:,i),'filled')
    ax = gca;
    plot(ax.XLim,[0,0],'LineStyle','--','Color',[.5,.5,.5])
    plot([0,0],ax.YLim,'LineStyle','--','Color',[.5,.5,.5])
    xlabel('\Re','Fontsize',18)
    ylabel('\Im','Fontsize',18)
    hold off
    if bsave == 1
        ssave = strcat('../Figures/Julia',num2str(i),'.png');
        saveas(gcf,ssave)
    end
end

%% Part 4: Computing the Fractal Dimension

% FractalDimension(M{end},.02)

%% Part 5: Connectivity of the Julia Set
clearvars -except c xRange yRange pts bsave
psi = @(z) z^2 + 3;
z = 0;
max_iter = 1000;
for i = 1:max_iter
    z = psi(z);
    if abs(z) > 100
        fprintf('The orbit diverged after %i iterations, the set is not connected\n',i)
        break
    end
end
if abs(z) < 100
    fprintf('The set did not diverge after %i iterations\n',max_iter)
    fprintf('It is reasonable to assume the Julia set is connected\n')
end
    

%% Part 6: Coloring Divergent Orbits
clearvars -except c xRange yRange pts bsave
phi = @(z,c) z^2 - c;
rl = -1; ru = -rl; %1.6
il = -1;  iu = -il; %.7
a = linspace(rl,ru,100);
b = linspace(il,iu,100);
M = cell(length(c),1);
for k = 1:length(c)
    M{k} = ones(length(a),length(b));
    for r = 1:length(a)
        for i = 1:length(b)
            clear z;
            z = a(r) + 1i*b(i);
            for j = 1:100
                z(j+1) = phi(z(j),c(k));
                if abs(z(j+1)) > 100
                    M{k}(r,i) = j;
                    break;
                end
            end
        end
    end
end

for i = 1:length(c)
    figure(); hold on
    image( [rl ru], [il iu], M{i}')
    axis xy
    axis equal
    ax = gca;
    ax.XLim = [rl,-rl]; ax.YLim = [il,-il];
    plot(ax.XLim,[0,0],'LineStyle','--','Color',[.5,.5,.5])
    plot([0,0],ax.YLim,'LineStyle','--','Color',[.5,.5,.5])
    xlabel('\Re','Fontsize',18)
    ylabel('\Im','Fontsize',18)
    colormap(jet(max(max(M{i}))))
    colorbar
    hold off
    if bsave == 1
        ssave = strcat('../Figures/ColoredJulia',num2str(i),'.png');
        saveas(gcf,ssave)
    end
end

%% Part 7: Newton's Method in the Complex Plane
clearvars -except c xRange yRange pts bsave
g = @(z,n) (z^n - 1)/(n*z^(n-1));
a = linspace(-5,5,500);
b = linspace(-5,5,500);
nmax = 5;
M = cell(nmax,1);
for n = 2:nmax
    M{n-1} = 100*ones(length(a),length(b));
    gn = @(z) (z^n - 1)/(n*z^(n-1));
    for r = 1:length(a)
        for i = 1:length(b)
            z = a(r) + 1i*b(i);
            for j = 1:100
                diff = z^n - 1;
                if abs(z^n-1) > 0.001
                    z = z - gn(z);
                else
                    if j < 3
                        fprintf('For n = %i, The root is near %2.4f, %2.4f\n',n,a(r),b(i));
                    end
                    M{n-1}(r,i) = j;
                    break;
                end
            end
        end
    end
end

stitle1 = 'Fixed points of';
for i = 1:nmax-1
    figure(); hold on
    stitle2 = strcat( ' $z^',num2str(i+1),' - 1$');
    image( [min(a) max(a)], [min(b) max(b)], M{i}')
    colormap(jet(max(max(M{i}))))
    colorbar
    axis xy
    ax = gca;
    axis equal
    ax.XLim = [min(a) max(a)]; ax.YLim = [min(b) max(b)];
    title(strcat(stitle1,stitle2),'Interpreter','Latex')
    xlabel('\Re','Fontsize',18)
    ylabel('\Im','Fontsize',18)
    hold off
    if bsave == 1
        ssave = strcat('../Figures/Newton',num2str(i),'.png');
        saveas(gcf,ssave)
    end
end
    

%% Part 8: Mandelbrot Set
clearvars -except c xRange yRange pts bsave
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
            if abs(z(j+1)) > 100
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
image( [-1 1], [-1 1], M')
colorbar
axis xy
axis('equal')
axis([ -1 1 -1 1])
if bsave == 1
    saveas(gcf,'../Figures/Madelbrot.png')
end
%% Part 8 (cont) Zoom in
% zoom in on a fractal by changing limits
clear M;
phi = @(z,c) z^2 + c;
a = linspace(-.3,0.05,500);
b = linspace(0.6,1,500);
M = ones(length(a),length(b));

for r = 1:length(a)
    for i = 1:length(b)
        clear z;
        z = 0;
        c = a(r) + b(i)*1i;
        for j = 1:100
            z(j+1) = phi(z(j),c);
            if abs(z(j+1)) > 100
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
image( [min(a) max(a)], [min(b) max(b)], M')
colorbar
axis xy
axis('equal')

%% Functions
function [result] = R(x,y)
    result = sqrt(x^2 + y^2);
end

function [result] = T(x,y)
    if x > 0
        t = atan(y/x);
    elseif x == 0
        if y > 0
            t = pi/2;
        else
            t = 3*pi/2;
        end
    else
        t = atan(y/x) + pi;
    end
    result = t;
end

function FractalDimension(M,r)
    n = length(M);
    N = 0; % Initialize value of N for box-counting
    for i = 1:n
        for j = 1:n
            if M(i,j) == 1
                N = N+1;
            end
        end
    end
    d = log(N)/log(1/r);
    fprintf('The fractal dimension is %5.4f\n', d);
end

% This function determines if a point is in a Julia Set
function [ M ] = FilledJuliaSet(phi, xrange, yrange, pts, c)
    a = linspace(-xrange, xrange, pts);
    b = linspace(-yrange, yrange, pts);
    for k = 1:length(c)
        M = ones(length(a), length(b));
        for r = 1:length(a)
            for i = 1:length(b)
                clear z;
                z = a(r) + 1i*b(i);
                for j = 1:100
                    z(j+1) = phi( z(j), c(k) );
                    if abs( z(j+1) ) > 2
                        M(r,i) = 2;
                        break;
                    end
                end
            end
        end
    end
end

% This function determines if a point is in a Julia Set
function [ M ] = FilledJuliaSet2(phi, xrange, yrange, pts, c, colored, maxvalue)
    a = linspace(-xrange, xrange, pts);
    b = linspace(-yrange, yrange, pts);
    for k = 1:length(c)
        M = ones(length(a), length(b));
        for r = 1:length(a)
            for i = 1:length(b)
                clear z;
                z = a(r) + 1i*b(i);
                for j = 1:100
                    z(j+1) = phi( z(j), c(k) );
                    if abs( z(j+1) ) > maxvalue
                        if (colored == 1)
                            % If the colored option is true, color the
                            % julia set by the number of iterations to
                            % converge
                            M(r,i) = j;
                        else
                            % if colored option is used, then give a single
                            % value based on convergence, independent of
                            % the number of iterations required to converge
                            M(r,i) = 2;
                        end
                        break
                    end
                end
            end
        end
    end
end