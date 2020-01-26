M = importdata('matlab.dat');
colormap([1 0 0; 1 1 1]);
%image( [-1 1], [-1 1], M)
image( [-1.8 1.8], [-.7 .7], M)
axis xy
axis('equal')