syms r c n

A = symsum(symsum(1+2*(n-c),r,c+1,n),c,1,n-1);
collect(A,n)