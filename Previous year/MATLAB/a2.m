fileid = fopen(s, 'r');
I = fscanf(fileid, '%f', inf);
fclose(fileid);
m = size(I);
n = I(1);
x = double.empty(n,0); 
A = zeros(n,n+1,'double');
i = 1;
while(i<=n)
    j = 1;
    while(j<=n+1)
        A(i,j) = I((i-1)*n + (i+j));
        j = j + 1;
    end
    i = i + 1;
end
r = 1;
while(r<=n)
    tmp = A(r,r);
    tmp1 = r + 1;
    while(tmp1<=n)
        if(A(tmp1,r)~=0)
            A(tmp1,:) = A(tmp1,:) - A(tmp1,r)*A(r,:)/tmp;
        end
        tmp1 = tmp1 + 1;
    end
    r = r + 1;
end
x(n) = A(n,n+1)/A(n,n);
j = n-1;
while(j>0)
y = n;
while(y>j)
    A(j,n+1) = A(j,n+1) - x(y)*A(j,y);
    y = y - 1;
end
x(j) = A(j,n+1)/A(j,j);
j = j - 1;
end
fprintf('\t\tX\n');
disp(x);


