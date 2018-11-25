function[] = cholesky()
fileid = fopen(s, 'r');
I = fscanf(fileid, '%f', inf);
fclose(fileid);
m = size(I);
n = I(1);
x = double.empty(n,0);
z = double.empty(n,0); 
A = zeros(n,n+1,'double');
h = 1;
while(h<=n)
    j = 1;
    while(j<=n+1)
        A(h,j) = I((h-1)*n + (h+j));
        j = j + 1;
    end
    h = h + 1;
end
L = zeros(n,n,'double');
u = zeros(n,n,'double');
i = 1;
while(i<=n)
    j = 1;
    while(j<=n)
        if(i<j)
            L(i,j) = 0;
        end
        if(i>j)
            u(i,j) = 0;
        end
        j = j + 1;
    end
    i = i + 1;
end
k = 1;
while(k<=n)
    p = 1;
    tmp1 = A(k,k);
    while(p<k)
        tmp1 = tmp1 - L(k,p)*L(k,p);
        p = p + 1;
    end
    L(k,k) = sqrt(tmp1);
    I = k+1;
    while(I<=n)
        m = 1;
        tmp2 = A(I,k);
        while(m < k)
            tmp2 = tmp2 - L(I,m)*L(k,m);
            m = m +1;
        end
        L(I,k) = tmp2/L(k,k);
        I = I+1;
    end
    k = k + 1;
end
u = trans_pose(L);
z(1) = double(A(1,n+1)/L(1,1));
j = 2;
while(j<=n)
    y = 1;
    while(y<j)
        A(j,n+1) = A(j,n+1) - z(y)*L(j,y);
        y = y + 1;
    end
    z(j) = double(A(j,n+1)/L(j,j));
    j = j + 1;
end
x(n) = double(z(n)/u(n,n));
j = n-1;
while(j>0)
y = n;
while(y>j)
    z(j) = z(j) - x(y)*u(j,y);
    y = y - 1;
end
x(j) = double(z(j)/u(j,j));
j = j - 1;
end
fprintf('\t\tL\n');
disp(L);
fprintf('\n');
fprintf('\t\tU\n');
disp(u);
fprintf('\n');
fprintf('\t\tX\n');
disp(x);
end
function[c] = trans_pose(a,m,n)
c = zeros(n,m,'double');
for i = 1:m
    for j = 1:n
        if(i~=j)
            c(j,i) = a(i,j);
        end
    end
end
end