fileid = fopen('c:\\t.txt', 'r');
I = fscanf(fileid, '%f', inf);
fclose(fileid);
m = size(I);
n = I(1);

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
    if(tmp~=1)
        A(r,:) = A(r,:)/tmp;
    end
    tmp1 = r + 1;
    while(tmp1<=n)
        if(A(tmp1,r)~=0)
            A(tmp1,:) = A(tmp1,:) - A(tmp1,r)*A(r,:);
        end
        tmp1 = tmp1 + 1;
    end
    r = r + 1;
end
    disp(A);
    