function[] = q1_pa2()
[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
n = input(['Choose method.',... 
    '\n1. Gauss Elimination without pivoting',... 
    '\n2. Gauss Elimination with partial pivoting',... 
    '\n3. Doolittle Method',... 
    '\n4. Crout Method',... 
    '\n5. Cholesky Decomposition: ']);
switch n
    case 1
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
        %fprintf('\t\tX\n');
        %disp(x);
        [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        fprintf(fileid1,'Gauss Elimination without pivoting \r\n');
        fprintf(fileid1,'X \r\n');
        formatspec1 = '%0.5f\r\n';
        fprintf(fileid1,formatspec1,x);
    case 2
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
            tmpr = r;
            z = r+1;
            while(z<=n)
                if(abs(A(z,r))>=abs(tmp))
                    tmp = A(z,r);
                    tmpr = z;
                end
                z = z+1;
            end
            if(tmpr~=r)
                q = 1;
                while(q<=n+1)
                    tmp2 = A(r,q);
                    A(r,q) = A(tmpr,q);
                    A(tmpr,q) = tmp2;
                    q = q + 1;
                end
            end
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
        [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        fprintf(fileid1,'Gauss Elimination with partial pivoting \r\n');
        fprintf(fileid1,'X \r\n');
        formatspec1 = '%0.5f\r\n';
        fprintf(fileid1,formatspec1,x);
    case 3
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
                if(i==j)
                    L(i,j) = 1;
                end
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
            J = k;
            while(J<=n)
                m = 1;
                while(m<k)
                    A(k,J) = A(k,J) - L(k,m)*u(m,J);
                    m = m + 1;
                end
                u(k,J) = A(k,J);
                J = J+1;
            end
            I = k + 1;
            while(I<=n)
                m = 1;
                while(m<k)
                    A(I,k) = A(I,k) - L(I,m)*u(m,k);
                    m = m + 1;
                end
                L(I,k) = A(I,k)/u(k,k);
                I = I + 1;
            end
            k = k + 1;
        end
        tmp1 = A(1,n+1);
        tmp2 = L(1,1);
        z(1) = tmp1/tmp2;
        j = 2;
        while(j<=n)
        y = 1;
        while(y<j)
            A(j,n+1) = A(j,n+1) - z(y)*L(j,y);
            y = y + 1;
        end
        z(j) = A(j,n+1)/L(j,j);
        j = j + 1;
        end
        x(n) = z(n)/u(n,n);
        j = n-1;
        while(j>0)
        y = n;
        while(y>j)
            z(j) = z(j) - x(y)*u(j,y);
            y = y - 1;
        end
        x(j) = z(j)/u(j,j);
        j = j - 1;
        end
        [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        fprintf(fileid1,'Doolittle Method \r\n');
        fprintf(fileid1,'X \r\n');
        formatspec1 = '%0.5f\r\n';
        fprintf(fileid1,formatspec1,x);
        fprintf(fileid1,'L \r\n');
        formatspec2 = '%0.5f  ';
        for g = 1:n
        fprintf(fileid1,formatspec2,L(g,:));
        fprintf(fileid1,'\r\n');
        end
        fprintf(fileid1,'U \r\n');
        formatspec2 = '%0.5f  ';
        for g = 1:n
        fprintf(fileid1,formatspec2,u(g,:));
        fprintf(fileid1,'\r\n');
        end
    case 4
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
        l = zeros(n,n,'double');
        u = zeros(n,n,'double');
        i = 1;
        while(i<=n)
            j = 1;
            while(j<=n)
                if(i==j)
                    u(i,j) = 1;
                end
                if(i<j)
                    l(i,j) = 0;
                end
                if(i>j)
                    u(i,j) = 0;
                end
                j = j + 1;
            end
            i = i + 1;
        end
        k = 1;
        while(k<=n);
             I = k;
            while(I<=n)
                m = 1;
                while(m<k)
                    A(I,k) = A(I,k) - l(I,m)*u(m,k);
                    m = m + 1;
                end
                l(I,k) = A(I,k);
                I = I + 1;
            end
             J = k + 1;
            while(J<=n)
                m = 1;
                while(m<k)
                    A(k,J) = A(k,J) - l(k,m)*u(m,J);
                    m = m + 1;
                end
                u(k,J) = A(k,J)/l(k,k);
                J = J+1;
            end
            k = k + 1;
        end
        z(1) = double(A(1,n+1)/l(1,1));
        j = 2;
        while(j<=n)
            y = 1;
            while(y<j)
                A(j,n+1) = A(j,n+1) - z(y)*l(j,y);
                y = y + 1;
            end
            z(j) = double(A(j,n+1)/l(j,j));
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
        [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        fprintf(fileid1,'Crout Method \r\n');
        fprintf(fileid1,'X \r\n');
        formatspec1 = '%0.5f\r\n';
        fprintf(fileid1,formatspec1,x);
        fprintf(fileid1,'L \r\n');
        formatspec2 = '%0.5f  ';
        for g = 1:n
        fprintf(fileid1,formatspec2,l(g,:));
        fprintf(fileid1,'\r\n');
        end
        fprintf(fileid1,'U \r\n');
        formatspec2 = '%0.5f  ';
        for g = 1:n
        fprintf(fileid1,formatspec2,u(g,:));
        fprintf(fileid1,'\r\n');
        end
    case 5
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
        u = trans_pose(L,n,n);
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
        [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        fprintf(fileid1,'Cholesky Method \r\n');
        fprintf(fileid1,'X \r\n');
        formatspec1 = '%0.5f\r\n';
        fprintf(fileid1,formatspec1,x);
        fprintf(fileid1,'L \r\n');
        formatspec2 = '%0.5f  ';
        for g = 1:n
        fprintf(fileid1,formatspec2,L(g,:));
        fprintf(fileid1,'\r\n');
        end
end
end
function[c] = trans_pose(a,m,n)
c = zeros(n,m,'double');
for i = 1:m
    for j = 1:n
            c(j,i) = a(i,j);
        
    end
end
end
%x = -1.5:0.1:1.5


