function[] = q2_pa2()
[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
n = input(['Choose method.',... 
    '\n1. Power Mehtod',... 
    '\n2. Inverse Power Method',... 
    '\n3. Inverse Power Method with Shift',... 
    '\n4. QR Decomposition: ',]);
switch n
    case 1
        fileid = fopen(s , 'r');
        I = fscanf(fileid, '%f', inf);
        fclose(fileid);
        n = I(1);
        m = size(I);
        max_iter = I(m(1) - 2);
        rel_error = I(m(1) - 1);
        close_guess = I(m);
         A = zeros(n,n,'double');
                i = 1;
                while(i<=n)
                    j = 1;
                    while(j<=n)
                        A(i,j) = I((i-1)*(n-1) + (i+j));
                        j = j + 1;
                    end
                    i = i + 1;
                end      
        z = double.empty(n,0);
        i1 = 2;
        z(1) = 1;
        while(i1<= n)
            z(i1) = 0;
            i1 = i1 + 1;
        end
        e = 100;
        y1 = double.empty(n,0);
        y = trans_pose(z,1,n);
        %disp(z_tr);
        lambda = 0;
        i = 1;
        while((e>rel_error)&&(i<=max_iter))
            z = A * y;
            y1 = trans_pose(z,n,1);
            tmp = y1(1);
            iter = 2;
            while(iter<=n)
                tmp = sqrt((tmp(1))^2 + (y1(iter))^2);
                iter = iter + 1;
            end
            y_tr = trans_pose(y,n,1);
            %disp(y);
            lambda1 = (mul(mul(y_tr , A,1,n,n), y,1,n,1))/(mul(y_tr, y,1,n,1));
            e = abs((lambda1 - lambda)/(lambda1));
            %fprintf('Lamda\t%.6f %d %.9f\n',lambda1,i,e);
            lambda = lambda1;
            y = z / abs(tmp);
            %disp(lambda);
            %disp(i);
            i = i+1;
        end
        [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        fprintf(fileid1,'Power Method- \r\n');
        fprintf(fileid1,'%0.5f \r\n',lambda);
        fprintf(fileid1,'eigen vector \r\n');
        formatspec1 = '%0.5f\r\n';
        fprintf(fileid1,formatspec1,y);
        fprintf(fileid1,'Iterations = %d',i);
    case 2
        fileid = fopen(s , 'r');
        I = fscanf(fileid, '%f', inf);
        fclose(fileid);
        n = I(1);
        m = size(I);
        max_iter = I(m(1)-2);
        rel_error = I(m(1)-1);
        close_guess = I(m);
         A = zeros(n,n,'double');
                i = 1;
                while(i<=n)
                    j = 1;
                    while(j<=n)
                        A(i,j) = I((i-1)*(n-1) + (i+j));
                        j = j + 1;
                    end
                    i = i + 1;
                end    
        z = double.empty(n,0);
        i1 = 2;
        z(1) = 1;
        while(i1<= n)
            z(i1) = 0;
            i1 = i1 + 1;
        end
        e = 100;
        y1 = double.empty(n,0);
        y = trans_pose(z,1,n);
        %disp(z_tr);
        lambda = 0;
        i = 1;
        while((e>rel_error)&&(i<=max_iter))
            z = LU_crout(A , y,n);
            y1 = trans_pose(z,n,1);
            tmp = y1(1);
            iter = 2;
            while(iter<=n)
                tmp = sqrt((tmp(1))^2 + (y1(iter))^2);
                iter = iter + 1;
            end
            y_tr = trans_pose(y,n,1);
            %disp(y);
            lambda1 = (mul(mul(y_tr , A,1,n,n), y,1,n,1))/(mul(y_tr, y,1,n,1));
            e = abs((lambda1 - lambda)/(lambda1));
            %fprintf('Lamda\t%.6f %d %.9f\n',lambda1,i,e);
            lambda = lambda1;
            y = z / abs(tmp);
            %disp(lambda);
            %disp(i);
            i = i+1;
        end
        [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        fprintf(fileid1,'Inverse Power Method- \r\n');
        fprintf(fileid1,'%0.5f \r\n',lambda);
        fprintf(fileid1,'eigen vector \r\n');
        formatspec1 = '%0.5f\r\n';
        fprintf(fileid1,formatspec1,y);
        fprintf(fileid1,'Iterations = %d',i+1);
    case 3
        fileid = fopen(s , 'r');
        I = fscanf(fileid, '%f', inf);
        fclose(fileid);
        n = I(1);
        m = size(I);
        max_iter = I(m(1)-2);
        rel_error = I(m(1)-1);
        close_guess = I(m(1));
         B = zeros(n,n,'double');
                i = 1;
                while(i<=n)
                    j = 1;
                    while(j<=n)
                        B(i,j) = I((i-1)*(n-1) + (i+j));
                        j = j + 1;
                    end
                    i = i + 1;
                end
         A = zeros(n,n,'double');
         for j1 = 1:n
             for j2 = 1:n
                 if(j1==j2)
                     A(j1,j2) = (B(j1,j2) - close_guess);
                     if(A(j1,j2) ==0)
                         A(j1,j2) = 0.00000000001;
                     end
                 end
                 if(j1~=j2)
                     A(j1,j2) = B(j1,j2);
                 end
             end
         end
        z = double.empty(n,0);
        i1 = 2;
        z(1) = 1;
        while(i1<= n)
            z(i1) = 0;
            i1 = i1 + 1;
        end
        e = 100;
        y1 = double.empty(n,0);
        y = trans_pose(z,1,n);
        %disp(z_tr);
        lambda = 0;
        i = 1;
        while((e>rel_error)&&(i<=max_iter))
            z = LU_crout(A , y,n);
            y1 = trans_pose(z,n,1);
            tmp = y1(1);
            iter = 2;
            while(iter<=n)
                tmp = sqrt((tmp(1))^2 + (y1(iter))^2);
                iter = iter + 1;
            end
            y_tr = trans_pose(y,n,1);
            %disp(y);
            lambda1 = (mul(mul(y_tr , A,1,n,n), y,1,n,1))/(mul(y_tr, y,1,n,1));
            e = abs((lambda1 - lambda)/(lambda1));
            %fprintf('Lamda\t%.6f %d %.9f\n',lambda1,i,e);
            lambda = lambda1;
            y = z / abs(tmp);
            %disp(lambda);
            %disp(i);
            i = i+1;
        end
        eigen = lambda + close_guess;
        [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        fprintf(fileid1,'Inverse Power Method with Shift- \r\n');
        fprintf(fileid1,'%0.5f \r\n',eigen);
        fprintf(fileid1,'eigen vector \r\n');
        formatspec1 = '%0.5f\r\n';
        fprintf(fileid1,formatspec1,y);
        fprintf(fileid1,'Iterations = %d',i+1);
    case 4
        fileid = fopen(s , 'r');
        I = fscanf(fileid, '%f', inf);
        fclose(fileid);
        n = I(1,1);
        m = size(I);
        max_error = I(m(1) - 1);
        A = zeros(n,n,'double');
                i = 1;
                while(i<=n)
                    j = 1;
                    while(j<=n)
                        A(i,j) = I((i-1)*(n-1) + (i+j));
                        j = j + 1;
                    end
                    i = i + 1;
                end
                Q = q_r_dec(A,n);
                %disp(Q);
                R = q_r_dec1(A,Q,n);
                %disp(R);
                tmp_e = 100;
                iter2 = 1;
                while(tmp_e>max_error)
                    A = mul(R,Q,n,n,n);
                    Q = q_r_dec(A,n);
                    R = q_r_dec1(A,Q,n);
                    B = mul(R,Q,n,n,n);
                    tmp_e = abs((B(1,1) - A(1,1))/B(1,1));
                    for x = 2:n
                        tmp_e1 = abs((B(x,x) - A(x,x))/B(x,x));
                        if(tmp_e1>tmp_e)
                            tmp_e = tmp_e1;
                        end
                    end
                    iter2 = iter2 + 1;
                end
                eigen = double.empty(n,0);
                for u = 1:n
                   eigen(u) =  B(u,u);
                end
                [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        format_spec = '%0.4f\r\n ';
        fprintf(fileid1,'QR Decomposition Method\r\n');
        fprintf(fileid1,format_spec,eigen);
        fprintf(fileid1,'Iterations = %d\r\n',iter2);
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
function[c] = mul(a,b,m,n,p)
c = zeros(m,p,'double');
for i = 1:m
    for j = 1:p
        sum = 0;
        for k = 1:n
            sum = sum + a(i,k)*b(k,j);
        end
        c(i,j) = sum;
    end
end
end
function[c] = LU_crout(a,y,n)
A = zeros(n,n+1,'double');
for b = 1:n
    for r = 1:n
        A(b,r) = a(b,r);
    end
end
for w = 1:n
    A(w,n+1) = y(w);
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
c = trans_pose(x,1,n);
end
function[sum] = tra_ce(b,n)
sum = 0;
for i = 1:n
    sum = sum + b(i,i);
end
end
function[q] = q_r_dec(a,n)
q = zeros(n,n,'double');
q(:,1) = a(:,1)/no_rm(a(:,1),n);
for k = 1:(n-1)
    z1 = double.empty(n,0);
    sum = 0;
    for i = 1:k
        sum = sum + (mul(trans_pose(a(:,k+1),n,1),q(:,i),1,n,1))*q(:,i);
    end
    z1 = a(:,k+1) - sum;
    q(:,k+1) = z1/no_rm(z1,n);
end
end
function[r] = q_r_dec1(a,q,n)
q_t = trans_pose(q,n,n);
r = mul(q_t,a,n,n,n);
end
function[tmp] = no_rm(y1,n)
 tmp = y1(1);
    iter = 2;
    while(iter<=n)
        tmp = sqrt((tmp(1))^2 + (y1(iter))^2);
        iter = iter + 1;
    end
end
    