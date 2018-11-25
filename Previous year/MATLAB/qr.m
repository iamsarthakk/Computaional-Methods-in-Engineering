function[] = qrd()
[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
fileid = fopen(s , 'r');
I = fscanf(fileid, '%f', inf);
fclose(fileid);
n = I(1,1);
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
        disp(A(:,1));
end
%function[q,r] = qr_dec(a,n)
%q = zeros(n,n,'double');
%r = zeros(n,n,'double');
