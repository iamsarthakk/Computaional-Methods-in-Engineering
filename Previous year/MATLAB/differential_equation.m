[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
fileid = fopen(s, 'r');
tline = fgetl(fileid);
ip = cell(4,1);
tmp = 1;
while ischar(tline)
    ip(tmp) = {tline};
    tline = fgetl(fileid);
    tmp = tmp + 1;
end
fclose(fileid);
str = char(ip(1));
f = inline(str);
str1 = char(ip(2));
str2 = strsplit(str1,',');
t_not = str2double(str2(1));
y_not = str2double(str2(2));
str3 = char(ip(3));
t_f = str2double(str3);
str4 = char(ip(4));
h = str2double(str4);
Z = input(['Choose method.',...
    '\n1. Forward Euler Method',... 
    '\n2. 2nd order RK Method',... 
    '\n3. 4th order RK Method: ']);
switch Z
    case 1
        n_p = double((t_f - t_not)/h);
        t = double.empty(n_p+1,0);
        y = double.empty(n_p+1,0);
        for k = 1:(n_p+1)
            t(k) = t_not + (k-1)*h;
            if k == 1
                y(k) = y_not;
            end
            if k ~= 1
                sum = h * f(t(k-1),y(k-1));
                y(k) = y(k-1)+ sum;
            end
        end
    case 2
        n_p = double((t_f - t_not)/h);
        t = double.empty(n_p+1,0);
        y = double.empty(n_p+1,0);
        for k = 1:(n_p+1)
            t(k) = t_not + (k-1)*h;
            if k == 1
                y(k) = y_not;
            end
            if k ~= 1
                sum = h/3*(f(t(k-1),y(k-1)) + 2*f(t(k-1) + 0.75*h,y(k-1) + 0.75*h*f(t(k-1),y(k-1))));
                y(k) = y(k-1)+ sum;
            end
        end
    case 3
        n_p = double((t_f - t_not)/h);
        t = double.empty(n_p+1,0);
        y = double.empty(n_p+1,0);
        for k = 1:(n_p+1)
            t(k) = t_not + (k-1)*h;
            if k == 1
                y(k) = y_not;
            end
            if k ~= 1
                k1 = f(t(k-1),y(k-1));
                k2 = f(t(k-1) + 0.5*h, y(k-1) + 0.5*h*k1);
                k3 = f(t(k-1) + 0.5*h, y(k-1) + 0.5*h*k2);
                k4 = f(t(k-1) + h, y(k-1) + h*k1);
                sum = h/6*(k1 + 2*k2 + 2*k3 + k4);
                y(k) = y(k-1)+ sum;
            end
        end
end
        [filename,filepath] = uigetfile('.txt','File Selector');
        s1= fullfile(filepath,filename);
        fileid1 = fopen(s1 , 'w');
        fprintf(fileid1,'t\ty \r\n');
        formatspec1 = '%0.4f\t%0.8f\r\n';
        for ind = 1:(n_p+1)
        fprintf(fileid1,formatspec1,t(ind),y(ind));
        end
        plot(t,y,'.','MarkerSize',25)
        xlabel('t')
        ylabel('y')
        grid on