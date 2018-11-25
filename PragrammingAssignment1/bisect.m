%% How To Convert A User Input String To A Function
%% Wrote By Masoud Ghanbari
%% Visit Us On Our Website: www.mghanbari.ir
%% Contact Me There Is If Any Question: ghanbari.masoud7@gmail.com
%% Primaries
clc;clear all;
%% Function Acquisition
syms x
z = input(' Enter The Equation ( Like sin(x) ): ','s');
S = vectorize(char(z)); % For Signs Like %,^,&,.
f = str2func(['@(x) ', z]);
%% Diff & Int Of Function
d1f=diff(f,x); % First Order Diff
%fprintf(dlf);
intf=int(f,x);   % Integral
%% PLotting Section
subplot(131);ezplot(intf);grid on;title('The Functions Integral');
subplot(132);ezplot(f);grid on;title('The Function');
subplot(133);ezplot(d1f);grid on;title('The Functions Derivative');
%% End