function[u]=Mackey_glass()
u(1:17) = 0.5;
a = 0.9;
b = 0.2;
c = 17;

for i=c+1:110020
   u(i) = a*u(i-1) + b*u(i-c)/(1+u(i-c)^10);  
end
u = u(18:end);
end