function [ y ] = Euler(y0, t0, t_end, dt, f);
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N=ceil((t_end-t0)/dt);

y=zeros(N+1, 1);
y(1)=y0;

for k = 1:1:N;

    t=t0+(k-1)*dt;
    y(k+1)=y(k)+dt*f(t, y(k));

end


end
