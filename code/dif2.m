function y=dif2(x,t)
    x_diff=diff(x)./diff(t);
    y = interp1((t(1:end-1)+t(2:end))/2,x_diff,t,'spline',0);
end