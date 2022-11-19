function wout=w(t)
if (0.25<=t&&t<=0.4)
    wout=0.05*sin(pi/0.15*(t-0.25))+0.01*rand(1);
elseif (0.65<=t&&t<=0.8)
    wout=0.05*sin(pi/0.15*(t-0.65))+0.01*rand(1);
else
    wout=0.01*rand(1);
end
end