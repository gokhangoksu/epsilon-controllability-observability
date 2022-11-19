function i=sigma(t)
if (t<2)||(2.1<=t&&t<4)||(4.2<=t)
    i=1;
elseif (2<=t&&t<2.1)
    i=2;
elseif (4<=t&&t<4.2)
    i=3;
end
end