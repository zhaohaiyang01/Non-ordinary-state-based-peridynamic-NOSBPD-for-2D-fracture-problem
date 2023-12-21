function wf = HigherF(x,wf) 
length = unique(single(x));
quado = size(length,1);
if quado == 3 % 2 <= m <= sqrt(5)
%     w1 = 1;
%     w2 = 0;
%     w3 = -1/16;
    w = [1,0,-1/16];
elseif quado == 4 % sqrt(5) <= m <= sqrt(8)
%     w1 = 1;
%     w2 = -8/46;
%     w3 = -4/46;
%     w4 = 1/46;
    w = [1,-8/46,-4/46,1/46];
elseif quado == 5 % sqrt(8) <= m <= 3
%     w1 = 1;
%     w2 = 8/14;
%     w3 = -2/14;
%     w4 = 1/14;
%     w5 = -1/14;
    w = [1,8/14,-2/14,1/14,-1/14];
elseif quado == 6 % 3 <= m <= sqrt(10)
    w = [1,32/87,-15/87,-8/87,2/87,1/87];
else 
    error('check quado in HigherF function')
end
for i = 1:quado
    index = single(x)==length(i);
    wf(index) = w(i);
end