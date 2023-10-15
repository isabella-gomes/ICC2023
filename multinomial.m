function ind = multinomial(p)
r = rand(1)*sum(p);
total = 0;
i = 1;
while total <= r
total = total + p(i);
i=i+1;
end
ind = i-1;
end
