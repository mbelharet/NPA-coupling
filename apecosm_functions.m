clearvars;

CRISTAL_SLOPE = 3;
CRISTAL_CRIT = 10:200;

oope = 0:1e3;

for i = 1:length(CRISTAL_CRIT)
    
    school_tmp = (oope ).^ CRISTAL_SLOPE;

    school(i,:) = school_tmp ./ (school_tmp + CRISTAL_CRIT(i) .^ CRISTAL_SLOPE);

end

%%
figure;

plot(oope , school([1 50 100 150],:))