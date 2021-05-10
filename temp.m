clear all
clc

gr = 10;
i=10;
for m = gr:gr:100
    passo = m/gr;
    disp("=======")
    j = m;
    for k = passo: passo: m
     disp('i')
        i
        disp('j')
        j = j - 1
    end
    i = i + 1;
disp("---")
end