PT = ;
writetable(table(sol(1:5000, :, 1)), '../new5000e.txt','Delimiter',',')
writetable(table(sol(1:5000, :, 2)), '../new5000s.txt','Delimiter',',')
writetable(table(sol(1:5000, :, 3)), '../new5000p.txt','Delimiter',',')