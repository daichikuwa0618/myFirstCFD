if(exist("n")==0 || n<0) n = 0

filename = sprintf("output/time%04dd-3",n)
plot filename using 1:2 with lines

n = n + 1
if(n < 300) reread
undefine n