A=[31 -13 0 0 0 -10 0 0 0;
-13 35 -9 0 -11 0 0 0 0;
0 -9 31 -10 0 0 0 0 0;
0 0 -10 79 -30 0 0 0 -9;
0 0 0 -30 57 -7 0 -5 0;
0 0 0 0 -7 47 -30 0 0;
0 0 0 0 0 -30 41 0 0;
0 0 0 0 -5 0 0 27 -2;
0 0 0 -9 0 0 0 -2 29];
b=[-15;27;-23;0;-20;12;-7;7;10];
fprintf("Guass_seidel\n");
y = guass_seidel(A,b,100,1e-7);
fprintf("SOR\n");
for i=1:99
    w = i/50;
    y = SOR(A,b,w,2000,1e-7);
end