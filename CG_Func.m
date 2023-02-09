function dF = CG_Func(N,Y,Ca,Cc,m,DOB,S,acrit,acheck)

[Mkma,Mkba,Mkmc,Mkbc,Yma,Yba,Ymc,Ybc] = geom_magtub(Y(1),Y(2),acrit,acheck);

dF = zeros(2,1);    % a column vector
% dF(1)= Ca*(S*(abs(Yma)*abs(Mkma)*(1-DOB) + abs(Yba)*abs(Mkba)*DOB)*sqrt(pi*Y(1)))^m;
% dF(2)= Cc*(S*(abs(Ymc)*abs(Mkmc)*(1-DOB) + abs(Ybc)*abs(Mkbc)*DOB)*sqrt(pi*Y(1)))^m;
dF(1)= Ca*(S*(Yma*Mkma*(1-DOB) + Yba*Mkba*DOB)*sqrt(pi*Y(1)))^m;
dF(2)= Cc*(S*(Ymc*Mkmc*(1-DOB) + Ybc*Mkbc*DOB)*sqrt(pi*Y(1)))^m;
end