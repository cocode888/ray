function [Tau1,Tau2]=equation_Tau(a,p)
deltz=eye(3,'like',complex(1,0));

%---------P_wave-----------
p1=p;

G1(1,1)=a(1)*p1(1)*p1(1)+a(4)*p1(3)*p1(3);
G1(1,3)=(a(2)+a(4))*p1(1)*p1(3);
G1(2,2)=a(5)*p1(1)*p1(1)+a(4)*p1(3)*p1(3);
G1(3,3)=a(4)*p1(1)*p1(1)+a(3)*p1(3)*p1(3);
GM1=G1-deltz;

% ---------- 9 possiable rank(GM1)=2 --------------
B1(1)=GM1(1,1)*GM1(2,2);
B1(2)=-GM1(1,3)*GM1(2,2);
B1(3)=GM1(1,1)*GM1(3,3)-GM1(1,3)*GM1(1,3);
B1(4)=GM1(2,2)*GM1(3,3);

Tau1=abs(B1(1));
Tau2=abs(B1(4));

end