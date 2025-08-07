function [c,p]=phase_velocity(a,k,zheta,zt)

n=[sin(zheta+1i*zt),0,cos(zheta+1i*zt)];

G11=a(1)*n(1)^2+a(4)*n(3)^2;
G22=a(5)*n(1)^2+a(4)*n(3)^2;
G33=a(4)*n(1)^2+a(3)*n(3)^2;
G13=(a(2)+a(4))*n(1)*n(3);
c(1)=sqrt(0.5*(G11+G33+sqrt((G11+G33)^2-4*(G11*G33-G13^2))));
c(2)=sqrt(0.5*(G11+G33-sqrt((G11+G33)^2-4*(G11*G33-G13^2))));
c(3)=sqrt(G22);

switch k
    case 1
        p=n/c(1);
    case 2
        p=n/c(2);
    case 3
        p=n/c(3);
end


end