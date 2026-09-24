function dGT=odefun_alphadeltaRIO(t,GT)
global  lam1 lam2 lam3 lam4

dGT=zeros(4,1);
dGT(1)=GT(3);
dGT(2)=GT(4);
dGT(3)=-2*GT(3)*GT(4)/tan(GT(2))+(lam4*cos(GT(1))-lam3*sin(GT(1)))/sin(GT(2))-GT(4)*(lam1*sin(GT(1))-lam2*cos(GT(1)));
dGT(4)=GT(3)^2*sin(GT(2))*cos(GT(2))+cos(GT(2))*(lam3*cos(GT(1))+lam4*sin(GT(1)))...
    +GT(3)*sin(GT(2))^2*(lam1*sin(GT(1))-lam2*cos(GT(1)));

