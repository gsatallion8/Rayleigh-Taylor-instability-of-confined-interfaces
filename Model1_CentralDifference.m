N=20;
h1=1;
h2=1;
a=1;
ca=-9.8696;
r1=10;
r2=19.8696;
A=zeros(2*N^2+N,2*N^2+N);
B=zeros(2*N^2+N,2*N^2+N);
%MATRIX A:
%Differential equation:____________________________________________________
i=N+1;
while i<=N^2
    A(i,i-1)=1;
    A(i,i)=-(2+8/h2^2);
    A(i,i+1)=1;
    A(i,i+N)=4/h2^2;
    A(i,i-N)=4/h2^2;
    i=i+1;    
end
%__________________________________________________________________________
i=N^2+1;
while i<=N^2+N
    A(i,i+1)=1;
    A(i,i)=-(2+4*ca/N^2);
    A(i,i-1)=1;    
    i=i+1;
end
%__________________________________________________________________________
i=N^2+N+1;
while i<=2*N^2
    A(i,i-1)=1;
    A(i,i)=-(2+8/h1^2);
    A(i,i+1)=1;
    A(i,i+N)=4/h1^2;
    A(i,i-N)=4/h1^2;
    i=i+1;    
end
%Boundary conditions:______________________________________________________
i=1;
while i<=N
    j=1;
    while j<=2*N^2+N
        A(i,j)=0;
        j=j+1;
    end
    A(i,i+N)=1;
    A(i,i)=-1;
    i=i+1;
end
i=N^2-N+1;
while i<=N^2+1
    j=1;
    while j<=2*N^2+N
        A(i,j)=0;
        j=j+1;
    end
    A(i,i-N)=-1;
    A(i,i)=1;
    i=i+1;
end
i=1;
while i<=N-2
    j=1;
    while j<=2*N^2+N
        A(i*N+1,j)=0;
        j=j+1;
    end
    A(i*N+1,i*N+1)=1;
    A(i*N+1,i*N+2)=-1;
    i=i+1;
end
i=1;
while i<=N-2
    j=1;
    while j<=2*N^2+N
        A(i*N+N,j)=0;
        j=j+1;
    end
    A(i*N+N,i*N+N-1)=-1;
    A(i*N+N,i*N+N)=1;
    i=i+1;
end
%__________________________________________________________________________
i=1+N^2;
j=1;
while j<=2*N^2+N
    A(i,j)=0;
    j=j+1;
end
A(i,i)=1;
A(i,i+1)=-1;
i=N+N^2;
j=1;
while j<=2*N^2+N
    A(i,j)=0;
    j=j+1;
end
A(i,i)=1;
A(i,i-1)=-1;
%__________________________________________________________________________
i=N^2+N+1;
while i<=N^2+2*N
    j=1;
    while j<=2*N^2+N
        A(i,j)=0;
        j=j+1;
    end
    A(i,i+N)=1;
    A(i,i)=-1;
    i=i+1;
end
i=2*N^2+1;
while i<=2*N^2+N
    j=1;
    while j<=2*N^2+N
        A(i,j)=0;
        j=j+1;
    end
    A(i,i-N)=-1;
    A(i,i)=1;
    i=i+1;
end
i=1;
while i<=N-2
    j=1;
    while j<=2*N^2+N
        A((i+1)*N+N^2+1,j)=0;
        j=j+1;
    end
    A((i+1)*N+N^2+1,(i+1)*N+N^2+2)=-1;
    A((i+1)*N+N^2+1,(i+1)*N+N^2+1)=1;
    i=i+1;
end
i=1;
while i<=N-2
    j=1;
    while j<=2*N^2+N
        A((i+1)*N+N^2+N,j)=0;
        j=j+1;
    end
    A((i+1)*N+N^2+N,(i+1)*N+N^2+N-1)=-1;
    A((i+1)*N+N^2+N,(i+1)*N+N^2+N)=1;
    i=i+1;
end
%__________________________________________________________________________
%MATRIX B:
i=1;
while i<=N
    B(i,i+N^2)=h2*a*1i/N;
    B(i+N^2,i+N+N^2)=4*r1*1i/N^2;
    B(i+N^2,i)=-4*r2*1i/N^2;
    B(i+N+N^2,i+N^2)=h1*a*1i/N;
    i=i+1;
end
%__________________________________________________________________________
%EIGEN VALUES:
[V,D]=eig(A,B);
e=diag(D);
plot(e,'.');

% e=eig(A,B);
% plot(e,'.');

% %EIGEN VECTORS:
% k=8;
% plot(real(V(N^2+1:N^2+N,k)));