%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%PML_BPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%format long;
lam = 1.55;
w0=10*lam; %beam width in um, FWHM=2ln(2)w
dx=lam;
dy=dx; 
Nx=100;
Ny=Nx;
Nz=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Random generator%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
A=rand(Nx,Ny);
for ii=1:Nx 
   for jj=1:Ny 
     if ((ii-Nx/2).^2+(jj-Ny/2).^2<10^2)
         if A(ii,jj)<0.5
          A(ii,jj)=1;
        else
        A(ii,jj)=0;
         end
     else 
       A(ii,jj)=0; 
     end
   end
end 
imshow(A); 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
x0=Nx/2*dx;
y0=Ny/2*dy; 
clear E0;
for ii=1:Nx
  for jj=1:Ny
	E0(ii,jj)=exp(-((ii.*dx-x0).^2+(jj.*dy-y0).^2)/(w0).^2);
  end
end
k0=2*pi/lam;
n=1;
n0=1.5;
ng=1.5; %permittivity of waveguide
n=0.5*(n0+ng);
dz=0.5.*k0*dx^2*n0;
E=zeros(Nx,Ny,2);
k1=zeros(Nx,Ny,1);
k2=zeros(Nx,Ny,1);
k3=zeros(Nx,Ny,1);
k4=zeros(Nx,Ny,1);
E(1:Nx,1:Ny,1)=E0; 
for mm=1:Nz
display(mm*dz);
  for ii=2:Nx-1
     for jj=2:Ny-1
%%%%%%%%%%%%%
%{
if ((ii-Nx/2).^2+(jj-Ny/2).^2<10^2)
        if A(ii,jj)<0.5
          n=ng;
        else
          n=n0;
         end
    else
      n=n0;
end
%}
%%%%%%%%%%%%%
k1(ii,jj)=(-j/2/n0/k0)*((E(ii+1,jj,1)-2*E(ii,jj,1)+E(ii-1,jj,1))./dx^2+(E(ii,jj+1, 1)-2*E(ii,jj,1)+E(ii,jj-1,1))./dy^2+(n.^2-n0.^2).*k0^2.*E(ii,jj,1));
end
end
E(:,:,2)=E(:,:,1)+0.5*dz.*k1;
for ii=2:Nx-1
 for jj=2:Ny-1
%%%%%%%%%%%%%
%{
if ((ii-Nx/2).^2+(jj-Ny/2).^2<10^2)
        if A(ii,jj)<0.5
          n=ng;
        else
          n=n0;
         end
    else
      n=n0;
end
%}
%%%%%%%%%%%%%

k2(ii,jj)=(-j/2/n0/k0)*((E(ii+1,jj,2)-2*E(ii,jj,2)+E(ii-1,jj,2))./dx^2+(E(ii,jj+1,2)-2*E(ii,jj,2)+E(ii,jj-1,2))./dy^2+(n.^2-n0.^2).*k0.^2.*E(ii,jj,2));
end
end
E(:,:,2)=E(:,:,1)+0.5*k2.*dz;
for ii=2:Nx-1
 for jj=2:Ny-1
%%%%%%%%%%%%%
%{
if ((ii-Nx/2).^2+(jj-Ny/2).^2<10^2)
        if A(ii,jj)<0.5
          n=ng;
        else
          n=n0;
         end
    else
      n=n0;
end
%}
%%%%%%%%%%%%%%%
k3(ii,jj)=(-j/2/n0/k0)*((E(ii+1,jj,2)-2*E(ii,jj,2)+E(ii-1,jj,2))./dx^2+(E(ii,jj+1,2)-2*E(ii,jj,2)+E(ii,jj-1,2))./dy^2+(n.^2-n0.^2).*k0.^2.*E(ii,jj,2));
end
end
E(:,:,2)=E(:,:,1)+k3.*dz;
for ii=2:Nx-1
 for jj=2:Ny-1
%%%%%%%%%%%%%
%{
if ((ii-Nx/2).^2+(jj-Ny/2).^2<10^2)
        if A(ii,jj)<0.5
          n=ng;
        else
          n=n0;
         end
    else
      n=n0;
end
%}
%%%%%%%%%%%%%%%
k4(ii,jj)=(-j/2/n0/k0)*((E(ii+1,jj,2)-2*E(ii,jj,2)+E(ii-1,jj,2))./dx^2+(E(ii,jj+1,2)-2*E(ii,jj,2)+E(ii,jj-1,2))./dy^2+(n.^2-n0.^2).*k0.^2.*E(ii,jj,2));
end
end
E(:,:,2)=E(:,:,1)+(k1+2*k2+2*k3+k4)*dz/6;
A=abs(E(:,:,2));
%display(abs(E(:,:,mm)))
%%%%%%%%%%%%%%%%TBC Right boundary%%%%%%%%%%%%%%%%%%%%%
kr1=j/sqrt(2)./dx.*log(E(Nx-1,:,2)./E(Nx-2,:,2)); 
E(Nx,:,2)=E(Nx-1,:,2).*exp(-j*kr1*dx*sqrt(2));

%%%%%%%%%%%%%%%%TBC Left boundary%%%%%%%%%%%%%%%%%%%%%
kr2=j/sqrt(2)./dx.*log(E(2,:,2)./E(3,:,2));
E(1,:,2)=E(2,:,2).*exp(-j*kr2*dx*sqrt(2));

%%%%%%%%%%%%%%%%%TBC Upper boundary%%%%%%%%%%%%%%%%%%%%%
kr3=j./dx/sqrt(2).*log(E(:,2,2)./E(:,3,2));
%display(abs(E(:,2,mm)));
E(:,1,2)=E(:,2,2).*exp(-j*kr3*sqrt(2)*dx);


%%%%%%%%%%%%%%%%%%TBC Lower boundary %%%%%%%%%%%%%%%%%%
kr4=j./dx/sqrt(2).*log(E(:,Ny-1,2)./E(:,Ny-2,2));
E(:,Ny,2)=E(:,Ny-1,2).*exp(-j*kr4*sqrt(2)*dx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,:,1)=E(:,:,2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surf(abs(E(:,:,2))); 
shading interp; 
