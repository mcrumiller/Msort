function [x,y,z]=ellipsoid_cov(C,R,c,n)
%--------------------------------------------------------------------------
% ellipsoid_cov.m - Generates an ellipoid of (n+1) x (n+1) x (n+1) points
% with center at [xc,yc,yz] and with radii [xr,yr,zr], oriented with
% covariance matric c.
%
% Note: that the radii are in units of standard deviations, or
% sqrt(eigenvalues)
%
% Note 2: This will also work for a 2-dimensional ellipse
%
% Usage: ellipsoid_cov();
%
% Input:  
% Output: 
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
[U,L]=eig(c);
radii=R .* sqrt(diag(L))';

% ellipsoid
if(length(C)==3)
    % generate unrotated ellipsoid
    [xe,ye,ze]=ellipsoid(0,0,0,radii(1),radii(2),radii(3),n);
    a=kron(U(:,1),xe); b=kron(U(:,2),ye);c=kron(U(:,3),ze);
    data=a+b+c;
    n=n+1;
    x = data(1:n,:)+C(1); y = data(n+1:2*n,:)+C(2); z = data(2*n+1:end,:)+C(3);
% ellipse
else
    dt=2*pi/n; t=0:dt:(dt*n);
    xe=radii(1)*cos(t); ye=radii(2)*sin(t);
    
    % determine rotation angle
    alpha=angle(complex(U(1,1),U(2,1)));
    x = xe*cos(alpha) - ye*sin(alpha)+C(1);
    y = xe*sin(alpha) + ye*cos(alpha)+C(2);
end