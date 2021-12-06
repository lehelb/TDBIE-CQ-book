function U = RKCQoperator(F,g,k,d1,flag,varargin)

% U = RKCQoperator(F,g,k,d1,flag)
% U = RKCQoperator(F,g,k,d1,flag,A)
%
% Input:
%   F(z,h)  : d1 x 1 valued function handle of ...
%             z : complex number, h : d2 x 1 complex valued vector
%   g       : S*d2 x (N+1) matrix (S=2 if RadauIIa is used)
%   d1      : Dimension of function space of the output    
%   k       : Real number, length of the time step
%   flag    : 1, return all stages, 0 return last stage (step) only
%   A       : S x S real matrix for RK methods. If not given 2-stage 
%             RadauIIa is implemented
%
% Output:
%   U    : S*d1 x (N+1) matrix (flag=1)
%   U    : d1 x (N+1) matrix (flag=0)
%
% Last modified: September 19, 2018

% If there is no A given, RadauIIa method is implemented
if nargin == 5
    A=[5/12 -1/12; 3/4 1/4];
else
    A=varargin{1};
end
S   = size(A,1); % number of stages
bT = A(end,:);
one = ones(S,1);

N  = size(g,2)-1; % number of time steps
d2 = size(g,1)/S; % dimension of the input space
omega= exp(2*pi*1i/(N+1));
R = eps^(0.5/(N+1));

g=bsxfun(@times,g,R.^(0:N));
g=fft(g,[],2); % hhat

idx=@(s) (s-1)*d1+1:(s*d1);
idy=@(s) (s-1)*d2+1:(s*d2); 

U=zeros(S*d1,N+1);
% Compute half the hermitian sequence
parfor l=0:floor((N+1)/2) 
    z = R*omega^(-l);
    Delta = inv(A + z/(1-z) * one * bT);
    [P,Lambda]=eig(Delta);
    Lambda=diag(Lambda)/k;
    gl=kron(inv(P),speye(d2))*g(:,l+1);
    ul=zeros(S*d1,1);
    for s=1:S
        ul(idx(s))=F(Lambda(s), gl(idy(s))); %#ok
    end
    U(:,l+1)=kron(P,speye(d1))*ul;
end
% Mirror the hermitian sequence
U(:,N+2-(1:floor(N/2)))=conj(U(:,2:floor(N/2)+1));
U=ifft(U,[],2);

U=bsxfun(@times,U,R.^(-(0:N)));
U=real(U);
% If all stages are asked, all rows are returned
if flag
    return
end
% If only the last stage was asked, the last block of the rows are returned
U=U((S-1)*d1+1:end,:);
return


