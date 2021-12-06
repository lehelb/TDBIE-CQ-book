  function [x,w] = gauss(Q)
  %function [x,w] = gauss(Q)
  %Gauss nodes and weights

  beta = .5./sqrt(1-(2*(1:Q-1)).^(-2));
  T = diag(beta,1) + diag(beta,-1);
  [V,D] = eig(T);
  x = diag(D); [x,i] = sort(x);
  w = 2*V(1,i).^2;