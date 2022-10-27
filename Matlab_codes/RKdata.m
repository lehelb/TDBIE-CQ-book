function [A,b,c] = RKdata(RK)
%function [A,b,c] = RKdata(RK)
%returns coefficients in Butcher notation
%RK = 0 for 2 stage Radau IIA
%RK = 1 for 3 stage Radau IIA
%RK = 2 for 3 stage Lobatto IIIC
%RK = 3 for 2 stage Gauss method
%RK = 4 for 3 stage Gauss

if (RK==0)
    %radau IIa of order 3 stage order 2
    A = [5/12 -1/12; 3/4 1/4]; c = [1/3; 1]; 
    b = [3/4; 1/4];
elseif (RK==1)
%radau IIa order 5, stage order 3
    A = [(88-7*sqrt(6))/360 (296-169*sqrt(6))/1800 (-2+3*sqrt(6))/225;
         (296+169*sqrt(6))/1800 (88+7*sqrt(6))/360 (-2-3*sqrt(6))/225;
        (16-sqrt(6))/36 (16+sqrt(6))/36 1/9];
    b = [(16-sqrt(6))/36 (16+sqrt(6))/36 1/9].';
    c = [(4-sqrt(6))/10 ; (4+sqrt(6))/10 ; 1];
elseif(RK==2)
%Lobatto IIIC, order 6, stage order 3
    a1 = sqrt(5);
    A = [1/12 -a1/12 a1/12 -1/12; 1/12 1/4 (10-7*a1)/60 a1/60;
         1/12 (10+7*a1)/60 1/4  -a1/60; 1/12 5/12 5/12 1/12];
    b = [1/12; 5/12; 5/12; 1/12];
    c = [0 ;(5-a1)/10; (5+a1)/10; 1];
elseif (RK==3)
    %Gauss method, order 4
    a = sqrt(3)/6; c = [1/2-a; 1/2+a];
    A = [1/4 1/4-a; 1/4+a 1/4];
    b = [1/2; 1/2];
elseif (RK==4)
    a = sqrt(15)/5; c = [1/2-a/2; 1/2; 1/2+a/2];
    A = [5/36 2/9-a/3 5/36-a/6; 5/36+a*5/24 2/9 5/36-a*5/24; 5/36+a/6 2/9+a/3 5/36];
    b = [5/18; 4/9; 5/18];
else
    fprintf('only RK types 0, 1,...,4  possible, see help\n');
    return;
end
