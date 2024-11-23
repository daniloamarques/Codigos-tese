function [flag,x] = EncontraCruzamento(p1,p2,p,d)
    
    tol = 1e-10;  
    flag = 0; % 0 se não se cruzam; 1 caso contrario.
    
    M = [p2(1)-p1(1) -d(1); p2(2)-p1(2) -d(2); p2(3)-p1(3) -d(3)];
    b = [p(1)-p1(1); p(2)-p1(2); p(3)-p1(3)];
    x = (M\b)';
    if(x(1) > 0+tol && x(1) < 1-tol && x(2) > 0+tol && x(2) < 1-tol )
    flag = 1;
    end
  
endfunction

 %Função que encontra o cruzamento entre dois segmentos