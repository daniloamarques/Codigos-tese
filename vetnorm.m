 function [vn] = vetnorm(A,B,C)
    vn=cross((B-A),(C-A));
    vn=vn/norm(vn);  
  endfunction 