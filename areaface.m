function [af] = areaface(A,B,C)
    af=1/2*norm(cross(B-A,C-A));  
  endfunction 