function [ang] = angulo(A,B,C)
  vet = dot(B-A,C-A);
  co = vet/(norm(B-A)*norm(C-A));
  ang = acos(co);
endfunction