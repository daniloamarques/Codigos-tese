function [A] = deletarLinha(Matrix,linha)
  
  if (linha == 1)
    A = Matrix(2:size(Matrix)(1),:);
    return;
  end
  
  if (linha == size(Matrix)(1))
    A = Matrix(1:size(Matrix)(1)-1,:);
    return;
  end
  
  A1 = Matrix(1:linha-1,:);
  A2 = Matrix(linha+1:size(Matrix)(1),:);

  A = [A1 ; A2];
    
end