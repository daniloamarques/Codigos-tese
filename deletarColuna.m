function [A] = deletarColuna(Matrix,coluna)
  
  if (coluna == 1)
    A = Matrix(:,2:size(Matrix)(2));
    return;
  end
  
  if (coluna == size(Matrix)(2))
    A = Matrix(:,1:size(Matrix)(2)-1);
    return;
  end
  
  A1 = Matrix(:,1:coluna-1);
  A2 = Matrix(:,coluna+1:size(Matrix)(2));

  A = [A1  A2];
    
end