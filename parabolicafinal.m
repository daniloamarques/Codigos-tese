close all;
clear;
clc;
tic

  %Leitura do Arquivo '.obj'
  [V,F] = leitura_obj('sela2.obj');
  
  %Numero de vértices e de faces

  nverts = size(V,1);
  nfaces = size(F,1);

##%------------------------------ Contas Discretas -------------------------------
##
##
##
## 
##
  S=[];
  Face=[];
  L=[];
  
  %Definindo a estrela, a quantidade de ligação de cada ponto e os triângulos que 
  %o ponto pertence

  for i=1:nverts
      count3=1; %conta a quantidade de vertices na estrela
      count3a=1;%conta a quantidade de faces na estrela
      for j=1:nfaces
        if (i == F(j,1))
          S(i,count3)=F(j,2);
          Face(i,count3a)=j;
          count3=count3+1;
          S(i,count3)=F(j,3);
          count3=count3+1;
          count3a=count3a+1;
        elseif (i == F(j,2))
          S(i,count3)=F(j,3);
          Face(i,count3a)=j;
          count3=count3+1;
          S(i,count3)=F(j,1);
          count3=count3+1;
          count3a=count3a+1;
        elseif (i == F(j,3))
          S(i,count3)=F(j,1);
          Face(i,count3a)=j;
          count3=count3+1;  
          S(i,count3)=F(j,2);
          count3=count3+1; 
          count3a=count3a+1;
        end
        L(i,1)=count3-1; %numero de pontos na estrela
        LL(i,1)=count3a-1; %numero de triangulos na estrela
      end  
  end
  
  
  Sc=[];
  for i=1:size(S,1)
    X=S(i,1:L(i))';
    Xc=unique(X,"rows");
    for j=1:size(Xc,1)
      Sc(i,j)=Xc(j,1);
    endfor
  endfor


for i=1:nverts
  count=0;  
  for j=1:size(Sc,2)
    if(Sc(i,j)!=0)
    count=count+1;
    endif
  endfor
  L(i,1)=count;
endfor
  
  printf('Parte 1');
  
% Calculando o vetor normal no vértice 
  
 for i=1:nverts
    for j=1:LL(i)
      areafaces1(i,j)=areaface(V(F(Face(i,j),1),:),V(F(Face(i,j),2),:),V(F(Face(i,j),3),:));
    endfor
  endfor
  
  b=zeros(nverts,3);
  c=zeros(nverts,3);
  Nvt=zeros(nverts,3);
  
  for i=1:nverts
    sum33=[0 0 0];
    for j=1:LL(i)
      sum33 = sum33 + areafaces1(i,j) * vetnorm(V(F(Face(i,j),1),:),V(F(Face(i,j),2),:),V(F(Face(i,j),3),:));
    end  
    Nvt(i,:)=sum33/norm(sum33); 
    
    xn=Nvt(i,:)';
    [Q R]=qr(xn);
    x1=Q(:,1);
    if (x1'*xn<0)
      Q=-Q;
    end
    x1=Q(:,1);  %vetor normal no vértice
    b(i,:)=Q(:,2)'; %vetor1 da base do plano tangente
    c(i,:)=Q(:,3)'; %vetor2 da base do plano tangente
  end 
  
% Cálculo da Área de Célula (utilizando Área Baricêntrica) 

Amixed=zeros(nverts,1);

for i=1:nverts
  for j=1:LL(i)
    Amixed(i)= Amixed(i) + 1/3 * areafaces1(i,j);
  endfor  
endfor

%Curvatura Gaussiana


Kd1=[];
Kd2=[];

for i=1:nverts
  sumK = 0.0;
    for j = 1:LL(i)
    if (F(Face(i,j),1) == i)
      sumK = sumK + angulo(V(F(Face(i,j),1),:),V(F(Face(i,j),2),:),V(F(Face(i,j),3),:));  
    endif
    if (F(Face(i,j),2) == i)
      sumK = sumK + angulo(V(F(Face(i,j),2),:),V(F(Face(i,j),1),:),V(F(Face(i,j),3),:));  
    endif
    if (F(Face(i,j),3) == i)
      sumK = sumK + angulo(V(F(Face(i,j),3),:),V(F(Face(i,j),2),:),V(F(Face(i,j),1),:));  
    endif    
  endfor
  
  % Curvatura Gaussiana 1: Sulivan
    
  Kd1(i,:)=2*pi-sumK;
  
  %Curvatura Gaussiana 2: K= (2pi-soma angulos)/A_mixed : Meyer
  
  Kd2(i,:)=Kd1(i,:)/Amixed(i,:);
endfor


  printf('Parte 2');
  
  %Cálculo dos Pontos Parabólicos
  
 Positivo=[];
 Negativo=[];
 count55=1;
 count56=1;
  
  for i=1:nverts
    if(Kd2(i,1)>1e-10)
      Positivo(count55,:)=V(i,:);
      count55=count55+1;
    end
    if (Kd2(i,1)<-1e-10)
      Negativo(count56,:)=V(i,:);
      count56=count56+1;
    end
  endfor
  
    count57=1;
    
    for i=1:nfaces
    if ((Kd2(F(i,1),1)>1e-10 && Kd2(F(i,2),1)<-1e-10))
      t=Kd2(F(i,1),1)/(Kd2(F(i,1),1)-Kd2(F(i,2),1));
      MM(count57,:)=V(F(i,1),:)+t*(V(F(i,2),:)-V(F(i,1),:));
      PMM(count57,1)=F(i,1);
      PMM(count57,2)=F(i,2);
      count57=count57+1;
    end
    if ((Kd2(F(i,2),1)>1e-10 && Kd2(F(i,3),1)<-1e-10))
      t=Kd2(F(i,2),1)/(Kd2(F(i,2),1)-Kd2(F(i,3),1));
      MM(count57,:)=V(F(i,2),:)+t*(V(F(i,3),:)-V(F(i,2),:));
      PMM(count57,1)=F(i,2);
      PMM(count57,2)=F(i,3);
      count57=count57+1;
    end
    if ((Kd2(F(i,3),1)>1e-10 && Kd2(F(i,1),1)<-1e-10))
      t=Kd2(F(i,3),1)/(Kd2(F(i,3),1)-Kd2(F(i,1),1));
      MM(count57,:)=V(F(i,3),:)+t*(V(F(i,1),:)-V(F(i,3),:));
      PMM(count57,1)=F(i,3);
      PMM(count57,2)=F(i,1);
      count57=count57+1;
    endif
  end
  
  
  printf('Parte 3');
  
  % Determinação de cada curva parabólica (juntando os pontos para
  % formação das curvas)
  
  MM3=[MM PMM];

MM2=[];
jafoi=zeros(size(MM3)(1),1);
kk=1;
for ii=1:size(MM3)(1)
  if (jafoi(ii)==0)
    MM2(1,3*kk-2:3*kk)=MM3(ii,1:3);
    PMM2(1,2*kk-1:2*kk)=MM3(ii,4:5);
cc=1;
ll=ii;
jjj=0;

jafoi(ii)=1;
while (ll != jjj)
  jjj=ll;
  for i=1:nfaces
    if (MM3(ll,4)==F(i,1) && MM3(ll,5)==F(i,2))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3(ll,4)==F(i,1) && MM3(ll,5)==F(i,2))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end

      if (MM3(ll,4)==F(i,2) && MM3(ll,5)==F(i,3))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)))
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3(ll,4)==F(i,2) && MM3(ll,5)==F(i,3))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)))
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3(ll,4)==F(i,3) && MM3(ll,5)==F(i,2))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)))
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3(ll,4)==F(i,2) && MM3(ll,5)==F(i,1))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
  end

      if (MM3(ll,4)==F(i,3) && MM3(ll,5)==F(i,1))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)))
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3(ll,4)==F(i,3) && MM3(ll,5)==F(i,1))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3(ll,4)==F(i,1) && MM3(ll,5)==F(i,3))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)))
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3(ll,4)==F(i,1) && MM3(ll,5)==F(i,3))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end

endfor
end

% Fechando as curvas fechadas

for i=1:nfaces
if((PMM2(1,2*kk-1)==F(i,1) && PMM2(1,2*kk)==F(i,2) && PMM2(cc,2*kk-1)==F(i,3)) || (PMM2(1,2*kk-1)==F(i,1) && PMM2(1,2*kk)==F(i,2) && PMM2(cc,2*kk)==F(i,3)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,2) && PMM2(1,2*kk)==F(i,3) && PMM2(cc,2*kk-1)==F(i,1)) || (PMM2(1,2*kk-1)==F(i,2) && PMM2(1,2*kk)==F(i,3) && PMM2(cc,2*kk)==F(i,1)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,3) && PMM2(1,2*kk)==F(i,1) && PMM2(cc,2*kk-1)==F(i,2)) || (PMM2(1,2*kk-1)==F(i,3) && PMM2(1,2*kk)==F(i,1) && PMM2(cc,2*kk)==F(i,2)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,2) && PMM2(1,2*kk)==F(i,1) && PMM2(cc,2*kk-1)==F(i,3)) || (PMM2(1,2*kk-1)==F(i,2) && PMM2(1,2*kk)==F(i,1) && PMM2(cc,2*kk)==F(i,3)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,3) && PMM2(1,2*kk)==F(i,2) && PMM2(cc,2*kk-1)==F(i,1)) || (PMM2(1,2*kk-1)==F(i,3) && PMM2(1,2*kk)==F(i,2) && PMM2(cc,2*kk)==F(i,1)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,1) && PMM2(1,2*kk)==F(i,3) && PMM2(cc,2*kk-1)==F(i,2)) || (PMM2(1,2*kk-1)==F(i,1) && PMM2(1,2*kk)==F(i,3) && PMM2(cc,2*kk)==F(i,2)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
endif
endfor
pontocurvas(kk,1)=cc;
kk=kk+1;
endif
endfor

printf('Parte 4');

%Contando curvas parabólicas fechadas e abertas

con6=0;
ncon6=0;
for tt=1:3:size(MM2,2)
  if (pontocurvas((tt+2)/3,1)>1 && MM2(1,tt)==MM2(pontocurvas((tt+2)/3,1),tt) && MM2(1,tt+1)==MM2(pontocurvas((tt+2)/3,1),tt+1))
    circulo(con6+1,1)=tt;
    con6=con6+1; %Quantidade de curvas parabólicas fechada
  else
    ncirculo(ncon6+1,1)=tt;
    ncon6=ncon6+1; %Quantidade de curvas parabólicas aberta
  endif
endfor



   %Plotagem dos pontos para triangulacao
x = V(:,1);
y = V(:,2);
z = V(:,3);

##  Resultado da implementacao octave
hold on;

##  trimesh (F,x,y,z);   %Ative caso queira a imagem com a malha

  for j=1:size(pontocurvas)(1)
    plot3(MM2(1:pontocurvas(j,1),3*j-2), MM2(1:pontocurvas(j,1),3*j-1), MM2(1:pontocurvas(j,1),3*j),'b-','LineWidth',[3]);
  endfor
  
axis square
hold off;
  
tempo=toc
