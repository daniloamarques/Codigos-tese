close all;
clear;
clc;
tic

  %Leitura do Arquivo '.obj'
  [V,F] = leitura_obj('Helencrianca.obj');
  
  %Decisão de fazer apenas as curvas ridges ou 
  %curvas ridges e curvas subparabolicas
  
  resposta = input('Digite 1 para fazer apenas curvas ridges ou digite 2 para fazer curvas ridges e curvas subparabolicas: ');
  while (resposta!=1 && resposta!=2)
    clc
    printf('Resposta Inválida! Tente novamente! \n');
    resposta = input('Digite 1 para fazer apenas curvas ridges ou digite 2 para fazer curvas ridges e curvas subparabolicas: ');
  endwhile

  
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
  
  %Definindo a estrela, a quantidade de ligação de cada ponto 
  %e os triângulos que o ponto pertence

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

SS=zeros(nverts,max(L));
idx=zeros(nverts,1);

  for i=1:nverts
      temp=1;
      
      while(idx(i)==0)
      count33=1;
      m=Sc(i,temp); 
      for k=1:L(i)
        for j=1:nfaces
          if (i == F(j,1) && m == F(j,2))
            SS(i,count33)=F(j,2);
            if (F(j,3) == SS(i,1))
              break
            else
              SS(i,count33+1)=F(j,3);
              m=SS(i,count33+1);
              count33=count33+1;
              break
            endif
            end
          if (i == F(j,2)  && m == F(j,3))
            SS(i,count33)=F(j,3);
            if (F(j,1) == SS(i,1))
              break
            else
              SS(i,count33+1)=F(j,1);
              m=SS(i,count33+1);
              count33=count33+1;
              break
            endif
            end
          if (i == F(j,3)  && m == F(j,1))
            SS(i,count33)=F(j,1);
            if (F(j,2) == SS(i,1))
              break
            else
              SS(i,count33+1)=F(j,2);
              m=SS(i,count33+1);
              count33=count33+1;
              break  
            endif
          end
        end 
       if (SS(i,1)==0)
          m=Sc(i,2);
       endif 
      Lx(i)=count33;
    end
    
    if(L(i)==Lx(i))
      idx(i)=1;
    else
      temp=temp+1;
    end
    endwhile       
  end
  
  printf('Parte 1');

  %Calculando o vetor normal no vértice 
  
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
    x1=Q(:,1); %vetor normal no vértice
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

%Cálculo da Curvatura Gaussiana


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


%Cálculo da Curvatura Média

sumHn3=zeros(nverts,1);

for i=1:nverts
  vp=zeros(2,L(i));
  countp=1;
  vp(1,countp)= 0;
  vp(2,countp)= 0;
  for j=1:LL(i)
    
    if (F(Face(i,j),1) == i)
      vj=F(Face(i,j),2);
      vp(2,countp)=F(Face(i,j),3);
      
      for k=1:nverts
        for l=1:LL(k)
          if ((F(Face(k,l),2) == i) && (F(Face(k,l),1) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),3);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),3) == i) && (F(Face(k,l),2) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),1);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),1) == i) && (F(Face(k,l),3) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),2);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          else
           continue;
          endif
        endfor
      endfor
      countp=countp+1;
    endif
      
     if (F(Face(i,j),2) == i)
      vj=F(Face(i,j),3);
      vp(2,countp)=F(Face(i,j),1);
      
      for k=1:nverts
        for l=1:LL(k)
          if ((F(Face(k,l),2) == i) && (F(Face(k,l),1) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),3);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),3) == i) && (F(Face(k,l),2) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),1);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),1) == i) && (F(Face(k,l),3) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),2);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          else
            continue;
          endif
        endfor
      endfor
      countp=countp+1;
    endif
    if (F(Face(i,j),3) == i)
      vj=F(Face(i,j),1);
      vp(2,countp)=F(Face(i,j),2);
      
      for k=1:nverts
        for l=1:LL(k)
          if ((F(Face(k,l),2) == i) && (F(Face(k,l),1) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),3);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),3) == i) && (F(Face(k,l),2) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),1);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),1) == i) && (F(Face(k,l),3) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),2);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          else
            continue;
          endif
        endfor
      endfor
      countp=countp+1;
    endif
  endfor
  
  
  % Curvatura Media 2: Crane
  
  Hd3(i,:) = (1/4*sumHn3(i,:))/Amixed(i,:);

endfor

printf('Parte 3');

% Cálculo das Curvaturas Principais

for i=1:nverts
  k14(i,:)=Hd3(i,:)+sqrt((Hd3(i,:))^2-(Kd2(i,:)));
  k24(i,:)=Hd3(i,:)-sqrt((Hd3(i,:))^2-(Kd2(i,:)));
  if ((Hd3(i,:))^2-(Kd2(i,:))>=1e-10)
    k13(i,:)=k14(i,:);
    k23(i,:)=k24(i,:);
  else
    k13(i,:)=Hd3(i,:);
    k23(i,:)=Hd3(i,:);
  endif
  
endfor

%Cálculo dos Possíveis Pontos Umbílicos


  Positivo3=[];
  Negativo3=[];
  count5551=1;
  count5651=1;
  for i=1:nverts  
    if((Hd3(i,:))^2-(Kd2(i,:))>1e-10)
      Positivo3(count5551,:)=V(i,:);
      count5551=count5551+1;
    end
    if ((Hd3(i,:))^2-(Kd2(i,:))<-1e-10)
      Negativo3(count5651,:)=V(i,:);
      count5651=count5651+1;
    end
  endfor
 
    count5721=1;
  for i=1:nfaces
    if ((((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1)))>1e-10 && ((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1)))<-1e-10))
      t2=((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1)))/(((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1)))-((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1))));
      MMw1(count5721,:)=V(F(i,1),:)+t2*(V(F(i,2),:)-V(F(i,1),:));
      PMM1(count5721,1)=F(i,1);
      PMM1(count5721,2)=F(i,2);
      count5721=count5721+1;
    end
    if ((((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1)))>1e-10 && ((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1)))<-1e-10))
      t2=((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1)))/(((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1)))-((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1))));
      MMw1(count5721,:)=V(F(i,2),:)+t2*(V(F(i,3),:)-V(F(i,2),:));
      PMM1(count5721,1)=F(i,2);
      PMM1(count5721,2)=F(i,3);
      count5721=count5721+1;
    end
    if ((((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1)))>1e-10 && ((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1)))<-1e-10) )
      t2=((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1)))/(((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1)))-((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1))));
      MMw1(count5721,:)=V(F(i,3),:)+t2*(V(F(i,1),:)-V(F(i,3),:));
      PMM1(count5721,1)=F(i,3);
      PMM1(count5721,2)=F(i,1);
      count5721=count5721+1;
    endif
  end
  

  
  printf('Parte 4');
  
%Cálculo da Curvatura Normal na direção da aresta vv_j
%Cálculo do vetor unitário no plano tangente da aresta vv_j em R^3

kn=[];
d=[];
for i=1:nverts
  for j=L(i):-1:1
    kn(i,j)=2*dot(V(i,:)-V(SS(i,j),:),Nvt(i,:))/(norm(V(i,:)-V(SS(i,j),:)))^2;
    d(i,3*(j-1)+1:3*(j-1)+3)=((V(SS(i,j),:)-V(i,:))-((dot(V(SS(i,j),:)-V(i,:),Nvt(i,:)))*Nvt(i,:)))/norm((V(SS(i,j),:)-V(i,:))-((dot(V(SS(i,j),:)-V(i,:),Nvt(i,:)))*Nvt(i,:)));
  endfor
endfor

%Cálculo do vetor unitário no plano tangente da aresta vv_j em R^2
ujvj=zeros(nverts,2*max(L)+2);
  for i=1:nverts
    count5=1;
    for j=1:3:3*L(i)
      ujvj(i,count5)=d(i,j)*b(i,1)+d(i,j+1)*b(i,2)+d(i,j+2)*b(i,3);
      count5=count5+1;
      ujvj(i,count5)=d(i,j)*c(i,1)+d(i,j+1)*c(i,2)+d(i,j+2)*c(i,3);
      count5=count5+1;
    end   
  end
  
% Cálculo da matriz tensor de curvatura

  P=zeros(nverts,3);
  for i=1:nverts
    count6=1;
    A=zeros(L(i)+1,3);
    B=zeros(L(i)+1,1);
    for j=1:2:2*L(i)
      A(count6,:)=[(ujvj(i,j))^2 ujvj(i,j)*ujvj(i,j+1) (ujvj(i,j+1))^2];
      B(count6,:)=[kn(i,j-(count6-1))];
      count6=count6+1;
    end
    A(count6,:)=[1 0 1];
    B(count6,:)=[-2*Hd3(i)];
    P(i,:)=(A\B)';  %matriz tensor de curvatura
  end

  %Obtenção das direções principais em R^2
  for i=1:nverts
    BB(1,1)=P(i,1);
    BB(1,2)=P(i,2);
    BB(2,1)=P(i,2);
    BB(2,2)=P(i,3);
    [Vb,D] = eig(BB); %Obtenção das direções principais em R^2
    AV(i,1)=Vb(1,1);
    AV(i,2)=Vb(2,1);
    AV(i,3)=Vb(1,2);
    AV(i,4)=Vb(2,2);
    k1D(i,1)=D(1,1);
    k2D(i,1)=D(2,2);
  end

%Obtenção das direções principais em R^3 
   for i=1:nverts
    ttt=AV(i,1)*b(i,:)+AV(i,2)*c(i,:);
    qqq=AV(i,3)*b(i,:)+AV(i,4)*c(i,:);
    AV2(i,1)=ttt(1,1);
    AV2(i,2)=ttt(1,2);
    AV2(i,3)=ttt(1,3);
    AV2(i,4)=qqq(1,1);
    AV2(i,5)=qqq(1,2);
    AV2(i,6)=qqq(1,3);
  end
  
  
%Encontrando o cruzamento das direções principais com a estrela  
    
cruzat11=zeros(1,nverts); %u1 do cruzamento de t1 com a estrela
cruzat12=zeros(1,nverts); %u2 do cruzamento de t1 com a estrela
cruzamt11=zeros(1,nverts); %u1 do cruzamento de t2 com a estrela
cruzamt12=zeros(1,nverts); %u2 do cruzamento de t2 com a estrela

cruzat21=zeros(1,nverts); %u1 do cruzamento de -t1 com a estrela
cruzat22=zeros(1,nverts); %u2 do cruzamento de -t1 com a estrela
cruzamt21=zeros(1,nverts); %u1 do cruzamento de -t2 com a estrela
cruzamt22=zeros(1,nverts); %u2 do cruzamento de -t2 com a estrela

xcruz1=zeros(nverts,3);
xcruz2=zeros(nverts,3);

for i=1:nverts
 for jj=L(i)-1:-1:1
      [cruz1,xc1]=EncontraCruzamento(V(SS(i,jj),:),V(SS(i,jj+1),:),V(i,:),AV2(i,1:3));
      if(cruz1==1)
        cruzat11(i)=SS(i,jj);
        cruzat12(i)=SS(i,jj+1);
        tt1(i)=xc1(1,1); %t>0  do cruzamento de t1 com a estrela
        %ponto do cruzamento
        xcruz1(i,:)=V(cruzat11(i),:)+tt1(i)*(V(cruzat12(i),:)-V(cruzat11(i),:));
      end
      [cruz2,xc2]=EncontraCruzamento(V(SS(i,jj),:),V(SS(i,jj+1),:),V(i,:),AV2(i,4:6));
      if(cruz2==1)
        cruzamt11(i)=SS(i,jj);
        cruzamt12(i)=SS(i,jj+1);
        tt2(i)=xc2(1,1);  %t>0  do cruzamento de t2 com a estrela
        %ponto do cruzamento
        xcruz2(i,:)=V(cruzamt11(i),:)+tt2(i)*(V(cruzamt12(i),:)-V(cruzamt11(i),:));
      end
    end
    
    %verificação se o cruzamento for entre u6 e u1
    if(cruzat11(i)==0)
      [cruz1,xc1]=EncontraCruzamento(V(SS(i,L(i)),:),V(SS(i,1),:),V(i,:),AV2(i,1:3));
      if(cruz1==1)
        cruzat11(i)=SS(i,L(i));
        cruzat12(i)=SS(i,1);
        tt1(i)=xc1(1,1); %t>0  do cruzamento de t1 com a estrela
        %ponto do cruzamento
        xcruz1(i,:)=V(cruzat11(i),:)+tt1(i)*(V(cruzat12(i),:)-V(cruzat11(i),:));
      end
    end
    if(cruzamt11(i)==0)
      [cruz2,xc2]=EncontraCruzamento(V(SS(i,L(i)),:),V(SS(i,1),:),V(i,:),AV2(i,4:6));
      if(cruz2==1)
        cruzamt11(i)=SS(i,L(i));
        cruzamt12(i)=SS(i,1);
        tt2(i)=xc2(1,1);
        xcruz2(i,:)=V(cruzamt11(i),:)+tt2(i)*(V(cruzamt12(i),:)-V(cruzamt11(i),:));
      end
    end
  end
  
  xcruz21=zeros(nverts,3);
  xcruz22=zeros(nverts,3);
  for i=1:nverts
 for jj=L(i)-1:-1:1
      [cruz1,xc1]=EncontraCruzamento(V(SS(i,jj),:),V(SS(i,jj+1),:),V(i,:),-AV2(i,1:3));
      if(cruz1==1)
        cruzat21(i)=SS(i,jj);
        cruzat22(i)=SS(i,jj+1);
        tt21(i)=xc1(1,1); %t>0  do cruzamento de -t1 com a estrela
        %ponto do cruzamento
        xcruz21(i,:)=V(cruzat21(i),:)+tt21(i)*(V(cruzat22(i),:)-V(cruzat21(i),:));
      end
      [cruz2,xc2]=EncontraCruzamento(V(SS(i,jj),:),V(SS(i,jj+1),:),V(i,:),-AV2(i,4:6));
      if(cruz2==1)
        cruzamt21(i)=SS(i,jj);
        cruzamt22(i)=SS(i,jj+1);
        tt22(i)=xc2(1,1);  %t>0  do cruzamento de -t2 com a estrela
        %ponto do cruzamento
        xcruz223(i,:)=[ujvj(i,jj) ujvj(i,jj+1)]+tt22(i)*([ujvj(i,jj+2) ujvj(i,jj+3)]-[ujvj(i,jj) ujvj(i,jj+1)]);
        xcruz22(i,:)=V(cruzamt21(i),:)+tt22(i)*(V(cruzamt22(i),:)-V(cruzamt21(i),:));
      end
    end
    
    %verificação se o cruzamento for entre u6 e u1
    if(cruzat21(i)==0)
      [cruz1,xc1]=EncontraCruzamento(V(SS(i,L(i)),:),V(SS(i,1),:),V(i,:),-AV2(i,1:3));
      if(cruz1==1)
        cruzat21(i)=SS(i,L(i));
        cruzat22(i)=SS(i,1);
        tt21(i)=xc1(1,1);
        xcruz21(i,:)=V(cruzat21(i),:)+tt21(i)*(V(cruzat22(i),:)-V(cruzat21(i),:));
      end
    end
    if(cruzamt21(i)==0)
      [cruz2,xc2]=EncontraCruzamento(V(SS(i,L(i)),:),V(SS(i,1),:),V(i,:),-AV2(i,4:6));
      if(cruz2==1)
        cruzamt21(i)=SS(i,L(i));
        cruzamt22(i)=SS(i,1);
        tt22(i)=xc2(1,1);
        xcruz22(i,:)=V(cruzamt21(i),:)+tt22(i)*(V(cruzamt22(i),:)-V(cruzamt21(i),:));
      end
    end
  end


  
  printf('Parte 5');
  
%Cálculo do valor da curvatura principal nos pontos de cruzamento
  
  kxcruz1=zeros(nverts,1);
  kxcruz2=zeros(nverts,1);
  kxcruz21=zeros(nverts,1);
  kxcruz22=zeros(nverts,1);
  
  for i=1:nverts
    if(cruzat11(i) == 0 || cruzat12(i) == 0)
      continue
    else
      kxcruz1(i,1)=k13(cruzat11(i))+tt1(i)*(k13(cruzat12(i))-k13(cruzat11(i)));
    end
    if(cruzamt11(i) == 0 || cruzamt12(i) == 0)
      continue
    else  
      kxcruz2(i,1)=k23(cruzamt11(i))+tt2(i)*(k23(cruzamt12(i))-k23(cruzamt11(i)));
    end
    if(cruzat21(i) == 0 || cruzat22(i) == 0)
      continue
    else      
      kxcruz21(i,1)=k13(cruzat21(i))+tt21(i)*(k13(cruzat22(i))-k13(cruzat21(i)));
    end
    if(cruzamt21(i) == 0 || cruzamt22(i) == 0)
      continue
    else
      kxcruz22(i,1)=k23(cruzamt21(i))+tt22(i)*(k23(cruzamt22(i))-k23(cruzamt21(i)));
    endif
  endfor
  
% Cálculo da Derivada Direcional de k1 na direção e1
% Cálculo da Derivada Direcional de k2 na direção e2
  
  Ddk1e1=zeros(nverts,1);
  Ddk2e2=zeros(nverts,1);
  Positivo=[];
  Negativo=[];
  Positivo2=[];
  Negativo2=[];
  count55=1;
  count56=1;
  count555=1;
  count565=1;
  MM=[];
  MMw=[];
  for i=1:nverts
    if(abs(xcruz1(i,:)-xcruz21(i,:)) < 1e-10 || abs(xcruz2(i,:)-xcruz22(i,:)) < 1e-10)
      continue
    else
      Ddk1e1(i,1)=(kxcruz1(i,1)-kxcruz21(i,1))/(dot(xcruz1(i,:)-xcruz21(i,:),[AV2(i,1) AV2(i,2) AV2(i,3)]/norm([AV2(i,1) AV2(i,2) AV2(i,3)])));
      Ddk2e2(i,1)=(kxcruz2(i,1)-kxcruz22(i,1))/(dot(xcruz2(i,:)-xcruz22(i,:),[AV2(i,4) AV2(i,5) AV2(i,6)]/norm([AV2(i,4) AV2(i,5) AV2(i,6)])));
    end
    if(Ddk1e1(i,1)>1e-10)
      Positivo(count55,:)=V(i,:);
      count55=count55+1;
    end
    if (Ddk1e1(i,1)<-1e-10)
      Negativo(count56,:)=V(i,:);
      count56=count56+1;
    end
    if(Ddk2e2(i,1)>1e-10)
      Positivo2(count555,:)=V(i,:);
      count555=count555+1;
    end
    if (Ddk2e2(i,1)<-1e-10)
      Negativo2(count565,:)=V(i,:);
      count565=count565+1;
    end
  endfor
  
% Cálculo dos pontos ridges azuis

    count57=1;
    count572=1;
  for i=1:nfaces
    if ((Ddk1e1(F(i,1),1)>1e-10 && Ddk1e1(F(i,2),1)<-1e-10))
      t=Ddk1e1(F(i,1),1)/(Ddk1e1(F(i,1),1)-Ddk1e1(F(i,2),1));
      MM(count57,:)=V(F(i,1),:)+t*(V(F(i,2),:)-V(F(i,1),:));
      PMM(count57,1)=F(i,1);
      PMM(count57,2)=F(i,2);
      count57=count57+1;
      end
    if ((Ddk1e1(F(i,2),1)>1e-10 && Ddk1e1(F(i,3),1)<-1e-10))
      t=Ddk1e1(F(i,2),1)/(Ddk1e1(F(i,2),1)-Ddk1e1(F(i,3),1));
      MM(count57,:)=V(F(i,2),:)+t*(V(F(i,3),:)-V(F(i,2),:));
      PMM(count57,1)=F(i,2);
      PMM(count57,2)=F(i,3);
      count57=count57+1;
    end
    if ((Ddk1e1(F(i,3),1)>1e-10 && Ddk1e1(F(i,1),1)<-1e-10))
      t=Ddk1e1(F(i,3),1)/(Ddk1e1(F(i,3),1)-Ddk1e1(F(i,1),1));
      MM(count57,:)=V(F(i,3),:)+t*(V(F(i,1),:)-V(F(i,3),:));
      PMM(count57,1)=F(i,3);
      PMM(count57,2)=F(i,1);
      count57=count57+1;
    end
  if ((Ddk1e1(F(i,2),1)>1e-10 && Ddk1e1(F(i,1),1)<-1e-10))
      t=Ddk1e1(F(i,1),1)/(Ddk1e1(F(i,1),1)-Ddk1e1(F(i,2),1));
      MM(count57,:)=V(F(i,1),:)+t*(V(F(i,2),:)-V(F(i,1),:));
      PMM(count57,1)=F(i,1);
      PMM(count57,2)=F(i,2);
      count57=count57+1;
    end
  if ((Ddk1e1(F(i,3),1)>1e-10 && Ddk1e1(F(i,2),1)<-1e-10))
      t=Ddk1e1(F(i,2),1)/(Ddk1e1(F(i,2),1)-Ddk1e1(F(i,3),1));
      MM(count57,:)=V(F(i,2),:)+t*(V(F(i,3),:)-V(F(i,2),:));
      PMM(count57,1)=F(i,2);
      PMM(count57,2)=F(i,3);
      count57=count57+1;
    end
  if ((Ddk1e1(F(i,1),1)>1e-10 && Ddk1e1(F(i,3),1)<-1e-10))
      t=Ddk1e1(F(i,3),1)/(Ddk1e1(F(i,3),1)-Ddk1e1(F(i,1),1));
      MM(count57,:)=V(F(i,3),:)+t*(V(F(i,1),:)-V(F(i,3),:));
      PMM(count57,1)=F(i,3);
      PMM(count57,2)=F(i,1);
      count57=count57+1;
    endif
  end
  
% Cálculo dos pontos ridges vermelhos
  
  for i=1:nfaces
    if ((Ddk2e2(F(i,1),1)>1e-10 && Ddk2e2(F(i,2),1)<-1e-10))
      t2=Ddk2e2(F(i,1),1)/(Ddk2e2(F(i,1),1)-Ddk2e2(F(i,2),1));
      MMw(count572,:)=V(F(i,1),:)+t2*(V(F(i,2),:)-V(F(i,1),:));
      PMMw(count572,1)=F(i,1);
      PMMw(count572,2)=F(i,2);
      count572=count572+1;
    end
    if ((Ddk2e2(F(i,2),1)>1e-10 && Ddk2e2(F(i,3),1)<-1e-10))
      t2=Ddk2e2(F(i,2),1)/(Ddk2e2(F(i,2),1)-Ddk2e2(F(i,3),1));
      MMw(count572,:)=V(F(i,2),:)+t2*(V(F(i,3),:)-V(F(i,2),:));
      PMMw(count572,1)=F(i,2);
      PMMw(count572,2)=F(i,3);
      count572=count572+1;
    end
    if ((Ddk2e2(F(i,3),1)>1e-10 && Ddk2e2(F(i,1),1)<-1e-10))
      t2=Ddk2e2(F(i,3),1)/(Ddk2e2(F(i,3),1)-Ddk2e2(F(i,1),1));
      MMw(count572,:)=V(F(i,3),:)+t2*(V(F(i,1),:)-V(F(i,3),:));
      PMMw(count572,1)=F(i,3);
      PMMw(count572,2)=F(i,1);
      count572=count572+1;
    endif
    if ((Ddk2e2(F(i,2),1)>1e-10 && Ddk2e2(F(i,1),1)<-1e-10))
      t2=Ddk2e2(F(i,1),1)/(Ddk2e2(F(i,1),1)-Ddk2e2(F(i,2),1));
      MMw(count572,:)=V(F(i,1),:)+t2*(V(F(i,2),:)-V(F(i,1),:));
      PMMw(count572,1)=F(i,1);
      PMMw(count572,2)=F(i,2);
      count572=count572+1;
    end
    if ((Ddk2e2(F(i,3),1)>1e-10 && Ddk2e2(F(i,2),1)<-1e-10))
      t2=Ddk2e2(F(i,2),1)/(Ddk2e2(F(i,2),1)-Ddk2e2(F(i,3),1));
      MMw(count572,:)=V(F(i,2),:)+t2*(V(F(i,3),:)-V(F(i,2),:));
      PMMw(count572,1)=F(i,2);
      PMMw(count572,2)=F(i,3);
      count572=count572+1;
    end
    if ((Ddk2e2(F(i,1),1)>1e-10 && Ddk2e2(F(i,3),1)<-1e-10))
      t2=Ddk2e2(F(i,3),1)/(Ddk2e2(F(i,3),1)-Ddk2e2(F(i,1),1));
      MMw(count572,:)=V(F(i,3),:)+t2*(V(F(i,1),:)-V(F(i,3),:));
      PMMw(count572,1)=F(i,3);
      PMMw(count572,2)=F(i,1);
      count572=count572+1;
    endif
  end
  
if (size(MM,1)==0 || size(MMw,1)==0)
  fprintf('\n NÃO HÁ CURVAS RIDGES!! \n');
else
  
%Limpando pontos ridges repetidos

MM3=[MM PMM];
MM3w=[MMw PMMw];

MM3=unique(MM3,'rows');
MM3w=unique(MM3w,'rows');

for i=size(MM3,1):-1:2
  if (MM3(i,5)==MM3(i-1,4) && MM3(i,4)==MM3(i-1,5))
    MM3=deletarLinha(MM3,i);
  endif
endfor

for i=size(MM3w,1):-1:2
  if (MM3w(i,5)==MM3w(i-1,4) && MM3w(i,4)==MM3w(i-1,5))
    MM3w=deletarLinha(MM3w,i);
  endif
endfor

%%%%Continuação da obtenção dos pontos umbílicos

MM=MM3(:,1:3);
PMM=MM3(:,4:5);
MMw=MM3w(:,1:3);
PMMw=MM3w(:,4:5);

  
  for i=1:size(PMM1,1)
    if (norm(MMw1(i,:)-V(PMM1(i,1),:))<norm(MMw1(i,:)-V(PMM1(i,2),:)))
      MMw1(i,:)=V(PMM1(i,1),:);
    else
      MMw1(i,:)=V(PMM1(i,2),:);
    endif
  endfor  
  

MMcopia=MM;
MMwcopia=MMw;
  
  for i=1:size(MM)
    for j=1:size(MMw1)
      if((PMM(i,1)==PMM1(j,1) & PMM(i,2)==PMM1(j,2)) || (PMM(i,1)==PMM1(j,2) & PMM(i,2)==PMM1(j,1)))
        MMcopia(i,:)=MMw1(j,:);
      end
    endfor
  endfor
  
  for i=1:size(MMw)
    for j=1:size(MMw1)
      if((PMMw(i,1)==PMM1(j,1) & PMMw(i,2)==PMM1(j,2)) || (PMMw(i,1)==PMM1(j,2) & PMMw(i,2)==PMM1(j,1)))
        MMwcopia(i,:)=MMw1(j,:);
      end
    endfor
  endfor
  
  
uuw=0;
uuwx=0;
for i=1:size(MM,1)
    for j=1:size(MMw,1)
        if((MMcopia(i,1)==MMwcopia(j,1) && MMcopia(i,2)==MMwcopia(j,2)))
          posumbw(uuw+1,:)=MMcopia(i,:);
          uuw=uuw+1;
          for k=1:size(MMw1)
            if(posumbw(uuw,:)==MMw1(k,:))
              umb(uuwx+1,:)=posumbw(uuw,:);
              uuwx=uuwx+1;
            end
          endfor
        end
    end
endfor

  umb=unique(umb,"rows");
  totalumb=0;
  for i=1:nverts
    for j=1:size(umb,1)
      if (V(i,:)==umb(j,:))
        Vumb(totalumb+1,:)=i;
        totalumb=totalumb+1;
      endif
    endfor
  endfor
 

%Encontrando os "pontos umbílicos" que se conectam por alguma aresta  
 
 
ligado=[];
pppp=1;    
for i=1:size(Vumb,1)
  ppp=1;
  ligado(i,ppp)=Vumb(i);
    for k=1:size(Vumb,1)
      for j=1:size(SS,2)
          if (Vumb(k)==SS(Vumb(i),j))
            ligado(i,ppp+1)=Vumb(k);
            ppp=ppp+1;
          endif          
      endfor
    endfor
    ligadoind(i,1)=ppp;
endfor   

ppppp=0;
for i=1:size(ligado,1)
  if (ligado(i,1)!=0)
        ppppp=ppppp+1;
        for k=1:ligadoind(ppppp,1);
          for l=1:size(Vumb,1)
            if (Vumb(l,1)==ligado(ppppp,k))
              xt=l;
              break;
            endif
          endfor
          for m=1:ligadoind(xt,1)
            kkkkk=0;
          for j=1:ligadoind(i,1);
            if (ligado(xt,m)==ligado(i,j))
              continue;
            else
              kkkkk=kkkkk+1;
            endif
          endfor
          kkkkk;
          if (kkkkk==ligadoind(i,1))
            ligado(i,ligadoind(i,1)+1)=ligado(xt,m);
            ligadoind(i,1)=ligadoind(i,1)+1;

          endif
          end
      end
        
    endif    
end

ppppp=size(ligado,1)+1;
for i=size(ligado,1):-1:1
  if (ligado(i,1)!=0)
        ppppp=ppppp-1;
        for k=1:ligadoind(ppppp,1);
          for l=1:size(Vumb,1)
            if (Vumb(l,1)==ligado(ppppp,k))
              xt=l;
              break;
            endif
          endfor
          for m=1:ligadoind(xt,1)
            kkkkk=0;
          for j=1:ligadoind(i,1);
            if (ligado(xt,m)==ligado(i,j))
              continue;
            else
              kkkkk=kkkkk+1;
            endif
          endfor
          kkkkk;
          if (kkkkk==ligadoind(i,1))
            ligado(i,ligadoind(i,1)+1)=ligado(xt,m);
            ligadoind(i,1)=ligadoind(i,1)+1;
          endif
          end
      end
    endif    
end   
 
 ligadosr=[];
  for i=1:size(ligado,1)
    C=ligado(i,1:ligadoind(i,1))';
    Cc=unique(C,"rows");
    for j=1:size(Cc,1)
      ligadosr(i,j)=Cc(j,1);
    endfor
  endfor
ligadosr=unique(ligadosr,"rows");
numlig=zeros(size(ligadosr,1),1);
for i=1:size(ligadosr,1)
  for j=1:size(ligadosr,2)
    if (ligadosr(i,j) != 0)
    numlig(i,1)=numlig(i,1)+1;
    endif
  endfor
endfor

%Encontrando o baricentro dos pontos umbílicos conectados

for i=1:size(ligadosr,1)
  if (numlig(i)==1)
    Pumb(i,:)=V(ligadosr(i,1),:);
  else
  soma=[0 0 0];
  for j=1:numlig(i)
    soma=soma+V(ligadosr(i,j),:);
    Pumb(i,:)=soma/numlig(i);
  endfor
  endif
endfor

%Encontrando a face que o baricentro pertence

for i=1:size(ligadosr,1)
  if (numlig(i)==1)
    aa(i)=ligadosr(i,1);
  endif
  if (numlig(i)==2)
    aa(i)=ligadosr(i,1);
    aaa(i)=ligadosr(i,2);
  endif
endfor


for i=1:size(ligadosr,1)
  if (numlig(i)>=3)
    for k=2:numlig(i)
    for j=1:nfaces
      if ((ligadosr(i,1)==F(j,1) & ligadosr(i,k)==F(j,2)) || (ligadosr(i,1)==F(j,2) & ligadosr(i,k)==F(j,3)) || (ligadosr(i,1)==F(j,3) & ligadosr(i,k)==F(j,1)))
        aa(i)=ligadosr(i,1);
        aaa(i)=ligadosr(i,k);
        aainicial(i)=aa(i);
        aaainicial(i)=aaa(i);
        break      
      endif
      endfor
    endfor
  endif
endfor

Borda = [11, 339, 298, 333, 285, 252, 390, 357, 455, 324, 362, 289, 398, 366, 380, 379, 401, 378, 153, 149, 177, 150, 151, 137, 173, 59, 133, 94, 235, 128, 163, 22, 55, 104, 68, 110];

for i=1:size(ligadosr,1)
  if (numlig(i)>=3)
    for j=1:nfaces
      if (aa(i)==F(j,1) && aaa(i)==F(j,2))
        aaaa(i)=F(j,3);
        break;
      elseif (aa(i)==F(j,2) && aaa(i)==F(j,3))
        aaaa(i)=F(j,1);
        break;
      elseif (aa(i)==F(j,3) && aaa(i)==F(j,1))
        aaaa(i)=F(j,2);
        break;
      endif
    endfor
    P1=[V(aa(i),1) V(aaa(i),1) V(aaaa(i),1) ;  V(aa(i),2) V(aaa(i),2) V(aaaa(i),2) ; 1 1 1];
    VV=[Pumb(i,1); Pumb(i,2); 1];
    XX=P1\VV;
    tol=-1e-15;
    contador=0;
    aaant2(i)=0;
    aaaant2(i)=0;
    aaaaant2(i)=0;
    while (XX(1)<0+tol || XX(2)<0+tol || XX(3)<0+tol)
      conf=0;
      contador=contador+1;
      XX=P1\VV;
      
      if (contador>=2)
        aaant2(i)=aaant(i);
        aaaant2(i)=aaaant(i);
        aaaaant2(i)=aaaaant(i);
      endif
      
      aaant(i)=aa(i);
      aaaant(i)=aaa(i);
      aaaaant(i)=aaaa(i);
      
      
      
      if (min(XX)==XX(1))
        for j=1:nfaces
          if (aaa(i)==F(j,1) && aaaa(i)==F(j,2) && aa(i)!=F(j,3))
            aa(i)=F(j,3);
            break;
          elseif (aaa(i)==F(j,2) && aaaa(i)==F(j,3) && aa(i)!=F(j,1))
            aa(i)=F(j,1);
            break;
          elseif (aaa(i)==F(j,3) && aaaa(i)==F(j,1) && aa(i)!=F(j,2))
            aa(i)=F(j,2);
            break;
          elseif (aaa(i)==F(j,2) && aaaa(i)==F(j,1) && aa(i)!=F(j,3))
            aa(i)=F(j,3);
            break;
          elseif (aaa(i)==F(j,3) && aaaa(i)==F(j,2) && aa(i)!=F(j,1))
            aa(i)=F(j,1);
            break;
          elseif (aaa(i)==F(j,1) && aaaa(i)==F(j,3) && aa(i)!=F(j,2))
            aa(i)=F(j,2);
            break;
          endif
        endfor
      elseif (min(XX)==XX(2))
        for j=1:nfaces
          if (aa(i)==F(j,1) && aaaa(i)==F(j,2) && aaa(i)!=F(j,3))
            aaa(i)=F(j,3);
            break;
          elseif (aa(i)==F(j,2) && aaaa(i)==F(j,3) && aaa(i)!=F(j,1))
            aaa(i)=F(j,1);
            break;
          elseif (aa(i)==F(j,3) && aaaa(i)==F(j,1) && aaa(i)!=F(j,2))
            aaa(i)=F(j,2);
            break;
          elseif (aa(i)==F(j,2) && aaaa(i)==F(j,1) && aaa(i)!=F(j,3))
            aaa(i)=F(j,3);
            break;
          elseif (aa(i)==F(j,3) && aaaa(i)==F(j,2) && aaa(i)!=F(j,1))
            aaa(i)=F(j,1);
            break;
          elseif (aa(i)==F(j,1) && aaaa(i)==F(j,3) && aaa(i)!=F(j,2))
            aaa(i)=F(j,2);
            break;
          endif
        endfor
      elseif (min(XX)==XX(3))
        for j=1:nfaces
          if (aa(i)==F(j,1) && aaa(i)==F(j,2) && aaaa(i)!=F(j,3))
            aaaa(i)=F(j,3);
            break;
          elseif (aa(i)==F(j,2) && aaa(i)==F(j,3) && aaaa(i)!=F(j,1))
            aaaa(i)=F(j,1);
            break;
          elseif (aa(i)==F(j,3) && aaa(i)==F(j,1) && aaaa(i)!=F(j,2))
            aaaa(i)=F(j,2);
            break;
          elseif (aa(i)==F(j,2) && aaa(i)==F(j,1) && aaaa(i)!=F(j,3))
            aaaa(i)=F(j,3);
            break;
          elseif (aa(i)==F(j,3) && aaa(i)==F(j,2) && aaaa(i)!=F(j,1))
            aaaa(i)=F(j,1);
            break;
          elseif (aa(i)==F(j,1) && aaa(i)==F(j,3) && aaaa(i)!=F(j,2))
            aaaa(i)=F(j,2);
            break;
          endif
        endfor
      endif
     
      
      for p=1:size(Borda,2)
        for r=1:size(Borda,2)
          if (((aa(i)==Borda(1,r) && aaa(i)==Borda(1,p)) || (aa(i)==Borda(1,r) && aaaa(i)==Borda(1,p)) || (aaa(i)==Borda(1,r) && aaaa(i)==Borda(1,p))) && (aaant(i)==aa(i) && aaaant(i)==aaa(i) && aaaaant(i)==aaaa(i)))
            conf=1;
          endif
        endfor
      endfor
      
      if (conf==0)
        P1=[V(aa(i),1) V(aaa(i),1) V(aaaa(i),1) ;  V(aa(i),2) V(aaa(i),2) V(aaaa(i),2) ; 1 1 1];
        VV=[Pumb(i,1); Pumb(i,2); 1];
        XX=P1\VV;
      else
        aa(i)=aaainicial(i);
        aaa(i)=aainicial(i);
        for j=1:nfaces
      if (aa(i)==F(j,1) && aaa(i)==F(j,2))
        aaaa(i)=F(j,3);
        break;
      elseif (aa(i)==F(j,2) && aaa(i)==F(j,3))
        aaaa(i)=F(j,1);
        break;
      elseif (aa(i)==F(j,3) && aaa(i)==F(j,1))
        aaaa(i)=F(j,2);
        break;
      endif
    endfor
    P1=[V(aa(i),1) V(aaa(i),1) V(aaaa(i),1) ;  V(aa(i),2) V(aaa(i),2) V(aaaa(i),2) ; 1 1 1];
    VV=[Pumb(i,1); Pumb(i,2); 1];
    XX=P1\VV;
    endif
        
    if (aaant2(i)==aa(i) && aaaant2(i)==aaa(i) && aaaaant2(i)==aaaa(i))
      for k=2:numlig(i)
      for j=1:nfaces
        if (((ligadosr(i,1)==F(j,1) & ligadosr(i,k)==F(j,2)) || (ligadosr(i,1)==F(j,2) & ligadosr(i,k)==F(j,3)) || (ligadosr(i,1)==F(j,3) & ligadosr(i,k)==F(j,1))) && (aainicial(i)!=ligadosr(i,1) || aaainicial(i) !=ligadosr(i,k)))
          aa(i)=ligadosr(i,1);
          aaa(i)=ligadosr(i,k);
          break      
        endif
      endfor
      endfor
      for j=1:nfaces
      if (aa(i)==F(j,1) && aaa(i)==F(j,2))
        aaaa(i)=F(j,3);
        break;
      elseif (aa(i)==F(j,2) && aaa(i)==F(j,3))
        aaaa(i)=F(j,1);
        break;
      elseif (aa(i)==F(j,3) && aaa(i)==F(j,1))
        aaaa(i)=F(j,2);
        break;
      endif
      endfor
      P1=[V(aa(i),1) V(aaa(i),1) V(aaaa(i),1) ;  V(aa(i),2) V(aaa(i),2) V(aaaa(i),2) ; 1 1 1];
      VV=[Pumb(i,1); Pumb(i,2); 1];
      XX=P1\VV
    endif  
    endwhile
  endif
  if (numlig(i)>=3)
    AA(1)=XX(1)*V(aa(i),1)+XX(2)*V(aaa(i),1)+XX(3)*V(aaaa(i),1);
    AA(2)=XX(1)*V(aa(i),2)+XX(2)*V(aaa(i),2)+XX(3)*V(aaaa(i),2);
    AA(3)=XX(1)*V(aa(i),3)+XX(2)*V(aaa(i),3)+XX(3)*V(aaaa(i),3);
    Pumb(i,:)=[AA(1) AA(2) AA(3)];
  endif
endfor

%Determinação dos pontos umbílicos FINAL

Pumbc=Pumb;

%Ponto umbílico sendo o vértice mais próximo do baricentro

for i=1:size(ligadosr,1)
  dist=0;
  if (numlig(i)>=2)
    for j=1:numlig(i)
      dist(1,j)=norm(Pumb(i,:)-V(ligadosr(i,j),:));
    endfor
    [Vumbmin(i),iVumbmin(i)]=min(dist);
    Pumbc(i,:)=V(ligadosr(i,iVumbmin(i)),:);
  endif
endfor


%Fazendo com que as curvas ridges passem pelos pontos umbílicos


for k=1:size(ligadosr,1)
  if (numlig(k)==1)
    for i=1:size(MM)
        if (PMM(i,1)==aa(1,k) || PMM(i,2)==aa(1,k))
          MM(i,:)=Pumb(k,:);
        endif
    endfor
  endif
endfor

for k=1:size(ligadosr,1)
  if (numlig(k)==1)
    for i=1:size(MMw)
        if (PMMw(i,1)==aa(1,k) || PMMw(i,2)==aa(1,k))
          MMw(i,:)=Pumb(k,:);
        endif
    endfor
  endif
endfor



for k=1:size(ligadosr,1)
  if (numlig(k)>=2)
    for i=1:size(MM)
        if (PMM(i,1)==ligadosr(k,iVumbmin(k)) || PMM(i,2)==ligadosr(k,iVumbmin(k)))
          MM(i,:)=Pumbc(k,:);
        endif
    endfor
  endif
endfor

for k=1:size(ligadosr,1)
  if (numlig(k)>=2)
    for i=1:size(MMw)
        if (PMMw(i,1)==ligadosr(k,iVumbmin(k)) || PMMw(i,2)==ligadosr(k,iVumbmin(k)))
          MMw(i,:)=Pumbc(k,:);
        endif
    endfor
  endif
endfor

MM3=[MM PMM];
MM3w=[MMw PMMw];

%Encontrando os vértices que possuem derivada direcional nula
  
ct=1;
ct2=1;
zero=[];
zero2=[];


for i=1:nverts
  if(abs(Ddk1e1(i))<=1e-10)
    zero(ct,1)=i;
    ct=ct+1;
  end
  if(abs(Ddk2e2(i))<=1e-10)
    zero2(ct2,1)=i;
    ct2=ct2+1;
  end
endfor

guardar=[];
for i=1:size(zero,1)
  ccccc=1;
  guardar(i,ccccc)=zero(i,1);
  for j=1:L(i)
    for k=1:size(zero,1)
      if (SS(zero(i),j)==zero(k))
        guardar(i,ccccc+1)=zero(k);
        ccccc=ccccc+1;
        break
        endif
    endfor
    if (ccccc>1)
      break
    endif
  endfor
  idx=0;
  while (idx!=1)
    ccc=ccccc; 
   if (j+1<=L(i))  
    for k=1:size(zero,1)
      if (SS(zero(i),j+1)==zero(k))
      guardar(i,ccccc+1)=zero(k);
      ccccc=ccccc+1;
      j=j+1;
      break
      endif
    endfor
   endif
   if(ccc==ccccc)
      idx=1;
    else
      ccc=ccccc;
    endif
  endwhile
  nguardar(i,1)=ccccc;
end

guardar2=[];
for i=1:size(zero2,1)
  ccccc=1;
  guardar2(i,ccccc)=zero2(i,1);
  for j=1:L(i)
    for k=1:size(zero2,1)
      if (SS(zero2(i),j)==zero2(k))
        guardar2(i,ccccc+1)=zero2(k);
        ccccc=ccccc+1;
        break
      endif
    endfor
    if (ccccc>1)
      break
    endif
  endfor
  idx=0;
  while (idx!=1)
    ccc=ccccc;
    if (j+1<=L(i))  
     for k=1:size(zero2,1)
      if (SS(zero2(i),j+1)==zero2(k))
      guardar2(i,ccccc+1)=zero2(k);
      ccccc=ccccc+1;
      j=j+1;
      break
      endif
     endfor
    endif
    if(ccc==ccccc)
      idx=1;
    else
      ccc=ccccc;
    endif
  endwhile
  nguardar2(i,1)=ccccc;
end



printf('Parte 6');


% Determinação de cada curva ridge (juntando os pontos para
  % formação das curvas)

%1ª Formação das curvas ridge vermelha  
  
MM2w=[];
jafoiw=zeros(size(MM3w)(1),1);
kkw=1;
for iiw=1:size(MM3w)(1)
  if (jafoiw(iiw)==0)
    MM2w(1,3*kkw-2:3*kkw)=MM3w(iiw,1:3);
    PMM2w(1,2*kkw-1:2*kkw)=MM3w(iiw,4:5);
ccw=1;
llw=iiw;
jjjw=0;

jafoiw(iiw)=1;
while (llw != jjjw)
  jjjw=llw;
  for i=1:nfaces
    if (MM3w(llw,4)==F(i,1) && MM3w(llw,5)==F(i,2))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,1) && MM3w(llw,5)==F(i,2))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,2) && MM3w(llw,5)==F(i,3))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)))
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,2) && MM3w(llw,5)==F(i,3))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)))
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3w(llw,4)==F(i,3) && MM3w(llw,5)==F(i,2))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)))
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3w(llw,4)==F(i,2) && MM3w(llw,5)==F(i,1))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
  end

      if (MM3w(llw,4)==F(i,3) && MM3w(llw,5)==F(i,1))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)))
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,3) && MM3w(llw,5)==F(i,1))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,1) && MM3w(llw,5)==F(i,3))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)))
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3w(llw,4)==F(i,1) && MM3w(llw,5)==F(i,3))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end

endfor
end
pontocurvasw(kkw,1)=ccw;
kkw=kkw+1;
endif
endfor

%1ª Formação das curvas ridge azul 

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
pontocurvas(kk,1)=cc;
kk=kk+1;
endif
endfor






######## Azul ######

%Juntando 2 curvas ridges azuis que o primeiro ponto da curva 1 
%e o primeiro ponto da curva 2 estão na mesma face

idx2=0;
idx3=0;
while(idx2!=1)
pontocurvasax=pontocurvas;
#1
idx=0;
while(idx!=1)
pontocurvasa=pontocurvas;
for tt=1:2:size(PMM2,2)
  ppt=size(PMM2,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMM2,2)
  pptz=size(PMM2,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp)==F(i,2) && PMM2(1,pp+1)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp+1)==F(i,2) && PMM2(1,pp)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp)==F(i,1) && PMM2(1,pp+1)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp+1)==F(i,1) && PMM2(1,pp)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp)==F(i,3) && PMM2(1,pp+1)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp+1)==F(i,3) && PMM2(1,pp)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp)==F(i,2) && PMM2(1,pp+1)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp+1)==F(i,2) && PMM2(1,pp)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp)==F(i,1) && PMM2(1,pp+1)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp+1)==F(i,1) && PMM2(1,pp)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp)==F(i,3) && PMM2(1,pp+1)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp+1)==F(i,3) && PMM2(1,pp)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp)==F(i,1) && PMM2(1,pp+1)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp+1)==F(i,1) && PMM2(1,pp)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp)==F(i,2) && PMM2(1,pp+1)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp+1)==F(i,2) && PMM2(1,pp)==F(i,3)))
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp)==F(i,2) && PMM2(1,pp+1)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp+1)==F(i,2) && PMM2(1,pp)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp)==F(i,3) && PMM2(1,pp+1)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp+1)==F(i,3) && PMM2(1,pp)==F(i,1)))
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp)==F(i,3) && PMM2(1,pp+1)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp+1)==F(i,3) && PMM2(1,pp)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp)==F(i,1) && PMM2(1,pp+1)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp+1)==F(i,1) && PMM2(1,pp)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasa,1)==size(pontocurvas,1))
idx=1;
end
endwhile

printf('Azul 1');

%Juntando 2 curvas ridges azuis que o último ponto da curva 2 
%e o primeiro ponto da curva 1 estão na mesma face
#2
idx=0;
while(idx!=1)
pontocurvasa=pontocurvas;
for tt=1:2:size(PMM2,2)
  ppt=size(PMM2,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMM2,2)
  pptz=size(PMM2,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasa,1)==size(pontocurvas,1))
idx=1;
end
endwhile

printf('Azul 2');

%Juntando 2 curvas ridges azuis que o último ponto da curva 1 
%e o primeiro ponto da curva 2 estão na mesma face
#3
idx=0;
while(idx!=1)
pontocurvasa=pontocurvas;
for pp=1:2:size(PMM2,2)
  ppt=size(PMM2,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMM2,2)
  pptz=size(PMM2,2);
  if (tt>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasa,1)==size(pontocurvas,1))
idx=1;
end
endwhile

printf('Azul 3');

%Fazendo a curva ridge azul, que o último ponto da curva está na 
%mesma face de um vértice que é ridge azul, passar por esse vértice
#4
idx=0;
while(idx!=1)
mx=0;
pontocurvasa=pontocurvas;
for tt=1:2:size(PMM2,2)
for pp=1:size(zero,1)
for i=1:nfaces
  if((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && zero(pp,1)==F(i,3)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && zero(pp,1)==F(i,1)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && zero(pp,1)==F(i,2)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && zero(pp,1)==F(i,3)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && zero(pp,1)==F(i,1)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && zero(pp,1)==F(i,2)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
endif
endfor
endfor
endfor
if(mx==0)
idx=1;
end
endwhile

printf('Azul 4');

%Fazendo com que a curva ridge azul passe pelo ponto umbílico
% Caso em que o ponto umbílico se originou de 2 pontos umbílicos ligados
abc=size(PMM2,2);
for j=1:size(ligadosr,1)
  if (numlig(j) ==2)
    for tt=1:2:abc
      idx=0;
      for i=1:pontocurvas((tt+1)/2)
        if (idx == 1)
          break;
        endif
        if (PMM2(i,tt)==aa(1,j) ||PMM2(i,tt)==aaa(1,j) || PMM2(i,tt+1)==aa(1,j) ||PMM2(i,tt+1)==aaa(1,j) )
          if   ((PMM2(i,tt)==aa(1,j) && PMM2(i,tt+1)==aaa(1,j)) || (PMM2(i,tt+1)==aa(1,j) && PMM2(i,tt)==aaa(1,j)) )
            break
          else
        for k=1:nfaces
          if ((PMM2(i,tt)==F(k,1) && PMM2(i,tt+1)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3)) || (PMM2(i,tt)==F(k,1) && PMM2(i,tt+1)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt)==F(k,2) && PMM2(i,tt+1)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1)) || (PMM2(i,tt)==F(k,2) && PMM2(i,tt+1)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt)==F(k,3) && PMM2(i,tt+1)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2)) || (PMM2(i,tt)==F(k,3) && PMM2(i,tt+1)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt)==F(k,2) && PMM2(i,tt+1)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3)) || (PMM2(i,tt)==F(k,2) && PMM2(i,tt+1)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt)==F(k,1) && PMM2(i,tt+1)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2)) || (PMM2(i,tt)==F(k,1) && PMM2(i,tt+1)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt)==F(k,3) && PMM2(i,tt+1)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1)) || (PMM2(i,tt)==F(k,3) && PMM2(i,tt+1)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
            elseif ((PMM2(i,tt+1)==F(k,1) && PMM2(i,tt)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3)) || (PMM2(i,tt+1)==F(k,1) && PMM2(i,tt)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt+1)==F(k,2) && PMM2(i,tt)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1)) || (PMM2(i,tt+1)==F(k,2) && PMM2(i,tt)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt+1)==F(k,3) && PMM2(i,tt)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2)) || (PMM2(i,tt+1)==F(k,3) && PMM2(i,tt)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt+1)==F(k,2) && PMM2(i,tt)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3)) || (PMM2(i,tt+1)==F(k,2) && PMM2(i,tt)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt+1)==F(k,1) && PMM2(i,tt)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2)) || (PMM2(i,tt+1)==F(k,1) && PMM2(i,tt)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2(i,tt+1)==F(k,3) && PMM2(i,tt)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1)) || (PMM2(i,tt+1)==F(k,3) && PMM2(i,tt)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1)))
            MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvas((tt+1)/2))
            MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
            PMM2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvas((tt+1)/2))
            PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
            end
            pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            MM2=deletarColuna(MM2,3*(tt+1)/2-2);
            PMM2=deletarColuna(PMM2,tt);
            PMM2=deletarColuna(PMM2,tt);
            pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
            idx=1;
            break
          endif
    endfor
  endif
    endif
  endfor
endfor
endif
endfor
% Caso em que o ponto umbílico se originou de 3 ou mais pontos umbílicos ligados

for j=1:size(ligadosr,1)
  if (numlig(j) >=3)
    for tt=1:2:size(PMM2,2)
      idx=0;
      for i=1:pontocurvas((tt+1)/2)-1
        if (idx==1)
            break
          endif
         if (PMM2(i,tt)==aa(1,j) && PMM2(i,tt+1)==aaa(1,j) && ((PMM2(i+1,tt)==aaa(1,j) && PMM2(i+1,tt+1)==aaaa(1,j)) | (PMM2(i+1,tt)==aaaa(1,j) && PMM2(i+1,tt+1)==aaa(1,j)) | (PMM2(i+1,tt)==aa(1,j) && PMM2(i+1,tt+1)==aaaa(1,j)) | (PMM2(i+1,tt)==aaaa(1,j) && PMM2(i+1,tt+1)==aa(1,j))))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          elseif (PMM2(i,tt)==aaa(1,j) && PMM2(i,tt+1)==aa(1,j)&& ((PMM2(i+1,tt)==aaa(1,j) && PMM2(i+1,tt+1)==aaaa(1,j)) | (PMM2(i+1,tt)==aaaa(1,j) && PMM2(i+1,tt+1)==aaa(1,j)) | (PMM2(i+1,tt)==aa(1,j) && PMM2(i+1,tt+1)==aaaa(1,j)) | (PMM2(i+1,tt)==aaaa(1,j) && PMM2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(i,tt)==aa(1,j) && PMM2(i,tt+1)==aaaa(1,j)&& ((PMM2(i+1,tt)==aaa(1,j) && PMM2(i+1,tt+1)==aaaa(1,j)) | (PMM2(i+1,tt)==aaaa(1,j) && PMM2(i+1,tt+1)==aaa(1,j)) | (PMM2(i+1,tt)==aa(1,j) && PMM2(i+1,tt+1)==aaa(1,j)) | (PMM2(i+1,tt)==aaa(1,j) && PMM2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(i,tt)==aaaa(1,j) && PMM2(i,tt+1)==aa(1,j) && ((PMM2(i+1,tt)==aaa(1,j) && PMM2(i+1,tt+1)==aaaa(1,j)) | (PMM2(i+1,tt)==aaaa(1,j) && PMM2(i+1,tt+1)==aaa(1,j)) | (PMM2(i+1,tt)==aa(1,j) && PMM2(i+1,tt+1)==aaa(1,j)) | (PMM2(i+1,tt)==aaa(1,j) && PMM2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(i,tt)==aaaa(1,j) && PMM2(i,tt+1)==aaa(1,j) && ((PMM2(i+1,tt)==aa(1,j) && PMM2(i+1,tt+1)==aaaa(1,j)) | (PMM2(i+1,tt)==aaaa(1,j) && PMM2(i+1,tt+1)==aa(1,j)) | (PMM2(i+1,tt)==aa(1,j) && PMM2(i+1,tt+1)==aaa(1,j)) | (PMM2(i+1,tt)==aaa(1,j) && PMM2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(i,tt)==aaa(1,j) && PMM2(i,tt+1)==aaaa(1,j) && ((PMM2(i+1,tt)==aa(1,j) && PMM2(i+1,tt+1)==aaaa(1,j)) | (PMM2(i+1,tt)==aaaa(1,j) && PMM2(i+1,tt+1)==aa(1,j)) | (PMM2(i+1,tt)==aa(1,j) && PMM2(i+1,tt+1)==aaa(1,j)) | (PMM2(i+1,tt)==aaa(1,j) && PMM2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2(1:i,end+1:end+3)=MM2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvas((tt+1)/2))
              MM2(i+2:pontocurvas((tt+1)/2)+1,end-2:end)=MM2(i+1:pontocurvas((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2(1:i,end+1:end+2)=PMM2(1:i,tt:(tt+1));
              PMM2(i+1,end-1)=nverts+10;
              PMM2(i+1,end)=nverts+10;
              if (i+1<=pontocurvas((tt+1)/2))
              PMM2(i+2:pontocurvas((tt+1)/2)+1,end-1:end)=PMM2(i+1:pontocurvas((tt+1)/2),tt:(tt+1));
              end
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          elseif (PMM2(1,tt)==aa(1,j) && PMM2(1,tt+1)==aaa(1,j) && (PMM2(2,tt)!=aaaa(1,j) | PMM2(2,tt+1)!=aaaa(1,j)))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          elseif (PMM2(1,tt)==aaa(1,j) && PMM2(1,tt+1)==aa(1,j) && (PMM2(2,tt)!=aaaa(1,j) | PMM2(2,tt+1)!=aaaa(1,j)))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(1,tt)==aa(1,j) && PMM2(1,tt+1)==aaaa(1,j) && (PMM2(2,tt)!=aaa(1,j) | PMM2(2,tt+1)!=aaa(1,j)))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(1,tt)==aaaa(1,j) && PMM2(1,tt+1)==aa(1,j) && (PMM2(2,tt)!=aaa(1,j) | PMM2(2,tt+1)!=aaa(1,j)))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(1,tt)==aaaa(1,j) && PMM2(1,tt+1)==aaa(1,j) && (PMM2(2,tt)!=aa(1,j) | PMM2(2,tt+1)!=aa(1,j)))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(1,tt)==aaa(1,j) && PMM2(1,tt+1)==aaaa(1,j) && (PMM2(2,tt)!=aa(1,j) | PMM2(2,tt+1)!=aa(1,j)))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2(1,end+1:end+3)=Pumbc(j,:);
              MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2(1,end+1)=nverts+10;
              PMM2(1,end+1)=nverts+10;
              PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          elseif (PMM2(pontocurvas((tt+1)/2,1),tt)==aa(1,j) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==aaa(1,j))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          elseif (PMM2(pontocurvas((tt+1)/2,1),tt)==aaa(1,j) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==aa(1,j))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(pontocurvas((tt+1)/2,1),tt)==aa(1,j) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==aaaa(1,j))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(pontocurvas((tt+1)/2,1),tt)==aaaa(1,j) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==aa(1,j))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(pontocurvas((tt+1)/2,1),tt)==aaaa(1,j) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==aaa(1,j))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2(pontocurvas((tt+1)/2,1),tt)==aaa(1,j) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==aaaa(1,j))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
          elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2(pontocurvas((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
              PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2(pontocurvas((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              MM2=deletarColuna(MM2,3*(tt+1)/2-2);
              PMM2=deletarColuna(PMM2,tt);
              PMM2=deletarColuna(PMM2,tt);
              pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          endif
      endfor
    endfor
  endif
endfor

printf('Azul 5');

%Fazendo a curva ridge azul, que o primeiro ponto da curva está na 
%mesma face de um vértice que é ridge azul, passar por esse vértice
idx=0;
while(idx!=1)
mx=0;
pontocurvasa=pontocurvas;
for tt=1:2:size(PMM2,2)
if (pontocurvas((tt+1)/2,1)>=2)
for pp=1:size(zero,1)
for i=1:nfaces
  if((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && zero(pp,1)==F(i,3) && zero(pp,1)!=PMM2(2,tt)) )
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && zero(pp,1)==F(i,1) && zero(pp,1)!=PMM2(2,tt)))
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && zero(pp,1)==F(i,2) && zero(pp,1)!=PMM2(2,tt)))
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && zero(pp,1)==F(i,3) && zero(pp,1)!=PMM2(2,tt)))
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && zero(pp,1)==F(i,1) && zero(pp,1)!=PMM2(2,tt)))
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && zero(pp,1)==F(i,2) && zero(pp,1)!=PMM2(2,tt)) )
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
endif
endfor
endfor
endif
endfor
if(mx==0)
idx=1;
endif
endwhile

printf('Azul 6');

%Juntando 2 curvas ridges azuis que o último ponto da curva 1 
%e o primeiro ponto da curva 2 são os mesmos pontos

idx=0;
while(idx!=1)
pontocurvasa=pontocurvas;
for pp=1:2:size(PMM2,2)
  ppt=size(PMM2,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMM2,2)
  pptz=size(PMM2,2);
  if (tt>pptz)
    break
  endif
  if(MM2(pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2)==MM2(1,3*(tt+1)/2-2:3*(tt+1)/2))
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1)-1,end-2:end)=MM2(2:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1)-1,end-1:end)=PMM2(2:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1)-1;
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
endif
endfor
endfor
if(size(pontocurvasa,1)==size(pontocurvas,1))
idx=1;
end
endwhile

%Juntando 2 curvas ridges azuis que o primeiro ponto da curva 1 
%e o primeiro ponto da curva 2 são os mesmos pontos

idx=0;
while(idx!=1)
pontocurvasa=pontocurvas;
for pp=1:2:size(PMM2,2)
  ppt=size(PMM2,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMM2,2)
  pptz=size(PMM2,2);
  if (tt>pptz)
    break
  endif
  if (MM2(1,3*(pp+1)/2-2:3*(pp+1)/2)==MM2(1,3*(tt+1)/2-2:3*(tt+1)/2))
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1)-1,end-2:end)=MM2(2:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1)-1,end-1:end)=PMM2(2:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1)-1;
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
endif
endfor
endfor
if(size(pontocurvasa,1)==size(pontocurvas,1))
idx=1;
end
endwhile

%Verificando se em todos os pontos umbílicos há curvas ridges azuis
%passando por eles.

if (idx2==0 && idx3==0)
  nPumbsem=[];
  MM2c=[];
  pontodegenerado=[];
  count=0;
  ptd=0;
  for i=1:size(Pumbc)
    for j=1:size(pontocurvas)(1)
      MM2(1:pontocurvas(j,1),3*j-2:3*j)==Pumbc(i,:);
      til = ans ;
      if (max(til(:,1)) == 1 && max(til(:,2)) == 1 && max(til(:,3)) == 1)
        if (numlig(i)>=1 && MM2(1,3*j-2:3*j)== MM2(2:pontocurvas(j,1),3*j-2:3*j))
          ptd=ptd+1;
          pontodegenerado(ptd,:)=i;
        endif
        count=count+1;
        nPumbsem(count,:)=i;
        break
      endif
    endfor
  endfor

  Pumbsem=Pumbc;
  uu=0;
  uu(1:size(Pumbc,1),1)=[1:size(Pumbc,1)];
  
  for i=size(nPumbsem,1):-1:1
    Pumbsem=deletarLinha(Pumbsem,nPumbsem(i));
    uu=deletarLinha(uu,nPumbsem(i));
  endfor
       
  if (size(Pumbsem)!=0)
  for i=1:size(Pumbsem,1)
      if (numlig(uu(i,1))==2 && iVumbmin(uu(i,1))==1);
        Pumbc(uu(i,1),:)=V(ligadosr(uu(i,1),2),:);  
          for l=1:2:size(PMM2,2)
            for k=1:pontocurvas((l+1)/2,1)
              if ((PMM2(k,l)==ligadosr(uu(i,1),iVumbmin(uu(i,1)))) || (PMM2(k,l+1)==ligadosr(uu(i,1),iVumbmin(uu(i,1)))))
                for j=1:size(MM)
                  if ((PMM2(k,l)==PMM(j,1) & PMM2(k,l+1)==PMM(j,2)) || (PMM2(k,l)==PMM(j,2) & PMM2(k,l+1)==PMM(j,1)))
                    MM2(k,3*(l+1)/2-2:3*(l+1)/2)=MMcopia(j,:);
                  endif
                endfor
              endif
            endfor
          endfor
        for l=1:2:size(PMM2,2)
          for k=1:pontocurvas((l+1)/2,1)
            if (PMM2(k,l)==ligadosr(uu(i,1),2) || PMM2(k,l+1)==ligadosr(uu(i,1),2))
              MM2(k,3*(l+1)/2-2:3*(l+1)/2)=Pumbc(uu(i,1),:);
              idx2=0;
              idx3=1;
            endif
          endfor
        endfor
      elseif (numlig(uu(i,1))==2 && iVumbmin(uu(i,1))==2);
        Pumbc(uu(i,1),:)=V(ligadosr(uu(i,1),1),:);
        for l=1:2:size(PMM2,2)
            for k=1:pontocurvas((l+1)/2,1)
              if ((PMM2(k,l)==ligadosr(uu(i,1),iVumbmin(uu(i,1)))) || (PMM2(k,l+1)==ligadosr(uu(i,1),iVumbmin(uu(i,1)))))
                for j=1:size(MM)
                  if ((PMM2(k,l)==PMM(j,1) & PMM2(k,l+1)==PMM(j,2)) || (PMM2(k,l)==PMM(j,2) & PMM2(k,l+1)==PMM(j,1)))
                    MM2(k,3*(l+1)/2-2:3*(l+1)/2)=MMcopia(j,:);
                  endif
                endfor
              endif
            endfor
          endfor
        for l=1:2:size(PMM2,2)
          for k=1:pontocurvas((l+1)/2,1)
            if (PMM2(k,l)==ligadosr(uu(i,1),1) || PMM2(k,l+1)==ligadosr(uu(i,1),1))
              MM2(k,3*(l+1)/2-2:3*(l+1)/2)=Pumbc(uu(i,1),:);
              idx2=0;
              idx3=1;
            endif
          endfor
        endfor
      endif        
  endfor
  endif
endif

if(size(pontocurvasax,1)==size(pontocurvas,1))
idx2=1;
endif


endwhile

% Fechando as curvas ridges azuis fechadas

for tt=1:2:size(PMM2,2)
for i=1:nfaces
if((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
endif
endfor
endfor


printf('Azul fecho');



####### Vermelha ########

%Juntando 2 curvas ridges vermelhas que o primeiro ponto da curva 1 
%e o primeiro ponto da curva 2 estão na mesma face
idx2=0;
idx3=0;
while(idx2!=1)
pontocurvasaxw=pontocurvasw;
idx=0;
while(idx!=1)
pontocurvasaw=pontocurvasw;
for tt=1:2:size(PMM2w,2)
  ppt=size(PMM2w,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMM2w,2)
  pptz=size(PMM2w,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp)==F(i,2) && PMM2w(1,pp+1)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp+1)==F(i,2) && PMM2w(1,pp)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp)==F(i,1) && PMM2w(1,pp+1)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp+1)==F(i,1) && PMM2w(1,pp)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp)==F(i,3) && PMM2w(1,pp+1)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp+1)==F(i,3) && PMM2w(1,pp)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp)==F(i,2) && PMM2w(1,pp+1)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp+1)==F(i,2) && PMM2w(1,pp)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp)==F(i,1) && PMM2w(1,pp+1)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp+1)==F(i,1) && PMM2w(1,pp)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp)==F(i,3) && PMM2w(1,pp+1)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp+1)==F(i,3) && PMM2w(1,pp)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp)==F(i,1) && PMM2w(1,pp+1)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp+1)==F(i,1) && PMM2w(1,pp)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp)==F(i,2) && PMM2w(1,pp+1)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp+1)==F(i,2) && PMM2w(1,pp)==F(i,3)))
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp)==F(i,2) && PMM2w(1,pp+1)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp+1)==F(i,2) && PMM2w(1,pp)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp)==F(i,3) && PMM2w(1,pp+1)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp+1)==F(i,3) && PMM2w(1,pp)==F(i,1)))
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp)==F(i,3) && PMM2w(1,pp+1)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp+1)==F(i,3) && PMM2w(1,pp)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp)==F(i,1) && PMM2w(1,pp+1)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp+1)==F(i,1) && PMM2w(1,pp)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasaw,1)==size(pontocurvasw,1))
idx=1;
end
endwhile

printf('Vermelha 1');

%Juntando 2 curvas ridges vermelhas que o último ponto da curva 2 
%e o primeiro ponto da curva 1 estão na mesma face

#2
idx=0;
while(idx!=1)
pontocurvasaw=pontocurvasw;
for tt=1:2:size(PMM2w,2)
  ppt=size(PMM2w,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMM2w,2)
  pptz=size(PMM2w,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasaw,1)==size(pontocurvasw,1))
idx=1;
end
endwhile

printf('Vermelha 2');

%Juntando 2 curvas ridges vermelhas que o último ponto da curva 1 
%e o primeiro ponto da curva 2 estão na mesma face

idx=0;
while(idx!=1)
pontocurvasaw=pontocurvasw;
for pp=1:2:size(PMM2w,2)
  ppt=size(PMM2w,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMM2w,2)
  pptz=size(PMM2w,2);
  if (tt>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasaw,1)==size(pontocurvasw,1))
idx=1;
end
endwhile

printf('Vermelha 3');

%Fazendo a curva ridge vermelha, que o último ponto da curva está na 
%mesma face de um vértice que é ridge vermelha, passar por esse vértice
#4
idx=0;
while(idx!=1)
mx=0;
for tt=1:2:size(PMM2w,2)
for pp=1:size(zero2,1)
for i=1:nfaces
  if((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && zero2(pp,1)==F(i,3)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && zero2(pp,1)==F(i,1)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && zero2(pp,1)==F(i,2)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && zero2(pp,1)==F(i,3)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && zero2(pp,1)==F(i,1)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && zero2(pp,1)==F(i,2)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
endif
endfor
endfor
endfor
if(mx==0)
idx=1;
end
endwhile

printf('Vermelha 4');

%Fazendo com que a curva ridge vermelha passe pelo ponto umbílico
% Caso em que o ponto umbílico se originou de 2 pontos umbílicos ligados


abc=size(PMM2w,2);
for j=1:size(ligadosr,1)
  if (numlig(j) ==2)
    for tt=1:2:abc
      idx=0;
      for i=1:pontocurvasw((tt+1)/2)
        if (idx == 1)
          break;
        endif
        if (PMM2w(i,tt)==aa(1,j) ||PMM2w(i,tt)==aaa(1,j) || PMM2w(i,tt+1)==aa(1,j) ||PMM2w(i,tt+1)==aaa(1,j) )
          if   ((PMM2w(i,tt)==aa(1,j) && PMM2w(i,tt+1)==aaa(1,j)) || (PMM2w(i,tt+1)==aa(1,j) && PMM2w(i,tt)==aaa(1,j)) )
            break
          else
        for k=1:nfaces
          if ((PMM2w(i,tt)==F(k,1) && PMM2w(i,tt+1)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3)) || (PMM2w(i,tt)==F(k,1) && PMM2w(i,tt+1)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt)==F(k,2) && PMM2w(i,tt+1)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1)) || (PMM2w(i,tt)==F(k,2) && PMM2w(i,tt+1)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt)==F(k,3) && PMM2w(i,tt+1)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2)) || (PMM2w(i,tt)==F(k,3) && PMM2w(i,tt+1)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt)==F(k,2) && PMM2w(i,tt+1)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3)) || (PMM2w(i,tt)==F(k,2) && PMM2w(i,tt+1)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt)==F(k,1) && PMM2w(i,tt+1)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2)) || (PMM2w(i,tt)==F(k,1) && PMM2w(i,tt+1)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt)==F(k,3) && PMM2w(i,tt+1)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1)) || (PMM2w(i,tt)==F(k,3) && PMM2w(i,tt+1)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
            elseif ((PMM2w(i,tt+1)==F(k,1) && PMM2w(i,tt)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3)) || (PMM2w(i,tt+1)==F(k,1) && PMM2w(i,tt)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt+1)==F(k,2) && PMM2w(i,tt)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1)) || (PMM2w(i,tt+1)==F(k,2) && PMM2w(i,tt)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt+1)==F(k,3) && PMM2w(i,tt)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2)) || (PMM2w(i,tt+1)==F(k,3) && PMM2w(i,tt)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt+1)==F(k,2) && PMM2w(i,tt)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3)) || (PMM2w(i,tt+1)==F(k,2) && PMM2w(i,tt)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt+1)==F(k,1) && PMM2w(i,tt)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2)) || (PMM2w(i,tt+1)==F(k,1) && PMM2w(i,tt)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          elseif ((PMM2w(i,tt+1)==F(k,3) && PMM2w(i,tt)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1)) || (PMM2w(i,tt+1)==F(k,3) && PMM2w(i,tt)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1)))
            MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MM2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvasw((tt+1)/2))
            MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
            PMM2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvasw((tt+1)/2))
            PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
            end
            pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
            PMM2w=deletarColuna(PMM2w,tt);
            PMM2w=deletarColuna(PMM2w,tt);
            pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
            idx=1;
            break
          endif
    endfor
  endif
    endif
  endfor
endfor
endif
endfor
% Caso em que o ponto umbílico se originou de 3 ou mais pontos umbílicos ligados

for j=1:size(ligadosr,1)
  if (numlig(j) >=3)
    for tt=1:2:size(PMM2w,2)
      idx=0;
      for i=1:pontocurvasw((tt+1)/2)
        if (idx==1)
            break
          endif
         if (PMM2w(i,tt)==aa(1,j) && PMM2w(i,tt+1)==aaa(1,j) && ((PMM2w(i+1,tt)==aaa(1,j) && PMM2w(i+1,tt+1)==aaaa(1,j)) | (PMM2w(i+1,tt)==aaaa(1,j) && PMM2w(i+1,tt+1)==aaa(1,j)) | (PMM2w(i+1,tt)==aa(1,j) && PMM2w(i+1,tt+1)==aaaa(1,j)) | (PMM2w(i+1,tt)==aaaa(1,j) && PMM2w(i+1,tt+1)==aa(1,j))))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          elseif (PMM2w(i,tt)==aaa(1,j) && PMM2w(i,tt+1)==aa(1,j) && ((PMM2w(i+1,tt)==aaa(1,j) && PMM2w(i+1,tt+1)==aaaa(1,j)) | (PMM2w(i+1,tt)==aaaa(1,j) && PMM2w(i+1,tt+1)==aaa(1,j)) | (PMM2w(i+1,tt)==aa(1,j) && PMM2w(i+1,tt+1)==aaaa(1,j)) | (PMM2w(i+1,tt)==aaaa(1,j) && PMM2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(i,tt)==aa(1,j) && PMM2w(i,tt+1)==aaaa(1,j) && ((PMM2w(i+1,tt)==aaa(1,j) && PMM2w(i+1,tt+1)==aaaa(1,j)) | (PMM2w(i+1,tt)==aaaa(1,j) && PMM2w(i+1,tt+1)==aaa(1,j)) | (PMM2w(i+1,tt)==aa(1,j) && PMM2w(i+1,tt+1)==aaa(1,j)) | (PMM2w(i+1,tt)==aaa(1,j) && PMM2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(i,tt)==aaaa(1,j) && PMM2w(i,tt+1)==aa(1,j) && ((PMM2w(i+1,tt)==aaa(1,j) && PMM2w(i+1,tt+1)==aaaa(1,j)) | (PMM2w(i+1,tt)==aaaa(1,j) && PMM2w(i+1,tt+1)==aaa(1,j)) | (PMM2w(i+1,tt)==aa(1,j) && PMM2w(i+1,tt+1)==aaa(1,j)) | (PMM2w(i+1,tt)==aaa(1,j) && PMM2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(i,tt)==aaaa(1,j) && PMM2w(i,tt+1)==aaa(1,j) && ((PMM2w(i+1,tt)==aa(1,j) && PMM2w(i+1,tt+1)==aaaa(1,j)) | (PMM2w(i+1,tt)==aaaa(1,j) && PMM2w(i+1,tt+1)==aa(1,j)) | (PMM2w(i+1,tt)==aa(1,j) && PMM2w(i+1,tt+1)==aaa(1,j)) | (PMM2w(i+1,tt)==aaa(1,j) && PMM2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(i,tt)==aaa(1,j) && PMM2w(i,tt+1)==aaaa(1,j) && ((PMM2w(i+1,tt)==aa(1,j) && PMM2w(i+1,tt+1)==aaaa(1,j)) | (PMM2w(i+1,tt)==aaaa(1,j) && PMM2w(i+1,tt+1)==aa(1,j)) | (PMM2w(i+1,tt)==aa(1,j) && PMM2w(i+1,tt+1)==aaa(1,j)) | (PMM2w(i+1,tt)==aaa(1,j) && PMM2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2w(1:i,end+1:end+3)=MM2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvasw((tt+1)/2))
              MM2w(i+2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(i+1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMM2w(1:i,end+1:end+2)=PMM2w(1:i,tt:(tt+1));
              PMM2w(i+1,end-1)=nverts+10;
              PMM2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvasw((tt+1)/2))
              PMM2w(i+2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(i+1:pontocurvasw((tt+1)/2),tt:(tt+1));
              end
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          elseif (PMM2w(1,tt)==aa(1,j) && PMM2w(1,tt+1)==aaa(1,j) && (PMM2w(2,tt)!=aaaa(1,j) | PMM2w(2,tt+1)!=aaaa(1,j)))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          elseif (PMM2w(1,tt)==aaa(1,j) && PMM2w(1,tt+1)==aa(1,j) && (PMM2w(2,tt)!=aaaa(1,j) | PMM2w(2,tt+1)!=aaaa(1,j)))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(1,tt)==aa(1,j) && PMM2w(1,tt+1)==aaaa(1,j) && (PMM2w(2,tt)!=aaa(1,j) | PMM2w(2,tt+1)!=aaa(1,j)))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
          elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(1,tt)==aaaa(1,j) && PMM2w(1,tt+1)==aa(1,j) && (PMM2w(2,tt)!=aaa(1,j) | PMM2w(2,tt+1)!=aaa(1,j)))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(1,tt)==aaaa(1,j) && PMM2w(1,tt+1)==aaa(1,j) && (PMM2w(2,tt)!=aa(1,j) | PMM2w(2,tt+1)!=aa(1,j)))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(1,tt)==aaa(1,j) && PMM2w(1,tt+1)==aaaa(1,j) && (PMM2w(2,tt)!=aa(1,j) | PMM2w(2,tt+1)!=aa(1,j)))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
          elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
          elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2w(1,end+1:end+3)=Pumbc(j,:);
              MM2w(2:pontocurvasw((tt+1)/2)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMM2w(1,end+1)=nverts+10;
              PMM2w(1,end+1)=nverts+10;
              PMM2w(2:pontocurvasw((tt+1)/2)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2),tt:(tt+1));
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(pontocurvasw((tt+1)/2,1),tt)==aa(1,j) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==aaa(1,j))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          elseif (PMM2w(pontocurvasw((tt+1)/2,1),tt)==aaa(1,j) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==aa(1,j))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(pontocurvasw((tt+1)/2,1),tt)==aa(1,j) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==aaaa(1,j))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(pontocurvasw((tt+1)/2,1),tt)==aaaa(1,j) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==aa(1,j))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(pontocurvasw((tt+1)/2,1),tt)==aaaa(1,j) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==aaa(1,j))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
            endfor
          elseif (PMM2w(pontocurvasw((tt+1)/2,1),tt)==aaa(1,j) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==aaaa(1,j))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
              elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MM2w(pontocurvasw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
              PMM2w=deletarColuna(PMM2w,tt);
              PMM2w=deletarColuna(PMM2w,tt);
              pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
              idx=1;
              break
            endif
          endfor
          endif
      endfor
    endfor
  endif
endfor
printf('Vermelha 5');

%Fazendo a curva ridge vermelha, que o primeiro ponto da curva está na 
%mesma face de um vértice que é ridge vermelha, passar por esse vértice

idx=0;
while(idx!=1)
mx=0;
for tt=1:2:size(PMM2w,2)
if (pontocurvasw((tt+1)/2,1)>=2)
for pp=1:size(zero2,1)
for i=1:nfaces
  if((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && zero2(pp,1)==F(i,3) && zero2(pp,1)!=PMM2w(2,tt)) )
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && zero2(pp,1)==F(i,1) && zero2(pp,1)!=PMM2w(2,tt)))
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && zero2(pp,1)==F(i,2) && zero2(pp,1)!=PMM2w(2,tt)))
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && zero2(pp,1)==F(i,3) && zero2(pp,1)!=PMM2w(2,tt)))
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && zero2(pp,1)==F(i,1) && zero2(pp,1)!=PMM2w(2,tt)))
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && zero2(pp,1)==F(i,2) && zero2(pp,1)!=PMM2w(2,tt)) )
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
endif
endfor
endfor
endif
endfor
if(mx==0)
idx=1;
end
endwhile

printf('Vermelha 6');

%Juntando 2 curvas ridges vermelhas que o último ponto da curva 1 
%e o primeiro ponto da curva 2 são os mesmos pontos


idx=0;
while(idx!=1)
pontocurvasaw=pontocurvasw;
for pp=1:2:size(PMM2w,2)
  ppt=size(PMM2w,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMM2w,2)
  pptz=size(PMM2w,2);
  if (tt>pptz)
    break
  endif
  if(MM2w(pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2)==MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2))
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1)-1,end-2:end)=MM2w(2:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1)-1,end-1:end)=PMM2w(2:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1)-1;
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
endif
endfor
endfor
if(size(pontocurvasaw,1)==size(pontocurvasw,1))
idx=1;
end
endwhile

%Juntando 2 curvas ridges vermelhas que o primeiro ponto da curva 1 
%e o primeiro ponto da curva 2 são os mesmos pontos

idx=0;
while(idx!=1)
pontocurvasaw=pontocurvasw;
for pp=1:2:size(PMM2w,2)
  ppt=size(PMM2w,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMM2w,2)
  pptz=size(PMM2w,2);
  if (tt>pptz)
    break
  endif
  if (MM2w(1,3*(pp+1)/2-2:3*(pp+1)/2)==MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2))
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1)-1,end-2:end)=MM2w(2:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1)-1,end-1:end)=PMM2w(2:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1)-1;
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
endif
endfor
endfor
if(size(pontocurvasaw,1)==size(pontocurvasw,1))
idx=1;
end
endwhile

if(size(pontocurvasaxw,1)==size(pontocurvasw,1))
idx2=1;
end
endwhile

%Verificando se há curvas ridges vermelhas degeneradas

nPumbsemw=[];
  MM2wc=[];
  pontodegeneradow=[];
  countw=0;
  ptdw=0;
  for i=1:size(Pumbc)
    for j=1:size(pontocurvasw)(1)
      MM2w(1:pontocurvasw(j,1),3*j-2:3*j)==Pumbc(i,:);
      tilw = ans ;
      if (max(tilw(:,1)) == 1 && max(tilw(:,2)) == 1 && max(tilw(:,3)) == 1)
        if (numlig(i)>=1 && MM2w(1,3*j-2:3*j)== MM2w(2:pontocurvasw(j,1),3*j-2:3*j))
          ptdw=ptdw+1;
          pontodegeneradow(ptdw,:)=i;
        endif
        countw=countw+1;
        nPumbsemw(countw,:)=i;
        break
      endif
    endfor
  endfor

  Pumbsemw=Pumbc;
  uuw=0;
  uuw(1:size(Pumbc,1),1)=[1:size(Pumbc,1)];
  
  for i=size(nPumbsemw,1):-1:1
    Pumbsemw=deletarLinha(Pumbsemw,nPumbsemw(i));
    uuw=deletarLinha(uuw,nPumbsemw(i));
  endfor

% Fechando as curvas ridges vermelhas fechadas
for tt=1:2:size(PMM2w,2)
for i=1:nfaces
if((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
endif
endfor
endfor

printf('Vermelha fecho');


printf('Parte 7');  

%Quantidade de curvas ridges abertas e fechadas
con6=0;
ncon6=0;
for tt=1:3:size(MM2,2)
  if (pontocurvas((tt+2)/3,1)>1 && MM2(1,tt)==MM2(pontocurvas((tt+2)/3,1),tt) && MM2(1,tt+1)==MM2(pontocurvas((tt+2)/3,1),tt+1))
    circulo(con6+1,1)=tt;
    con6=con6+1; %Quantidade de curvas ridge azul fechada
  else
    ncirculo(ncon6+1,1)=tt;
    ncon6=ncon6+1; %Quantidade de curvas ridge azul aberta
  endif
endfor

con6w=0;
ncon6w=0;
for tt=1:3:size(MM2w,2)
  if (pontocurvasw((tt+2)/3,1)>1 && MM2w(1,tt)==MM2w(pontocurvasw((tt+2)/3,1),tt) && MM2w(1,tt+1)==MM2w(pontocurvasw((tt+2)/3,1),tt+1))
    circulow(con6w+1,1)=tt;
    con6w=con6w+1; %Quantidade de curvas ridge vermelha fechada
  else
    ncirculow(ncon6w+1,1)=tt;
    ncon6w=ncon6w+1; %Quantidade de curvas ridge vermelha fechada
  endif
endfor

totalazul = con6+ncon6; %Total de curvas ridge azul
totalvermelha= con6w+ncon6w; %Total de curvas ridge vermelha

%Determinação do número do vértice do ponto umbílico
for i=1:nverts
  for j=1:size(Pumbc,1)
    if (Pumbc(j,:)==V(i,:))
      VPumbc(j,1)=i;
    endif
  endfor
endfor

%%%%%%%%%%%% Ativar caso queira encontrar! %%%%%%%%%%%%%%%%

%%%%%Uso exclusivo para faces obtidas pelo MediaPipe

% Número dos vértices de cada região da face
% Obtenção da região de cada ponto umbílico
% Obtenção do número total de umbílico interior e de bordo



##

          % Número dos vértices de cada região da face
##Borda = [11, 339, 298, 333, 285, 252, 390, 357, 455, 324, 362, 289, 398, 366, 380, 379, 401, 378, 153, 149, 177, 150, 151, 137, 173, 59, 133, 94, 235, 128, 163, 22, 55, 104, 68, 110];
##
##Nariz = [2, 3, 4, 5, 6, 7, 20, 21, 45, 46, 49, 50, 52, 60, 61, 65, 76, 80, 95, 98, 99,  100, 103, 116, 123, 126, 130, 132, 135, 142, 167, 169, 175, 189, 194, 196, 197, 198, 199, 210, 218, 219, 220, 221, 236, 237, 238, 239, 240, 241, 242, 243, 249, 251, 275, 276, 279, 280, 282, 290, 291, 295, 306, 310,  327, 328, 329, 332, 345, 352,  355, 359, 361, 364, 371, 393, 400, 413, 418, 420, 421, 430, 438, 439, 440, 441, 456, 457, 458, 459, 460, 461, 462, 463];  
##
##Boca = [1, 12, 13, 14, 15, 16, 17, 18, 38, 39, 40, 41, 42, 43, 62, 63, 73, 74, 75, 77, 78, 79, 81, 82, 83, 85, 86, 87, 88, 89, 90, 91, 92, 96, 97, 147, 179, 180, 181, 182, 184, 185, 186, 192, 268, 269, 270, 271, 272, 273, 292, 293, 303, 304, 305, 307, 308, 309, 311, 312, 313, 315, 316, 317, 318, 319, 320, 321, 322, 325, 326, 376, 403, 404, 405, 406, 408, 409, 410, 416];
##
##Testa = [9, 10, 47, 53, 54, 56, 64, 66, 67, 69, 70, 71, 72, 105, 106, 108, 109, 125, 140, 152, 157, 277, 283, 284, 286, 294, 296, 297, 299, 300, 301, 302, 334, 335, 337, 338, 354, 369, 384];
##
##Olho_Direito = [8, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 57, 111, 113, 114, 131, 134, 145, 146, 154, 155, 156, 158, 159, 160, 161, 162, 164, 174, 190, 191, 222, 223, 224, 225, 226, 227, 229, 230, 231, 232, 233, 234, 244, 245, 246, 247, 248];
##
##Bochecha_Direita = [35, 36, 37, 48, 51, 58, 93, 101, 102, 112, 115, 117, 118, 119, 120, 121, 122, 124, 127, 129, 138, 143, 144, 148, 178, 187, 188, 193, 204, 206, 207, 208, 213, 214, 215, 216, 217, 228]; 
##
##Queixo = [19, 33, 44, 84, 107, 136, 139, 141, 170, 171, 172, 176, 183, 195, 200, 201, 202, 203, 205, 209, 211, 212, 263, 274, 314, 336, 365, 368, 370, 395, 396, 397, 407, 419, 422, 423, 425, 429, 431, 432];
##
##Bigode = [165, 166, 168, 392, 394]; 
##
##Olho_Esquerdo = [250, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 264,  287, 340, 342, 343,  360, 363, 374, 375, 381, 382, 383, 385, 386, 387, 388, 389, 391, 399, 414, 415, 442, 443, 444, 445, 446, 447, 449, 450, 451, 452, 453, 454, 464, 465, 466, 467, 468];
##
##Bochecha_Esquerda = [265, 266, 267, 278, 281, 288,  323, 330, 331, 341, 344, 346, 347, 348, 349, 350, 351, 353, 356, 358, 367, 372, 373, 377, 402, 411, 412, 417, 424, 426, 427, 428, 433, 434, 435, 436, 437, 448];
##



##totalumbfinal=size(Pumbc,1);
##totalumbfront=0;
##Vumbint=VPumbc;
##
##for i=1:size(Pumbc,1)
##    for j=1:size(Borda,2)
##        if (VPumbc(i)==Borda(1,j))
##            Vumbint = deletarLinha(Vumbint,i-totalumbfront);
##            totalumbfront = totalumbfront+1;
##       endif
##    endfor
##endfor
##
##totalumbint=totalumbfinal-totalumbfront;   %número total de umbílico interior
##
##umbnariz=0;
##VNariz=[];
##for i=1:size(Vumbint,1)
##  for j=1:size(Nariz,2)
##    if (Vumbint(i,1) == Nariz(1,j))
##      umbnariz=umbnariz+1;
##      VNariz(umbnariz,1)=Vumbint(i,1);  %Número dos vértices dos umbílicos no nariz
##      break;
##    endif
##  endfor
##endfor
##
##umbboca=0;
##VBoca=[];
##for i=1:size(Vumbint,1)
##  for j=1:size(Boca,2)
##    if (Vumbint(i,1) == Boca(1,j))
##      umbboca=umbboca+1;
##      VBoca(umbboca,1)=Vumbint(i,1);  %Número dos vértices dos umbílicos na boca
##      break;
##    endif
##  endfor
##endfor
##
##umbtesta=0;
##VTesta=[];
##for i=1:size(Vumbint,1)
##  for j=1:size(Testa,2)
##    if (Vumbint(i,1) == Testa(1,j))
##      umbtesta=umbtesta+1;
##      VTesta(umbtesta,1)=Vumbint(i,1);  %Número dos vértices dos umbílicos na testa
##      break;
##    endif
##  endfor
##endfor
##
##umbolhoD=0;
##VOlhoD=[];
##for i=1:size(Vumbint,1)
##  for j=1:size(Olho_Direito,2)
##    if (Vumbint(i,1) == Olho_Direito(1,j))
##      umbolhoD=umbolhoD+1;
##      VOlhoD(umbolhoD,1)=Vumbint(i,1); %Número dos vértices dos umbílicos no olho direito
##      break;
##    endif
##  endfor
##endfor
##
##umbbochechaD=0;
##VBochechaD=[];
##for i=1:size(Vumbint,1)
##  for j=1:size(Bochecha_Direita,2)
##    if (Vumbint(i,1) == Bochecha_Direita(1,j))
##      umbbochechaD=umbbochechaD+1;
##      VBochechaD(umbbochechaD,1)=Vumbint(i,1); %Número dos vértices dos umbílicos na buchecha direita
##      break;
##    endif
##  endfor
##endfor
##
##umbqueixo=0;
##VQueixo=[];
##for i=1:size(Vumbint,1)
##  for j=1:size(Queixo,2)
##    if (Vumbint(i,1) == Queixo(1,j))
##      umbqueixo=umbqueixo+1;
##      VQueixo(umbqueixo,1)=Vumbint(i,1);  %Número dos vértices dos umbílicos no queixo
##      break;
##    endif
##  endfor
##endfor
##
##umbbigode=0;
##VBigode=[];
##for i=1:size(Vumbint,1)
##  for j=1:size(Bigode,2)
##    if (Vumbint(i,1) == Bigode(1,j))
##      umbbigode=umbbigode+1;
##      VBigode(umbbigode,1)=Vumbint(i,1); %Número dos vértices dos umbílicos na buço
##      break;
##    endif
##  endfor
##endfor
##
##umbolhoE=0;
##VOlhoE=[];
##for i=1:size(Vumbint,1)
##  for j=1:size(Olho_Esquerdo,2)
##    if (Vumbint(i,1) == Olho_Esquerdo(1,j))
##      umbolhoE=umbolhoE+1;
##      VOlhoE(umbolhoE,1)=Vumbint(i,1); %Número dos vértices dos umbílicos no olho esquerdo
##      break;
##    endif
##  endfor
##endfor
##
##umbbochechaE=0;
##VBochechaE=[];
##for i=1:size(Vumbint,1)
##  for j=1:size(Bochecha_Esquerda,2)
##    if (Vumbint(i,1) == Bochecha_Esquerda(1,j))
##      umbbochechaE=umbbochechaE+1;
##      VBochechaE(umbbochechaE,1)=Vumbint(i,1); %Número dos vértices dos umbílicos na bochecha esquerda
##      break;
##    endif
##  endfor
##endfor
##
##totalumbint;     %Número de pontos umbílicos interiores 
##umbnariz;        %Número de pontos umbílicos no nariz
##umbboca;         %Número de pontos umbílicos na boca
##umbtesta;        %Número de pontos umbílicos na testa
##umbolhoD;        %Número de pontos umbílicos no olho direito
##umbbochechaD;    %Número de pontos umbílicos na buchecha direita
##umbqueixo;       %Número de pontos umbílicos no queixo
##umbbigode;       %Número de pontos umbílicos no buço
##umbolhoE;        %Número de pontos umbílicos no olho esquerdo
##umbbochechaE;    %Número de pontos umbílicos na bochecha esquerda
##
endif


if(resposta==2)




%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPARABÓLICA$$$$$$$$$$$$$$$$$$$$$$

%Cálculo do valor da curvatura principal nos pontos de cruzamento
  kxcruzs1=zeros(nverts,1);
  kxcruzs2=zeros(nverts,1);
  kxcruzs21=zeros(nverts,1);
  kxcruzs22=zeros(nverts,1);
  
  for i=1:nverts
    if(cruzat11(i) == 0 || cruzat12(i) == 0)
      continue
    else
      kxcruzs1(i,1)=k23(cruzat11(i))+tt1(i)*(k23(cruzat12(i))-k23(cruzat11(i)));
    end
    if(cruzamt11(i) == 0 || cruzamt12(i) == 0)
      continue
    else  
      kxcruzs2(i,1)=k13(cruzamt11(i))+tt2(i)*(k13(cruzamt12(i))-k13(cruzamt11(i)));
    end
    if(cruzat21(i) == 0 || cruzat22(i) == 0)
      continue
    else      
      kxcruzs21(i,1)=k23(cruzat21(i))+tt21(i)*(k23(cruzat22(i))-k23(cruzat21(i)));
    end
    if(cruzamt21(i) == 0 || cruzamt22(i) == 0)
      continue
    else
      kxcruzs22(i,1)=k13(cruzamt21(i))+tt22(i)*(k13(cruzamt22(i))-k13(cruzamt21(i)));
    endif
  endfor
  
% Cálculo da Derivada Direcional de k2 na direção e1
% Cálculo da Derivada Direcional de k1 na direção e2
  
  Ddsk2e1=zeros(nverts,1);
  Ddsk1e2=zeros(nverts,1);
  Positivos=[];
  Negativos=[];
  Positivos2=[];
  Negativos2=[];
  counts55=1;
  counts56=1;
  counts555=1;
  counts565=1;
  MMs=[];
  MMsw=[];
  
  for i=1:nverts
    if(abs(xcruz1(i,:)-xcruz21(i,:)) < 1e-10 || abs(xcruz2(i,:)-xcruz22(i,:)) < 1e-10)
      continue
    else
      Ddsk2e1(i,1)=(kxcruzs1(i,1)-kxcruzs21(i,1))/(dot(xcruz1(i,:)-xcruz21(i,:),[AV2(i,1) AV2(i,2) AV2(i,3)]/norm([AV2(i,1) AV2(i,2) AV2(i,3)])));
      Ddsk1e2(i,1)=(kxcruzs2(i,1)-kxcruzs22(i,1))/(dot(xcruz2(i,:)-xcruz22(i,:),[AV2(i,4) AV2(i,5) AV2(i,6)]/norm([AV2(i,4) AV2(i,5) AV2(i,6)])));
    end
    
    if(Ddsk2e1(i,1)>1e-10)
      Positivos(counts55,:)=V(i,:);
      counts55=counts55+1;
    end
    if (Ddsk2e1(i,1)<-1e-10)
      Negativos(counts56,:)=V(i,:);
      counts56=counts56+1;
    end
    if(Ddsk1e2(i,1)>1e-10)
      Positivos2(counts555,:)=V(i,:);
      counts555=counts555+1;
    end
    if (Ddsk1e2(i,1)<-1e-10)
      Negativos2(counts565,:)=V(i,:);
      counts565=counts565+1;
    end
  endfor

% Cálculo dos pontos subparabolicos vermelhos

    counts57=1;
    counts572=1;
  for i=1:nfaces
    if ((Ddsk2e1(F(i,1),1)>1e-10 && Ddsk2e1(F(i,2),1)<-1e-10))
      t=Ddsk2e1(F(i,1),1)/(Ddsk2e1(F(i,1),1)-Ddsk2e1(F(i,2),1));
      MMs(counts57,:)=V(F(i,1),:)+t*(V(F(i,2),:)-V(F(i,1),:));
      PMMs(counts57,1)=F(i,1);
      PMMs(counts57,2)=F(i,2);
      counts57=counts57+1;
      end
    if ((Ddsk2e1(F(i,2),1)>1e-10 && Ddsk2e1(F(i,3),1)<-1e-10))
      t=Ddsk2e1(F(i,2),1)/(Ddsk2e1(F(i,2),1)-Ddsk2e1(F(i,3),1));
      MMs(counts57,:)=V(F(i,2),:)+t*(V(F(i,3),:)-V(F(i,2),:));
      PMMs(counts57,1)=F(i,2);
      PMMs(counts57,2)=F(i,3);
      counts57=counts57+1;
    end
    if ((Ddsk2e1(F(i,3),1)>1e-10 && Ddsk2e1(F(i,1),1)<-1e-10))
      t=Ddsk2e1(F(i,3),1)/(Ddsk2e1(F(i,3),1)-Ddsk2e1(F(i,1),1));
      MMs(counts57,:)=V(F(i,3),:)+t*(V(F(i,1),:)-V(F(i,3),:));
      PMMs(counts57,1)=F(i,3);
      PMMs(counts57,2)=F(i,1);
      counts57=counts57+1;
    end
  if ((Ddsk2e1(F(i,2),1)>1e-10 && Ddsk2e1(F(i,1),1)<-1e-10))
      t=Ddsk2e1(F(i,1),1)/(Ddsk2e1(F(i,1),1)-Ddsk2e1(F(i,2),1));
      MMs(counts57,:)=V(F(i,1),:)+t*(V(F(i,2),:)-V(F(i,1),:));
      PMMs(counts57,1)=F(i,1);
      PMMs(counts57,2)=F(i,2);
      counts57=counts57+1;
    end
  if ((Ddsk2e1(F(i,3),1)>1e-10 && Ddsk2e1(F(i,2),1)<-1e-10))
      t=Ddsk2e1(F(i,2),1)/(Ddsk2e1(F(i,2),1)-Ddsk2e1(F(i,3),1));
      MMs(counts57,:)=V(F(i,2),:)+t*(V(F(i,3),:)-V(F(i,2),:));
      PMMs(counts57,1)=F(i,2);
      PMMs(counts57,2)=F(i,3);
      counts57=counts57+1;
    end
  if ((Ddsk2e1(F(i,1),1)>1e-10 && Ddsk2e1(F(i,3),1)<-1e-10))
      t=Ddsk2e1(F(i,3),1)/(Ddsk2e1(F(i,3),1)-Ddsk2e1(F(i,1),1));
      MMs(counts57,:)=V(F(i,3),:)+t*(V(F(i,1),:)-V(F(i,3),:));
      PMMs(counts57,1)=F(i,3);
      PMMs(counts57,2)=F(i,1);
      counts57=counts57+1;
    endif
  end

% Cálculo dos pontos subparabolicos azuis

  for i=1:nfaces
    if ((Ddsk1e2(F(i,1),1)>1e-10 && Ddsk1e2(F(i,2),1)<-1e-10))
      t2=Ddsk1e2(F(i,1),1)/(Ddsk1e2(F(i,1),1)-Ddsk1e2(F(i,2),1));
      MMsw(counts572,:)=V(F(i,1),:)+t2*(V(F(i,2),:)-V(F(i,1),:));
      PMMsw(counts572,1)=F(i,1);
      PMMsw(counts572,2)=F(i,2);
      counts572=counts572+1;
    end
    if ((Ddsk1e2(F(i,2),1)>1e-10 && Ddsk1e2(F(i,3),1)<-1e-10))
      t2=Ddsk1e2(F(i,2),1)/(Ddsk1e2(F(i,2),1)-Ddsk1e2(F(i,3),1));
      MMsw(counts572,:)=V(F(i,2),:)+t2*(V(F(i,3),:)-V(F(i,2),:));
      PMMsw(counts572,1)=F(i,2);
      PMMsw(counts572,2)=F(i,3);
      counts572=counts572+1;
    end
    if ((Ddsk1e2(F(i,3),1)>1e-10 && Ddsk1e2(F(i,1),1)<-1e-10))
      t2=Ddsk1e2(F(i,3),1)/(Ddsk1e2(F(i,3),1)-Ddsk1e2(F(i,1),1));
      MMsw(counts572,:)=V(F(i,3),:)+t2*(V(F(i,1),:)-V(F(i,3),:));
      PMMsw(counts572,1)=F(i,3);
      PMMsw(counts572,2)=F(i,1);
      counts572=counts572+1;
    endif
    if ((Ddsk1e2(F(i,2),1)>1e-10 && Ddsk1e2(F(i,1),1)<-1e-10))
      t2=Ddsk1e2(F(i,1),1)/(Ddsk1e2(F(i,1),1)-Ddsk1e2(F(i,2),1));
      MMsw(counts572,:)=V(F(i,1),:)+t2*(V(F(i,2),:)-V(F(i,1),:));
      PMMsw(counts572,1)=F(i,1);
      PMMsw(counts572,2)=F(i,2);
      counts572=counts572+1;
    end
    if ((Ddsk1e2(F(i,3),1)>1e-10 && Ddsk1e2(F(i,2),1)<-1e-10))
      t2=Ddsk1e2(F(i,2),1)/(Ddsk1e2(F(i,2),1)-Ddsk1e2(F(i,3),1));
      MMsw(counts572,:)=V(F(i,2),:)+t2*(V(F(i,3),:)-V(F(i,2),:));
      PMMsw(counts572,1)=F(i,2);
      PMMsw(counts572,2)=F(i,3);
      counts572=counts572+1;
    end
    if ((Ddsk1e2(F(i,1),1)>1e-10 && Ddsk1e2(F(i,3),1)<-1e-10))
      t2=Ddsk1e2(F(i,3),1)/(Ddsk1e2(F(i,3),1)-Ddsk1e2(F(i,1),1));
      MMsw(counts572,:)=V(F(i,3),:)+t2*(V(F(i,1),:)-V(F(i,3),:));
      PMMsw(counts572,1)=F(i,3);
      PMMsw(counts572,2)=F(i,1);
      counts572=counts572+1;
    endif
  end
  
if (size(MMs,1)==0 || size(MMsw,1)==0)
  fprintf('\n NÃO HÁ CURVAS SUBPARABÓLICAS!! \n');
else


%Limpando pontos subparabolicos repetidos


MMs3=[MMs PMMs];
MMs3w=[MMsw PMMsw];

MMs3=unique(MMs3,'rows');
MMs3w=unique(MMs3w,'rows');

for i=size(MMs3,1):-1:2
  if (MMs3(i,5)==MMs3(i-1,4) && MMs3(i,4)==MMs3(i-1,5))
    MMs3=deletarLinha(MMs3,i);
  endif
endfor

for i=size(MMs3w,1):-1:2
  if (MMs3w(i,5)==MMs3w(i-1,4) && MMs3w(i,4)==MMs3w(i-1,5))
    MMs3w=deletarLinha(MMs3w,i);
  endif
endfor


MMs=MMs3(:,1:3);
PMMs=MMs3(:,4:5);
MMsw=MMs3w(:,1:3);
PMMsw=MMs3w(:,4:5);
  

MMscopia=MMs;
MMswcopia=MMsw;
  
  for i=1:size(MMs)
    for j=1:size(MMw1)
      if((PMMs(i,1)==PMM1(j,1) & PMMs(i,2)==PMM1(j,2)) || (PMMs(i,1)==PMM1(j,2) & PMMs(i,2)==PMM1(j,1)))
        MMscopia(i,:)=MMw1(j,:);
      end
    endfor
  endfor
  
  for i=1:size(MMsw)
    for j=1:size(MMw1)
      if((PMMsw(i,1)==PMM1(j,1) & PMMsw(i,2)==PMM1(j,2)) || (PMMsw(i,1)==PMM1(j,2) & PMMsw(i,2)==PMM1(j,1)))
        MMswcopia(i,:)=MMw1(j,:);
      end
    endfor
  endfor
  
%Fazendo com que as curvas subparabolicas passem pelos pontos umbílicos
  
for k=1:size(ligadosr,1)
  if (numlig(k)==1)
    for i=1:size(MMs)
        if (PMMs(i,1)==aa(1,k) || PMMs(i,2)==aa(1,k))
          MMs(i,:)=Pumb(k,:);
        endif
    endfor
  endif
endfor

for k=1:size(ligadosr,1)
  if (numlig(k)==1)
    for i=1:size(MMsw)
        if (PMMsw(i,1)==aa(1,k) || PMMsw(i,2)==aa(1,k))
          MMsw(i,:)=Pumb(k,:);
        endif
    endfor
  endif
endfor


for k=1:size(ligadosr,1)
  if (numlig(k)>=2)
    for i=1:size(MMs)
      for hh=1:size(VPumbc,1)
        if (PMMs(i,1)==VPumbc(hh,1) || PMMs(i,2)==VPumbc(hh,1))
          MMs(i,:)=Pumbc(hh,:);
        endif
      endfor
    endfor
  endif
endfor

for k=1:size(ligadosr,1)
  if (numlig(k)>=2)
    for i=1:size(MMsw)
      for hh=1:size(VPumbc,1)
        if (PMMsw(i,1)==VPumbc(hh,1) || PMMsw(i,2)==VPumbc(hh,1))
          MMsw(i,:)=Pumbc(hh,:);
        endif
      endfor
    endfor
  endif
endfor


MMs3=[MMs PMMs];
MMs3w=[MMsw PMMsw];

%Encontrando os vértices que possuem derivada direcional nula
  
cts=1;
cts2=1;
zerosub=[];
zerosub2=[];


for i=1:nverts
  if(abs(Ddsk2e1(i))<=1e-10)
    zerosub(cts,1)=i;
    cts=cts+1;
  end
  if(abs(Ddsk1e2(i))<=1e-10)
    zerosub2(cts2,1)=i;
    cts2=cts2+1;
  end
endfor

guardars=[];
for i=1:size(zerosub,1)
  ccsccsc=1;
  guardars(i,ccsccsc)=zerosub(i,1);
  for j=1:L(i)
    for k=1:size(zerosub,1)
      if (SS(zerosub(i),j)==zerosub(k))
        guardars(i,ccsccsc+1)=zerosub(k);
        ccsccsc=ccsccsc+1;
        break
        endif
    endfor
    if (ccsccsc>1)
      break
    endif
  endfor
  idxs=0;
  while (idxs!=1)
    ccsc=ccsccsc; 
   if (j+1<=L(i))  
    for k=1:size(zerosub,1)
      if (SS(zerosub(i),j+1)==zerosub(k))
      guardars(i,ccsccsc+1)=zerosub(k);
      ccsccsc=ccsccsc+1;
      j=j+1;
      break
      endif
    endfor
   endif
   if(ccsc==ccsccsc)
      idxs=1;
    else
      ccsc=ccsccsc;
    endif
  endwhile
  nguardars(i,1)=ccsccsc;
end

guardars2=[];
for i=1:size(zerosub2,1)
  ccsccsc=1;
  guardars2(i,ccsccsc)=zerosub2(i,1);
  for j=1:L(i)
    for k=1:size(zerosub2,1)
      if (SS(zerosub2(i),j)==zerosub2(k))
        guardars2(i,ccsccsc+1)=zerosub2(k);
        ccsccsc=ccsccsc+1;
        break
      endif
    endfor
    if (ccsccsc>1)
      break
    endif
  endfor
  idxs=0;
  while (idxs!=1)
    ccsc=ccsccsc;
    if (j+1<=L(i))  
     for k=1:size(zerosub2,1)
      if (SS(zerosub2(i),j+1)==zerosub2(k))
      guardars2(i,ccsccsc+1)=zerosub2(k);
      ccsccsc=ccsccsc+1;
      j=j+1;
      break
      endif
     endfor
    endif
    if(ccsc==ccsccsc)
      idxs=1;
    else
      ccsc=ccsccsc;
    endif
  endwhile
  nguardars2(i,1)=ccsccsc;
end






  
  printf('Parte 8');


% Determinação de cada curva subparabolica (juntando os pontos para
  % formação das curvas)

%subparabolica azul  

MMs2w=[];
jafoisw=zeros(size(MMs3w)(1),1);
kksw=1;
for iisw=1:size(MMs3w)(1)
  if (jafoisw(iisw)==0)
    MMs2w(1,3*kksw-2:3*kksw)=MMs3w(iisw,1:3);
    PMMs2w(1,2*kksw-1:2*kksw)=MMs3w(iisw,4:5);
ccsw=1;
llsw=iisw;
jjjsw=0;

jafoisw(iisw)=1;
while (llsw != jjjsw)
  jjjsw=llsw;
  for i=1:nfaces
    if (MMs3w(llsw,4)==F(i,1) && MMs3w(llsw,5)==F(i,2))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3)))
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MMs3w(llsw,4)==F(i,1) && MMs3w(llsw,5)==F(i,2))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3)))
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MMs3w(llsw,4)==F(i,2) && MMs3w(llsw,5)==F(i,3))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3)) || (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1)))
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MMs3w(llsw,4)==F(i,2) && MMs3w(llsw,5)==F(i,3))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1)))
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3w(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MMs3w(llsw,4)==F(i,3) && MMs3w(llsw,5)==F(i,2))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3)) || (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1)))
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3w(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MMs3w(llsw,4)==F(i,2) && MMs3w(llsw,5)==F(i,1))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3)))
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3w(k,1:3);
          break
      endif
      end
    endfor
  end

      if (MMs3w(llsw,4)==F(i,3) && MMs3w(llsw,5)==F(i,1))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3)) || (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2)))
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MMs3w(llsw,4)==F(i,3) && MMs3w(llsw,5)==F(i,1))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3)))
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MMs3w(llsw,4)==F(i,1) && MMs3w(llsw,5)==F(i,3))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3)) || (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2)))
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3w(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MMs3w(llsw,4)==F(i,1) && MMs3w(llsw,5)==F(i,3))
      for k=1:size(MMs3w)(1)
        if(jafoisw(k)==0)
        if ((MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2)) || (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1)) || (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3)))
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,2))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,2);
          endif
          if (MMs3w(k,4)==F(i,3) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,3);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,2) && MMs3w(k,5)==F(i,1))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,2);
            PMMs2w(ccsw+1,2*kksw)=F(i,1);
          endif
          if (MMs3w(k,4)==F(i,1) && MMs3w(k,5)==F(i,3))
            PMMs2w(ccsw+1,2*kksw-1)=F(i,1);
            PMMs2w(ccsw+1,2*kksw)=F(i,3);
          endif
          llsw=k;
          jafoisw(k)=1;
          ccsw=ccsw+1;
          MMs2w(ccsw,3*kksw-2:3*kksw)=MMs3w(k,1:3);
          break
      endif
      end
    endfor
    end

endfor
end
pontocurvassw(kksw,1)=ccsw;
kksw=kksw+1;
endif
endfor

%subparabolica vermelha
MMs2=[];
jafois=zeros(size(MMs3)(1),1);
kks=1;
for iis=1:size(MMs3)(1)
  if (jafois(iis)==0)
    MMs2(1,3*kks-2:3*kks)=MMs3(iis,1:3);
    PMMs2(1,2*kks-1:2*kks)=MMs3(iis,4:5);
ccs=1;
lls=iis;
jjjs=0;

jafois(iis)=1;
while (lls != jjjs)
  jjjs=lls;
  for i=1:nfaces
    if (MMs3(lls,4)==F(i,1) && MMs3(lls,5)==F(i,2))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3)))
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MMs3(lls,4)==F(i,1) && MMs3(lls,5)==F(i,2))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3)))
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
    end

      if (MMs3(lls,4)==F(i,2) && MMs3(lls,5)==F(i,3))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3)) || (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1)))
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MMs3(lls,4)==F(i,2) && MMs3(lls,5)==F(i,3))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1)))
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MMs3(lls,4)==F(i,3) && MMs3(lls,5)==F(i,2))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3)) || (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1)))
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MMs3(lls,4)==F(i,2) && MMs3(lls,5)==F(i,1))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3)))
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
  end

      if (MMs3(lls,4)==F(i,3) && MMs3(lls,5)==F(i,1))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3)) || (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2)))
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MMs3(lls,4)==F(i,3) && MMs3(lls,5)==F(i,1))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3)))
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MMs3(lls,4)==F(i,1) && MMs3(lls,5)==F(i,3))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3)) || (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2)))
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MMs3(lls,4)==F(i,1) && MMs3(lls,5)==F(i,3))
      for k=1:size(MMs3)(1)
        if(jafois(k)==0)
        if ((MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2)) || (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1)) || (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3)))
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,2))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,2);
          endif
          if (MMs3(k,4)==F(i,3) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,3);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,2) && MMs3(k,5)==F(i,1))
            PMMs2(ccs+1,2*kks-1)=F(i,2);
            PMMs2(ccs+1,2*kks)=F(i,1);
          endif
          if (MMs3(k,4)==F(i,1) && MMs3(k,5)==F(i,3))
            PMMs2(ccs+1,2*kks-1)=F(i,1);
            PMMs2(ccs+1,2*kks)=F(i,3);
          endif
          lls=k;
          jafois(k)=1;
          ccs=ccs+1;
          MMs2(ccs,3*kks-2:3*kks)=MMs3(k,1:3);
          break
      endif
      end
    endfor
    end

endfor
end
pontocurvass(kks,1)=ccs;
kks=kks+1;
endif
endfor






######## Vermelha #####
##

%Juntando 2 curvas subparabolicas vermelhas que o primeiro ponto da curva 1 
%e o primeiro ponto da curva 2 estão na mesma face
idxs2=0;
idxs3=0;
while(idxs2!=1)
pontocurvassax=pontocurvass;
#1
idxs=0;
while(idxs!=1)
pontocurvassa=pontocurvass;
for tt=1:2:size(PMMs2,2)
  ppt=size(PMMs2,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMMs2,2)
  pptz=size(PMMs2,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,pp)==F(i,2) && PMMs2(1,pp+1)==F(i,3)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,pp+1)==F(i,2) && PMMs2(1,pp)==F(i,3)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,pp)==F(i,1) && PMMs2(1,pp+1)==F(i,3)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,pp+1)==F(i,1) && PMMs2(1,pp)==F(i,3)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(pontocurvass((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(pontocurvass((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,pp)==F(i,3) && PMMs2(1,pp+1)==F(i,1)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,pp+1)==F(i,3) && PMMs2(1,pp)==F(i,1)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,pp)==F(i,2) && PMMs2(1,pp+1)==F(i,1)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,pp+1)==F(i,2) && PMMs2(1,pp)==F(i,1)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(pontocurvass((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(pontocurvass((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,pp)==F(i,1) && PMMs2(1,pp+1)==F(i,2)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,pp+1)==F(i,1) && PMMs2(1,pp)==F(i,2)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,pp)==F(i,3) && PMMs2(1,pp+1)==F(i,2)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,pp+1)==F(i,3) && PMMs2(1,pp)==F(i,2)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(pontocurvass((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(pontocurvass((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,pp)==F(i,1) && PMMs2(1,pp+1)==F(i,3)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,pp+1)==F(i,1) && PMMs2(1,pp)==F(i,3)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,pp)==F(i,2) && PMMs2(1,pp+1)==F(i,3)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,pp+1)==F(i,2) && PMMs2(1,pp)==F(i,3)))
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(pontocurvass((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(pontocurvass((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,pp)==F(i,2) && PMMs2(1,pp+1)==F(i,1)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,pp+1)==F(i,2) && PMMs2(1,pp)==F(i,1)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,pp)==F(i,3) && PMMs2(1,pp+1)==F(i,1)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,pp+1)==F(i,3) && PMMs2(1,pp)==F(i,1)))
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(pontocurvass((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(pontocurvass((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,pp)==F(i,3) && PMMs2(1,pp+1)==F(i,2)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,pp+1)==F(i,3) && PMMs2(1,pp)==F(i,2)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,pp)==F(i,1) && PMMs2(1,pp+1)==F(i,2)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,pp+1)==F(i,1) && PMMs2(1,pp)==F(i,2)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(pontocurvass((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(pontocurvass((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvassa,1)==size(pontocurvass,1))
idxs=1;
end
endwhile

printf('Vermelha 1');

%Juntando 2 curvas subparabolicas vermelhas que o último ponto da curva 2 
%e o primeiro ponto da curva 1 estão na mesma face

#2
idxs=0;
while(idxs!=1)
pontocurvassa=pontocurvass;
for tt=1:2:size(PMMs2,2)
  ppt=size(PMMs2,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMMs2,2)
  pptz=size(PMMs2,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,tt)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,tt)==F(i,3)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,tt)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,tt)==F(i,1)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,tt)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,tt)==F(i,2)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,tt)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,tt)==F(i,3)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,tt)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,tt)==F(i,1)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,tt)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,tt)==F(i,2)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-5);
  PMMs2=deletarColuna(PMMs2,pp-2);
  PMMs2=deletarColuna(PMMs2,pp-2);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvassa,1)==size(pontocurvass,1))
idxs=1;
end
endwhile

printf('Vermelha 2');

%Juntando 2 curvas subparabolicas vermelhas que o último ponto da curva 1 
%e o primeiro ponto da curva 2 estão na mesma face
#3
idxs=0;
while(idxs!=1)
pontocurvassa=pontocurvass;
for pp=1:2:size(PMMs2,2)
  ppt=size(PMMs2,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMMs2,2)
  pptz=size(PMMs2,2);
  if (tt>pptz)
    break
  endif
for i=1:nfaces
  if((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,tt)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,tt)==F(i,3)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  PMMs2=deletarColuna(PMMs2,pp);
  PMMs2=deletarColuna(PMMs2,pp);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  PMMs2=deletarColuna(PMMs2,tt-2);
  PMMs2=deletarColuna(PMMs2,tt-2);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,tt)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,tt)==F(i,1)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  PMMs2=deletarColuna(PMMs2,pp);
  PMMs2=deletarColuna(PMMs2,pp);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  PMMs2=deletarColuna(PMMs2,tt-2);
  PMMs2=deletarColuna(PMMs2,tt-2);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,tt)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,tt)==F(i,2)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  PMMs2=deletarColuna(PMMs2,pp);
  PMMs2=deletarColuna(PMMs2,pp);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  PMMs2=deletarColuna(PMMs2,tt-2);
  PMMs2=deletarColuna(PMMs2,tt-2);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,tt)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,2) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,tt)==F(i,3)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  PMMs2=deletarColuna(PMMs2,pp);
  PMMs2=deletarColuna(PMMs2,pp);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  PMMs2=deletarColuna(PMMs2,tt-2);
  PMMs2=deletarColuna(PMMs2,tt-2);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt+1)==F(i,2) && PMMs2(1,tt)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,3) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,tt)==F(i,1)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  PMMs2=deletarColuna(PMMs2,pp);
  PMMs2=deletarColuna(PMMs2,pp);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  PMMs2=deletarColuna(PMMs2,tt-2);
  PMMs2=deletarColuna(PMMs2,tt-2);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2-1);
  break
  elseif((PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt+1)==F(i,3) && PMMs2(1,tt)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2)) || (PMMs2(pontocurvass((pp+1)/2,1),pp)==F(i,1) && PMMs2(pontocurvass((pp+1)/2,1),pp+1)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(1,tt)==F(i,2)) )
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  PMMs2=deletarColuna(PMMs2,pp);
  PMMs2=deletarColuna(PMMs2,pp);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  PMMs2=deletarColuna(PMMs2,tt-2);
  PMMs2=deletarColuna(PMMs2,tt-2);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvassa,1)==size(pontocurvass,1))
idxs=1;
end
endwhile

printf('Vermelha 3');

%Fazendo a curva subparabolica vermelha, que o último ponto da curva está na 
%mesma face de um vértice que é subparabolica vermelha, passar por esse vértice
#4
idxs=0;
while(idxs!=1)
mxs=0;
pontocurvassa=pontocurvass;
for tt=1:2:size(PMMs2,2)
for pp=1:size(zerosub,1)
for i=1:nfaces
  if((PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2) && zerosub(pp,1)==F(i,3)))
  MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-2:end)=V(guardars(pp,1:nguardars(pp,1)),:);
  PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1)=guardars(pp,1:nguardars(pp,1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end)=guardars(pp,1:nguardars(pp,1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3) && zerosub(pp,1)==F(i,1)))
  MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-2:end)=V(guardars(pp,1:nguardars(pp,1)),:);
  PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1)=guardars(pp,1:nguardars(pp,1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end)=guardars(pp,1:nguardars(pp,1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1) && zerosub(pp,1)==F(i,2)))
  MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-2:end)=V(guardars(pp,1:nguardars(pp,1)),:);
  PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1)=guardars(pp,1:nguardars(pp,1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end)=guardars(pp,1:nguardars(pp,1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1) && zerosub(pp,1)==F(i,3)))
  MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-2:end)=V(guardars(pp,1:nguardars(pp,1)),:);
  PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1)=guardars(pp,1:nguardars(pp,1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end)=guardars(pp,1:nguardars(pp,1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2) && zerosub(pp,1)==F(i,1)))
  MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-2:end)=V(guardars(pp,1:nguardars(pp,1)),:);
  PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1)=guardars(pp,1:nguardars(pp,1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end)=guardars(pp,1:nguardars(pp,1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3) && zerosub(pp,1)==F(i,2)))
  MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-2:end)=V(guardars(pp,1:nguardars(pp,1)),:);
  PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1)=guardars(pp,1:nguardars(pp,1));
  PMMs2(pontocurvass((tt+1)/2,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end)=guardars(pp,1:nguardars(pp,1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
endif
endfor
endfor
endfor
if(mxs==0)
idxs=1;
end
endwhile

printf('Vermelha 4');

%Fazendo com que a curva subparabolica vermelha passe pelo ponto umbílico
% Caso em que o ponto umbílico se originou de 2 pontos umbílicos ligados


abcs=size(PMMs2,2);
for j=1:size(ligadosr,1)
  if (numlig(j) ==2)
    for tt=1:2:abcs
      idxs=0;
      for i=1:pontocurvass((tt+1)/2)
        if (idxs == 1)
          break;
        endif
        if (PMMs2(i,tt)==aa(1,j) ||PMMs2(i,tt)==aaa(1,j) || PMMs2(i,tt+1)==aa(1,j) ||PMMs2(i,tt+1)==aaa(1,j) )
          if   ((PMMs2(i,tt)==aa(1,j) && PMMs2(i,tt+1)==aaa(1,j)) || (PMMs2(i,tt+1)==aa(1,j) && PMMs2(i,tt)==aaa(1,j)) )
            break
          else
        for k=1:nfaces
          if ((PMMs2(i,tt)==F(k,1) && PMMs2(i,tt+1)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3)) || (PMMs2(i,tt)==F(k,1) && PMMs2(i,tt+1)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt)==F(k,2) && PMMs2(i,tt+1)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1)) || (PMMs2(i,tt)==F(k,2) && PMMs2(i,tt+1)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt)==F(k,3) && PMMs2(i,tt+1)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2)) || (PMMs2(i,tt)==F(k,3) && PMMs2(i,tt+1)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt)==F(k,2) && PMMs2(i,tt+1)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3)) || (PMMs2(i,tt)==F(k,2) && PMMs2(i,tt+1)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt)==F(k,1) && PMMs2(i,tt+1)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2)) || (PMMs2(i,tt)==F(k,1) && PMMs2(i,tt+1)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt)==F(k,3) && PMMs2(i,tt+1)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1)) || (PMMs2(i,tt)==F(k,3) && PMMs2(i,tt+1)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
            elseif ((PMMs2(i,tt+1)==F(k,1) && PMMs2(i,tt)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3)) || (PMMs2(i,tt+1)==F(k,1) && PMMs2(i,tt)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt+1)==F(k,2) && PMMs2(i,tt)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1)) || (PMMs2(i,tt+1)==F(k,2) && PMMs2(i,tt)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt+1)==F(k,3) && PMMs2(i,tt)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2)) || (PMMs2(i,tt+1)==F(k,3) && PMMs2(i,tt)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt+1)==F(k,2) && PMMs2(i,tt)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3)) || (PMMs2(i,tt+1)==F(k,2) && PMMs2(i,tt)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt+1)==F(k,1) && PMMs2(i,tt)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2)) || (PMMs2(i,tt+1)==F(k,1) && PMMs2(i,tt)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2(i,tt+1)==F(k,3) && PMMs2(i,tt)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1)) || (PMMs2(i,tt+1)==F(k,3) && PMMs2(i,tt)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1)))
            MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvass((tt+1)/2))
            MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
            PMMs2(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvass((tt+1)/2))
            PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
            end
            pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
            PMMs2=deletarColuna(PMMs2,tt);
            PMMs2=deletarColuna(PMMs2,tt);
            pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
            idxs=1;
            break
          endif
    endfor
  endif
    endif
  endfor
endfor
endif
endfor
% Caso em que o ponto umbílico se originou de 3 ou mais pontos umbílicos ligados

for j=1:size(ligadosr,1)
  if (numlig(j) >=3)
    for tt=1:2:size(PMMs2,2)
      idxs=0;
      for i=1:pontocurvass((tt+1)/2)-1
        if (idxs==1)
            break
          endif
         if (PMMs2(i,tt)==aa(1,j) && PMMs2(i,tt+1)==aaa(1,j) && ((PMMs2(i+1,tt)==aaa(1,j) && PMMs2(i+1,tt+1)==aaaa(1,j)) | (PMMs2(i+1,tt)==aaaa(1,j) && PMMs2(i+1,tt+1)==aaa(1,j)) | (PMMs2(i+1,tt)==aa(1,j) && PMMs2(i+1,tt+1)==aaaa(1,j)) | (PMMs2(i+1,tt)==aaaa(1,j) && PMMs2(i+1,tt+1)==aa(1,j))))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          elseif (PMMs2(i,tt)==aaa(1,j) && PMMs2(i,tt+1)==aa(1,j)&& ((PMMs2(i+1,tt)==aaa(1,j) && PMMs2(i+1,tt+1)==aaaa(1,j)) | (PMMs2(i+1,tt)==aaaa(1,j) && PMMs2(i+1,tt+1)==aaa(1,j)) | (PMMs2(i+1,tt)==aa(1,j) && PMMs2(i+1,tt+1)==aaaa(1,j)) | (PMMs2(i+1,tt)==aaaa(1,j) && PMMs2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(i,tt)==aa(1,j) && PMMs2(i,tt+1)==aaaa(1,j)&& ((PMMs2(i+1,tt)==aaa(1,j) && PMMs2(i+1,tt+1)==aaaa(1,j)) | (PMMs2(i+1,tt)==aaaa(1,j) && PMMs2(i+1,tt+1)==aaa(1,j)) | (PMMs2(i+1,tt)==aa(1,j) && PMMs2(i+1,tt+1)==aaa(1,j)) | (PMMs2(i+1,tt)==aaa(1,j) && PMMs2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(i,tt)==aaaa(1,j) && PMMs2(i,tt+1)==aa(1,j) && ((PMMs2(i+1,tt)==aaa(1,j) && PMMs2(i+1,tt+1)==aaaa(1,j)) | (PMMs2(i+1,tt)==aaaa(1,j) && PMMs2(i+1,tt+1)==aaa(1,j)) | (PMMs2(i+1,tt)==aa(1,j) && PMMs2(i+1,tt+1)==aaa(1,j)) | (PMMs2(i+1,tt)==aaa(1,j) && PMMs2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(i,tt)==aaaa(1,j) && PMMs2(i,tt+1)==aaa(1,j) && ((PMMs2(i+1,tt)==aa(1,j) && PMMs2(i+1,tt+1)==aaaa(1,j)) | (PMMs2(i+1,tt)==aaaa(1,j) && PMMs2(i+1,tt+1)==aa(1,j)) | (PMMs2(i+1,tt)==aa(1,j) && PMMs2(i+1,tt+1)==aaa(1,j)) | (PMMs2(i+1,tt)==aaa(1,j) && PMMs2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(i,tt)==aaa(1,j) && PMMs2(i,tt+1)==aaaa(1,j) && ((PMMs2(i+1,tt)==aa(1,j) && PMMs2(i+1,tt+1)==aaaa(1,j)) | (PMMs2(i+1,tt)==aaaa(1,j) && PMMs2(i+1,tt+1)==aa(1,j)) | (PMMs2(i+1,tt)==aa(1,j) && PMMs2(i+1,tt+1)==aaa(1,j)) | (PMMs2(i+1,tt)==aaa(1,j) && PMMs2(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2(1:i,end+1:end+3)=MMs2(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvass((tt+1)/2))
              MMs2(i+2:pontocurvass((tt+1)/2)+1,end-2:end)=MMs2(i+1:pontocurvass((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2(1:i,end+1:end+2)=PMMs2(1:i,tt:(tt+1));
              PMMs2(i+1,end-1)=nverts+10;
              PMMs2(i+1,end)=nverts+10;
              if (i+1<=pontocurvass((tt+1)/2))
              PMMs2(i+2:pontocurvass((tt+1)/2)+1,end-1:end)=PMMs2(i+1:pontocurvass((tt+1)/2),tt:(tt+1));
              end
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          elseif (PMMs2(1,tt)==aa(1,j) && PMMs2(1,tt+1)==aaa(1,j) && (PMMs2(2,tt)!=aaaa(1,j) | PMMs2(2,tt+1)!=aaaa(1,j)))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          elseif (PMMs2(1,tt)==aaa(1,j) && PMMs2(1,tt+1)==aa(1,j) && (PMMs2(2,tt)!=aaaa(1,j) | PMMs2(2,tt+1)!=aaaa(1,j)))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(1,tt)==aa(1,j) && PMMs2(1,tt+1)==aaaa(1,j) && (PMMs2(2,tt)!=aaa(1,j) | PMMs2(2,tt+1)!=aaa(1,j)))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(1,tt)==aaaa(1,j) && PMMs2(1,tt+1)==aa(1,j) && (PMMs2(2,tt)!=aaa(1,j) | PMMs2(2,tt+1)!=aaa(1,j)))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(1,tt)==aaaa(1,j) && PMMs2(1,tt+1)==aaa(1,j) && (PMMs2(2,tt)!=aa(1,j) | PMMs2(2,tt+1)!=aa(1,j)))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(1,tt)==aaa(1,j) && PMMs2(1,tt+1)==aaaa(1,j) && (PMMs2(2,tt)!=aa(1,j) | PMMs2(2,tt+1)!=aa(1,j)))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2(1,end+1:end+3)=Pumbc(j,:);
              MMs2(2:pontocurvass((tt+1)/2,1)+1,end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2(1,end+1)=nverts+10;
              PMMs2(1,end+1)=nverts+10;
              PMMs2(2:pontocurvass((tt+1)/2,1)+1,end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          elseif (PMMs2(pontocurvass((tt+1)/2,1),tt)==aa(1,j) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==aaa(1,j))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          elseif (PMMs2(pontocurvass((tt+1)/2,1),tt)==aaa(1,j) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==aa(1,j))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(pontocurvass((tt+1)/2,1),tt)==aa(1,j) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==aaaa(1,j))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(pontocurvass((tt+1)/2,1),tt)==aaaa(1,j) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==aa(1,j))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(pontocurvass((tt+1)/2,1),tt)==aaaa(1,j) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==aaa(1,j))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2(pontocurvass((tt+1)/2,1),tt)==aaa(1,j) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==aaaa(1,j))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
          elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2(1:pontocurvass((tt+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2(pontocurvass((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2(1:pontocurvass((tt+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
              PMMs2(pontocurvass((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2(pontocurvass((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+1;
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
              PMMs2=deletarColuna(PMMs2,tt);
              PMMs2=deletarColuna(PMMs2,tt);
              pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          endif
      endfor
    endfor
  endif
endfor

printf('Vermelha 5');

%Fazendo a curva subparabolica vermelha, que o primeiro ponto da curva está na 
%mesma face de um vértice que é subparabolica vermelha, passar por esse vértice
idxs=0;
while(idxs!=1)
mxs=0;
pontocurvassa=pontocurvass;
for tt=1:2:size(PMMs2,2)
if (pontocurvass((tt+1)/2,1)>=2)
for pp=1:size(zerosub,1)
for i=1:nfaces
  if((PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && zerosub(pp,1)==F(i,3) && zerosub(pp,1)!=PMMs2(2,tt)) )
  MMs2(1:nguardars(pp,1),end+1:end+3)=V(guardars(pp,nguardars(pp,1):-1:1),:);
  MMs2(nguardars(pp,1)+1:nguardars(pp,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(nguardars(pp,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && zerosub(pp,1)==F(i,1) && zerosub(pp,1)!=PMMs2(2,tt)))
  MMs2(1:nguardars(pp,1),end+1:end+3)=V(guardars(pp,nguardars(pp,1):-1:1),:);
  MMs2(nguardars(pp,1)+1:nguardars(pp,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(nguardars(pp,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && zerosub(pp,1)==F(i,2) && zerosub(pp,1)!=PMMs2(2,tt)))
  MMs2(1:nguardars(pp,1),end+1:end+3)=V(guardars(pp,nguardars(pp,1):-1:1),:);
  MMs2(nguardars(pp,1)+1:nguardars(pp,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(nguardars(pp,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && zerosub(pp,1)==F(i,3) && zerosub(pp,1)!=PMMs2(2,tt)))
  MMs2(1:nguardars(pp,1),end+1:end+3)=V(guardars(pp,nguardars(pp,1):-1:1),:);
  MMs2(nguardars(pp,1)+1:nguardars(pp,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(nguardars(pp,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && zerosub(pp,1)==F(i,1) && zerosub(pp,1)!=PMMs2(2,tt)))
  MMs2(1:nguardars(pp,1),end+1:end+3)=V(guardars(pp,nguardars(pp,1):-1:1),:);
  MMs2(nguardars(pp,1)+1:nguardars(pp,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(nguardars(pp,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && zerosub(pp,1)==F(i,2) && zerosub(pp,1)!=PMMs2(2,tt)) )
  MMs2(1:nguardars(pp,1),end+1:end+3)=V(guardars(pp,nguardars(pp,1):-1:1),:);
  MMs2(nguardars(pp,1)+1:nguardars(pp,1)+pontocurvass((tt+1)/2,1),end-2:end)=MMs2(1:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(1:nguardars(pp,1),end+1)=guardars(pp,nguardars(pp,1):-1:1);
  PMMs2(nguardars(pp,1)+1:pontocurvass((tt+1)/2,1)+nguardars(pp,1),end-1:end)=PMMs2(1:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+nguardars(pp,1);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-2);
  PMMs2=deletarColuna(PMMs2,tt);
  PMMs2=deletarColuna(PMMs2,tt);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2);
  mxs=1;
  break
endif
endfor
endfor
endif
endfor
if(mxs==0)
idxs=1;
endif
endwhile

printf('Vermelha 6');

%Juntando 2 curvas subparabolicas vermelhas que o último ponto da curva 1 
%e o primeiro ponto da curva 2 são os mesmos pontos

idxs=0;
while(idxs!=1)
pontocurvassa=pontocurvass;
for pp=1:2:size(PMMs2,2)
  ppt=size(PMMs2,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMMs2,2)
  pptz=size(PMMs2,2);
  if (tt>pptz)
    break
  endif
  if(MMs2(pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2)==MMs2(1,3*(tt+1)/2-2:3*(tt+1)/2))
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(1:pontocurvass((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1)-1,end-2:end)=MMs2(2:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(1:pontocurvass((pp+1)/2,1),pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((pp+1)/2,1)+pontocurvass((tt+1)/2,1)-1,end-1:end)=PMMs2(2:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1)-1;
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  PMMs2=deletarColuna(PMMs2,pp);
  PMMs2=deletarColuna(PMMs2,pp);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  PMMs2=deletarColuna(PMMs2,tt-2);
  PMMs2=deletarColuna(PMMs2,tt-2);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2-1);
endif
endfor
endfor
if(size(pontocurvassa,1)==size(pontocurvass,1))
idxs=1;
end
endwhile

%Juntando 2 curvas subparabolicas vermelhas que o primeiro ponto da curva 1 
%e o primeiro ponto da curva 2 são os mesmos pontos

idxs=0;
while(idxs!=1)
pontocurvassa=pontocurvass;
for pp=1:2:size(PMMs2,2)
  ppt=size(PMMs2,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMMs2,2)
  pptz=size(PMMs2,2);
  if (tt>pptz)
    break
  endif
  if (MMs2(1,3*(pp+1)/2-2:3*(pp+1)/2)==MMs2(1,3*(tt+1)/2-2:3*(tt+1)/2))
  MMs2(1:pontocurvass((pp+1)/2,1),end+1:end+3)=MMs2(pontocurvass((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1)-1,end-2:end)=MMs2(2:pontocurvass((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2(1:pontocurvass((pp+1)/2,1),end+1:end+2)=PMMs2(pontocurvass((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2(pontocurvass((pp+1)/2,1)+1:pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1)-1,end-1:end)=PMMs2(2:pontocurvass((tt+1)/2,1),tt:(tt+1));
  pontocurvass(end+1,1)=pontocurvass((tt+1)/2,1)+pontocurvass((pp+1)/2,1)-1;
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  MMs2=deletarColuna(MMs2,3*(pp+1)/2-2);
  PMMs2=deletarColuna(PMMs2,pp);
  PMMs2=deletarColuna(PMMs2,pp);
  pontocurvass=deletarLinha(pontocurvass,(pp+1)/2);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  MMs2=deletarColuna(MMs2,3*(tt+1)/2-5);
  PMMs2=deletarColuna(PMMs2,tt-2);
  PMMs2=deletarColuna(PMMs2,tt-2);
  pontocurvass=deletarLinha(pontocurvass,(tt+1)/2-1);
endif
endfor
endfor
if(size(pontocurvassa,1)==size(pontocurvass,1))
idxs=1;
end
endwhile

if(size(pontocurvassax,1)==size(pontocurvass,1))
idxs2=1;
end
endwhile

%Verificando se há curvas subparabolicas vermelhas degeneradas

nPumbsems=[];
  MMs2c=[];
  pontodegenerados=[];
  counts=0;
  ptds=0;
  for i=1:size(Pumbc)
    for j=1:size(pontocurvass)(1)
      MMs2(1:pontocurvass(j,1),3*j-2:3*j)==Pumbc(i,:);
      tils = ans ;
      if (max(tils(:,1)) == 1 && max(tils(:,2)) == 1 && max(tils(:,3)) == 1)
        if (numlig(i)>=1 && MMs2(1,3*j-2:3*j)== MMs2(2:pontocurvass(j,1),3*j-2:3*j))
          ptds=ptds+1;
          pontodegenerados(ptds,:)=i;
        endif
        counts=counts+1;
        nPumbsems(counts,:)=i;
        break
      endif
    endfor
  endfor

  Pumbsems=Pumbc;
  uus=0;
  uus(1:size(Pumbc,1),1)=[1:size(Pumbc,1)];
  
  for i=size(nPumbsems,1):-1:1
    Pumbsems=deletarLinha(Pumbsems,nPumbsems(i));
    uus=deletarLinha(uus,nPumbsems(i));
  endfor

% Fechando as curvas subparabolicas vermelhas fechadas

for tt=1:2:size(PMMs2,2)
for i=1:nfaces
if((PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3)) )
  MMs2(pontocurvass((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvass((tt+1)/2,1)=pontocurvass((tt+1)/2,1)+1;
  break
  elseif((PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1)) )
  MMs2(pontocurvass((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvass((tt+1)/2,1)=pontocurvass((tt+1)/2,1)+1;
  break
  elseif((PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2)) )
  MMs2(pontocurvass((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvass((tt+1)/2,1)=pontocurvass((tt+1)/2,1)+1;
  break
  elseif((PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3)) || (PMMs2(1,tt)==F(i,2) && PMMs2(1,tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3)) )
  MMs2(pontocurvass((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvass((tt+1)/2,1)=pontocurvass((tt+1)/2,1)+1;
  break
  elseif((PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1)) || (PMMs2(1,tt)==F(i,3) && PMMs2(1,tt+1)==F(i,2) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1)) )
  MMs2(pontocurvass((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvass((tt+1)/2,1)=pontocurvass((tt+1)/2,1)+1;
  break
  elseif((PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,2)) || (PMMs2(1,tt)==F(i,1) && PMMs2(1,tt+1)==F(i,3) && PMMs2(pontocurvass((tt+1)/2,1),tt+1)==F(i,1) && PMMs2(pontocurvass((tt+1)/2,1),tt)==F(i,2)) )
  MMs2(pontocurvass((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvass((tt+1)/2,1)=pontocurvass((tt+1)/2,1)+1;
  break
endif
endfor
endfor


printf('Vermelha fecho');


####### Azul ######
##

%Juntando 2 curvas subparabolicas azuis que o primeiro ponto da curva 1 
%e o primeiro ponto da curva 2 estão na mesma face

idxs2=0;
idxs3=0;
while(idxs2!=1)
pontocurvassaxw=pontocurvassw;
idxs=0;
while(idxs!=1)
pontocurvassaw=pontocurvassw;
for tt=1:2:size(PMMs2w,2)
  ppt=size(PMMs2w,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMMs2w,2)
  pptz=size(PMMs2w,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,pp)==F(i,2) && PMMs2w(1,pp+1)==F(i,3)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,pp+1)==F(i,2) && PMMs2w(1,pp)==F(i,3)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,pp)==F(i,1) && PMMs2w(1,pp+1)==F(i,3)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,pp+1)==F(i,1) && PMMs2w(1,pp)==F(i,3)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(pontocurvassw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(pontocurvassw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,pp)==F(i,3) && PMMs2w(1,pp+1)==F(i,1)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,pp+1)==F(i,3) && PMMs2w(1,pp)==F(i,1)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,pp)==F(i,2) && PMMs2w(1,pp+1)==F(i,1)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,pp+1)==F(i,2) && PMMs2w(1,pp)==F(i,1)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(pontocurvassw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(pontocurvassw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,pp)==F(i,1) && PMMs2w(1,pp+1)==F(i,2)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,pp+1)==F(i,1) && PMMs2w(1,pp)==F(i,2)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,pp)==F(i,3) && PMMs2w(1,pp+1)==F(i,2)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,pp+1)==F(i,3) && PMMs2w(1,pp)==F(i,2)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(pontocurvassw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(pontocurvassw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,pp)==F(i,1) && PMMs2w(1,pp+1)==F(i,3)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,pp+1)==F(i,1) && PMMs2w(1,pp)==F(i,3)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,pp)==F(i,2) && PMMs2w(1,pp+1)==F(i,3)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,pp+1)==F(i,2) && PMMs2w(1,pp)==F(i,3)))
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(pontocurvassw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(pontocurvassw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,pp)==F(i,2) && PMMs2w(1,pp+1)==F(i,1)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,pp+1)==F(i,2) && PMMs2w(1,pp)==F(i,1)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,pp)==F(i,3) && PMMs2w(1,pp+1)==F(i,1)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,pp+1)==F(i,3) && PMMs2w(1,pp)==F(i,1)))
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(pontocurvassw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(pontocurvassw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,pp)==F(i,3) && PMMs2w(1,pp+1)==F(i,2)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,pp+1)==F(i,3) && PMMs2w(1,pp)==F(i,2)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,pp)==F(i,1) && PMMs2w(1,pp+1)==F(i,2)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,pp+1)==F(i,1) && PMMs2w(1,pp)==F(i,2)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(pontocurvassw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(pontocurvassw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvassaw,1)==size(pontocurvassw,1))
idxs=1;
end
endwhile

printf('Azul 1');

%Juntando 2 curvas subparabolicas azuis que o último ponto da curva 2 
%e o primeiro ponto da curva 1 estão na mesma face

#2
idxs=0;
while(idxs!=1)
pontocurvassaw=pontocurvassw;
for tt=1:2:size(PMMs2w,2)
  ppt=size(PMMs2w,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMMs2w,2)
  pptz=size(PMMs2w,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,tt)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,tt)==F(i,3)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,tt)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,tt)==F(i,1)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,tt)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,tt)==F(i,2)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,tt)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,tt)==F(i,3)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,tt)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,tt)==F(i,1)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,tt)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,tt)==F(i,2)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  PMMs2w=deletarColuna(PMMs2w,pp-2);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvassaw,1)==size(pontocurvassw,1))
idxs=1;
end
endwhile

printf('Azul 2');

%Juntando 2 curvas subparabolicas azuis que o último ponto da curva 1 
%e o primeiro ponto da curva 2 estão na mesma face
idxs=0;
while(idxs!=1)
pontocurvassaw=pontocurvassw;
for pp=1:2:size(PMMs2w,2)
  ppt=size(PMMs2w,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMMs2w,2)
  pptz=size(PMMs2w,2);
  if (tt>pptz)
    break
  endif
for i=1:nfaces
  if((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,tt)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,tt)==F(i,3)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,pp);
  PMMs2w=deletarColuna(PMMs2w,pp);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,tt)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,tt)==F(i,1)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,pp);
  PMMs2w=deletarColuna(PMMs2w,pp);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,tt)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,tt)==F(i,2)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,pp);
  PMMs2w=deletarColuna(PMMs2w,pp);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,tt)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,2) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,tt)==F(i,3)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,pp);
  PMMs2w=deletarColuna(PMMs2w,pp);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(1,tt)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,3) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,tt)==F(i,1)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,pp);
  PMMs2w=deletarColuna(PMMs2w,pp);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2-1);
  break
  elseif((PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(1,tt)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2)) || (PMMs2w(pontocurvassw((pp+1)/2,1),pp)==F(i,1) && PMMs2w(pontocurvassw((pp+1)/2,1),pp+1)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(1,tt)==F(i,2)) )
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,pp);
  PMMs2w=deletarColuna(PMMs2w,pp);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvassaw,1)==size(pontocurvassw,1))
idxs=1;
end
endwhile

printf('Azul 3');

%Fazendo a curva subparabolica azul, que o último ponto da curva está na 
%mesma face de um vértice que é subparabolica azul, passar por esse vértice
#4
idxs=0;
while(idxs!=1)
mxs=0;
for tt=1:2:size(PMMs2w,2)
for pp=1:size(zerosub2,1)
for i=1:nfaces
  if((PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2) && zerosub2(pp,1)==F(i,3)))
  MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-2:end)=V(guardars2(pp,1:nguardars2(pp,1)),:);
  PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1)=guardars2(pp,1:nguardars2(pp,1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end)=guardars2(pp,1:nguardars2(pp,1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3) && zerosub2(pp,1)==F(i,1)))
  MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-2:end)=V(guardars2(pp,1:nguardars2(pp,1)),:);
  PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1)=guardars2(pp,1:nguardars2(pp,1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end)=guardars2(pp,1:nguardars2(pp,1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1) && zerosub2(pp,1)==F(i,2)))
  MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-2:end)=V(guardars2(pp,1:nguardars2(pp,1)),:);
  PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1)=guardars2(pp,1:nguardars2(pp,1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end)=guardars2(pp,1:nguardars2(pp,1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1) && zerosub2(pp,1)==F(i,3)))
  MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-2:end)=V(guardars2(pp,1:nguardars2(pp,1)),:);
  PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1)=guardars2(pp,1:nguardars2(pp,1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end)=guardars2(pp,1:nguardars2(pp,1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2) && zerosub2(pp,1)==F(i,1)))
  MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-2:end)=V(guardars2(pp,1:nguardars2(pp,1)),:);
  PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1)=guardars2(pp,1:nguardars2(pp,1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end)=guardars2(pp,1:nguardars2(pp,1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3) && zerosub2(pp,1)==F(i,2)))
  MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-2:end)=V(guardars2(pp,1:nguardars2(pp,1)),:);
  PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1)=guardars2(pp,1:nguardars2(pp,1));
  PMMs2w(pontocurvassw((tt+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end)=guardars2(pp,1:nguardars2(pp,1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
endif
endfor
endfor
endfor
if(mxs==0)
idxs=1;
end
endwhile

printf('Azul 4');

%Fazendo com que a curva subparabolica azul passe pelo ponto umbílico
% Caso em que o ponto umbílico se originou de 2 pontos umbílicos ligados


abcs=size(PMMs2w,2);
for j=1:size(ligadosr,1)
  if (numlig(j) ==2)
    for tt=1:2:abcs
      idxs=0;
      for i=1:pontocurvassw((tt+1)/2)
        if (idxs == 1)
          break;
        endif
        if (PMMs2w(i,tt)==aa(1,j) ||PMMs2w(i,tt)==aaa(1,j) || PMMs2w(i,tt+1)==aa(1,j) ||PMMs2w(i,tt+1)==aaa(1,j) )
          if   ((PMMs2w(i,tt)==aa(1,j) && PMMs2w(i,tt+1)==aaa(1,j)) || (PMMs2w(i,tt+1)==aa(1,j) && PMMs2w(i,tt)==aaa(1,j)) )
            break
          else
        for k=1:nfaces
          if ((PMMs2w(i,tt)==F(k,1) && PMMs2w(i,tt+1)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3)) || (PMMs2w(i,tt)==F(k,1) && PMMs2w(i,tt+1)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt)==F(k,2) && PMMs2w(i,tt+1)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1)) || (PMMs2w(i,tt)==F(k,2) && PMMs2w(i,tt+1)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt)==F(k,3) && PMMs2w(i,tt+1)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2)) || (PMMs2w(i,tt)==F(k,3) && PMMs2w(i,tt+1)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt)==F(k,2) && PMMs2w(i,tt+1)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3)) || (PMMs2w(i,tt)==F(k,2) && PMMs2w(i,tt+1)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt)==F(k,1) && PMMs2w(i,tt+1)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2)) || (PMMs2w(i,tt)==F(k,1) && PMMs2w(i,tt+1)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt)==F(k,3) && PMMs2w(i,tt+1)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1)) || (PMMs2w(i,tt)==F(k,3) && PMMs2w(i,tt+1)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
            elseif ((PMMs2w(i,tt+1)==F(k,1) && PMMs2w(i,tt)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3)) || (PMMs2w(i,tt+1)==F(k,1) && PMMs2w(i,tt)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt+1)==F(k,2) && PMMs2w(i,tt)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1)) || (PMMs2w(i,tt+1)==F(k,2) && PMMs2w(i,tt)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt+1)==F(k,3) && PMMs2w(i,tt)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2)) || (PMMs2w(i,tt+1)==F(k,3) && PMMs2w(i,tt)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt+1)==F(k,2) && PMMs2w(i,tt)==F(k,1) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3)) || (PMMs2w(i,tt+1)==F(k,2) && PMMs2w(i,tt)==F(k,1) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt+1)==F(k,1) && PMMs2w(i,tt)==F(k,3) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2)) || (PMMs2w(i,tt+1)==F(k,1) && PMMs2w(i,tt)==F(k,3) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          elseif ((PMMs2w(i,tt+1)==F(k,3) && PMMs2w(i,tt)==F(k,2) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1)) || (PMMs2w(i,tt+1)==F(k,3) && PMMs2w(i,tt)==F(k,2) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1)))
            MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
            MMs2w(i+1,end-2:end)=Pumbc(j,:);
            if (i+1<=pontocurvassw((tt+1)/2))
            MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
            end
            PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
            PMMs2w(i+1,end-1:end)=nverts+10;
            if (i+1<=pontocurvassw((tt+1)/2))
            PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
            end
            pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
            PMMs2w=deletarColuna(PMMs2w,tt);
            PMMs2w=deletarColuna(PMMs2w,tt);
            pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
            idxs=1;
            break
          endif
    endfor
  endif
    endif
  endfor
endfor
endif
endfor

% Caso em que o ponto umbílico se originou de 3 ou mais pontos umbílicos ligados

for j=1:size(ligadosr,1)
  if (numlig(j) >=3)
    for tt=1:2:size(PMMs2w,2)
      idxs=0;
      for i=1:pontocurvassw((tt+1)/2)
        if (idxs==1)
            break
          endif
         if (PMMs2w(i,tt)==aa(1,j) && PMMs2w(i,tt+1)==aaa(1,j) && ((PMMs2w(i+1,tt)==aaa(1,j) && PMMs2w(i+1,tt+1)==aaaa(1,j)) | (PMMs2w(i+1,tt)==aaaa(1,j) && PMMs2w(i+1,tt+1)==aaa(1,j)) | (PMMs2w(i+1,tt)==aa(1,j) && PMMs2w(i+1,tt+1)==aaaa(1,j)) | (PMMs2w(i+1,tt)==aaaa(1,j) && PMMs2w(i+1,tt+1)==aa(1,j))))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          elseif (PMMs2w(i,tt)==aaa(1,j) && PMMs2w(i,tt+1)==aa(1,j) && ((PMMs2w(i+1,tt)==aaa(1,j) && PMMs2w(i+1,tt+1)==aaaa(1,j)) | (PMMs2w(i+1,tt)==aaaa(1,j) && PMMs2w(i+1,tt+1)==aaa(1,j)) | (PMMs2w(i+1,tt)==aa(1,j) && PMMs2w(i+1,tt+1)==aaaa(1,j)) | (PMMs2w(i+1,tt)==aaaa(1,j) && PMMs2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(i,tt)==aa(1,j) && PMMs2w(i,tt+1)==aaaa(1,j) && ((PMMs2w(i+1,tt)==aaa(1,j) && PMMs2w(i+1,tt+1)==aaaa(1,j)) | (PMMs2w(i+1,tt)==aaaa(1,j) && PMMs2w(i+1,tt+1)==aaa(1,j)) | (PMMs2w(i+1,tt)==aa(1,j) && PMMs2w(i+1,tt+1)==aaa(1,j)) | (PMMs2w(i+1,tt)==aaa(1,j) && PMMs2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(i,tt)==aaaa(1,j) && PMMs2w(i,tt+1)==aa(1,j) && ((PMMs2w(i+1,tt)==aaa(1,j) && PMMs2w(i+1,tt+1)==aaaa(1,j)) | (PMMs2w(i+1,tt)==aaaa(1,j) && PMMs2w(i+1,tt+1)==aaa(1,j)) | (PMMs2w(i+1,tt)==aa(1,j) && PMMs2w(i+1,tt+1)==aaa(1,j)) | (PMMs2w(i+1,tt)==aaa(1,j) && PMMs2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(i,tt)==aaaa(1,j) && PMMs2w(i,tt+1)==aaa(1,j) && ((PMMs2w(i+1,tt)==aa(1,j) && PMMs2w(i+1,tt+1)==aaaa(1,j)) | (PMMs2w(i+1,tt)==aaaa(1,j) && PMMs2w(i+1,tt+1)==aa(1,j)) | (PMMs2w(i+1,tt)==aa(1,j) && PMMs2w(i+1,tt+1)==aaa(1,j)) | (PMMs2w(i+1,tt)==aaa(1,j) && PMMs2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(i,tt)==aaa(1,j) && PMMs2w(i,tt+1)==aaaa(1,j) && ((PMMs2w(i+1,tt)==aa(1,j) && PMMs2w(i+1,tt+1)==aaaa(1,j)) | (PMMs2w(i+1,tt)==aaaa(1,j) && PMMs2w(i+1,tt+1)==aa(1,j)) | (PMMs2w(i+1,tt)==aa(1,j) && PMMs2w(i+1,tt+1)==aaa(1,j)) | (PMMs2w(i+1,tt)==aaa(1,j) && PMMs2w(i+1,tt+1)==aa(1,j))))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2w(1:i,end+1:end+3)=MMs2w(1:i,3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(i+1,end-2:end)=Pumbc(j,:);
              if (i+1<=pontocurvassw((tt+1)/2))
              MMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(i+1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              end
              PMMs2w(1:i,end+1:end+2)=PMMs2w(1:i,tt:(tt+1));
              PMMs2w(i+1,end-1)=nverts+10;
              PMMs2w(i+1,end)=nverts+10;
              if (i+1<=pontocurvassw((tt+1)/2))
              PMMs2w(i+2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(i+1:pontocurvassw((tt+1)/2),tt:(tt+1));
              end
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          elseif (PMMs2w(1,tt)==aa(1,j) && PMMs2w(1,tt+1)==aaa(1,j) && (PMMs2w(2,tt)!=aaaa(1,j) | PMMs2w(2,tt+1)!=aaaa(1,j)))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          elseif (PMMs2w(1,tt)==aaa(1,j) && PMMs2w(1,tt+1)==aa(1,j) && (PMMs2w(2,tt)!=aaaa(1,j) | PMMs2w(2,tt+1)!=aaaa(1,j)))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(1,tt)==aa(1,j) && PMMs2w(1,tt+1)==aaaa(1,j) && (PMMs2w(2,tt)!=aaa(1,j) | PMMs2w(2,tt+1)!=aaa(1,j)))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
          elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(1,tt)==aaaa(1,j) && PMMs2w(1,tt+1)==aa(1,j) && (PMMs2w(2,tt)!=aaa(1,j) | PMMs2w(2,tt+1)!=aaa(1,j)))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(1,tt)==aaaa(1,j) && PMMs2w(1,tt+1)==aaa(1,j) && (PMMs2w(2,tt)!=aa(1,j) | PMMs2w(2,tt+1)!=aa(1,j)))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
          elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(1,tt)==aaa(1,j) && PMMs2w(1,tt+1)==aaaa(1,j) && (PMMs2w(2,tt)!=aa(1,j) | PMMs2w(2,tt+1)!=aa(1,j)))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
          elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
          elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2w(1,end+1:end+3)=Pumbc(j,:);
              MMs2w(2:pontocurvassw((tt+1)/2)+1,end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2),3*(tt+1)/2-2:3*(tt+1)/2);
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(1,end+1)=nverts+10;
              PMMs2w(2:pontocurvassw((tt+1)/2)+1,end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2),tt:(tt+1));
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(pontocurvassw((tt+1)/2,1),tt)==aa(1,j) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==aaa(1,j))
           for k=1:nfaces
            if (aa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          elseif (PMMs2w(pontocurvassw((tt+1)/2,1),tt)==aaa(1,j) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==aa(1,j))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(pontocurvassw((tt+1)/2,1),tt)==aa(1,j) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==aaaa(1,j))
            for k=1:nfaces
          if (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(pontocurvassw((tt+1)/2,1),tt)==aaaa(1,j) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==aa(1,j))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3) && aaa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1) && aaa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2) && aaa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(pontocurvassw((tt+1)/2,1),tt)==aaaa(1,j) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==aaa(1,j))
            for k=1:nfaces
          if (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,1) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,2) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaaa(1,j)==F(k,3) && aaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
            endfor
          elseif (PMMs2w(pontocurvassw((tt+1)/2,1),tt)==aaa(1,j) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==aaaa(1,j))
            for k=1:nfaces
          if (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,1) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,3) && aa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,2) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,3))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,1) && aa(1,j)==F(k,2))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
              elseif (aaa(1,j)==F(k,3) && aaaa(1,j)==F(k,2) && aa(1,j)==F(k,1))
              MMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
              MMs2w(pontocurvassw((tt+1)/2,1)+1,end-2:end)=Pumbc(j,:);
              PMMs2w(1:pontocurvassw((tt+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end-1)=nverts+10;
              PMMs2w(pontocurvassw((tt+1)/2,1)+1,end)=nverts+10;
              pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+1;
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
              PMMs2w=deletarColuna(PMMs2w,tt);
              PMMs2w=deletarColuna(PMMs2w,tt);
              pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
              idxs=1;
              break
            endif
          endfor
          endif
      endfor
    endfor
  endif
endfor
printf('Azul 5');

%Fazendo a curva subparabolica azul, que o primeiro ponto da curva está na 
%mesma face de um vértice que é subparabolica azul, passar por esse vértice
idxs=0;
while(idxs!=1)
mxs=0;
for tt=1:2:size(PMMs2w,2)
if (pontocurvassw((tt+1)/2,1)>=2)
for pp=1:size(zerosub2,1)
for i=1:nfaces
  if((PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && zerosub2(pp,1)==F(i,3) && zerosub2(pp,1)!=PMMs2w(2,tt)) )
  MMs2w(1:nguardars2(pp,1),end+1:end+3)=V(guardars2(pp,nguardars2(pp,1):-1:1),:);
  MMs2w(nguardars2(pp,1)+1:nguardars2(pp,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(nguardars2(pp,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && zerosub2(pp,1)==F(i,1) && zerosub2(pp,1)!=PMMs2w(2,tt)))
  MMs2w(1:nguardars2(pp,1),end+1:end+3)=V(guardars2(pp,nguardars2(pp,1):-1:1),:);
  MMs2w(nguardars2(pp,1)+1:nguardars2(pp,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(nguardars2(pp,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && zerosub2(pp,1)==F(i,2) && zerosub2(pp,1)!=PMMs2w(2,tt)))
  MMs2w(1:nguardars2(pp,1),end+1:end+3)=V(guardars2(pp,nguardars2(pp,1):-1:1),:);
  MMs2w(nguardars2(pp,1)+1:nguardars2(pp,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(nguardars2(pp,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && zerosub2(pp,1)==F(i,3) && zerosub2(pp,1)!=PMMs2w(2,tt)))
  MMs2w(1:nguardars2(pp,1),end+1:end+3)=V(guardars2(pp,nguardars2(pp,1):-1:1),:);
  MMs2w(nguardars2(pp,1)+1:nguardars2(pp,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(nguardars2(pp,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && zerosub2(pp,1)==F(i,1) && zerosub2(pp,1)!=PMMs2w(2,tt)))
  MMs2w(1:nguardars2(pp,1),end+1:end+3)=V(guardars2(pp,nguardars2(pp,1):-1:1),:);
  MMs2w(nguardars2(pp,1)+1:nguardars2(pp,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(nguardars2(pp,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
  elseif((PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && zerosub2(pp,1)==F(i,2) && zerosub2(pp,1)!=PMMs2w(2,tt)) )
  MMs2w(1:nguardars2(pp,1),end+1:end+3)=V(guardars2(pp,nguardars2(pp,1):-1:1),:);
  MMs2w(nguardars2(pp,1)+1:nguardars2(pp,1)+pontocurvassw((tt+1)/2,1),end-2:end)=MMs2w(1:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(1:nguardars2(pp,1),end+1)=guardars2(pp,nguardars2(pp,1):-1:1);
  PMMs2w(nguardars2(pp,1)+1:pontocurvassw((tt+1)/2,1)+nguardars2(pp,1),end-1:end)=PMMs2w(1:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+nguardars2(pp,1);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,tt);
  PMMs2w=deletarColuna(PMMs2w,tt);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2);
  mxs=1;
  break
endif
endfor
endfor
endif
endfor
if(mxs==0)
idxs=1;
end
endwhile

printf('Azul 6');

%Juntando 2 curvas subparabolicas azuis que o último ponto da curva 1 
%e o primeiro ponto da curva 2 são os mesmos pontos

idxs=0;
while(idxs!=1)
pontocurvassaw=pontocurvassw;
for pp=1:2:size(PMMs2w,2)
  ppt=size(PMMs2w,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMMs2w,2)
  pptz=size(PMMs2w,2);
  if (tt>pptz)
    break
  endif
  if(MMs2w(pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2)==MMs2w(1,3*(tt+1)/2-2:3*(tt+1)/2))
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(1:pontocurvassw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1)-1,end-2:end)=MMs2w(2:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(1:pontocurvassw((pp+1)/2,1),pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((pp+1)/2,1)+pontocurvassw((tt+1)/2,1)-1,end-1:end)=PMMs2w(2:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1)-1;
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,pp);
  PMMs2w=deletarColuna(PMMs2w,pp);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2-1);
endif
endfor
endfor
if(size(pontocurvassaw,1)==size(pontocurvassw,1))
idxs=1;
end
endwhile

%Juntando 2 curvas subparabolicas azuis que o primeiro ponto da curva 1 
%e o primeiro ponto da curva 2 são os mesmos pontos

idxs=0;
while(idxs!=1)
pontocurvassaw=pontocurvassw;
for pp=1:2:size(PMMs2w,2)
  ppt=size(PMMs2w,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMMs2w,2)
  pptz=size(PMMs2w,2);
  if (tt>pptz)
    break
  endif
  if (MMs2w(1,3*(pp+1)/2-2:3*(pp+1)/2)==MMs2w(1,3*(tt+1)/2-2:3*(tt+1)/2))
  MMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+3)=MMs2w(pontocurvassw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1)-1,end-2:end)=MMs2w(2:pontocurvassw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMMs2w(1:pontocurvassw((pp+1)/2,1),end+1:end+2)=PMMs2w(pontocurvassw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMMs2w(pontocurvassw((pp+1)/2,1)+1:pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1)-1,end-1:end)=PMMs2w(2:pontocurvassw((tt+1)/2,1),tt:(tt+1));
  pontocurvassw(end+1,1)=pontocurvassw((tt+1)/2,1)+pontocurvassw((pp+1)/2,1)-1;
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  MMs2w=deletarColuna(MMs2w,3*(pp+1)/2-2);
  PMMs2w=deletarColuna(PMMs2w,pp);
  PMMs2w=deletarColuna(PMMs2w,pp);
  pontocurvassw=deletarLinha(pontocurvassw,(pp+1)/2);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  MMs2w=deletarColuna(MMs2w,3*(tt+1)/2-5);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  PMMs2w=deletarColuna(PMMs2w,tt-2);
  pontocurvassw=deletarLinha(pontocurvassw,(tt+1)/2-1);
endif
endfor
endfor

if(size(pontocurvassaw,1)==size(pontocurvassw,1))
idxs=1;
end
endwhile

if(size(pontocurvassaxw,1)==size(pontocurvassw,1))
idxs2=1;
end
endwhile

%Verificando se há curvas subparabolicas azuis degeneradas

nPumbsemsw=[];
  MMs2wc=[];
  pontodegeneradosw=[];
  countsw=0;
  ptdsw=0;
  for i=1:size(Pumbc)
    for j=1:size(pontocurvassw)(1)
      MMs2w(1:pontocurvassw(j,1),3*j-2:3*j)==Pumbc(i,:);
      tilsw = ans ;
      if (max(tilsw(:,1)) == 1 && max(tilsw(:,2)) == 1 && max(tilsw(:,3)) == 1)
        if (numlig(i)>=1 && MMs2w(1,3*j-2:3*j)== MMs2w(2:pontocurvassw(j,1),3*j-2:3*j))
          ptdsw=ptdsw+1;
          pontodegeneradosw(ptdsw,:)=i;
        endif
        countsw=countsw+1;
        nPumbsemsw(countsw,:)=i;
        break
      endif
    endfor
  endfor

  Pumbsemsw=Pumbc;
  uusw=0;
  uusw(1:size(Pumbc,1),1)=[1:size(Pumbc,1)];
  
  for i=size(nPumbsemsw,1):-1:1
    Pumbsemsw=deletarLinha(Pumbsemsw,nPumbsemsw(i));
    uusw=deletarLinha(uusw,nPumbsemsw(i));
  endfor



% Fechando as curvas subparabolicas azuis fechadas

for tt=1:2:size(PMMs2w,2)
for i=1:nfaces
if((PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3)) )
  MMs2w(pontocurvassw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvassw((tt+1)/2,1)=pontocurvassw((tt+1)/2,1)+1;
  break
  elseif((PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1)) )
  MMs2w(pontocurvassw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvassw((tt+1)/2,1)=pontocurvassw((tt+1)/2,1)+1;
  break
  elseif((PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2)) )
  MMs2w(pontocurvassw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvassw((tt+1)/2,1)=pontocurvassw((tt+1)/2,1)+1;
  break
  elseif((PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3)) || (PMMs2w(1,tt)==F(i,2) && PMMs2w(1,tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3)) )
  MMs2w(pontocurvassw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvassw((tt+1)/2,1)=pontocurvassw((tt+1)/2,1)+1;
  break
  elseif((PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1)) || (PMMs2w(1,tt)==F(i,3) && PMMs2w(1,tt+1)==F(i,2) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1)) )
  MMs2w(pontocurvassw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvassw((tt+1)/2,1)=pontocurvassw((tt+1)/2,1)+1;
  break
  elseif((PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,2)) || (PMMs2w(1,tt)==F(i,1) && PMMs2w(1,tt+1)==F(i,3) && PMMs2w(pontocurvassw((tt+1)/2,1),tt+1)==F(i,1) && PMMs2w(pontocurvassw((tt+1)/2,1),tt)==F(i,2)) )
  MMs2w(pontocurvassw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MMs2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvassw((tt+1)/2,1)=pontocurvassw((tt+1)/2,1)+1;
  break
endif
endfor
endfor

printf('Azul fecho');

%Quantidade de curvas subparabolicas abertas e fechadas

printf('Parte 9');  
con6s=0;
ncon6s=0;
for tt=1:3:size(MMs2,2)
  if (pontocurvass((tt+2)/3,1)>1 && MMs2(1,tt)==MMs2(pontocurvass((tt+2)/3,1),tt) && MMs2(1,tt+1)==MMs2(pontocurvass((tt+2)/3,1),tt+1))
    circulos(con6s+1,1)=tt;
    con6s=con6s+1; %Quantidade de curvas subparabolica vermelha fechada
  else
    ncirculos(ncon6s+1,1)=tt;
    ncon6s=ncon6s+1; %Quantidade de curvas subparabolica vermelha aberta
  endif
endfor
con6sw=0;
ncon6sw=0;
for tt=1:3:size(MMs2w,2)
  if (pontocurvassw((tt+2)/3,1)>1 && MMs2w(1,tt)==MMs2w(pontocurvassw((tt+2)/3,1),tt) && MMs2w(1,tt+1)==MMs2w(pontocurvassw((tt+2)/3,1),tt+1))
    circulosw(con6sw+1,1)=tt;
    con6sw=con6sw+1; %Quantidade de curvas subparabolica azul fechada
  else
    ncirculosw(ncon6sw+1,1)=tt; %Quantidade de curvas subparabolica azul aberta
    ncon6sw=ncon6sw+1;
  endif
endfor

totalvermelhas = con6s+ncon6s; %Total de curvas subparabolica vermelha
totalazuls= con6sw+ncon6sw;  %Total de curvas subparabolica azul
endif

endif

%Plotagem dos pontos para triangulacao
  x = V(:,1);
  y = V(:,2);
  z = V(:,3);
  
if(resposta==1 && size(MM!=0) && size(MMw!=0) )

##  Resultado da implementacao octave
  hold on;
##  trimesh (F,x,y,z);    %Ative caso queira a imagem com a malha
 
%Pontos Umbílicos - Ative caso queira os pontos umbílicos na imagem
 %plot3(Pumbc(:,1), Pumbc(:,2), Pumbc(:,3),'k*','LineWidth',[15]); 

%Curvas Ridges Azul - Ative caso queira na imagem
  for j=1:size(pontocurvas)(1)
    plot3(MM2(1:pontocurvas(j,1),3*j-2), MM2(1:pontocurvas(j,1),3*j-1), MM2(1:pontocurvas(j,1),3*j),'b-','LineWidth',[2]);
  endfor
  
%Curvas Ridges Vermelha - Ative caso queira na imagem

  for j=1:size(pontocurvasw)(1)
    plot3(MM2w(1:pontocurvasw(j,1),3*j-2), MM2w(1:pontocurvasw(j,1),3*j-1), MM2w(1:pontocurvasw(j,1),3*j),'r-','LineWidth',[2]);
endfor

axis square
hold off;

elseif(resposta==2 && size(MMs!=0) && size(MMsw!=0) )

##  Resultado da implementacao octave
  hold on;
##  trimesh (F,x,y,z);    %Ative caso queira a imagem com a malha
 
%Pontos Umbílicos - Ative caso queira os pontos umbílicos na imagem
 plot3(Pumbc(:,1), Pumbc(:,2), Pumbc(:,3),'k*','LineWidth',[15]); 

%Curvas Ridges Azul - Ative caso queira na imagem
  for j=1:size(pontocurvas)(1)
    plot3(MM2(1:pontocurvas(j,1),3*j-2), MM2(1:pontocurvas(j,1),3*j-1), MM2(1:pontocurvas(j,1),3*j),'b-','LineWidth',[2]);
  endfor
  
%Curvas Ridges Vermelha - Ative caso queira na imagem

  for j=1:size(pontocurvasw)(1)
    plot3(MM2w(1:pontocurvasw(j,1),3*j-2), MM2w(1:pontocurvasw(j,1),3*j-1), MM2w(1:pontocurvasw(j,1),3*j),'r-','LineWidth',[2]);
endfor

%Curvas Subparabólicas Vermelha - Ative caso queira os pontos umbílicos na imagem

  for j=1:size(pontocurvass)(1)
    plot3(MMs2(1:pontocurvass(j,1),3*j-2), MMs2(1:pontocurvass(j,1),3*j-1), MMs2(1:pontocurvass(j,1),3*j),'m-','LineWidth',[2]);
  endfor
  
%Curvas Subparabólicas Azul - Ative caso queira os pontos umbílicos na imagem

  for j=1:size(pontocurvassw)(1)
    plot3(MMs2w(1:pontocurvassw(j,1),3*j-2), MMs2w(1:pontocurvassw(j,1),3*j-1), MMs2w(1:pontocurvassw(j,1),3*j),'c-','LineWidth',[2]);
endfor


axis square
hold off;

endif

toc
 
 
