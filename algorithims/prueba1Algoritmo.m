clear all
clc
format shortG
%% Verificación algoritmo
% Parámetros base
vref = 48;
Pbase = 500;

rlineas = 0.05; 
rbase = vref^2/Pbase;
rpu = rlineas/rbase;

step = 1;
droop = 2;
potencia = 3;
%% Datos de Lineas y Nodos
%    Lineas   i j r
        A = [ 4 1 rpu
              1 2 rpu
              2 3 rpu
              3 4 rpu];
%   nodo tipo     Potgen/gdroop  PotCarga   
B = [1   droop      rbase/0.2           0
     2   potencia    0            466.25/Pbase
     3   droop      rbase/0.5           0
     4   potencia    0             697.5/Pbase  ];

NumN = max(max(A(:,1:2))); %Número de Nodos
NumL = length(A(:,1)); %Número de Lineas

%% Separación de variables

Nstep = B(find(B(:,2)==step),1); %Nodos step
Vstep = B(find(B(:,2)==step),2); %Variable step
Ndroop = B(find(B(:,2)==droop),1);%Nodos droop
Vdroop = B(find(B(:,2)==droop),2); %Variable droop
NP = B(find(B(:,2)==potencia),1);%Nodos de Pot const
Vpotencia = B(find(B(:,2)==potencia),2); %Variable Pot const

%Matriz clasificación de variables

Var = [[Nstep; Ndroop; NP] [Vstep;Vdroop;Vpotencia]];
NumVar = length(Var);

%% Matriz Ybus y declaración J y H 
%Cálculo de Ybus
Y = zeros(NumN,NumN);
for k = 1:NumL
    n1 = A(k,1);
    n2 = A(k,2);
    ykm = 1/A(k,3);
    Y(n1,n1) = Y(n1,n1) + ykm;
    Y(n1,n2) = Y(n1,n2) - ykm;
    Y(n2,n1) = Y(n2,n1) - ykm;
    Y(n2,n2) = Y(n2,n2) + ykm;
end

v = ones(NumVar,1); %Condiciones iniciales de tensión
H=zeros(NumVar,1); %Inicialización desajustes
J = zeros(NumVar,NumVar);%Declaración matriz jacobiana
deltaV=zeros(NumVar,1);%Inicialización de variaciones

%% Método Newton Raphson
tic%contador de tiempo
%inicio de iteraciones
for j=1:10
    
%Hacer un barrido por todas las variables    
for k=1:NumVar
 %Si es step, aplicar ecuaciones step 
  if Var(k,2)==step
      pstep = 0;
      pstepl = v(Var(k,1))*(Y(Var(k,1),:)*v);
      H(Var(k,1)) = pstep-pstepl;
      for i=1:NumVar
       if Var(k,1)==i
       J(Var(k,1),i) = -v(Var(k,1))*Y(Var(k,1),Var(k,1)) - Y(Var(k,1),:)*v;
       else
       J(Var(k,1),i) = -v(Var(k,1))*Y(Var(k,1),i);
       end
      end
  else
 %Si es droop, aplicar ecuaciones droop
  if Var(k,2)==droop
      %inyecciones de potencia al nodo droop
      pgdroop = B(Var(k,1),4);
      %expresión incluyendo ecuas droop
      pdroop = B(Var(k,1),3)*(1*v(Var(k,1))-v(Var(k,1))^2);
      %potencias por las lineas
      pdroopl = v(Var(k,1))*(Y(Var(k,1),:)*v);
      %ecuación de desajuste
      H(Var(k,1)) = pdroop-pgdroop-pdroopl;
      %Cálculo de las derivadas parciales para droop k específico
      for i=1:NumVar
        %si k=i
        if Var(k,1)==i
          Jdroop = B(Var(k,1),3)*(1 - 2*v(Var(k,1)));
          Jdroopl = v(Var(k,1))*Y(Var(k,1),Var(k,1)) - Y(Var(k,1),:)*v;
          J(Var(k,1),i) = Jdroop - Jdroopl;
        else
          %k != i
          J(Var(k,1),i) = -v(Var(k,1))*Y(Var(k,1),i);
        end
      end
  else
 %Si es potencia, aplicar ecuaciones potencia   
  if Var(k,2)==potencia
      ppot = B(Var(k,1),3)-B(Var(k,1),4); %PG-Pcarga
      ppotl = v(Var(k,1))*(Y(Var(k,1),:)*v);%Plineas
      H(Var(k,1)) = ppot-ppotl;%ecuación h(k) de potencia
      for i=1:NumVar
        %si k=i
        if Var(k,1)==i
      J(Var(k,1),i) = -v(Var(k,1))*Y(Var(k,1),Var(k,1)) - Y(Var(k,1),:)*v;
        else
        %k != i
      J(Var(k,1),i) = -v(Var(k,1))*Y(Var(k,1),i);
        end
     end
  end
  end
  end
    
end
deltaV=inv(J)*H;%vector de variacion
 v = v - deltaV;%Actalización de variables
 I = Y*v;%calculo de corrientes nodales
 P = diag(v)*I;%potencias nodales
 %iter(j) = j;
 errAbs = abs(max(H));
end
%fin de iteraciones
t=toc%fin contador de tiempo

%% Muestra de Resultados
disp('         nodo      v(pu)         P(pu)        v(V)        P(W)')
solucion = [B(:,1)    v        P     vref*v  Pbase*P];
disp(solucion)

disp('Maximo error abs')
disp(errAbs)