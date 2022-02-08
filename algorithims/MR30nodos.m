clear all
clc
format shortG

%% Verificación algoritmo
% Parámetros base
vref = 1000;
Pbase = 70000;
rbase = vref^2/Pbase;

step = 1;
droop = 2;
potencia = 3;
%% Datos de Lineas y Nodos
%           i   j    r(pu)
    A = [   1	2	0.01919232
            1	3	0.04518192
            2	4	0.0569772
            3	4	0.01319472
            2	5	0.04718112
            2	6	0.05807676
            4	6	0.01189524
            5	7	0.0459816
            6	7	0.02668932
            6	8	0.0119952
            6	9	0.2079168
            6	10	0.5557776
            9	11	0.2079168
            9	10	0.109956
            4	12	0.2558976
            12	13	0.139944
            12	14	0.12305076
            12	15	0.06617352
            12	16	0.0944622
            14	15	0.2209116
            16	17	0.08236704
            15	18	0.1069572
            18	19	0.06387444
            19	20	0.339864
            10	20	0.09356256
            10	17	0.03238704
            10	21	0.03478608
            10	22	0.07267092
            21	22	0.01159536
            15	23	0.09996
            22	24	0.114954
            23	24	0.1319472
            24	25	0.1884246
            25	26	0.25429824
            25	27	0.10925628
            28	27	0.3958416
            27	29	0.21971208
            27	30	0.32007192
            29	30	0.23980404
            8	28	0.06357456
            6	28	0.1689324 ]; 
            
%       nodo tipo       Potgen/gdroop (p.u)            PotCarga(pu)   
B = [   1	droop       1/0.0437491250174997            0
        2	droop       1/0.0874982500349993            0.31
        3	potencia	0                               0.0342857142857143
        4	potencia	0                               0.108571428571429
        5	droop       1/0.0945927027405398            0.271428571428571
        6	step        0                               0
        7	potencia	0                               0.325714285714286
        8	droop       1/0.092103421089473             0.214285714285714
        9	step        0                               0
        10	potencia	0                               0.0828571428571429
        11	droop       1/0.0777762222533327            0
        12	potencia	0                               0.16
        13	droop       1/0.0897417949076916            0
        14	potencia	0                               0.0885714285714286
        15	potencia	0                               0.117142857142857
        16	potencia	0                               0.05
        17	potencia	0                               0.128571428571429
        18	potencia	0                               0.0457142857142857
        19	potencia	0                               0.135714285714286
        20	potencia	0                               0.0314285714285714
        21	potencia	0                               0.25
        22	step        0                               0
        23	potencia	0                               0.0457142857142857
        24	potencia	0                               0.124285714285714
        25	step        0                               0
        26	potencia	0                               0.05
        27	step        0                               0
        28	step        0                               0
        29	potencia	0                               0.0342857142857143
        30	potencia	0                               0.151428571428571  ];

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
 I = Y*v;%calculo de corrientes nodales  j
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