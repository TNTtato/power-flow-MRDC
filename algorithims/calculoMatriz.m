%Datos lineas
%           i   j    r(pu)
    A = [   1	2	0.0192
            1	3	0.0452
            2	4	0.057
            3	4	0.0132
            2	5	0.0472
            2	6	0.0581
            4	6	0.0119
            5	7	0.046
            6	7	0.0267
            6	8	0.012
            6	9	0.208
            6	10	0.556
            9	11	0.208
            9	10	0.11
            4	12	0.256
            12	13	0.14
            12	14	0.1231
            12	15	0.0662
            12	16	0.0945
            14	15	0.221
            16	17	0.0824
            15	18	0.107
            18	19	0.0639
            19	20	0.34
            10	20	0.0936
            10	17	0.0324
            10	21	0.0348
            10	22	0.0727
            21	22	0.0116
            15	23	0.1
            22	24	0.115
            23	24	0.132
            24	25	0.1885
            25	26	0.2544
            25	27	0.1093
            28	27	0.396
            27	29	0.2198
            27	30	0.3202
            29	30	0.2399
            8	28	0.0636
            6	28	0.169   ]; 
   
NumN = max(max(A(:,1:2))); %Encuentra el numero de nodos
NumL = length(A(:,1)); %Encuantra el numero de lineas
Y = zeros(NumN,NumN); %Matriz de ceros de NxN
%C?lculo de la matriz YBUS
for k = 1:NumL
    n1 = A(k,1);
    n2 = A(k,2);
    ykm = 1/A(k,3);
    Y(n1,n1) = Y(n1,n1) + ykm;
    Y(n1,n2) = Y(n1,n2) - ykm;
    Y(n2,n1) = Y(n2,n1) - ykm;
    Y(n2,n2) = Y(n2,n2) + ykm;
end
