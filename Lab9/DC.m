 f = [-66.67; -33.33; 0;0;0;0;0;0;0;0;0;0];  
 A = [-1 -1 0 0 0 0 0 0 0 0 0 0; 1  1 0 0 0 0 0 0 0 0 0 0];  
 b = [-0.5; 1.5];  
 Aeq = [0 0  1 -1 0 0 0 -1 0 0 0 0;
0 0  1 0 -1 0 0 0 -1 0 0 0;
0 0  0 0 1 0 -1 0 0 -1 0 0;
0 0  0 1 0 -1 0 0 0 0 -1 0;
0 0  0 0 0 1 -1 0 0 0 0 -1;
0 0  1 0 0 0 0 0 0 0 0 0;
1 1 0 0 0 0 0 1 1 0 0 0;
0 0 0 0 0 0 0 1 0 1 0 0 ;
0 0 0 0 0 0 0 0 1 0 1 0;
1 0 0 0 0 0 0 0 0 -1 0 -1;
0 1 0 0 0 0 0 0 0 0 -1 -1];
beq=[0 0 0 0 0 0 1.5 -1.5 0 0 0];
 lb = [0; 0;-3.14;-3.14;-3.14;-3.14;-3.14;-1;-1;-1;-1;-1] ; 
 ub= [1;1;3.14;3.14;3.14;3.14;3.14;1;1;1;1;1];
[x,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub)  

