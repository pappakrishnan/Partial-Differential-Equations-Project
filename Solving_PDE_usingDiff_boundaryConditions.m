
function [] = Solving_PDE_usingDiff_boundaryConditions
%global values ar br cr Dr xref yref Xminn Xmaxx Yminn Ymaxx

% matrix defining triangles
load EP.txt;  
% Node matrix
load nodes.txt;
% removes the 1st column of matrix "nodes"
nodes(:,1) = [];

%Boundaries 		
%               BD(:,1) = 0 interior 
%				BD(:,1) = 1 Dirichlet 
%				BD(:,1) = 2 Neumann BC
%				BD(:,2) = value of Dirichlet BC

%loading the boundary matrix
load BD.txt;

BD(:,1) = []; % removes 1st column of boundary matrix

A        = zeros(size(nodes, 1), size(nodes, 1));
F        = zeros(size(nodes, 1), 1); 

     for j = 1:size(EP,1)
        A(EP(j,:), EP(j,:))  = A(EP(j,:),EP(j,:)) + Ks(nodes(EP(j,:), :));  
        F(EP(j,:), 1) = F(EP(j,:), 1)  - Qs(nodes(EP(j,:),:),@f0);
     end
  
      for m = 1:size(nodes, 1)
          if BD(m, 1) == 1
              % Dirichlet BC 
              for i = 1:size(nodes, 1)
                  F(i)    = F(i) - A(i,m)*BD(m,2); 
                  A(i, m) = 0.0; 
                  A(m, i) = 0.0; 
              end % i loop 
             A(m, m) = 1.0;
             F(m)    = BD(m, 2); 
          end 
      end 
%gauss-Siedel iterative method:

C=ones(length(A),1);
iteration1=0;
for i=1:10000
    count=0;
    iteration1=iteration1+1;
    for j=1:length(A)
        if(i>1)
            error=(C(j,i)-C(j,i-1));
            if(error>0.2)
                     B1=0;
                     B2=0;
                        if (j>1)
                            for h=1:j-1
                                B1=B1+A(j,h)* C(h,i+1);
                            end
                        end 
                            for h=j+1:length(A)
                                B2=B2+A(j,h)* C(h,i);
                            end               
                      C(j,i+1)= (1/A(j,j))* (F(j)-B1-B2);                
            else
                count=count+1;
                C(j,i+1)=C(j,i);
            end             
        else
            B1=0;
            B2=0;
                if (j>1)
                    for h=1:j-1
                        B1=B1+A(j,h)* C(h,i+1);
                    end
                end 
                    for h=j+1:length(A)
                        B2=B2+A(j,h)* C(h,i);
                    end               
                 C(j,i+1)= (1/A(j,j))* (F(j)-B1-B2);   
         end          
    end
                    if (count == length(A))                        
                        break;
                    else
                        continue;
                    end
                    
end


% To display the output 
disp('Full array A: ');
A
disp('Final right hand side*12:') 
F 
disp('Node Values');
C
iteration1
%Plotting the output
figure(1);
trisurf(EP, nodes(:,1)', nodes(:,2)',C(:,iteration1+1));
zlabel('U(x,y)'); xlabel('X'); ylabel('Y');title('Distribution of U for 10 nodes');
figure(2);
trimesh(EP, nodes(:,1)', nodes(:,2)');
xlabel('X'); ylabel('Y');title('Illustration of nodes (10)');

  function r = Ks(coordinates)  
  x    = coordinates(:,1);
  y    = coordinates(:,2);
  a1   = x(2)*y(3) - x(3)*y(2); b1 = y(2) - y(3); c1 = x(3) - x(2); 
  a2   = x(3)*y(1) - x(1)*y(3); b2 = y(3) - y(1); c2 = x(1) - x(3);
  a3   = x(1)*y(2) - x(2)*y(1); b3 = y(1) - y(2); c3 = x(2) - x(1);
  
  r    = [(b1*b1 + c1*c1) (b1*b2 + c1*c2) (b1*b3 + c1*c3); 
          (b2*b1 + c2*c1) (b2*b2 + c2*c2) (b2*b3 + c2*c3);
          (b3*b1 + c3*c1) (b3*b2 + c3*c2) (b3*b3 + c3*c3)];

 Delta = @(x, y) abs( (x(2)*y(3) - x(3)*y(2) + x(3)*y(1)- x(1)*y(3) + x(1)*y(2) - x(2)*y(1))/2.0 );
 uv = [1/6 1/12 1/12;
       1/12 1/6 1/12;
       1/12 1/12 1/6];
 uv = uv./Delta(x,y);
 r  = r./(4.0*Delta(x, y));%+uv;
 
 
function r = single_pt_rule(p, q, r, a,  b, c, D, f)
j         = abs( (q(1) - p(1))*(r(2) - p(2)) - (q(2) - p(2))*(r(1) - p(1)) ); 
one       = (p + q + r)./3.0;
% beta is a basis function 
beta     = @(x, a, b, c, D)((a + b*x(1) + c*x(2))/(2*D)); 
r        = (j/2.0)* ( (f(one)*beta(one, a, b, c, D)) );  

 function r = f0(x, y)
 r = 0;
 function  r = Qs(coordinates, f)

  x    = coordinates(:,1);
  y    = coordinates(:,2);
  a1   = x(2)*y(3) - x(3)*y(2); b1 = y(2) - y(3); c1 = x(3) - x(2); 
  a2   = x(3)*y(1) - x(1)*y(3); b2 = y(3) - y(1); c2 = x(1) - x(3);
  a3   = x(1)*y(2) - x(2)*y(1); b3 = y(1) - y(2); c3 = x(2) - x(1);
  
  Delta= @(x, y) (x(2)*y(3) - x(3)*y(2) + x(3)*y(1) - x(1)*y(3) + x(1)*y(2) - x(2)*y(1))/2.0;
  D    = Delta(x, y); 
  
  xt   = x; yt = y;
  r1   = single_pt_rule([x(1), y(1)], [x(2), y(2)], [x(3), y(3)], a1,  b1, c1, D, f); 
  
  x(1) = xt(2); x(2) = xt(3); x(3) = xt(1);
  y(1) = yt(2); y(2) = yt(3); y(53) = yt(1);
  r2   = single_pt_rule([x(1), y(1)], [x(2), y(2)], [x(3), y(3)], a2,  b2, c2, D, f); 

  xt   = x; yt = y;
  x(1) = xt(2); x(2) = xt(3); x(3) = xt(1);
  y(1) = yt(2); y(2) = yt(3); y(3) = yt(1);
  r3   = single_pt_rule([x(1), y(1)], [x(2), y(2)], [x(3), y(3)], a3,  b3, c3, D, f); 
  r    = [r1; r2; r3];