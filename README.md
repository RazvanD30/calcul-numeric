
A = diag(3 * ones(6,1),0) +  diag(-ones(5,1),1) +  diag(-ones(5,1),-1);
b = [2; 1; 1;1;1;2];
res = gauss(A,b,1e-5,b)


function x = lab11(A,b,er,xOld)
  n = length(b);
  xNew = zeros(n,1);
  for i = 1:n
   xNew(i) = (1 / A(i,i)) * (b(i) - A(i,1:i-1) * xNew(1:i-1) - A(i,i+1:n) * xOld(i+1:n)); 
  end
  # stop if norm(xNew - xOld) < er
  if(norm(xNew - xOld) < er)
    x = xNew;
  else
    x = lab11(A,b,er,xNew);
  end
end


-------------------------------------------------------------------------------


A = diag(3 * ones(6,1),0) +  diag(-ones(5,1),1) +  diag(-ones(5,1),-1)
b = [2; 1; 1;1;1;2];
M = diag(diag(A));
N = M - A;
T = M \ N;
c = M \ b;

function x = lab11b(T,c,er,xOld)
  xNew = T * xOld + c;
  normT = norm(T);
  if(norm(xOld - xNew) < er * ((1 - normT )/ normT))
    x = xNew;
  else
    x = lab11(T,c,er,xNew);
 end

