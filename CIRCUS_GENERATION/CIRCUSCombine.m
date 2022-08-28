function [py_unique,pz_unique,n_unique] = CIRCUSCombine(py1,pz1,n1,py2,pz2,n2)
% Combine two CIRCUS patterns, keeping all entries in the first pattern
% and only those in the second pattern which are not already in the first

py_unique = zeros(size(py1,1) + size(py2,1),1);
pz_unique = zeros(size(pz1,1) + size(pz2,1),1);
n_unique = n1;

for n=1:n1
  py_unique(n) = py1(n);
  pz_unique(n) = pz1(n);
end
for n=1:n2
     keepit = 1;
     % Iterate through phase table and keep only new values - not
     % efficient but should work
     for m = 1:n1
       if (py1(m)==py2(n)) && (pz1(m)==pz2(n))
          keepit = 0;
       end
     end
     if keepit
       n_unique = n_unique + 1;
       py_unique(n_unique) = py2(n);
       pz_unique(n_unique) = pz2(n);
     end
end   
py_unique = py_unique(1:n_unique);
pz_unique = pz_unique(1:n_unique);

end