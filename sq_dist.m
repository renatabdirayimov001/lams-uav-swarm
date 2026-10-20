function D2 = sq_dist(p)
%SQ_DIST  Pairwise squared Euclidean distances between the rows of p.
s = sum(p.^2, 2);
D2 = max(s + s' - 2*(p*p'), 0);
end
