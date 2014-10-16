function A = make_matrix_A(Xs, Z, order)
% A = make_matrix_A(Xs, Z, order)
% Construct the matrix A (the linear operator L) involved in the model Conductivity_R2 and ElectricFish
% Inputs:
% Xs: coordinates of sources/receivers
% Z: reference center
% order: highest order of CGPT

N = size(Xs,2);
A = zeros(N, 2*order);

for n=1:N
	toto = Xs(:,n)-Z(:);
	[T,R] = cart2pol(toto(1), toto(2));
	
	for m=1:order
		A(n, 2*m-1:2*m)=[cos(m*T),sin(m*T)] / (2*pi*m*R^m);
	end
end

end
