function superpixel_c_prob = calculateProbability(msprime_superpix,ck,T)
% based on formula (2) in paper
% INPUTS:
% - msprime_superpix : ms' value for every superpixel (Nx3)
% - ck : current color palette (Kx3)
% - T : current temperature
%
% OUTPUTS:
% - superpix_c_prob : matrix every row represents a superpixel, every
%                     column the conditional prob for every color in
%                     palette (NxK)

%P(ck) ook megeven, ma hoe berekenen???

N = height(msprime_superpix);
K = heigh(ck);
superpix_c_prob = zeros(N,K);

for i=1:N
    for j=1:K
        distance = norm(msprime_superpix(i,:)-ck(K,:));
        superpix_c_prob(i,j) = Pck*exp(-1*distance/T);
    end
end


end

