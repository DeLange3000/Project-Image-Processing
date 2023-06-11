function [normalisation_factor] = GetNormalisationFactor(p_ck, ms, palette, T)

[a, b] = size(palette);

normalisation_factor = 0;
for i = 1:b
    normalisation_factor = normalisation_factor + p_ck(i)*exp(-norm(ms - palette(:,i)')/T);
end