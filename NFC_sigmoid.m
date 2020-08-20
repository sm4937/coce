function [cost] = NFC_sigmoid(rawcosts,NFC)
% Run an inverse sigmoid on components of cost

if ~isnan(NFC)
    for i = 1:length(rawcosts)
       cost(i) = exp(NFC*rawcosts(i))/sum(exp(NFC.*[rawcosts(i) max(rawcosts(1:i))]));
    end
else
    cost = NaN(1,length(rawcosts));
end

end

