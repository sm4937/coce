% Write Theta for all the possible parameters across models with empirical
% prior P(Theta)
% 
% take subject s, with posterior model weightings rho^s_j for model j
% [from Payam]
% 
% say model j involves parameters Theta_j but not Theta_{j'}
% 
% and subject s has posterior distribution P(Theta^s_j|D^s) for model j
% where D^s is the data for subject s [from Payam]
% 
% then, we might be interested in
% 
% P(Theta | D^s) =
%   \sum_j P(Theta,j|D^s)
% = \sum_j [\rho^s_j P(Theta^s_j|D^s)*P(Theta_{j'})]
% 
% using the population prior P(Theta_{j'}) for those parameters not in
% model j and the posterior
% 
% so then, if you have subjects s_1...s_k in your low tertile group (say),
% you can work out the approximate posterior distribution
% 
% P(Theta ; low_tertile) =
%    \prod_{i=1}^k \sum_j [\rho^{s_i}_j P(Theta^{s_i}_j|D^{s_i})*P(Theta_{j'})]
% 
% in the low_tertile [in practice you should calculate it carefully using
% logs - in the three dimensions of maintenance, lure and false-alarm cast]
% 
% and you could then compare these posteriors for the three groups.
% 
% In practice, if a subject isn't well fit by a model (rho^s_j is low);
% then that subject won't pull parameters that are only in that model away
% from their population prior very much - which is just the property you
% want. 
