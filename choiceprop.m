function [proportionH] = choiceprop(choicehigh,choicelow)

choicehigh = nnz(neutraltrials.choice == 'H');
choicelow = nnz(neutraltrials.choice == 'L');

proportionH = choicehigh / (choicehigh + choicelow);

end

