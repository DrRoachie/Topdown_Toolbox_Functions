function [eID_a, eID_b] = getAreaLabel(coh)
%gets the label of area following coherence calculation

coh = coh.labelcmb;

column1 = coh(1, 1);
str1 = char(column1);
eID_a = str1(1:2);

column2 = coh(1, 2);
str2 = char(column2);
eID_b = str2(1:2);

end