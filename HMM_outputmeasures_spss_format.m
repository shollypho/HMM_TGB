cnt = 0;
K = 6;
clear Spss
for ss = 1:2:size(FO,1)
    cnt = cnt + 1;
    Spss{cnt,1} = all_data{data_incl(ss),1};
    Spss{cnt,3} = all_data{data_incl(ss),2};
    if strcmp(Spss{cnt,3}, 'Con')
        Spss{cnt,2} = 1;
    else
        Spss{cnt,2} = 2;
    end
    
    
    if all_data{data_incl(ss),7} ==1
        for k = 1:K
            Spss{cnt,3+k} = FO(ss,k);
            Spss{cnt,3+K+k} = FO(ss+1,k);
            
            Spss{cnt,3+(2*K)+k} = LTmerged(ss,k);
            Spss{cnt,3+(3*K)+k} = LTmerged(ss+1,k);
            
            Spss{cnt,3+(4*K)+k} = ITmerged(ss,k);
            Spss{cnt,3+(5*K)+k} = ITmerged(ss+1,k);
        end
    else
        for k = 1:K
            Spss{cnt,3+K+k} = FO(ss,k);
            Spss{cnt,3+k} = FO(ss+1,k);
            
            Spss{cnt,3+(3*K)+k} = LTmerged(ss,k);
            Spss{cnt,3+(2*K)+k} = LTmerged(ss+1,k);
                        
            Spss{cnt,3+(5*K)+k} = ITmerged(ss,k);
            Spss{cnt,3+(4*K)+k} = ITmerged(ss+1,k);
        end
    end
end