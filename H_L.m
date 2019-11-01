function H = H_L(x, T_stage)
    %methanol
    H(1) = x(1)*integral(@(T) 105800-362.23.*T+0.9379.*T.^2, 330, T_stage);
    %acetone
    H(2) = x(2)*integral(@(T) 135600-177.*T+0.2837.*T.^2+0.000689.*T.^3+0.000689.*T.^4, 330, T_stage);
    %methyl acetate
    H(3) = x(3)*integral(@(T) 61260+270.9.*T, 330, T_stage);
    %benzene
    H(4) = x(4)*integral(@(T) 162940-344.94.*T+0.85562.*T.^2, 330, T_stage);
    %chloroform
    H(5) = x(5)*integral(@(T) 124850-166.34.*T+0.43209.*T.^2, 330, T_stage);
    H = sum(H);
    H = H + dep_L(x, T_stage)*10^3;
end