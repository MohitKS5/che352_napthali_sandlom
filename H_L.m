function H = H_L(x, T_stage)
    %methanol
    h(1) = x(1)*integral(105800-362.23.*T+0.9379.*T.^2, 330, T_stage);
    %acetone
    h(2) = x(2)*integral(135600-177.*T+0.2837.*T.^2+0.000689.*T.^3+0.000689.*T.^4, 330, T_Stage);
    %methyl acetate
    h(3) = x(3)*integral(61260+270.9.*T, 330, T_stage);
    %benzene
    h(4) = x(4)*integral(162940-344.94.*T+0.85562.*T.^2, 330, T_stage);
    %chloroform
    h(5) = x(5)*integral(124850-166.34.*T+0.43209.*T.^2, 330, T_stage);
    
    H = sum(h(1:5));
    H = H + dep_Liq(x, T_stage)*10^3;
end