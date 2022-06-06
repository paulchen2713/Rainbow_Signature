% 
% Multivariate cryptography ---Rainbow Signature:
%         multi-layer unbalanced Oil-Vinegar Signature
% 
clear;
clc;
global m n alpha alpha1;
%
starting_time = cputime;
m  = 8;
nn = 37;
uu = 5;
VV = [10 20 24 27 37];
%
OO = zeros(1, uu - 1);
for iu = 1 : uu - 1
    OO(iu) = VV(iu + 1) - VV(iu);
end
%
oo = nn - VV(1); % 27
q  = 2^m;
% n = v + o;
one  = q - 1;
zero = q;
%
% generation of qq(x)
%
qq = primitive_polynomial(m);
qq_size = size(qq, 1); 
qq_size = floor(qq_size * rand(1) + 1);
qq = qq(qq_size, :);
%
alpha = zeros(q, m);
alpha(1, 2) = 1; % alpha^1
for i = 2 : one % i = 2 : 2^m - 1
    if alpha(i-1, m) == 1
        alpha(i, 2:m) = alpha(i-1, 1 : m-1);
        alpha(i, :) = bitxor(alpha(i, :), qq(1:m));
    else
        alpha(i, 2:m) = alpha(i-1, 1 : m-1);
    end
end
%
alphaA = char(alpha + 48);
alphaA = bin2dec(alphaA);
alpha1 = zeros(q, 1); % zeros(2^m, 1);
for i = 1 : one % i = 1 : 2^m - 1
    alpha1(alphaA(i) + 1) = i;
end
alpha1(1) = zero;
%
% random generation of L1, C1 and L2, C2, generation of fb1, ..., fbo
%
indexS = ones(1, uu - 1);
while indexS ~= 0
    indexS = ones(1, uu - 1);
    %
    % random generation of L1, C1
    %
    eyee = zero * ones(oo, oo);
    for in = 1 : oo
        eyee(in, in) = one;
    end
    %
    indexL1 = 0; % set L1 flag to 0
    while indexL1 == 0
        L1  = floor(zero * rand(oo, oo) + 1);
        LL1 = [L1, eyee]; %
        LL1 = reduced_row_echelon_form_power(LL1);
        L1I = LL1(:, oo+1 : 2*oo);
        L1_L1I = matrix_multiplication_power(L1, L1I);
        if any(any(L1_L1I - eyee)) == 0
            indexL1 = 1; % set L1 flag to 1
        end
    end
    C1 = floor(zero * rand(oo, 1) + 1);
    %
    % random generation of L2, C2
    %
    eyee = zero * ones(nn, nn);
    for in = 1 : nn
        eyee(in, in) = one;
    end
    indexL2 = 0; % L2 flag set to 0
    while indexL2 == 0
        L2  = floor(zero * rand(nn, nn) + 1);
        LL2 = [L2, eyee]; %
        LL2 = reduced_row_echelon_form_power(LL2);
        L2I = LL2(:, nn+1 : 2*nn);
        L2_L2I = matrix_multiplication_power(L2, L2I);
        if any(any(L2_L2I - eyee)) == 0
            indexL2 = 1; % set L2 flag to 1
        end
    end
    C2 = floor(zero * rand(nn, 1) + 1);
    % 
    % random generation of f1, f2, ..., fo
    %
    F  = zero * ones(nn, nn, oo); % 37 x 37 x 27
    F1 = zero * ones( 1, nn, oo); %  1 x 37 x 27
    for iu = 1 : uu-1
        for io = 1 : OO(iu)
            AA = floor(zero * rand(VV(iu), VV(iu + 1)) + 1);
            F(1:VV(iu), 1:VV(iu + 1), VV(iu) + io - VV(1)) = AA;
            %
            BB = floor(zero * rand(1, VV(iu + 1)) + 1);
            F1(:, 1:VV(iu + 1), VV(iu) + io - VV(1)) = BB; % F1(1, :, io)?
        end
    end
    Fc = floor(zero * rand(1, 1, oo) + 1);
    %
    % generation of fb1, fb2, ..., fbo
    %
    FBB  = zero * ones(nn, nn, oo); % 37 x 37 x 27
    FBB1 = zero * ones( 1, nn, oo); %  1 x 37 x 27
    FBBc = zero * ones( 1,  1, oo); %  1 x  1 x 27
    for io = 1 : oo
        FBB(:, :, io) = matrix_multiplication_power(matrix_multiplication_power(L2', F(:, :, io)), L2);
        FBB11 = matrix_multiplication_power(matrix_multiplication_power(C2', (F(:, :, io))'), L2);
        FBB12 = matrix_multiplication_power(matrix_multiplication_power(C2', F(:, :, io)), L2);
        FBB13 = matrix_multiplication_power(F1(:, :, io), L2); % F1(1, :, io)?
        %
        FBB1(:, :, io) = matrix_addition_power(matrix_addition_power(FBB11, FBB12), FBB13); % FBB1(1, :, io)?
        FBBc1 = matrix_multiplication_power(matrix_multiplication_power(C2', F(:, :, io)), C2);
        FBBc2 = matrix_multiplication_power(F1(:, :, io), C2); % F1(1, :, io)?
        FBBc(1, 1, io) = matrix_addition_power(matrix_addition_power(FBBc1, FBBc2), Fc(1, 1, io));
    end
    %
    FB  = zero * ones(nn, nn, oo); % 37 x 37 x 27
    FB1 = zero * ones( 1, nn, oo); %  1 x 37 x 27
    FBc = zero * ones( 1,  1, oo); %  1 x  1 x 27
    for io1 = 1 : oo
        for io2 = 1 : oo
            LFL = scalar_matrix_multiplication_power(L1(io1, io2), FBB(:, :, io2));
            LFL1 = scalar_matrix_multiplication_power(L1(io1, io2), FBB1(:, :, io2)); % FBB1(1, :, io2))?
            LFLc = multiplication_power(L1(io1, io2), FBBc(1, 1, io2)); % 
            FB(:, :, io1) = matrix_addition_power(FB(:, :, io1), LFL);
            FB1(1, :, io1) = matrix_addition_power(FB1(1, :, io1), LFL1); % FB1(1, :, io1)?
            FBc(1, io1) = addition_power(FBc(1, io1), LFLc); % FBc(1, 1, io1)?
        end
        FBc(1, io1) = addition_power(FBc(1, io1), C1(io1)); % FBc(1, 1, io1)?
    end
    %
    % public key: m, n, fb1, fb2, ..., fboo
    % private key: L1, L2, f1, f2, ..., foo
    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %
    % generation of signature
    %
    % document: Y
    %
    message = 'Has To Be Exactly 27 Words!';
    Y = (double(message))';
    %
    % compute YBB = L1^(-1)(y)
    %
    YBB = matrix_addition_power(Y, C1);
    YBB = matrix_multiplication_power(L1I, YBB);
    %
    % randomly initial Vinegar: V
    %
    V = floor(zero * rand(VV(1), 1) + 1);
    % V = floor(zero * rand(v, 1) + 1);
    %
    for iu = 1 : uu-1
        v = VV(iu);
        o = OO(iu);
        n = v + o;
        %
        YB = YBB(VV(iu) + 1 - VV(1) : VV(iu + 1) - VV(1));
        S  = zero * ones(o, o);
        Sr = zero * ones(o, 1);
        for io = 1 : o
            A = F(1 : VV(iu), 1 : VV(iu), VV(iu) + io - VV(1));
            B = F(1 : VV(iu), VV(iu)+1 : VV(iu + 1), VV(iu) + io - VV(1));
            c = F1(1, 1 : VV(iu), VV(iu) + io - VV(1));
            d = F1(1, VV(iu)+1 : VV(iu + 1),  VV(iu) + io - VV(1));
            S(io, :) = matrix_addition_power(matrix_multiplication_power(V', B), d);
            Sr1 = matrix_multiplication_power(matrix_multiplication_power(V', A), V);
            Sr2 = matrix_multiplication_power(c, V);
            Sr(io) = addition_power(addition_power(addition_power(Sr1, Sr2), YB(io)), Fc(1, 1, VV(iu) + io - VV(1)));
        end
        %
        % SI = S^(-1)
        %
        eyee = zero * ones(o, o);
        for in = 1 : o
            eyee(in, in) = one;
        end
        SS = [S, eyee];
        SS = reduced_row_echelon_form_power(SS);
        SI = SS(:, o+1 : 2*o);
        S_SI = matrix_multiplication_power(S, SI);
        if any(any(S_SI - eyee)) == 0
            indexS(iu) = 0; % set S flag to 0
        else
            fprintf('Warning for matrix inversion!!\n');
            break;
        end
        X = matrix_multiplication_power(SI, Sr);
        XB = [V; X];
        V = XB;
    end
end
%
%
XBc = matrix_addition_power(XB, C2);
z = matrix_multiplication_power(L2I, XBc);
%
% send the message:
%   (signature-message pair) = (z, Y)
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% signature verification
%
check = zeros(oo, 1);
for io = 1 : oo
    check1 = matrix_multiplication_power(matrix_multiplication_power(z', FB(:, :, io)), z);
    check2 = matrix_multiplication_power(FB1(:, :, io), z);
    check(io) = matrix_addition_power(matrix_addition_power(check1, check2), FBc(io)); %
end
%
message_R = char(check');
%
%
if any(message - message_R) == 0
    fprintf('Valid signature. \n\n');
else
    fprintf('Invalid signature. \n\n');
end
%
ending_time = cputime;
computation_time = ending_time - starting_time;
%
fprintf('the original message is: %s\n', message);
fprintf('the recovery message is: %s\n', message_R);
fprintf('the computation time is: %f sec\n', computation_time);
%
% Valid signature. 
% 
% the original message is: Has To Be Exactly 27 Words!
% the recovery message is: Has To Be Exactly 27 Words!
% the computation time is: 81.796875 sec
%

