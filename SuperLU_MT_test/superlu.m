% ����ļ�ͨ����ȡsuperLU�ֽ�����ݣ�Ȼ����֤LU-A
% ѹ���б��
oldcol = [   1   5  10  14  17  21  25  30  33  36  40  44  49  53  56  60  65  69  74  78  83  88  93  97 101 106];
% ����Ԫ���к�
row = [   1   2   6  24   1   2   7  21  25   3   5  22  23   4   5   6   3   4   5   7   1   4   6   7   2   5   6   7  22   ...
    8  13  23   9  11  19  10  11  13  20   9  10  11  18  12  13  20  22  23   8  10  12  13  14  17  24  15  16  ...
    17  19  15  16  18  21  25  14  15  17  25  11  16  18  19  20   9  15  18  19  10  12  18  20  21   2  16  ...
    20  21  22   3   7  12  21  22   3   8  12  23   1  14  24  25   2  16  17  24  25];
% ����Ԫ��
nzval = [
    4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 4.00000000E+00 ...
-1.00000000E+00 -1.00000000E+00 -1.00000000E+00 4.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 -1.00000000E+00 -1.00000000E+00 4.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 -1.00000000E+00 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 -1.00000000E+00 4.00000000E+00 -1.00000000E+00 4.00000000E+00 ...
-1.00000000E+00 -1.00000000E+00 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 4.00000000E+00 -1.00000000E+00 4.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 4.00000000E+00 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 ...
 4.00000000E+00 -1.00000000E+00 -1.00000000E+00 4.00000000E+00 -1.00000000E+00 ...
-1.00000000E+00 -1.00000000E+00 -1.00000000E+00 -1.00000000E+00 4.00000000E+00];
col = zeros(size(row));
% ����δѹ������
for i=1:length(oldcol)-1
   for j=oldcol(i):oldcol(i+1)-1
      col(j) = i; 
   end
end
% ����ϡ�����
s = sparse(row,col,nzval);
% ʹ��MATLAB����LU�ֽ�
[L,U,P,Q] = lu(s);
Lfull = full(L);
Ufull = full(U);
Pfull = full(P);
Qfull = full(Q);
sfull = full(s);
F = Lfull + Ufull - eye(25,25);
% superLU�����L
superLUL = [
0 	 0 	 1
0 	 1 	 -0.25
0 	 6 	 -0.25
1 	 1 	 1
1 	 6 	 -0.0666667
1 	 17 	 -0.266667
1 	 20 	 -0.266667
2 	 2 	 1
2 	 3 	 -0.25
2 	 24 	 -0.25
3 	 3 	 1
3 	 4 	 -0.266667
3 	 5 	 -0.266667
3 	 24 	 -0.0666667
4 	 4 	 1
4 	 5 	 -0.0714286
4 	 24 	 -0.285714
4 	 15 	 -0.267857
4 	 19 	 -0.267857
5 	 5 	 1
5 	 6 	 -0.269231
5 	 24 	 -0.0384615
5 	 15 	 -0.0192308
5 	 19 	 -0.288462
6 	 6 	 1
6 	 24 	 -0.0111029
6 	 17 	 -0.019245
6 	 20 	 -0.30792
6 	 15 	 -0.00555144
6 	 19 	 -0.0832717
7 	 7 	 1
7 	 8 	 -0.25
7 	 13 	 -0.25
8 	 8 	 1
8 	 13 	 -0.0666667
8 	 14 	 -0.266667
8 	 16 	 -0.266667
9 	 9 	 1
9 	 10 	 -0.25
9 	 18 	 -0.25
10 	 10 	 1
10 	 11 	 -0.266667
10 	 12 	 -0.266667
10 	 18 	 -0.0666667
11 	 11 	 1
11 	 12 	 -0.0714286
11 	 18 	 -0.285714
11 	 21 	 -0.267857
11 	 23 	 -0.267857
12 	 12 	 1
12 	 13 	 -0.269231
12 	 18 	 -0.0384615
12 	 21 	 -0.0192308
12 	 23 	 -0.288462
13 	 13 	 1
13 	 14 	 -0.019245
13 	 18 	 -0.0111029
13 	 16 	 -0.30792
13 	 21 	 -0.00555144
13 	 23 	 -0.0832717
14 	 14 	 1
14 	 15 	 -0.267949
14 	 18 	 -0.000198334
14 	 16 	 -0.0769536
14 	 21 	 -9.9167e-05
14 	 23 	 -0.0014875
14 	 24 	 -0.267949
15 	 15 	 1
15 	 16 	 -0.311014
15 	 17 	 -0.00010688
15 	 18 	 -5.7277e-05
15 	 19 	 -0.0837675
15 	 20 	 -0.00171009
15 	 21 	 -2.86385e-05
15 	 22 	 -0.288791
15 	 23 	 -0.000429578
15 	 24 	 -0.160748
16 	 16 	 1
16 	 17 	 -3.77662e-05
16 	 18 	 -0.00392466
16 	 19 	 -0.0295993
16 	 20 	 -0.000604259
16 	 21 	 -0.00196233
16 	 22 	 -0.102044
16 	 23 	 -0.357537
16 	 24 	 -0.0820489
17 	 17 	 1
17 	 18 	 -0.267949
17 	 19 	 -0.00149672
17 	 20 	 -0.0769538
17 	 21 	 -0.267949
17 	 22 	 -3.17858e-05
17 	 23 	 -1.10699e-05
17 	 24 	 -0.000216805
18 	 18 	 1
18 	 19 	 -0.00059195
18 	 20 	 -0.0243961
18 	 21 	 -0.176453
18 	 22 	 -0.000415161
18 	 23 	 -0.106007
18 	 24 	 -0.000452813
19 	 19 	 1
19 	 20 	 -0.323198
19 	 21 	 -0.000596537
19 	 22 	 -0.324125
19 	 23 	 -0.0096609
19 	 24 	 -0.113923
20 	 20 	 1
20 	 21 	 -0.360287
20 	 22 	 -0.117275
20 	 23 	 -0.00638948
20 	 24 	 -0.0453012
21 	 21 	 1
21 	 22 	 -0.38014
21 	 23 	 -0.120556
21 	 24 	 -0.0170751
22 	 22 	 1
22 	 23 	 -0.441726
22 	 24 	 -0.121281
23 	 23 	 1
23 	 24 	 -0.107396
24 	 24 	 1];
sluL = zeros(25,25);
for i=1:length(superLUL)
   sluL(int32(superLUL(i,2)+1),int32(superLUL(i,1)+1))  = superLUL(i,3);
end
% superLU�����U
superluU = [
    0 	 0 	 4
1 	 0 	 -1
1 	 1 	 3.75
2 	 2 	 4
3 	 2 	 -1
3 	 3 	 3.75
4 	 3 	 -1
4 	 4 	 3.73333
5 	 3 	 -1
5 	 4 	 -0.266667
5 	 5 	 3.71429
6 	 0 	 -1
6 	 1 	 -0.25
6 	 5 	 -1
6 	 6 	 3.4641
7 	 7 	 4
8 	 7 	 -1
8 	 8 	 3.75
9 	 9 	 4
10 	 9 	 -1
10 	 10 	 3.75
11 	 10 	 -1
11 	 11 	 3.73333
12 	 10 	 -1
12 	 11 	 -0.266667
12 	 12 	 3.71429
13 	 7 	 -1
13 	 8 	 -0.25
13 	 12 	 -1
13 	 13 	 3.4641
14 	 14 	 3.73205
15 	 15 	 3.46271
16 	 15 	 -1.07695
16 	 16 	 3.04784
17 	 15 	 -0.000370096
17 	 16 	 -0.000115105
17 	 17 	 3.73205
18 	 15 	 -0.000198334
18 	 16 	 -0.0119617
18 	 17 	 -1
18 	 18 	 3.15465
19 	 15 	 -0.290063
19 	 16 	 -0.0902137
19 	 17 	 -0.00558585
19 	 18 	 -0.0018674
19 	 19 	 3.37208
20 	 15 	 -0.00592154
20 	 16 	 -0.00184168
20 	 17 	 -0.287195
20 	 18 	 -0.0769614
20 	 19 	 -1.08985
20 	 20 	 3.02866
21 	 15 	 -9.9167e-05
21 	 16 	 -0.00598086
21 	 17 	 -1
21 	 18 	 -0.556648
21 	 19 	 -0.00201157
21 	 20 	 -1.09119
21 	 21 	 2.97134
22 	 15 	 -1
22 	 16 	 -0.311014
22 	 17 	 -0.000118626
22 	 18 	 -0.00130969
22 	 19 	 -1.09297
22 	 20 	 -0.355186
22 	 21 	 -1.12952
22 	 22 	 2.85418
23 	 15 	 -0.0014875
23 	 16 	 -1.08971
23 	 17 	 -4.13133e-05
23 	 18 	 -0.334415
23 	 19 	 -0.0325773
23 	 20 	 -0.0193516
23 	 21 	 -0.358213
23 	 22 	 -1.26077
23 	 23 	 2.37345
24 	 15 	 -0.556624
24 	 16 	 -0.250072
24 	 17 	 -0.000809129
24 	 18 	 -0.00142847
24 	 19 	 -0.384157
24 	 20 	 -0.137202
24 	 21 	 -0.050736
24 	 22 	 -0.346159
24 	 23 	 -0.2549
24 	 24 	 2.9245
14 	 8 	 -1
14 	 13 	 -0.0666667
15 	 14 	 -1
15 	 4 	 -1
15 	 5 	 -0.0714286
15 	 6 	 -0.0192308
16 	 14 	 -0.287195
16 	 8 	 -1
16 	 13 	 -1.06667
17 	 1 	 -1
17 	 6 	 -0.0666667
18 	 14 	 -0.000740192
18 	 9 	 -1
18 	 10 	 -0.25
18 	 11 	 -1.06667
18 	 12 	 -0.142857
18 	 13 	 -0.0384615
19 	 4 	 -1
19 	 5 	 -1.07143
19 	 6 	 -0.288462
20 	 1 	 -1
20 	 6 	 -1.06667
21 	 14 	 -0.000370096
21 	 11 	 -1
21 	 12 	 -0.0714286
21 	 13 	 -0.0192308
23 	 14 	 -0.00555144
23 	 11 	 -1
23 	 12 	 -1.07143
23 	 13 	 -0.288462
24 	 2 	 -1
24 	 3 	 -0.25
24 	 4 	 -1.06667
24 	 5 	 -0.142857
24 	 6 	 -0.0384615
24 	 14 	 -1
];
sluU = zeros(25,25);
for i=1:length(superluU)
   sluU(int32(superluU(i,2)+1),int32(superluU(i,1)+1))  = superluU(i,3);
end
A = sluL*sluU;
% �������飬������ӵ������������������һ����
% ������˼�ǣ������һ��Ԫ����14������˵��һ���Ǵ�
% ԭ���ĵ�14+1���ƹ�����
perm_r = [
14
15
12
7
13
8
16
9
0
17
1
11
18
2
5
19
3
20
6
21
22
23
10
24
4
];
perm_r = perm_r + 1;
Aold = A;
% ���н��н���
for i=1:length(perm_r)
    A(i,:) = Aold(perm_r(i),:);
end
Aold = A;
% ���н��н���
for i=1:length(perm_r)
        A(:,i) = Aold(:,perm_r(i));
end
err = A - sfull;
level = [
    0
    1987
3124
3423
3617
3766
3899
4024
4148
4252
4345
4428
4504
4570
4621
4664
4697
4723
4743
4758
4773
4786
4799
4811
4821
4830
4839
4847
4855
4864
4872
4880
4887
4892
4898
4903
4908
4912
4916
4920
4924
4928
4932
4937
4941
4943
4945
4947
4949
4951
4953
4955
4957
4959];
level = 4959 - level;