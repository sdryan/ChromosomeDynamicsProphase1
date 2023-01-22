clc
clear

for k = 1:13

    k
A = load(strcat('homolog_data_',num2str(k), '.dat'));

B(1,2:3) = A(1,2:3);
B(2,2:3) = A(6,2:3);
B(3,2:3) = A(3,2:3);
B(4,2:3) = A(9,2:3);
B(5,2:3) = A(8,2:3);
B(6,2:3) = A(5,2:3);
B(7,2:3) = A(11,2:3);
B(8,2:3) = A(10,2:3);
B(9,2:3) = A(14,2:3);
B(10,2:3) = A(2,2:3);
B(11,2:3) = A(13,2:3);
B(12,2:3) = A(16,2:3);
B(13,2:3) = A(12,2:3);
B(14,2:3) = A(7,2:3);
B(15,2:3) = A(15,2:3);
B(16,2:3) = A(4,2:3);

for i = 1:16
    B(i,1) = i;
end

dlmwrite([strcat('sort_homolog_data_',num2str(k), '.dat')],B,'delimiter','\t','precision',5)
%dlmwrite(['sort_nonhomolog_data_1.dat'],B,'delimiter','\t','precision',5)

end

for k = 1:13
    k

A = load(strcat('nonhomolog_data_',num2str(k), '.dat'));

B(1,2:3) = A(1,2:3);
B(2,2:3) = A(6,2:3);
B(3,2:3) = A(3,2:3);
B(4,2:3) = A(9,2:3);
B(5,2:3) = A(8,2:3);
B(6,2:3) = A(5,2:3);
B(7,2:3) = A(11,2:3);
B(8,2:3) = A(10,2:3);
B(9,2:3) = A(14,2:3);
B(10,2:3) = A(2,2:3);
B(11,2:3) = A(13,2:3);
B(12,2:3) = A(16,2:3);
B(13,2:3) = A(12,2:3);
B(14,2:3) = A(7,2:3);
B(15,2:3) = A(15,2:3);
B(16,2:3) = A(4,2:3);

for i = 1:16
    B(i,1) = i;
end

dlmwrite([strcat('sort_nonhomolog_data_',num2str(k), '.dat')],B,'delimiter','\t','precision',5)
%dlmwrite(['sort_nonhomolog_data_1.dat'],B,'delimiter','\t','precision',5)


end

