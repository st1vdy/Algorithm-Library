/*
 * gen.exe是数据生成器
 * a.exe 和 std.exe 是对拍程序和标程
 */
while (1) {
    system("gen.exe > in.txt");
    system("a.exe < in.txt > a.out");
    system("std.exe < in.txt > std.out");
    if (system("fc a.out std.out")) {
        break;
    }
}