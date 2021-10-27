int zeller(int y, int m, int d) { // 蔡勒公式 返回星期几
	if (m <= 2) y--, m += 12;
	int c = y / 100; y %= 100;
	int w = ((c >> 2) - (c << 1) + y + (y >> 2) + 
		(13 * (m + 1) / 5) + d - 1) % 7;
	if (w < 0) w += 7;
	return (w);
}
int getId(int y, int m, int d) { // 返回到公元1年1月1日的天数
	if (m < 3) { y--; m += 12; }
	return 365 * y + y / 4 - y / 100 + y / 400 +
		(153 * (m - 3) + 2) / 5 + d - 307;
}