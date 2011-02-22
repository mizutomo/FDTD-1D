// FDTDメインプログラム
// 10.01.06 by T.Mizukusa
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "utility.h"

// 二乗関数
double square(double x)
{
	return x * x;
}

// 最小のセルを見つける関数
double find_min_cell(double* ary, int leng)
{
	int i;
	double min;

	min = 1e+30;
	for (i = 0; i < leng; i++) {
		if (ary[i] <= min) {
			min = ary[i];
		}
	}

	return min;
}

// X,Y平面に対して、グリッド間隔を設定(X軸)
void set_delta(double* dx)
{
	int i;
	for (i = 0; i < NX; i++) dx[i] = DX;
}

// CFL値の算出
double calc_cfl_constant()
{
	double dt;

	dt = CFL * 1.0 / (LIGHT * (1.0/DX));

	return dt;
}

// Eyの更新計算
void calc_ey(double* ey, double* hz, double* dx, double dt)
{
	int i;

	for (i = 1; i < NX; i++) {
		ey[i] = ey[i] - dt / (PERMITTIVITY * dx[i]) * (hz[i] - hz[i-1]);
	}
}

// Hzの更新計算
void calc_hz(double* hz, double* ey, double* dx, double dt)
{
	int i;

	for (i = 0; i < NX; i++) {
		hz[i] = hz[i] - dt / (PERMEABILITY * dx[i]) * (ey[i+1] - ey[i]);
	}
}

// Mur1st(Ey計算)
void ey_boundary(double* ey, double* hist_ey, double dt)
{
	double c = LIGHT * dt;

	ey[0]  = hist_ey[1] + (c - DX) / (c + DX) * (ey[1] - hist_ey[0]);
	ey[NX] = hist_ey[2] + (c - DX) / (c + DX) * (ey[NX-1] - hist_ey[3]);
}

// Mur1st(アップデート)
void hist_update(double* ey, double* hist_ey)
{
	hist_ey[0] = ey[0];
	hist_ey[1] = ey[1];
	hist_ey[2] = ey[NX-1];
	hist_ey[3] = ey[NX];
}

// データ出力1(ヘッダ)
void print_point_value_header(FILE* fp)
{
	fprintf(fp, "#Time(1), Stimulus(2), hz[10](3), hz[20](4), hz[30](5), hz[40](6), hz[50](7), hz[60](8), hz[70](9), hz[80](10), hz[90](11)\n");
}

// データ出力1
void print_point_value(FILE* fp, double* data, double time, double stimulus)
{
	fprintf(fp, "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
					time, stimulus, 
					data[100], data[200], data[300], 
					data[400], data[500], data[600], 
					data[700], data[800], data[900]);
}

// データ出力2
void print_line_value(FILE* fp, double* data, int step)
{
	int i;

	fprintf(fp, "# STEP = %d\n", step);

	for (i = 0; i < NX; i++) {
		fprintf(fp, "%d, %g\n", i, data[i]);
	}
	fprintf(fp, "\n");
}

// 入力源(sine波)
/*
double calc_stimulus(int step, double dt)
{
	//	double TT = 9.4e-9;
	//	double TT = 10e-9;
	double time = step * dt;

	//return square(sin(2.0 * PI * time / TT)) / 2.0;
	return sin(2.0 * PI * time / TT) / 2.0;
}
*/

 // ガウシアンパルス
double calc_stimulus(int step, double dt)
{
	double tc = 10e-9; 
	double pw = 0.1e-9;
	double t_current = ((float)step - 0.5) * dt;
	double t_eff;

	t_eff = (t_current - tc) / pw;

	return 1.0 * exp(-1 * (t_eff * t_eff));
}

// FDTDメインルーチン
void calc_fdtd(double* ey, double* hz, double* dx, double dt, double* hist_ey)
{
	int step;
	int last_step = (int)(STOP_TIME/dt);
	double stimulus;
	FILE *fp1, *fp2;
  double tstart, tend;
  long int memsize;

	// Hz[50] ~ Hz[90]の点のデータを出力するファイル
	fp1 = fopen("fdtd_point.csv", "w");
	print_point_value_header(fp1);
	fp2 = fopen("fdtd_line.csv", "w");

	// メインループ
  get_current_time_by_sec(&tstart);
	for (step = 0; step <= last_step; step++) {
		if (step % 100 == 0) {
			printf("%d / %d (%.2f %%)\n", step, last_step, 
						 (double)step/last_step*100);
		}

		// 磁流源の設定
		stimulus = calc_stimulus(step, dt);
		// 波源を(0.01mm, 0.05mm)に設定
		hz[500] += stimulus * dt;

		// Ex&Eyの計算
		calc_ey(ey, hz, dx, dt);

		// Mur1次吸収境界条件
		ey_boundary(ey, hist_ey, dt);
		hist_update(ey, hist_ey);

		// Hzの計算
		calc_hz(hz, ey, dx, dt);

		// 波形出力
		/* print_point_value(fp1, hz, step*dt, stimulus); */
		/* if (step % 100 == 0) print_line_value(fp2, hz, step); */
	}
  get_current_time_by_sec(&tend);
  get_use_memory_size_from_mac(&memsize);
  
	fclose(fp1);
	fclose(fp2);

  printf("All User Time: %.2f [sec], %.2f [min], %.2f [hour]\n", 
         tend-tstart, (tend-tstart)/60, (tend-tstart)/3600);

  if (0 <= memsize  && memsize < 1024) {
    printf("Memory       : %.2f [B]\n" , (float)memsize);
  } else if (1024 <= memsize && memsize/1024 < 1024) {
    printf("Memory       : %.2f [KB]\n" , (float)memsize/1024);
  } else if (1024 <= memsize/1024 && memsize/1024/1024 < 1024) {
    printf("Memory       : %.2f [MB]\n" , (float)memsize/1024/1024);
  } else if (1024 <= memsize/1024/1024) {
    printf("Memory       : %.2f [GB]\n" , (float)memsize/1024/1024/1024);
  }
}

// 出力ファイルヘッダ
int main(int argc, char** argv)
{
	double* dx;

	double* ey;
	double* hz;
	double hist_ey[4] = {0.0, 0.0, 0.0, 0.0};

	double dt;

	// グリッド配列の確保&初期化
	dx = (double *)malloc(sizeof(double) * NX);
	// グリッド初期化
	set_delta(dx);

	// CFL条件の設定(BCK込み)
	dt = calc_cfl_constant();

	// 解析情報の出力
	printf("Input Freq : %g [Hz]\n", FREQ);
	printf("Wave Length: %f [m]\n", LIGHT/FREQ);
	printf("Area X : %f [m]\n", PX);
	printf("Grid X : %d\n", NX);
	printf("Min X  : %f [m]\n", DX);
	printf("CFL    : %f\n", CFL);
	printf("dt     : %g [sec]\n", dt);

	// 電磁界配列の確保
	printf("Allocating Memory...");
	ey = (double*)calloc(NX+1, sizeof(double));
	hz = (double*)calloc(NX,   sizeof(double));
	printf("Done\n");

	// FDTDメインルーチン
	printf("FDTD Calculating...\n");
	calc_fdtd(ey, hz, dx, dt, hist_ey);
	printf("Finished FDTD Calculating.\n");

	free(dx);
	free(ey); free(hz);

	return 0;
}
