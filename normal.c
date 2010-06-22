// Normal(Starndard)-FDTD
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "common.h"

extern double** ex;
extern double** ey;
extern double** hz;
extern double** cex;
extern double** cexly;
extern double** cey;
extern double** ceylx;
extern double** chzlx;
extern double** chzly;
extern double* dx;
extern double* dy;
extern int nx;
extern int ny;
extern double dt;
extern mat_t material[];
extern conform_t** mat;
extern int bck;

// 関数プロトタイピング
void init_co_ex();
void init_co_ey();
void init_co_hz();

// Standard FDTD係数計算
void initialize_coefficient()
{
	init_co_ex();
	init_co_ey();
	init_co_hz();
}

void init_co_ex()
{
	int i, j;
	double eps, sig;
	double denom;

	for (i = L; i < NX+L; i++) {
		for (j = L; j < NY+L; j++) {
			eps = PERMITTIVITY * material[mat[i-L][j-L].master_material].eps;
			sig = material[mat[i-L][j-L].master_material].sigma;

			denom = eps + sig * dt;
			cex[i][j] = eps / denom;

			if (bck == 1) {
				cexly[i][j] = dt / denom / dy[j] * mat[i-L][j-L].lx_ratio;
			} else {
				cexly[i][j] = dt / denom / dy[j];
			}
		}
	}
}

void init_co_ey()
{
	int i, j;
	double eps, sig;
	double denom;

	for (i = L; i < NX+L; i++) {
		for (j = L; j < NY+L; j++) {
			eps = PERMITTIVITY * material[mat[i-L][j-L].master_material].eps;
			sig = material[mat[i-L][j-L].master_material].sigma;

			denom = eps + sig * dt;
			cey[i][j] = eps / denom;

			if (bck == 1) {
				ceylx[i][j] = dt / denom / dx[i] * mat[i-L][j-L].ly_ratio;
			} else {
				ceylx[i][j] = dt / denom / dx[i];
			}
		}
	}
}

void init_co_hz()
{
	int i, j;
	double mu;
	double common;

	for (i = L; i < NX+L; i++) {
		for (j = L; j < NY+L; j++) {
			mu = PERMEABILITY * material[mat[i-L][j-L].master_material].mu;
			if (bck == 1) {
				if (mat[i-L][j-L].a_ratio != 0.0) {
					common = dt / (mu * mat[i-L][j-L].a_ratio);
				} else {
					common = dt / mu;
				}
			} else {
				common = dt / mu;
			}

			chzlx[i][j] = common / dx[i];
			chzly[i][j] = common / dy[j];
		}
	}
}

// 解析エリアExの計算
void std_calc_ex()
{
	int i, j;

#pragma omp parallel for private(i,j)	
	for (i = L; i < NX+L; i++) {
		for (j = L; j < NY+L; j++) {
			ex[i][j] = cex[i][j] * ex[i][j] 
				+ cexly[i][j] * (hz[i][j] - hz[i][j-1]);
		}
	}
}

// 解析エリアEyの計算
void std_calc_ey()
{
	int i, j;
	
#pragma omp parallel for private(i,j)	
	for (i = L; i < NX+L; i++) {
		for (j = L; j < NY+L; j++) {
			ey[i][j] = cey[i][j] * ey[i][j]
				- ceylx[i][j] * (hz[i][j] - hz[i-1][j]);
		}
	}
}

// 解析エリアHzの計算
void std_calc_hz()
{
	int i, j;
	
#pragma omp parallel for private(i,j)	
	for (i = L; i < NX+L; i++) {
		for (j = L; j < NY+L; j++) {
			hz[i][j] = hz[i][j]
				- chzlx[i][j] * (ey[i+1][j] - ey[i][j])
				+ chzly[i][j] * (ex[i][j+1] - ex[i][j]);
		}
	}
}

