#ifndef __FDTD_H__
#define __FDTD_H__

// アルゴリズムを現すEnumerate Type
enum eALG {
	NORMAL,   // Normal-FDTD
	ADE,      // ADE(平均化する)-FDTD
	ADE2      // ADE(平均化しない)-FDTD
};

// プロトタイプ宣言
double find_min_cell(double* ary, int leng);
void set_delta_map(double* d, int leng);
double calc_cfl_constant(float cfl, double* dx, double* dy, 
												 int nx, int ny, 
												 double* min_dx, double* min_dy);
double** alloc_2d_array(int nx, int ny);

void calc_fdtd(double** ex, double** ey, double** hz,
							 double* dx, double* dy,
							 double* hzx, double* hzy,
							 int nx, int ny,
							 double dt, double stop_time, enum eALG alg);

void calc_hz(double** hz,
						 double* dx, double* dy,
						 double** ex, double** ey,
						 int nx, int ny, double dt);

double calc_stimulus(int step, double dt);
void print_point_value(FILE* fp, double** hz, double time, double stimulus);
void print_line_value(FILE* fp, double** hz, int step, int mode, int nx, int ny);
void print_map_value(FILE* fp, double** data, int step, int nx, int ny);
void initialize_2d_array(double** ary, int nx, int ny);

char* get_algorithm_name(enum eALG alg);

#endif //__FDTD_H__
