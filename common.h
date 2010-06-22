#ifndef __COMMON_H__
#define __COMMON_H__

#define PI            (3.14159)
#define LIGHT         (2.99792458e8)
#define PERMITTIVITY  (8.85418782e-12)
#define PERMEABILITY  (4.0 * PI * 1.0e-7)

// PML吸収境界条件
#define L             (10)
#define M             (4)
#define R0            (1.0e-6)

// 解析エリアの物理サイズ
#define PX            (0.100)     // 0.1m
#define PY            (0.100)     // 0.1m

// 解析領域の配列サイズ
#define K             (4)
#define NX            (400/K)
#define NY            (400/K)

// 解析領域のグリッドサイズ
#define DX            (PX/NX)
#define DY            (PY/NY)

// 入射波の周波数
#define FREQ          (100e+6)
#define TT            (1/FREQ)

// 過渡解析ストップタイム
#define STOP_TIME     (50e-9)

// 共通データ構造
// 材料定数の構造体
typedef struct _mat_t {
	float eps;
	float mu;
	float sigma;
} mat_t;

// Conformalタイプ
typedef struct _conform_t {
	int master_material;
	int second_material;
	float a_ratio;
	float lx_ratio;
	float ly_ratio;
	char conform_flag;
} conform_t;

#endif // __COMMON_H__
