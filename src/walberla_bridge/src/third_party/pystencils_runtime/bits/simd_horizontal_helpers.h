#pragma once

#include <algorithm>

#if defined(__SSE3__)
#include <immintrin.h>

inline double _mm_horizontal_add_pd(double dst, __m128d src) { 
	__m128d _v = src;
	return dst + _mm_cvtsd_f64(_mm_hadd_pd(_v, _v));
}

inline float _mm_horizontal_add_ps(float dst, __m128 src) { 
	__m128 _v = src;
	__m128 _h = _mm_hadd_ps(_v, _v);
	return dst + _mm_cvtss_f32(_mm_add_ps(_h, _mm_movehdup_ps(_h)));
}

inline double _mm_horizontal_mul_pd(double dst, __m128d src) { 
	__m128d _v = src;
	double _r = _mm_cvtsd_f64(_mm_mul_pd(_v, _mm_shuffle_pd(_v, _v, 1)));
	return dst * _r;
}

inline float _mm_horizontal_mul_ps(float dst, __m128 src) { 
	__m128 _v = src;
	__m128 _h = _mm_mul_ps(_v, _mm_shuffle_ps(_v, _v, 177));
	float _r = _mm_cvtss_f32(_mm_mul_ps(_h, _mm_shuffle_ps(_h, _h, 10)));
	return dst * _r;
}

inline double _mm_horizontal_min_pd(double dst, __m128d src) { 
	__m128d _v = src;
	double _r = _mm_cvtsd_f64(_mm_min_pd(_v, _mm_shuffle_pd(_v, _v, 1)));
	return std::min(_r, dst);
}

inline float _mm_horizontal_min_ps(float dst, __m128 src) { 
	__m128 _v = src;
	__m128 _h = _mm_min_ps(_v, _mm_shuffle_ps(_v, _v, 177));
	float _r = _mm_cvtss_f32(_mm_min_ps(_h, _mm_shuffle_ps(_h, _h, 10)));
	return std::min(_r, dst);
}

inline double _mm_horizontal_max_pd(double dst, __m128d src) { 
	__m128d _v = src;
	double _r = _mm_cvtsd_f64(_mm_max_pd(_v, _mm_shuffle_pd(_v, _v, 1)));
	return std::max(_r, dst);
}

inline float _mm_horizontal_max_ps(float dst, __m128 src) { 
	__m128 _v = src;
	__m128 _h = _mm_max_ps(_v, _mm_shuffle_ps(_v, _v, 177));
	float _r = _mm_cvtss_f32(_mm_max_ps(_h, _mm_shuffle_ps(_h, _h, 10)));
	return std::max(_r, dst);
}

#endif

#if defined(__AVX__)
#include <immintrin.h>

inline double _mm256_horizontal_add_pd(double dst, __m256d src) { 
	__m256d _v = src;
	__m256d _h = _mm256_hadd_pd(_v, _v);
	return dst + _mm_cvtsd_f64(_mm_add_pd(_mm256_extractf128_pd(_h,1), _mm256_castpd256_pd128(_h)));
}

inline float _mm256_horizontal_add_ps(float dst, __m256 src) { 
	__m256 _v = src;
	__m256 _h = _mm256_hadd_ps(_v, _v);
	__m128  _i = _mm_add_ps(_mm256_extractf128_ps(_h,1), _mm256_castps256_ps128(_h));
	return dst + _mm_cvtss_f32(_mm_hadd_ps(_i,_i));
}

inline double _mm256_horizontal_mul_pd(double dst, __m256d src) { 
	__m256d _v = src;
	__m128d _w = _mm_mul_pd(_mm256_extractf128_pd(_v,1), _mm256_castpd256_pd128(_v));
	double _r = _mm_cvtsd_f64(_mm_mul_pd(_w, _mm_permute_pd(_w,1))); 
	return dst * _r;
}

inline float _mm256_horizontal_mul_ps(float dst, __m256 src) { 
	__m256 _v = src;
	__m128 _w = _mm_mul_ps(_mm256_extractf128_ps(_v,1), _mm256_castps256_ps128(_v));
	__m128 _h = _mm_mul_ps(_w, _mm_shuffle_ps(_w, _w, 177));
	float _r = _mm_cvtss_f32(_mm_mul_ps(_h, _mm_shuffle_ps(_h, _h, 10)));
	return dst * _r;
}

inline double _mm256_horizontal_min_pd(double dst, __m256d src) { 
	__m256d _v = src;
	__m128d _w = _mm_min_pd(_mm256_extractf128_pd(_v,1), _mm256_castpd256_pd128(_v));
	double _r = _mm_cvtsd_f64(_mm_min_pd(_w, _mm_permute_pd(_w,1))); 
	return std::min(_r, dst);
}

inline float _mm256_horizontal_min_ps(float dst, __m256 src) { 
	__m256 _v = src;
	__m128 _w = _mm_min_ps(_mm256_extractf128_ps(_v,1), _mm256_castps256_ps128(_v));
	__m128 _h = _mm_min_ps(_w, _mm_shuffle_ps(_w, _w, 177));
	float _r = _mm_cvtss_f32(_mm_min_ps(_h, _mm_shuffle_ps(_h, _h, 10)));
	return std::min(_r, dst);
}

inline double _mm256_horizontal_max_pd(double dst, __m256d src) { 
	__m256d _v = src;
	__m128d _w = _mm_max_pd(_mm256_extractf128_pd(_v,1), _mm256_castpd256_pd128(_v));
	double _r = _mm_cvtsd_f64(_mm_max_pd(_w, _mm_permute_pd(_w,1))); 
	return std::max(_r, dst);
}

inline float _mm256_horizontal_max_ps(float dst, __m256 src) { 
	__m256 _v = src;
	__m128 _w = _mm_max_ps(_mm256_extractf128_ps(_v,1), _mm256_castps256_ps128(_v));
	__m128 _h = _mm_max_ps(_w, _mm_shuffle_ps(_w, _w, 177));
	float _r = _mm_cvtss_f32(_mm_max_ps(_h, _mm_shuffle_ps(_h, _h, 10)));
	return std::max(_r, dst);
}

#endif

#if defined(__AVX512F__)
#include <immintrin.h>

inline double _mm512_horizontal_add_pd(double dst, __m512d src) { 
	double _r = _mm512_reduce_add_pd(src);
	return dst + _r;
}

inline float _mm512_horizontal_add_ps(float dst, __m512 src) { 
	float _r = _mm512_reduce_add_ps(src);
	return dst + _r;
}

inline double _mm512_horizontal_mul_pd(double dst, __m512d src) { 
	double _r = _mm512_reduce_mul_pd(src);
	return dst * _r;
}

inline float _mm512_horizontal_mul_ps(float dst, __m512 src) { 
	float _r = _mm512_reduce_mul_ps(src);
	return dst * _r;
}

inline double _mm512_horizontal_min_pd(double dst, __m512d src) { 
	double _r = _mm512_reduce_min_pd(src);
	return std::min(_r, dst);
}

inline float _mm512_horizontal_min_ps(float dst, __m512 src) { 
	float _r = _mm512_reduce_min_ps(src);
	return std::min(_r, dst);
}

inline double _mm512_horizontal_max_pd(double dst, __m512d src) { 
	double _r = _mm512_reduce_max_pd(src);
	return std::max(_r, dst);
}

inline float _mm512_horizontal_max_ps(float dst, __m512 src) { 
	float _r = _mm512_reduce_max_ps(src);
	return std::max(_r, dst);
}

#endif

#if defined(_M_ARM64)
#include <arm_neon.h>

inline double vgetq_horizontal_add_f64(double dst, float64x2_t src) { 
	float64x2_t _v = src;
	double _r = vgetq_lane_f64(_v,0);
	_r += vgetq_lane_f64(_v,1);
	return dst + _r;
}

inline float vget_horizontal_add_f32(float dst, float32x4_t src) { 
	float32x4_t _v = src;
	float32x2_t _w = vadd_f32(vget_high_f32(_v), vget_low_f32(_v));
	float _r = vgetq_lane_f32(_w,0);
	_r += vget_lane_f32(_w,1);
	return dst + _r;
}

inline double vgetq_horizontal_mul_f64(double dst, float64x2_t src) { 
	float64x2_t _v = src;
	double _r = vgetq_lane_f64(_v,0);
	_r *= vgetq_lane_f64(_v,1);
	return dst * _r;
}

inline float vget_horizontal_mul_f32(float dst, float32x4_t src) { 
	float32x4_t _v = src;
	float32x2_t _w = vmul_f32(vget_high_f32(_v), vget_low_f32(_v));
	float _r = vgetq_lane_f32(_w,0);
	_r *= vget_lane_f32(_w,1);
	return dst * _r;
}

inline double vgetq_horizontal_min_f64(double dst, float64x2_t src) { 
	float64x2_t _v = src;
	double _r = vgetq_lane_f64(_v,0);
	_r = std::min(_r, vgetq_lane_f64(_v,1));
	return std::min(_r, dst);
}

inline float vget_horizontal_min_f32(float dst, float32x4_t src) { 
	float32x4_t _v = src;
	float32x2_t _w = vmin_f32(vget_high_f32(_v), vget_low_f32(_v));
	float _r = vgetq_lane_f32(_w,0);
	_r = std::min(_r, vget_lane_f32(_w,1));
	return std::min(_r, dst);
}

inline double vgetq_horizontal_max_f64(double dst, float64x2_t src) { 
	float64x2_t _v = src;
	double _r = vgetq_lane_f64(_v,0);
	_r = std::max(_r, vgetq_lane_f64(_v,1));
	return std::max(_r, dst);
}

inline float vget_horizontal_max_f32(float dst, float32x4_t src) { 
	float32x4_t _v = src;
	float32x2_t _w = vmax_f32(vget_high_f32(_v), vget_low_f32(_v));
	float _r = vgetq_lane_f32(_w,0);
	_r = std::max(_r, vget_lane_f32(_w,1));
	return std::max(_r, dst);
}

#endif