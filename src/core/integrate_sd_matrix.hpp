#pragma once

/// stuct for sparse matrix: (blocked version off ELLPACK-R, see \link http://elib.uni-stuttgart.de/opus/volltexte/2010/5033/pdf/DIP_2938.pdf \endlink
/// for dense matrizes use only data

struct wavepart{
public:
  real * vecs;
  real * matrices;
  real * cosines;
  real * sines;
  int num;
  int max;
  int & ldd_short;
  wavepart(int _ldd_short);
  ~wavepart();
private:
  wavepart(wavepart & other);
public:
  void print( int N, int ldd) const;
  void hash_vecs();
  void hash_matrices();
  void hash_cosines();
  void hash_sines();
  void assert_all() const;
#ifdef SD_DEBUG
private:
  unsigned char * vecs_hash;
  unsigned char * matrices_hash;
  unsigned char * cosines_hash;
  unsigned char * sines_hash;
#endif
};



struct matrix{
  //public:
  real * data;
  int  * col_idx;
  int  * row_l;
  wavepart * wavespace;
  bool is_sparse;
  int size;
  int ldd;
  int ldd_short;
  void _init();
  matrix();
  ~matrix();
  void _free_wavespace();
  void _free();
  void printWavespace() const;
  void print() const;
  void printStats() const;
private:
  matrix(matrix & other); // disable copy operator
#ifdef SD_DEBUG
  unsigned char * data_hash;
  unsigned char * dense_hash;
#endif
public:
  real * dense;
  void hash_data();
  void hash_dense();
  void assert_all() const;
};
  
  
#ifdef SD_DEBUG
  void getAndHash(const real * data, const int size, unsigned char * result);
#endif
  
