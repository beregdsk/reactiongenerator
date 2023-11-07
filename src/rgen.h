#pragma once
#if defined(_WIN32)
#  define DLL00_EXPORT_API __declspec(dllexport)
#else
#  define DLL00_EXPORT_API
#endif

bool contains_zero(int* x, int n);
int* minimal_nullspace(int* A, int m, int n, bool* found);
DLL00_EXPORT_API void generate(int* S, int* len_S, int* ref, int dim, int n_vecs, int maxn);