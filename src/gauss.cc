// Iterates several steps of a Jacobi method:
//   p_i,j <= 1/4 * (div_i,j + p_i+2,j + p_i-2,j + p_i,j+2 + p_i,j-2)
// and generate the corresponding code. There's one accumulation for
// the pressure (in a diamond pattern), and another one for the divergence
// (more a Gaussian blur on top of the diamond pattern).

#include <cstdio>
#include <cmath>
#include <algorithm>

static const int MaxN = 30;

struct Kernel {
 public:
  Kernel(int size) : size_(size) {  // recursive initialization
    At(0, 0) = size ? 0 : 1;
    if (size > 0) {
      const Kernel C(size - 1);
      Add(C, 1, 0);
      Add(C,-1, 0);
      Add(C, 0, 1);
      Add(C, 0,-1);
    }
  }

  int64_t At(int i, int j) const { return B_[j * S + i]; }
  int64_t& At(int i, int j) { return B_[j * S + i]; }
  void Add(const Kernel& C, int X = 0, int Y = 0, int64_t Mult = 1) {
    assert(C.size_ <= size_);
    for (int j = -C.size_; j <= C.size_; ++j) {
      for (int i = -C.size_; i <= C.size_; ++i) At(i + X, j + Y) += Mult * C.At(i, j);
    }
  }

  // pretty-print
  void Print(const char str[]) const {
    printf("%s\n", str);
    for (int j = 0; j <= size_; ++j) {
      for (int i = 0; i <= size_; ++i) printf("%7lld ", At(i, j));
      printf("\n");
    }
  }
  void Print(const char var[], const char macro[], int64_t Div) const {
    printf("  float %s = 0.;\n", var);
    for (int j = -size_; j <= size_; ++j) {
      for (int i = -size_; i <= size_; ++i) {
        const int64_t V = At(i, j);
        if (V) printf("  %s += %12lld. * %s(%3d,%3d);\n", var, V, macro, i, j);
      }
    }
    printf("  %s /= %lld.;\n", var, Div);
  }

 protected:
  static const int S = 2 * MaxN + 1;
  int64_t A_[S][S] = { 0 }, * const B_ = &A_[MaxN][MaxN];
  int size_ = 0;
};

////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char* argv[]) {
  const int N = std::min((argc > 1) ? atoi(argv[1]) : 10, MaxN - 1);
  const bool info_only = (argc > 2);
  Kernel A(N), B(N + 1);
  for (int d = 4, i = N - 1; i >= 0; --i, d *= 4) A.Add(Kernel(i), 0, 0, d);

  int64_t Div = 4;
  for (int i = 0; i < N; ++i) Div *= 4;

  if (info_only) {
    A.Print("====== DIV part =====");
    B.Print("======= P part ======");
    printf("Div = %llu\n", Div);
    printf("\n===========\n");
    return 0;
  }
  A.Print("div", "DIV", Div);
  B.Print("p", "P", Div);
  return 0;
}
