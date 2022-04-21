// Fluid flow simulation proto
//
// Copyright (c) 2022 Pascal Massimino (skal -> "pascal.massimino@gmail.com")
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <cstdio>
#include <cmath>
#include <cassert>

#include <algorithm>
#include <string>
#include <vector>

using std::string;
using std::vector;

#define DBG 1

#if defined(SIM_HAVE_SDL)
#include <SDL.h>
#endif  // SIM_HAVE_SDL

////////////////////////////////////////////////////////////////////////////////

typedef float f1;
struct f2 {
  f1 x = 0.f, y = 0.f;
  f2(float x = 0.f, float y = 0.f) : x(x), y(y) {}
  f2(int x, int y) : x((float)x), y((float)y) {}
  f2 operator*(float amp) const { return {x * amp, y * amp}; }
  f2 operator*(const f2& xy) const { return {x * xy.x, y * xy.y}; }
  f2 operator+(const f2& xy) const { return {x + xy.x, y + xy.y}; }
  f2 operator-(const f2& xy) const { return {x - xy.x, y - xy.y}; }
  float norm() const { return sqrt(x * x + y * y); }
};
f2 operator*(float amp, const f2& xy) { return xy * amp; }

////////////////////////////////////////////////////////////////////////////////

template<typename T> struct Field {
  typedef vector<vector<T>> Array;

  void Init(int W, int H) {
    W_ = W - 1;
    H_ = H - 1;
    data_.resize(H);
    for (auto& v : data_) v.resize(W);
  }
  void Reset() { for(auto& v : data_) std::fill(v.begin(), v.end(), T()); }

  T get(f2 xy) const {  // bilinear interpolation
    float x = std::max(0.f, std::min(xy.x, 1.f * W_));
    float y = std::max(0.f, std::min(xy.y, 1.f * H_));
    const int X0 = (int)x, Y0 = (int)y;
    const int X1 = std::min((int)(x + 1.), W_);
    const int Y1 = std::min((int)(y + 1.), H_);
    x -= X0;
    y -= Y0;
    const float e00 = (1.f - x) * (1.f - y);
    const float e01 = x * (1.f - y);
    const float e10 = (1.f - x) * y;
    const float e11 = x * y;
    const T* const p0 = data_[Y0].data();
    const T* const p1 = data_[Y1].data();
    return e00 * p0[X0] + e01 * p0[X1] + e10 * p1[X0] + e11 * p1[X1];
  }
  T at(int x, int y) const {  // clamped
    x = std::max(0, std::min(x, W_));
    y = std::max(0, std::min(y, H_));
    return data_[y][x];
  }

  template<typename FN> void Process(FN fn) {
    int y = 0;
    for (auto& v : data_) {
      int x = 0;
      for (auto& cell : v) fn(x++, y, cell);
      ++y;
    }
  }
  template<typename P, typename FN> void Emit(FN fn, P& p) const {
    int y = 0;
    for (auto& v : data_) {
      for (auto& cell : v) fn(p, cell);
    }
  }

  int W() const { return W_ + 1; }
  int H() const { return H_ + 1; }

  Array& data() const { return data_; }
  const vector<T>& operator[](int y) const { return data_[y]; }
  vector<T>& operator[](int y) { return data_[y]; }

  void GetDisplacement(const Field<T>& src, float dt) {

    Process([&](int x, int y, f2& out) {
      const f2 xy(x, y);
      const f2 k1 = dt * src.at(x, y);
      const f2 k2 = dt * src.get(xy - 0.5 * k1);
#if DBG     // simple Mid-Point step
      out = k2;
#else       // 4th-order Runge-Kutta
      const f2 k3 = dt * src.get(xy - 0.5f * k2);
      const f2 k4 = dt * src.get(xy -        k3);
      // final displacement
      out = (0.5f * (k1 + k4) + k2 + k3) * (1.f / 3.0f);
#endif
    });
  }
  void Apply(const Field<f2>& dsp) {   // apply the advection field 'v' to oneself
    const Field<T> tmp = *this;
    Process([&](int x, int y, T& out) {
      const f2 xy(x, y);
      out = tmp.get(xy - dsp.get(xy));
    });
  }   

 protected:
  friend void Swap(Field& A, Field& B);
  int W_, H_;
  Array data_;
};

template<typename T>
static void Swap(Field<T>& A, Field<T>& B) {
  assert(A.W_ == B.W_ && A.H_ == B.H_);
  std::swap(A.data_, B.data_);
}

typedef Field<f1> Vec1;
typedef Field<f2> Vec2;

////////////////////////////////////////////////////////////////////////////////
// Pressure solver

void Divergence(const Vec2& v, Vec1& div) {
  div.Process([&](int x, int y, float& f) {
    const float dv_dx = v.at(x + 1, y).x - v.at(x - 1, y).x;
    const float dv_dy = v.at(x, y + 1).y - v.at(x, y - 1).y;
    f = -0.5 * (dv_dx + dv_dy);
  });
}

void SolvePressure(const Vec1& div, Vec1& p, int nb_iters) {
  double diff = 0.;
  for (int n = 0; n < nb_iters; ++n) {
    Vec1 p1 = p;  // TODO: swap instead
    diff = 0.;
    p.Process([&](int x, int y, float& out) {
      const float f = 0.25 * (div.at(x, y)
        + p.at(x + 1, y) + p.at(x - 1, y)
        + p.at(x, y + 1) + p.at(x, y - 1));
      diff += fabs(out - f);
      out = f;
    });
    diff /= p.W() * p.H();
    if (n > 3 && diff < 1.e-8) break;
  }
  // printf("diff: %f\n", diff);
}

// unrolled 9x9 version
void SolvePressure9(const Vec1& div, Vec1& pressure, int nb_iters) {
#define D(DX, DY) div.at(x + (DX), y + (DY))
#define P(DX, DY) pressure.at(x + (DX), y + (DY))
  pressure.Process([&](int x, int y, float& out) {
    float d = 0.;
    d +=   1. * D(  0, -3);
    d +=   3. * D( -1, -2);
    d +=   4. * D(  0, -2);
    d +=   3. * D(  1, -2);
    d +=   3. * D( -2, -1);
    d +=   8. * D( -1, -1);
    d +=  25. * D(  0, -1);
    d +=   8. * D(  1, -1);
    d +=   3. * D(  2, -1);
    d +=   1. * D( -3,  0);
    d +=   4. * D( -2,  0);
    d +=  25. * D( -1,  0);
    d +=  80. * D(  0,  0);
    d +=  25. * D(  1,  0);
    d +=   4. * D(  2,  0);
    d +=   1. * D(  3,  0);
    d +=   3. * D( -2,  1);
    d +=   8. * D( -1,  1);
    d +=  25. * D(  0,  1);
    d +=   8. * D(  1,  1);
    d +=   3. * D(  2,  1);
    d +=   3. * D( -1,  2);
    d +=   4. * D(  0,  2);
    d +=   3. * D(  1,  2);
    d +=   1. * D(  0,  3);
    float p = 0.;
    p +=   1. * P(  0, -4);
    p +=   4. * P( -1, -3);
    p +=   4. * P(  1, -3);
    p +=   6. * P( -2, -2);
    p +=  16. * P(  0, -2);
    p +=   6. * P(  2, -2);
    p +=   4. * P( -3, -1);
    p +=  24. * P( -1, -1);
    p +=  24. * P(  1, -1);
    p +=   4. * P(  3, -1);
    p +=   1. * P( -4,  0);
    p +=  16. * P( -2,  0);
    p +=  36. * P(  0,  0);
    p +=  16. * P(  2,  0);
    p +=   1. * P(  4,  0);
    p +=   4. * P( -3,  1);
    p +=  24. * P( -1,  1);
    p +=  24. * P(  1,  1);
    p +=   4. * P(  3,  1);
    p +=   6. * P( -2,  2);
    p +=  16. * P(  0,  2);
    p +=   6. * P(  2,  2);
    p +=   4. * P( -1,  3);
    p +=   4. * P(  1,  3);
    p +=   1. * P(  0,  4);
    out = (d + p) / 256.;
  });
#undef P
#undef D
}

void ApplyGradient(const Vec1& p, Vec2& v) {
  v.Process([&](int x, int y, f2& out) {
    const f2 grad{
      p.at(x + 1, y) - p.at(x - 1, y),
      p.at(x, y + 1) - p.at(x, y - 1)
    };
    out = out - 0.5f * grad;
  });
}

////////////////////////////////////////////////////////////////////////////////
// I/O

static const auto Dump1 = [](FILE* f, float c) { fputc((int)(c * 255), f); };
static const auto Dump2 = [](FILE* f, f2 c) { fputc((int)(c.norm() * 255), f); };
static const auto DumpX = [](FILE* f, f2 c) { fputc((int)(c.x * 255), f); };
static const auto DumpY = [](FILE* f, f2 c) { fputc((int)(c.y * 255), f); };

template<typename T, typename FN>
bool Save(int n, const char* out, T& field, FN fn) {
  char tmp[100];
  snprintf(tmp, sizeof(tmp), "%s.%.3d.ppm", out, n);
  FILE* const file = fopen(tmp, "wb");
  if (file == NULL) {
    fprintf(stderr, "could not open '%s' for writing\n", tmp);
    return false;
  }
  fprintf(file, "P5\n%d %d\n255\n", field.W(), field.H());
  field.Emit(fn, file);
  fprintf(file, "\n");
  fclose(file);
  printf("saved output '%s'\n", tmp);
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// collective object, that's more practical

struct Simulation {
  // data fields
  Vec2 v, dsp;
  Vec1 p, c, div;

  // params
  float dt = 1.;
  float rho = 0.05;
  int W = 600;
  int H = 300;
  int N = 20000;
  int period = 0;
  int nb_iterations = 8;
  float Cx = 0.25, Cy = 0.52, R = 0.09;
  const char* out = "toto";
  bool use_display = true;

  bool ParseArgs(int argc, char* argv[]) {
    for (int c = 1; c < argc; ++c) {
      if (!strcmp(argv[c], "-dt")) {
        if (++c < argc) dt = atof(argv[c]);
      } else if (!strcmp(argv[c], "-W")) {
        if (++c < argc) W = std::max(10, atoi(argv[c]));
      } else if (!strcmp(argv[c], "-H")) {
        if (++c < argc) H = std::max(10, atoi(argv[c]));
      } else if (!strcmp(argv[c], "-N")) {
        if (++c < argc) N = std::max(10, atoi(argv[c]));
      } else if (!strcmp(argv[c], "-n")) {
        if (++c < argc) nb_iterations = std::min(20, atoi(argv[c]));
      } else if (!strcmp(argv[c], "-p")) {
        if (++c < argc) period = std::max(0, atoi(argv[c]));
      } else if (!strcmp(argv[c], "-R")) {
        if (++c < argc) R = atof(argv[c]);
      } else if (!strcmp(argv[c], "-rho")) {
        if (++c < argc) rho = std::max(0.00001f, (float)atof(argv[c]));
      } else if (!strcmp(argv[c], "-o")) {
        if (++c < argc) out = argv[c];
      } else if (!strcmp(argv[c], "-no_dsp")) {
        use_display = false;
      } else if (!strcmp(argv[c], "-h")) {
        printf("sim [options]\n");
        printf("  -dt <float> ........ time step\n");
        printf("  -W <int> ........... width\n");
        printf("  -H <int> ........... height\n");
        printf("  -n <int> ........... number of Jacobi solve iterations\n");
        printf("  -N <int> ........... number of steps\n");
        printf("  -p <int> ........... dumping period (0=off)\n");
        printf("  -R <float> ......... ball radius\n");
        printf("  -rho <float> ....... density\n");
        printf("  -o <string> ........ prefix for file dumps\n");
        return false;
      }
    }
    printf("Simulating %d x %d: dt=%.2f rho=%.4f  R=%.3f\n", W, H, dt, rho, R);
    return true;
  }

  bool InitAll() {
    v.Init(W, H);
    p.Init(W, H);
    c.Init(W, H);
    div.Init(W, H);
    dsp.Init(W, H);
    c.Reset();
    p.Reset();  // zero initial pressure
    // Init velocity
    v.Process([&](int x, int y, f2& out) {
      out = { 1.0f, 0.0f };
      const f2 rnd = {(float)drand48() * 0.1f, (float)(drand48() - .5) * 0.1f};
      out = out + rnd;
    });
    return true;
  }

  void OneStep() {
    // Boundary conditions
    const auto is_ball = [&](int x, int y) -> bool {
      const f2 d = {Cx - 1.f * x / H, Cy - 1.f * y / H};
      return (d.norm() < R);
    };
    const auto BoundaryV = [&](int x, int y, f2& f) {
      if (y < 2 || y >= H - 2) f = { 1.0f, 0.0f };
      if (x < 2) f = { 1.0f, 0.0f };
      if (x >= W - 2) f = { 1.0f, 0.0f };
      if (is_ball(x, y)) f = { 0.0f, 0.0f };
    };
    const auto BoundaryC = [&](int x, int y, float& c) {
      if (x < 5) c = cos(y * 80. / H) > .6 ? 1. : 0.;
      if (is_ball(x, y)) c = 0.;
    };

    // go
    dsp.GetDisplacement(v, dt);     // compute advection field    
    v.Apply(dsp);  
    v.Process(BoundaryV);
    if (true) {
      Divergence(v, div);
      SolvePressure(div, p, nb_iterations);
      ApplyGradient(p, v);
    }
    c.Process(BoundaryC);
    c.Apply(v);
  }

#if defined(SIM_HAVE_SDL)
  SDL_Surface* screen = NULL;
  SDL_Surface* surface = NULL;
  int show = 1;  // 1: trace, 2: divergence, 3:pressure, 4: v.norm()
  bool InitDisplay() {
    if (use_display) {
      SDL_Init(SDL_INIT_VIDEO);
      screen = SDL_SetVideoMode(W, H, 32, SDL_SWSURFACE);
      if (screen == NULL) {
        fprintf(stderr, "Unable to set video mode (32bpp %dx%d)!\n", W, H);
        return false;
      }
      surface = SDL_CreateRGBSurface(SDL_SWSURFACE, W, H, 32,
                                     0x000000ffu,   // R mask
                                     0x0000ff00u,   // G mask
                                     0x00ff0000u,   // B mask
                                     0xff000000u);  // A mask
      if (surface == NULL) {
        fprintf(stderr, "Unable to create %dx%d RGBA surface!\n", W, H);
        SDL_FreeSurface(screen);
        screen = NULL;
        return false;
      }
    }
    return true;
  }
  bool Show() {
    if (!use_display) return false;
    if (surface != NULL) {
      if (SDL_MUSTLOCK(surface)) SDL_LockSurface(surface);
      uint32_t* const dst = (uint32_t*)surface->pixels;
      const int stride = surface->pitch / sizeof(uint32_t);
      const auto display_1d = [&](int x, int y, float& v) {
          const uint8_t C = std::min(255.f, 128.f + v * 128.f);
          dst[x + y * stride] = (C * 0x010101u) | 0xff000000u;
      };
      const auto display_2d = [&](int x, int y, f2& v) {
          const uint8_t C = std::min(255.f, v.norm() * 255.f);
          dst[x + y * stride] = (C * 0x010101u) | 0xff000000u;
      };
      if (show == 1) c.Process(display_1d);
      if (show == 2) div.Process(display_1d);
      if (show == 3) p.Process(display_1d);
      if (show == 4) v.Process(display_2d);
      if (show == 5) dsp.Process(display_2d);
      if (SDL_MUSTLOCK(surface)) SDL_UnlockSurface(surface);
      if (SDL_BlitSurface(surface, NULL, screen, NULL) ||
          SDL_Flip(screen)) {
        fprintf(stderr, "Error blitting surface!\n");
        use_display = false;
        return true;
      } 
      SDL_Event event;
      while (SDL_PollEvent(&event)) {
        // printf("#event: %d\n", (int)event.type);
        switch (event.type) {
          case SDL_KEYUP:
            switch (event.key.keysym.sym) {
              case SDLK_q: return true; break;
              case SDLK_f: {
                uint32_t flags = surface->flags ^ SDL_FULLSCREEN;
                surface = SDL_SetVideoMode(surface->w, surface->h, surface->format->BitsPerPixel, flags);
              }
              break;
              default: show = 1; break;
            }
          case SDL_KEYDOWN:
            switch (event.key.keysym.sym) {
              default:
              case SDLK_1: show = 1; break;
              case SDLK_2: show = 2; break;
              case SDLK_3: show = 3; break;
              case SDLK_4: show = 4; break;
              case SDLK_5: show = 5; break;
            }
            break;
          case SDL_QUIT: return true; break;
          default: break;
        }
      }
    }
    return false;
  }
  void CloseDisplay() {
    if (use_display) {
      if (surface != NULL) SDL_FreeSurface(surface);
      if (screen != NULL) SDL_FreeSurface(screen);
      SDL_Quit();
    }
  }
#else
  bool InitDisplay() { return true; }
  bool Show() { return false; }
  void CloseDisplay() {}
#endif

  ~Simulation() { CloseDisplay(); }
};

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {  // *not* 'const char*', because SDL !
  Simulation S;

  if (!S.ParseArgs(argc, argv)) return 1;
  if (!S.InitDisplay()) return 2;
  if (!S.InitAll()) return 3;

  // go!
  bool stopped = false;  // for user-abort (or end of loop)
  for (int n = 0; !stopped; ++n) {  // simulation steps
    S.OneStep();

    // Show something
    stopped = S.Show() || (n == S.N - 1);  // in case

    // Save regularly
    if (stopped || (S.period > 0 && (n % S.period) == 0)) {
      if (S.out != nullptr) Save(n, S.out, S.c, Dump1);
    }
  }

  return 0;  
}
