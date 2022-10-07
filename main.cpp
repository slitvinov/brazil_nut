#include <GL/glut.h>
#include <assert.h>
#include <cstring>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;
#pragma pack(1)
struct TGAHEADER {
  GLbyte identsize;    // Size of ID field that follows header (0)
  GLbyte colorMapType; // 0 = None, 1 = paletted
  GLbyte imageType;    // 0 = none, 1 = indexed, 2 = rgb, 3 = grey, +8=rle
  unsigned short colorMapStart;  // First colour map entry
  unsigned short colorMapLength; // Number of colors
  unsigned char colorMapBits;    // bits per palette entry
  unsigned short xstart;         // image x origin
  unsigned short ystart;         // image y origin
  unsigned short width;          // width in pixels
  unsigned short height;         // height in pixels
  GLbyte bits;                   // bits per pixel (8 16, 24, 32)
  GLbyte descriptor;             // image descriptor
};
#pragma pack(8)
static int domain_size = 2;
static double rmin = domain_size / 40.0;
static double cutoff = 2.5 * rmin / pow(2, 1. / 6.);
static double cutoff2 = cutoff * cutoff;
static vector<string> runtime_inputs;
static int bStoreImages;
static double dt = 1e-4;
enum { n = 600 };
static double x[n];
static double y[n];
static double vx[n];
static double vy[n];
static double ax[n];
static double ay[n];
static double om[n];
static double to[n];
static double color[n][3];
enum { sX = 75, sY = 75 };

static void paintSphere(double x, double y, double r, double g, double b,
                        double radius) {
  glPushAttrib(GL_ENABLE_BIT);
  glEnable(GL_LIGHTING);
  GLfloat lightColor[] = {r * 1.2, g * 1.2, b * 1.2, 1};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
  glColor3f(r, g, b);
  glPushMatrix();
  glTranslated(x, y, 0);
  glutSolidSphere(radius, 20 * radius / 0.02, 20 * radius / 0.02);
  glPopMatrix();
  glPopAttrib();
};

struct {
  vector<int> data[sX * sY];
  vector<int> &operator()(int ix, int iy) {
    assert(ix >= 0 && ix < sX);
    assert(iy >= 0 && iy < sY);
    return data[sX * iy + ix];
  }
  void insert(int ix, int iy, int particle_id) {
    assert(ix >= 0 && ix < sX);
    assert(iy >= 0 && iy < sY);
    data[sX * iy + ix].push_back(particle_id);
  }
  void clear() {
    for (int i = 0; i < sX * sY; i++)
      data[i].clear();
  }
} cells;

struct Nut {
  double r;
  double x, y;
  double u, v;
  double ax, ay;
  double omega, domegadt;

  Nut()
      : r(0.20), x(0.0), y(-0.7), u(0), v(0), ax(0), ay(0), omega(0),
        domegadt(0) {}

  void update() {
    ay += -1.0;
    u += dt * ax;
    v += dt * ay;

    x += dt * u;
    y += dt * v;

    omega += dt * domegadt;

    ax = 0;
    ay = 0;
    domegadt = 0;
  }

  void draw() { paintSphere(x, y, 1.3 * 179. / 256, 1.3 * 89. / 256, 0, r); }
} nut;

struct Box {
  struct {
    double a, a_desired, r1, r2;
  } tinfo;
  double half_width;
  double aspect_ratio;
  double a2;

  double planes[4][3];
  double center[2];
  double angle;
  double t;
  double angular_speed;
  void update() {
    t += dt;

    static bool bGo = false;
    if (!bGo && t < 4) {
      tinfo.a = 0;
      angular_speed = 0;
      angle = 0;
    } else if (bGo == false && t >= 4) {
      bGo = true;
      t = 0;
    } else if (bGo) {
      tinfo.a = tinfo.a_desired;
      double maxangle = 2. / 180. * M_PI;
      angular_speed = a2 * maxangle * cos(a2 * t);
      angle = maxangle * sin(a2 * t);
    }

    center[0] = tinfo.r1 * cos(tinfo.a * t + M_PI / 2);
    center[1] = tinfo.r2 * sin(tinfo.a * t + M_PI / 2);

    double D[4] = {aspect_ratio, 1, aspect_ratio, 1};
    for (int i = 0; i < 4; i++) {
      planes[i][0] = cos(i * M_PI / 2 + M_PI + angle);
      planes[i][1] = sin(i * M_PI / 2 + M_PI + angle);
      planes[i][2] = -(half_width * D[i] + center[0] * planes[i][0] +
                       center[1] * planes[i][1]);
    }
  }

  double distance_to_plane(double x, double y, int p) {
    assert(p >= 0 && p < 4);
    return planes[p][0] * x + planes[p][1] * y + planes[p][2];
  }

  void velocity(double x, double y, double &vx, double &vy) {
    vx = -angular_speed * y - tinfo.r1 * tinfo.a * sin(tinfo.a * t);
    vy = +angular_speed * x + tinfo.r2 * tinfo.a * sin(tinfo.a * t);
  }

  double normal(int component, int plane) {
    assert(component >= 0 && component < 2);
    assert(plane >= 0 && plane < 4);

    return planes[plane][component];
  }

  void view() {
    glLineWidth(2.);
    glColor3f(1, 1, 1);
    glBegin(GL_LINE_LOOP);
    double p[4][2] = {-1.01, -1.01, +1.01, -1.01, +1.01, +1.01, -1.01, +1.01};
    double R[2][2] = {cos(angle), -sin(angle), sin(angle), cos(angle)};
    for (int i = 0; i < 4; i++) {
      double q[2] = {R[0][0] * p[i][0] * half_width * aspect_ratio +
                         R[0][1] * p[i][1] * half_width,
                     R[1][0] * p[i][0] * half_width * aspect_ratio +
                         R[1][1] * p[i][1] * half_width};

      glVertex2f(q[0] + center[0], q[1] + center[1]);
    }
    glEnd();
  }
} box;

void box_ini(double r1_over_r2, double a1, double a2) {
  box.angular_speed = 0;
  box.half_width = 1.0;
  box.aspect_ratio = 0.4;
  box.a2 = a2;
  box.angle = 0;
  box.t = 0;
  box.tinfo.a = a1;
  box.tinfo.a_desired = a1;
  box.tinfo.r1 = 0.05;
  box.tinfo.r2 = box.tinfo.r1 / r1_over_r2;
  box.update();
}

struct Collision {
  double ux;
  double uy;
  Collision() { ux = uy = 0; }
  // CORE OF GRANULAR (CHUTE) FLOW
  // this method solves a collision step between two particles
  // INPUT: dt, radius (constant for all the particles), relative position r,
  // v1, v2, respectively omega1, omega2 are the velocities respectively angular
  // velocities of the two colliding particles, OUTPUT: force1, force2: forces
  // felt by the 2 particles torque1, torque2: torques felt by the 2 particles
  void compute(double dt, double radius1, double radius2, double r[2],
               double v1[2], double v2[2], double omega1, double omega2,
               double force1[2], double force2[2], double &torque1,
               double &torque2) {
    double kn = 1e5;
    double kt = 2. / 70 * kn;
    double gn = 5e1;
    double gt = 0;

    double IrI = sqrt(pow(r[0], 2) + pow(r[1], 2));
    double invIrI = 1. / sqrt(pow(r[0], 2) + pow(r[1], 2));
    double delta = max(0., radius1 + radius2 - IrI);

    double v[2] = {v1[0] - v2[0], v1[1] - v2[1]};
    double n[2] = {invIrI * r[0], invIrI * r[1]};

    double vDOTn = v[0] * n[0] + v[1] * n[1];
    double vn[2] = {vDOTn * n[0], vDOTn * n[1]};

    double average_omega = 0.5 * (omega1 + omega2);

    double vt[2] = {v[0] - vn[0] + average_omega * r[1],
                    v[1] - vn[1] - average_omega * r[0]};

    ux += dt * v[0];
    uy += dt * v[1];

    double force_factor = sqrt(delta / (radius1 + radius2));

    double Fn[2] = {
        kn * delta * n[0] - gn * 0.5 * vn[0],
        kn * delta * n[1] - gn * 0.5 * vn[1],
    };

    double Ft[2] = {
        -kt * ux - gt * 0.5 * vt[0],
        -kt * uy - gt * 0.5 * vt[1],
    };

    double f[2] = {force_factor * (Fn[0] + Ft[0]),
                   force_factor * (Fn[1] + Ft[1])};

    force1[0] += f[0];
    force1[1] += f[1];

    force2[0] -= f[0];
    force2[1] -= f[1];

    torque1 -= -1. / 2 * (r[0] * Ft[1] - r[1] * Ft[0]);
    torque2 += -1. / 2 * (r[0] * Ft[1] - r[1] * Ft[0]);
  }
};

// table to check whether collisions were existing or not
struct BitTable {
  int n;              // number of particles
  int term;           // 2*n-1
  int words;          // number of integers
  unsigned int *data; // integers

  // compute the id of the collision of particle x with particle y
  int _bit_id(int x, int y) {
    int m = min(x, y);
#ifndef NDEBUG
    {
      int result = ((m * (term - m)) >> 1) + (x + y - 2 * m - 1);
      assert(result >= 0);
      assert(result < n * (n - 1) / 2);
    }
#endif

    return ((m * (term - m)) >> 1) + (x + y - 2 * m - 1);
  }

  // allocate the integers
  BitTable(int n)
      : n(n), words(ceil(n * (n - 1) / 2 / 32.)), data(NULL), term(2 * n - 1) {
    assert(n >= 0);

    data = new unsigned int[words];
    memset(data, 0, sizeof(unsigned int) * words);
  }

  ~BitTable() { delete[] data; }

  // find out if a collision (particle x, particle y) already exists (in that
  // case it returns 1)
  int operator()(int x, int y) {
    assert(x >= 0 && x < n);
    assert(y >= 0 && y < n);

    int bit_id = _bit_id(x, y);

    return (data[bit_id >> 5] >> (bit_id & 0x1f)) & 0x1;
  }

  // set the bit of the collision (particle x, particle y)
  void set(int x, int y) {
    int bit_id = _bit_id(x, y);

    data[bit_id >> 5] |= (1 << (bit_id & 0x1f));
  }

  // unset the bit of the collision (particle x, particle y)
  void clear(int x, int y) {
    int bit_id = _bit_id(x, y);
    data[bit_id >> 5] &= ~(1 << (bit_id & 0x1f));
  }
};

// this processing elements compute the force felt by particles due to
// particle-particle and particle-boundary collisions
struct GranularFlowCollisionProcessing {
  // global data
  double radius, dt;
  double threshold_collision_p2p;

  double t;

  // containers for collisions
  map<int, Collision> collisions; // particle-particle collision container
  map<int, Collision>
      boundary_collisions; // particle-boundary collision container
  map<int, Collision> nut_c2p;
  map<int, Collision> nut_c2b;

  BitTable bit_table;

  // called in _add_new_collisions
  void _buildCellList() {}

  // called in _add_new_collisions
  bool _add_collision(int a, int b) {
    if (pow(x[a] - x[b], 2) + pow(y[a] - y[b], 2) <= threshold_collision_p2p)
      if (!bit_table(a, b)) {
        int minab = min(a, b);
        collisions[minab + n * (a + b - minab)] = Collision();
        bit_table.set(a, b);

        return true;
      }

    return false;
  }

  void _remove_old_collisions() {
    {
      vector<int> to_remove;

      for (map<int, Collision>::iterator it = collisions.begin();
           it != collisions.end(); ++it) {
        int a = it->first % n;
        int b = it->first / n;

        if (pow(x[a] - x[b], 2) + pow(y[a] - y[b], 2) >
            threshold_collision_p2p) {
          to_remove.push_back(it->first);
          bit_table.clear(a, b);
        }
      }

      for (vector<int>::iterator it = to_remove.begin(); it != to_remove.end();
           it++)
        collisions.erase(*it);
    }

    {
      vector<int> to_remove;

      for (map<int, Collision>::iterator it = nut_c2p.begin();
           it != nut_c2p.end(); ++it) {
        int a = it->first;

        if (pow(x[a] - nut.x, 2) + pow(y[a] - nut.y, 2) >
            pow(radius + nut.r, 2))
          to_remove.push_back(it->first);
      }

      for (vector<int>::iterator it = to_remove.begin(); it != to_remove.end();
           it++)
        nut_c2p.erase(*it);
    }
  }

  void _add_new_collisions() {
    int counter = 0;
    int all = 0;

    for (int iy = 0; iy < sX; iy++)
      for (int ix = 0; ix < sY; ix++) {
        vector<int> &cellA = cells(ix, iy);
        vector<int>::const_iterator itA;
        vector<int>::const_iterator itAS = cellA.begin();
        vector<int>::const_iterator itAE = cellA.end();
        int sizeA = cellA.size();

        if (sizeA == 0)
          continue;

        for (int code = 0; code < 9; code++) {
          int isx = (ix + (code % 3) - 1);
          int isy = (iy + (code / 3) - 1);

          if (isx < 0 || isx >= sX || isy < 0 || isy >= sY)
            continue;

          vector<int> &cellB = cells(isx, isy);

          if (cellB.size() == 0)
            continue;

          if (code != 1 + 3)
            for (itA = itAS; itA != itAE; ++itA)
              for (vector<int>::const_iterator itB = cellB.begin();
                   itB != cellB.end(); ++itB, all++)
                counter += (int)_add_collision(*itA, *itB);
          else {
            assert(cellA == cellB);
            int m = sizeA / 2 + 1;
            for (int i = 0; i < m; i++)
              for (int j = i + 1; j < sizeA; j++, all++)
                counter += (int)_add_collision(cellA[i], cellA[j]);
          }
        }
      }

    for (int i = 0; i < n; i++)
      if (pow(x[i] - nut.x, 2) + pow(y[i] - nut.y, 2) <= pow(radius + nut.r, 2))
        if (nut_c2p.find(i) == nut_c2p.end())
          nut_c2p[i] = Collision();
  }

  void _update_collisions() {
    for (map<int, Collision>::iterator it = collisions.begin();
         it != collisions.end(); ++it) {
      int a = it->first % n;
      int b = it->first / n;

      double x1[2] = {x[a], y[a]};
      double x2[2] = {x[b], y[b]};
      double r[2] = {x1[0] - x2[0], x1[1] - x2[1]};

      double v1[2] = {vx[a], vy[a]};
      double v2[2] = {vx[b], vy[b]};

      double f1[2] = {0, 0};
      double f2[2] = {0, 0};

      it->second.compute(dt, radius, radius, r, v1, v2, om[a], om[b], f1, f2,
                         to[a], to[b]);

      ax[a] += f1[0];
      ay[a] += f1[1];

      ax[b] += f2[0];
      ay[b] += f2[1];
    }

    // m = r^2 P  = 1, M = k^2 r^2 P = k^2 , k = R/r
    double factor = pow(radius / nut.r, 2);
    for (map<int, Collision>::iterator it = nut_c2p.begin();
         it != nut_c2p.end(); ++it) {
      int a = it->first;

      double x1[2] = {x[a], y[a]};
      double r[2] = {x1[0] - nut.x, x1[1] - nut.y};

      assert(pow(x[a] - nut.x, 2) + pow(y[a] - nut.y, 2) <=
             pow(radius + nut.r, 2));

      double v1[2] = {vx[a], vy[a]};
      double v2[2] = {nut.u, nut.v};

      double f1[2] = {0, 0};
      double f2[2] = {0, 0};

      it->second.compute(dt, nut.r, radius, r, v1, v2, om[a], nut.omega, f1, f2,
                         to[a], nut.domegadt);

      assert(!isnan(f1[0]));
      assert(!isnan(f1[1]));

      assert(!isnan(f2[0]));
      assert(!isnan(f2[1]));

      ax[a] += f1[0];
      ay[a] += f1[1];

      nut.ax += factor * f2[0];
      nut.ay += factor * f2[1];
    }
  }

  // called by _add_new_bcollisions
  void _add_bcollision(int a) {
    for (int p = 0; p < 4; p++) {
      //	 double x[2] = {x[a], y[a]};
      double d = box.distance_to_plane(
          x[a], y[a],
          p); // box_planes[p][0]*x[0] +
              // box_planes[p][1]*x[1] + box_planes[p][2] ;

      if (fabs(d) <= radius) {
        int collision_id = a + n * p;

        if (boundary_collisions.find(collision_id) == boundary_collisions.end())
          boundary_collisions[collision_id] = Collision();
      }
    }
  }

  // detect and add new particle-boundary collisions to the particle-boundary
  // collision container
  void _add_new_bcollisions() {
    for (int i = 0; i < n; i++)
      _add_bcollision(i);

    for (int p = 0; p < 4; p++) {
      double d = box.distance_to_plane(nut.x, nut.y, p);

      if (fabs(d) <= nut.r) {
        int collision_id = p;

        if (nut_c2b.find(collision_id) == nut_c2b.end())
          nut_c2b[collision_id] = Collision();
      }
    }
  }

  // remove obsolete particle-boundary collisions from the particle-boundary
  // collision container
  void _remove_old_bcollisions() {
    {
      vector<int> to_remove;

      for (map<int, Collision>::iterator it = boundary_collisions.begin();
           it != boundary_collisions.end(); ++it) {
        int a = it->first % n;
        int p = it->first / n;

        double d = box.distance_to_plane(x[a], y[a], p);
        // box_planes[p][0]*particles.read(a, 0) +
        // box_planes[p][1]*particles.read(a,1) + box_planes[p][2] ;

        if (fabs(d) > radius)
          to_remove.push_back(it->first);
      }

      for (vector<int>::iterator it = to_remove.begin(); it != to_remove.end();
           it++)
        boundary_collisions.erase(*it);
    }

    // nut
    {
      vector<int> to_remove;

      for (int p = 0; p < 4; p++) {
        double d = box.distance_to_plane(nut.x, nut.y, p);

        if (fabs(d) > nut.r)
          to_remove.push_back(p);
      }

      for (vector<int>::iterator it = to_remove.begin(); it != to_remove.end();
           it++)
        nut_c2b.erase(*it);
    }
  }

  // perform 1 step for all particle-boundary collisions
  void _update_bcollision(int p, double radius, double x, double y, double u,
                          double v, double omega, Collision &collision,
                          double &fx, double &fy, double &domegadt) {
    // 1. create a mirrored particle at the other side of the boundary (ghost
    // particle)
    // 2. solve for the collision
    // 3. store the force&torque felt by the particle, discard the ghost infos

    double d = min(-1e-5, box.distance_to_plane(x, y, p));
    double x1[2] = {x, y};
    double v1[2] = {u, v};

    // 1.
    double ghost[2] = {x1[0] + (-2 * d - 1e-2) * box.normal(0, p),
                       x1[1] + (-2 * d - 1e-2) * box.normal(1, p)};
    double v1DOTn = v1[0] * box.normal(0, p) + v1[1] * box.normal(1, p);

    double vbox[2] = {0, 0};
    box.velocity(ghost[0], ghost[1], vbox[0], vbox[1]);

    double vghost[2] = {vbox[0] - 2 * v1DOTn * box.normal(0, p),
                        vbox[1] - 2 * v1DOTn * box.normal(1, p)};

    double ghost_omega = -omega;

    double r[2] = {x1[0] - ghost[0], x1[1] - ghost[1]};

    double f1[2] = {0, 0};

    {
      double f_dummy[2] = {0, 0};
      double dummy_domegadt = 0;

      // 2.
      collision.compute(dt, radius, radius, r, v1, vghost, omega, ghost_omega,
                        f1, f_dummy, domegadt, dummy_domegadt);

      assert(!isnan(f1[0]));
      assert(!isnan(f1[1]));
    }

    // 3.
    fx += f1[0];
    fy += f1[1];
  }

  void _update_bcollisions() {
    for (map<int, Collision>::iterator it = boundary_collisions.begin();
         it != boundary_collisions.end(); ++it) {
      // identify the particle and the boundary that are involved in the
      // particle-boundary collision
      int a = it->first % n;
      int p = it->first / n;

      _update_bcollision(p, radius, x[a], y[a], vx[a], vy[a], om[a], it->second,
                         ax[a], ay[a], to[a]);
    }

    for (map<int, Collision>::iterator it = nut_c2b.begin();
         it != nut_c2b.end(); ++it) {
      // identify the particle and the boundary that are involved in the
      // particle-boundary collision
      int p = it->first;

      _update_bcollision(p, nut.r, nut.x, nut.y, nut.u, nut.v, nut.omega,
                         it->second, nut.ax, nut.ay, nut.domegadt);
    }
  }

  // allocate cell lists covering from [-2, 2]x[-2, 2]
  GranularFlowCollisionProcessing()
      : collisions(), radius(0.02), t(0), bit_table(n),
        threshold_collision_p2p(4 * 0.02 * 0.02) {}

  void operator()() {
    // this is clear, right?
    _remove_old_collisions();
    cells.clear();
    double hx = 4. / sX;
    double hy = 4. / sY;
    for (int i = 0; i < n; i++) {
      int idx = max(0, min(sX - 1, (int)floor((x[i] + 2) / hx)));
      int idy = max(0, min(sY - 1, (int)floor((y[i] + 2) / hy)));
      cells.insert(idx, idy, i);
    }
    _add_new_collisions();
    _update_collisions();

    _remove_old_bcollisions();
    _add_new_bcollisions();
    _update_bcollisions();

    box.update();
    nut.update();
  }

  // retrieve Nut height
  double getNutHeight() { return -box.distance_to_plane(nut.x, nut.y, 1); }
} gfcp;

GLint gltWriteTGA(char *szFileName, int nSizeX, int nSizeY) {
  FILE *pFile;
  struct TGAHEADER tgaHeader;
  unsigned long lImageSize;
  GLbyte *pBits = NULL;
  GLint iViewport[4];
  GLint nImageSize[2];
  bool bUseViewport = (nSizeX == 0 && nSizeY == 0);

  if (bUseViewport) {
    glGetIntegerv(GL_VIEWPORT, iViewport);
    nImageSize[0] = iViewport[2];
    nImageSize[1] = iViewport[3];
  } else {
    nImageSize[0] = nSizeX;
    nImageSize[1] = nSizeY;
  }

  lImageSize = nImageSize[0] * nImageSize[1] * 4;

  // Allocate block. If this doesn't work, go home
  pBits = (GLbyte *)malloc(lImageSize);
  if (pBits == NULL)
    return 0;

  // Read bits from color buffer
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
  glReadPixels(0, 0, nImageSize[0], nImageSize[1], GL_BGRA_EXT,
               GL_UNSIGNED_BYTE, pBits);

  // Initialize the Targa header
  tgaHeader.identsize = 0;
  tgaHeader.colorMapType = 0;
  tgaHeader.imageType = 2;
  tgaHeader.colorMapStart = 0;
  tgaHeader.colorMapLength = 0;
  tgaHeader.colorMapBits = 0;
  tgaHeader.xstart = 0;
  tgaHeader.ystart = 0;
  tgaHeader.width = nImageSize[0];
  tgaHeader.height = nImageSize[1];
  tgaHeader.bits = 32;
  tgaHeader.descriptor = 8;

  // Attempt to open the file
  pFile = fopen(szFileName, "wb");
  if (pFile == NULL) {
    free(pBits); // Free buffer and return error
    return 0;
  }

  // Write the header
  fwrite(&tgaHeader, sizeof(TGAHEADER), 1, pFile);

  // Write the image data
  fwrite(pBits, lImageSize, 1, pFile);

  // Free temporary buffer and close the file
  free(pBits);
  fclose(pFile);

  // Success!
  return 1;
}

void launch_fitnessFunction() {
  double params[3];

  if (runtime_inputs.size() < 4) {
    printf("Aborting. It is missing some input.\n");
    abort();
  }

  params[0] = atof(runtime_inputs[1].c_str());
  params[1] = atof(runtime_inputs[2].c_str());
  params[2] = atof(runtime_inputs[3].c_str());
  box_ini(params[0], params[1], params[2]);

  bStoreImages = 0;
  if (runtime_inputs.size() > 4 && runtime_inputs[4] == "saveimages")
    bStoreImages = 1;

  for (int i = 0; i < n; i++) {
    double h = 0.04;
    int D = 0.8 / h;
    double ix = i % D;
    double iy = i / D;
    x[i] = -0.4 + (ix + 0.5) * h;
    y[i] = -0.5 + (iy + 0.5) * h;
  }
  //  gfcp = new GranularFlowCollisionProcessing();
  for (int i = 0; i < n; i++)
    if (x[i] > 0) {
      color[i][0] = 210. / 256;
      color[i][1] = 170. / 256;
      color[i][2] = 58. / 256;
    } else {
      color[i][0] = 121. / 256;
      color[i][1] = 61. / 256;
      color[i][2] = 0. / 256;
    }
  int steps_per_frame = 4 * 200;
  double final_time = 50;

  double time = 0;
  int iframe = 0;

  // collect the samples about the "altitude" of the nut with respect to the
  // box's bottom
  double sum = 0;
  int nsample = 0;
  for (int step = 0; time < final_time; step++) {
    for (int i = 0; i < n; i++)
      ay[i] -= 1;
    gfcp();
    for (int i = 0; i < n; i++) {
      vx[i] += dt * ax[i];
      vy[i] += dt * ay[i];
    }
    memset(ax, 0, sizeof(double) * n);
    memset(ay, 0, sizeof(double) * n);
    double radius = 0.02;
    double inv_momOFinertia = 1. / (radius * radius);
    double factor = dt * inv_momOFinertia;
    for (int i = 0; i < n; i++) {
      om[i] += factor * to[i];
      to[i] = 0;
    }
    for (int i = 0; i < n; i++) {
      x[i] += dt * vx[i];
      y[i] += dt * vy[i];
    }

    // if we are about to end, start measure the nuts y-position
    if (time > 0.95 * final_time) {
      sum += box.distance_to_plane(nut.x, nut.y, 1);
      nsample++;
    }

    // do some frames
    if (step % steps_per_frame == 0) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      box.view();
      nut.draw();

      for (int i = 0; i < n; i++)
        paintSphere(x[i], y[i], color[i][0], color[i][1], color[i][2], 0.02);

      // OPEN GL STUFF
      glPopAttrib();
      glutSwapBuffers();

      // write to file
      if (bStoreImages) {
        char frame_name[300];
        sprintf(frame_name, "%05d.tga", iframe++);
        gltWriteTGA(frame_name, 0, 0);
      }

      if (bStoreImages)
        printf("T=%2.2f\n", time);
    }

    time += dt;
  }

  printf("%e\n", sum / nsample);

  exit(0);
}

int main(int argc, char **argv) {
  assert(argc >= 4);
  for (int i = 0; i < argc; i++)
    runtime_inputs.push_back(string(argv[i]));

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(512, 512);
  glutCreateWindow("FitnessFunction");
  glutIdleFunc(launch_fitnessFunction);
  glPointSize(4);
  glClearColor(0, 0, 0, 0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1.5, 1.5, -1.5, 1.5, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  GLfloat lightPos[] = {0.0f, 0.0, -4., 2.0f};
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
  glutMainLoop();
}
