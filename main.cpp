#include <GL/glut.h>
#include <assert.h>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#pragma pack(1)
static struct {
  GLbyte identsize;
  GLbyte colorMapType;
  GLbyte imageType;
  unsigned short colorMapStart;
  unsigned short colorMapLength;
  unsigned char colorMapBits;
  unsigned short xstart;
  unsigned short ystart;
  unsigned short width;
  unsigned short height;
  GLbyte bits;
  GLbyte descriptor;
} tgaHeader;
#pragma pack(8)
enum { sX = 75, sY = 75, n = 600, term = 2 * n - 1 };
static int domain_size = 2;
static double rmin = domain_size / 40.0;
static double cutoff = 2.5 * rmin / pow(2, 1. / 6.);
static double cutoff2 = cutoff * cutoff;
static char **argv;
static int bStoreImages;
static double dt = 1e-4;
static const int words = ceil(n * (n - 1) / 2 / 32.);
static unsigned int table[words];
static double x[n];
static double y[n];
static double vx[n];
static double vy[n];
static double ax[n];
static double ay[n];
static double om[n];
static double to[n];
static double color[n][3];
static const double radius = 0.02;
static const double threshold_collision_p2p = 4 * 0.02 * 0.02;
struct Collision {
  double ux;
  double uy;
};
static std::map<int, Collision> collisions;
static std::map<int, Collision> boundary_collisions;
static std::map<int, Collision> nut_c2p;
static std::map<int, Collision> nut_c2b;

static int imin(int x, int y) { return x < y ? x : y; }
static int imax(int x, int y) { return x > y ? x : y; }
static double dmin(double x, double y) { return x < y ? x : y; }
static double dmax(double x, double y) { return x > y ? x : y; }

static void paintSphere(double x, double y, GLfloat r, GLfloat g, GLfloat b,
                        GLdouble radius) {
  glPushAttrib(GL_ENABLE_BIT);
  glEnable(GL_LIGHTING);
  GLfloat lightColor[] = {r * 1.2f, g * 1.2f, b * 1.2f, 1};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
  glColor3f(r, g, b);
  glPushMatrix();
  glTranslated(x, y, 0);
  glutSolidSphere(radius, 20 * radius / 0.02, 20 * radius / 0.02);
  glPopMatrix();
  glPopAttrib();
};

static struct {
  int *data[sX * sY];
  int n[sX * sY];
  int cap[sX * sY];
} cells;

static struct {
  double r;
  double x, y;
  double u, v;
  double ax, ay;
  double omega, domegadt;
} nut;

static void nut_update(void) {
  nut.ay += -1.0;
  nut.u += dt * nut.ax;
  nut.v += dt * nut.ay;

  nut.x += dt * nut.u;
  nut.y += dt * nut.v;

  nut.omega += dt * nut.domegadt;

  nut.ax = 0;
  nut.ay = 0;
  nut.domegadt = 0;
}

static struct Box {
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
} box;

static double box_distance_to_plane(double x, double y, int p) {
  return box.planes[p][0] * x + box.planes[p][1] * y + box.planes[p][2];
}

static void box_update(void) {
  box.t += dt;

  static bool bGo = false;
  if (!bGo && box.t < 4) {
    box.tinfo.a = 0;
    box.angular_speed = 0;
    box.angle = 0;
  } else if (bGo == false && box.t >= 4) {
    bGo = true;
    box.t = 0;
  } else if (bGo) {
    box.tinfo.a = box.tinfo.a_desired;
    double maxangle = 2. / 180. * M_PI;
    box.angular_speed = box.a2 * maxangle * cos(box.a2 * box.t);
    box.angle = maxangle * sin(box.a2 * box.t);
  }

  box.center[0] = box.tinfo.r1 * cos(box.tinfo.a * box.t + M_PI / 2);
  box.center[1] = box.tinfo.r2 * sin(box.tinfo.a * box.t + M_PI / 2);

  double D[4] = {box.aspect_ratio, 1, box.aspect_ratio, 1};
  for (int i = 0; i < 4; i++) {
    box.planes[i][0] = cos(i * M_PI / 2 + M_PI + box.angle);
    box.planes[i][1] = sin(i * M_PI / 2 + M_PI + box.angle);
    box.planes[i][2] =
        -(box.half_width * D[i] + box.center[0] * box.planes[i][0] +
          box.center[1] * box.planes[i][1]);
  }
}

static void collision_compute(double ux, double uy, double dt, double radius1,
                              double radius2, double r[2], double v1[2],
                              double v2[2], double omega1, double omega2,
                              double force1[2], double force2[2],
                              double &torque1, double &torque2) {
  double kn = 1e5;
  double kt = 2. / 70 * kn;
  double gn = 5e1;
  double gt = 0;
  double IrI = sqrt(pow(r[0], 2) + pow(r[1], 2));
  double invIrI = 1. / sqrt(pow(r[0], 2) + pow(r[1], 2));
  double delta = dmax(0., radius1 + radius2 - IrI);
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
};

static int table_id(int x, int y) {
  int m = imin(x, y);
  return ((m * (term - m)) >> 1) + (x + y - 2 * m - 1);
}

static int table_get(int x, int y) {
  int bit_id = table_id(x, y);
  return (table[bit_id >> 5] >> (bit_id & 0x1f)) & 0x1;
}

static void table_clear(int x, int y) {
  int bit_id = table_id(x, y);
  table[bit_id >> 5] &= ~(1 << (bit_id & 0x1f));
}

static void table_set(int x, int y) {
  int bit_id = table_id(x, y);
  table[bit_id >> 5] |= (1 << (bit_id & 0x1f));
}

int _add_collision(int a, int b) {
  if (pow(x[a] - x[b], 2) + pow(y[a] - y[b], 2) <= threshold_collision_p2p)
    if (!table_get(a, b)) {
      int minab = imin(a, b);
      collisions[minab + n * (a + b - minab)].ux = 0;
      collisions[minab + n * (a + b - minab)].uy = 0;
      table_set(a, b);
      return 1;
    }
  return 0;
}

void _remove_old_collisions() {}

void _add_new_collisions() {
  int all;
  int *cellA;
  int *cellB;
  int code;
  int counter;
  int i;
  int isx;
  int isy;
  int itA;
  int itB;
  int ix;
  int iy;
  int j;
  int m;
  int sizeA;
  int sizeB;

  all = 0;
  counter = 0;
  for (iy = 0; iy < sX; iy++)
    for (ix = 0; ix < sY; ix++) {
      cellA = cells.data[sX * iy + ix];
      sizeA = cells.n[sX * iy + ix];
      if (sizeA == 0)
        continue;

      for (code = 0; code < 9; code++) {
        isx = (ix + (code % 3) - 1);
        isy = (iy + (code / 3) - 1);

        if (isx < 0 || isx >= sX || isy < 0 || isy >= sY)
          continue;

        cellB = cells.data[sX * isy + isx];
        sizeB = cells.n[sX * isy + isx];
        if (sizeB == 0)
          continue;

        if (code != 1 + 3)
          for (itA = 0; itA != sizeA; ++itA)
            for (itB = 0; itB != sizeB; ++itB, all++)
              counter += (int)_add_collision(cellA[itA], cellB[itB]);
        else {
          m = sizeA / 2 + 1;
          for (i = 0; i < m; i++)
            for (j = i + 1; j < sizeA; j++, all++)
              counter += (int)_add_collision(cellA[i], cellA[j]);
        }
      }
    }

  for (i = 0; i < n; i++)
    if (pow(x[i] - nut.x, 2) + pow(y[i] - nut.y, 2) <= pow(radius + nut.r, 2))
      if (nut_c2p.find(i) == nut_c2p.end()) {
        nut_c2p[i].ux = 0;
        nut_c2p[i].uy = 0;
      }
}

void _update_collisions() {
  for (std::map<int, Collision>::iterator it = collisions.begin();
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

    collision_compute(it->second.ux, it->second.uy, dt, radius, radius, r, v1,
                      v2, om[a], om[b], f1, f2, to[a], to[b]);

    ax[a] += f1[0];
    ay[a] += f1[1];

    ax[b] += f2[0];
    ay[b] += f2[1];
  }

  double factor = pow(radius / nut.r, 2);
  for (std::map<int, Collision>::iterator it = nut_c2p.begin();
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

    collision_compute(it->second.ux, it->second.uy, dt, nut.r, radius, r, v1,
                      v2, om[a], nut.omega, f1, f2, to[a], nut.domegadt);

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

void _add_bcollision(int a) {
  for (int p = 0; p < 4; p++) {
    double d = box_distance_to_plane(x[a], y[a], p);

    if (fabs(d) <= radius) {
      int collision_id = a + n * p;

      if (boundary_collisions.find(collision_id) == boundary_collisions.end()) {
        boundary_collisions[collision_id].ux = 0;
        boundary_collisions[collision_id].uy = 0;
      }
    }
  }
}

// detect and add new particle-boundary collisions to the particle-boundary
// collision container
void _add_new_bcollisions() {
  for (int i = 0; i < n; i++)
    _add_bcollision(i);

  for (int p = 0; p < 4; p++) {
    double d = box_distance_to_plane(nut.x, nut.y, p);

    if (fabs(d) <= nut.r) {
      int collision_id = p;

      if (nut_c2b.find(collision_id) == nut_c2b.end()) {
        nut_c2b[collision_id].ux = 0;
        nut_c2b[collision_id].uy = 0;
      }
    }
  }
}

void _remove_old_bcollisions() {
  int *rem;
  int nrem;
  int crem;
  int i;
  int a;
  int p;
  double d;
  nrem = 0;
  crem = 0;
  rem = NULL;
  for (std::map<int, Collision>::iterator it = boundary_collisions.begin();
       it != boundary_collisions.end(); ++it) {
    a = it->first % n;
    p = it->first / n;
    d = box_distance_to_plane(x[a], y[a], p);
    if (fabs(d) > radius) {
      if (nrem >= crem) {
        crem = 2 * crem + 1;
        if ((rem = (int *)realloc(rem, crem * sizeof *rem)) == NULL) {
          fprintf(stderr, "realloc failed\n");
          exit(1);
        }
      }
      rem[nrem++] = it->first;
    }
  }
  for (i = 0; i < nrem; i++)
    boundary_collisions.erase(rem[i]);
  nrem = 0;
  for (p = 0; p < 4; p++) {
    d = box_distance_to_plane(nut.x, nut.y, p);
    if (fabs(d) > nut.r) {
      if (nrem >= crem) {
        crem = 2 * crem + 1;
        if ((rem = (int *)realloc(rem, crem * sizeof *rem)) == NULL) {
          fprintf(stderr, "realloc failed\n");
          exit(1);
        }
      }
      rem[nrem++] = p;
    }
  }
  for (i = 0; i < nrem; i++)
    nut_c2b.erase(rem[i]);
}

// perform 1 step for all particle-boundary collisions
void _update_bcollision(int p, double radius, double x, double y, double u,
                        double v, double omega, Collision &collision,
                        double &fx, double &fy, double &domegadt) {
  // 1. create a mirrored particle at the other side of the boundary (ghost
  // particle)
  // 2. solve for the collision
  // 3. store the force&torque felt by the particle, discard the ghost infos

  double d = dmin(-1e-5, box_distance_to_plane(x, y, p));
  double x1[2] = {x, y};
  double v1[2] = {u, v};

  // 1.
  double ghost[2] = {x1[0] + (-2 * d - 1e-2) * box.planes[p][0],
                     x1[1] + (-2 * d - 1e-2) * box.planes[p][1]};
  double v1DOTn = v1[0] * box.planes[p][0] + v1[1] * box.planes[p][1];

  double vbox[2] = {0, 0};
  vbox[0] = -box.angular_speed * ghost[1] -
            box.tinfo.r1 * box.tinfo.a * sin(box.tinfo.a * box.t);
  vbox[1] = +box.angular_speed * ghost[0] +
            box.tinfo.r2 * box.tinfo.a * sin(box.tinfo.a * box.t);
  double vghost[2] = {vbox[0] - 2 * v1DOTn * box.planes[p][0],
                      vbox[1] - 2 * v1DOTn * box.planes[p][1]};

  double ghost_omega = -omega;

  double r[2] = {x1[0] - ghost[0], x1[1] - ghost[1]};

  double f1[2] = {0, 0};

  {
    double f_dummy[2] = {0, 0};
    double dummy_domegadt = 0;

    // 2.
    collision_compute(collision.ux, collision.uy, dt, radius, radius, r, v1,
                      vghost, omega, ghost_omega, f1, f_dummy, domegadt,
                      dummy_domegadt);

    assert(!isnan(f1[0]));
    assert(!isnan(f1[1]));
  }

  // 3.
  fx += f1[0];
  fy += f1[1];
}

void _update_bcollisions() {
  for (std::map<int, Collision>::iterator it = boundary_collisions.begin();
       it != boundary_collisions.end(); ++it) {
    int a = it->first % n;
    int p = it->first / n;

    _update_bcollision(p, radius, x[a], y[a], vx[a], vy[a], om[a], it->second,
                       ax[a], ay[a], to[a]);
  }
  for (std::map<int, Collision>::iterator it = nut_c2b.begin();
       it != nut_c2b.end(); ++it) {
    int p = it->first;
    _update_bcollision(p, nut.r, nut.x, nut.y, nut.u, nut.v, nut.omega,
                       it->second, nut.ax, nut.ay, nut.domegadt);
  }
}

GLint gltWriteTGA(char *szFileName, int nSizeX, int nSizeY) {
  FILE *pFile;
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
  pBits = (GLbyte *)malloc(lImageSize);
  if (pBits == NULL)
    return 0;
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
  glReadPixels(0, 0, nImageSize[0], nImageSize[1], GL_BGRA_EXT,
               GL_UNSIGNED_BYTE, pBits);
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
  pFile = fopen(szFileName, "wb");
  if (pFile == NULL) {
    free(pBits);
    return 0;
  }
  fwrite(&tgaHeader, sizeof tgaHeader, 1, pFile);
  fwrite(pBits, lImageSize, 1, pFile);
  free(pBits);
  fclose(pFile);
  return 1;
}

static void loop() {
  double a1;
  double a2;
  double final_time = 50;
  double h;
  double hx = 4. / sX;
  double hy = 4. / sY;
  double ix;
  double iy;
  double r1_over_r2;
  double sum = 0;
  double time = 0;
  int a;
  int b = 0;
  int crem;
  int D;
  int i;
  int idx;
  int idy;
  int iframe = 0;
  int k;
  int nrem;
  int nsample = 0;
  int *rem;
  int step;
  int steps_per_frame = 4 * 200;

  argv++;
  r1_over_r2 = atof(*argv++);
  a1 = atof(*argv++);
  a2 = atof(*argv++);
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
  box_update();

  nut.r = 0.2;
  nut.x = 0.0;
  nut.y = -0.7;
  nut.u = 0;
  nut.v = 0;
  nut.ax = 0;
  nut.ay = 0;
  nut.omega = 0;
  nut.domegadt = 0;
  bStoreImages = *argv && strcmp("saveimages", *argv) == 0;

  for (i = 0; i < n; i++) {
    h = 0.04;
    D = 0.8 / h;
    ix = i % D;
    iy = i / D;
    x[i] = -0.4 + (ix + 0.5) * h;
    y[i] = -0.5 + (iy + 0.5) * h;
  }
  for (i = 0; i < n; i++)
    if (x[i] > 0) {
      color[i][0] = 210. / 256;
      color[i][1] = 170. / 256;
      color[i][2] = 58. / 256;
    } else {
      color[i][0] = 121. / 256;
      color[i][1] = 61. / 256;
      color[i][2] = 0. / 256;
    }

  for (step = 0; time < final_time; step++) {
    for (i = 0; i < n; i++)
      ay[i] -= 1;

    nrem = 0;
    crem = 0;
    rem = NULL;
    for (std::map<int, Collision>::iterator it = collisions.begin();
         it != collisions.end(); ++it) {
      a = it->first % n;
      b = it->first / n;
      if (pow(x[a] - x[b], 2) + pow(y[a] - y[b], 2) > threshold_collision_p2p) {
        if (nrem >= crem) {
          crem = 2 * crem + 1;
          if ((rem = (int *)realloc(rem, crem * sizeof *rem)) == NULL) {
            fprintf(stderr, "realloc failed\n");
            exit(1);
          }
        }
        rem[nrem++] = it->first;
        table_clear(a, b);
      }
    }
    for (i = 0; i < nrem; i++)
      collisions.erase(rem[i]);
    nrem = 0;
    for (std::map<int, Collision>::iterator it = nut_c2p.begin();
         it != nut_c2p.end(); ++it) {
      a = it->first;
      if (pow(x[a] - nut.x, 2) + pow(y[a] - nut.y, 2) >
          pow(radius + nut.r, 2)) {
        if (nrem >= crem) {
          crem = 2 * crem + 1;
          if ((rem = (int *)realloc(rem, crem * sizeof *rem)) == NULL) {
            fprintf(stderr, "realloc failed\n");
            exit(1);
          }
        }
        rem[nrem++] = it->first;
        table_clear(a, b);
      }
    }
    for (i = 0; i < nrem; i++)
      nut_c2p.erase(rem[i]);
    free(rem);
    memset(cells.n, 0, sizeof cells.n);
    for (i = 0; i < n; i++) {
      idx = imax(0, imin(sX - 1, (int)floor((x[i] + 2) / hx)));
      idy = imax(0, imin(sY - 1, (int)floor((y[i] + 2) / hy)));
      k = sX * idy + idx;
      if (cells.n[k] >= cells.cap[k]) {
        cells.cap[k] = 2 * cells.cap[k] + 1;
        if ((cells.data[k] = (int *)realloc(
                 cells.data[k], cells.cap[k] * sizeof *cells.data[k])) ==
            NULL) {
          fprintf(stderr, "realloc failed\n");
          exit(1);
        }
      }
      cells.data[k][cells.n[k]++] = i;
    }
    _add_new_collisions();
    _update_collisions();

    _remove_old_bcollisions();
    _add_new_bcollisions();
    _update_bcollisions();

    box_update();
    nut_update();

    for (int i = 0; i < n; i++) {
      vx[i] += dt * ax[i];
      vy[i] += dt * ay[i];
    }
    memset(ax, 0, sizeof ax);
    memset(ay, 0, sizeof ay);
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
      sum += box_distance_to_plane(nut.x, nut.y, 1);
      nsample++;
    }

    if (step % steps_per_frame == 0) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glLineWidth(2.);
      glColor3f(1, 1, 1);
      glBegin(GL_LINE_LOOP);
      double p[4][2] = {-1.01, -1.01, +1.01, -1.01, +1.01, +1.01, -1.01, +1.01};
      double R[2][2] = {cos(box.angle), -sin(box.angle), sin(box.angle),
                        cos(box.angle)};
      for (int i = 0; i < 4; i++) {
        double q[2] = {R[0][0] * p[i][0] * box.half_width * box.aspect_ratio +
                           R[0][1] * p[i][1] * box.half_width,
                       R[1][0] * p[i][0] * box.half_width * box.aspect_ratio +
                           R[1][1] * p[i][1] * box.half_width};
        glVertex2f(q[0] + box.center[0], q[1] + box.center[1]);
      }
      glEnd();

      paintSphere(nut.x, nut.y, 1.3 * 179. / 256, 1.3 * 89. / 256, 0, nut.r);

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

int main(int argc, char **argv0) {
  assert(argc >= 4);
  argv = argv0;
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(512, 512);
  glutCreateWindow("Brazil Nut");
  glutIdleFunc(loop);
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
