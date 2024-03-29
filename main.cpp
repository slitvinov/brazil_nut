#include <GL/glut.h>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { sX = 75, sY = 75, n = 600, term = 2 * n - 1 };
static char **argv;
static int bStoreImages;
static double x[n];
static double y[n];
static double vx[n];
static double vy[n];
static double ax[n];
static double ay[n];
static double om[n];
static double to[n];
static double color[n][3];
static const int words = ceil(n * (n - 1) / 2 / 32.);
static unsigned int table[words];
static const int xsize = 512;
static const int ysize = 512;
static const int domain_size = 2;
static const double rmin = domain_size / 40.0;
static const double radius = 0.02;
static const double cutoff = 2.5 * rmin / pow(2, 1. / 6.);
static const double cutoff2 = cutoff * cutoff;
static const double dt = 1e-4;
static const double inv_momOFinertia_dt = 1. / (radius * radius) * dt;
static const double threshold_collision_p2p = 4 * 0.02 * 0.02;
static const double maxangle = 2. / 180. * M_PI;
static std::map<int, double[2]> collisions;
static std::map<int, double[2]> boundary_collisions;
static std::map<int, double[2]> nut_c2p;
static std::map<int, double[2]> nut_c2b;
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

static int imin(int x, int y) { return x < y ? x : y; }
static int imax(int x, int y) { return x > y ? x : y; }
static double dmin(double x, double y) { return x < y ? x : y; }
static double dmax(double x, double y) { return x > y ? x : y; }
static void paintSphere(double x, double y, GLfloat r, GLfloat g, GLfloat b,
                        GLdouble radius) {
  glEnable(GL_LIGHTING);
  GLfloat lightColor[] = {r, g, b, 1};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
  glPushMatrix();
  glTranslated(x, y, 0);
  glutSolidSphere(radius, 20 * radius / 0.02, 20 * radius / 0.02);
  glPopMatrix();
};

static double box_distance_to_plane(double x, double y, int p) {
  return box.planes[p][0] * x + box.planes[p][1] * y + box.planes[p][2];
}

static void box_update(void) {
  box.t += dt;

  static int bGo = 0;
  if (!bGo && box.t < 4) {
    box.tinfo.a = 0;
    box.angular_speed = 0;
    box.angle = 0;
  } else if (bGo == 0 && box.t >= 4) {
    bGo = 1;
    box.t = 0;
  } else if (bGo) {
    box.tinfo.a = box.tinfo.a_desired;
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
                              double force1[2], double *torque1,
                              double *torque2) {
  double kn = 1e5;
  double kt = 2. / 70 * kn;
  double gn = 5e1;
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
  double Ft[2] = {-kt * ux, -kt * uy};
  force1[0] = force_factor * (Fn[0] + Ft[0]);
  force1[1] = force_factor * (Fn[1] + Ft[1]);
  *torque1 -= -1. / 2 * (r[0] * Ft[1] - r[1] * Ft[0]);
  *torque2 += -1. / 2 * (r[0] * Ft[1] - r[1] * Ft[0]);
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

int add_collision(int a, int b) {
  int minab;
  if (pow(x[a] - x[b], 2) + pow(y[a] - y[b], 2) <= threshold_collision_p2p)
    if (!table_get(a, b)) {
      minab = imin(a, b);
      collisions[minab + n * (a + b - minab)][0] = 0;
      collisions[minab + n * (a + b - minab)][1] = 0;
      table_set(a, b);
      return 1;
    }
  return 0;
}

void _update_bcollision(int p, double radius, double x, double y, double u,
                        double v, double omega, double collision[2], double *fx,
                        double *fy, double *domegadt) {
  double d = dmin(-1e-5, box_distance_to_plane(x, y, p));
  double v1[2] = {u, v};
  double ghost[2] = {x + (-2 * d - 1e-2) * box.planes[p][0],
                     y + (-2 * d - 1e-2) * box.planes[p][1]};
  double v1DOTn = v1[0] * box.planes[p][0] + v1[1] * box.planes[p][1];
  double vbx = -box.angular_speed * ghost[1] -
               box.tinfo.r1 * box.tinfo.a * sin(box.tinfo.a * box.t);
  double vby = +box.angular_speed * ghost[0] +
               box.tinfo.r2 * box.tinfo.a * sin(box.tinfo.a * box.t);
  double vghost[2] = {vbx - 2 * v1DOTn * box.planes[p][0],
                      vby - 2 * v1DOTn * box.planes[p][1]};
  double r[2] = {x - ghost[0], y - ghost[1]};
  double f1[2];
  double dummy_domegadt;
  collision_compute(collision[0], collision[0], dt, radius, radius, r, v1,
                    vghost, omega, -omega, f1, domegadt, &dummy_domegadt);
  *fx += f1[0];
  *fy += f1[1];
}

static void loop() {
  char path[FILENAME_MAX];
  double a1;
  double a2;
  double d;
  double final_time = 50;
  double h;
  double hx = 4. / sX;
  double hy = 4. / sY;
  double r1_over_r2;
  double sum;
  double time;
  FILE *file;
  GLbyte *pBits;
  int a;
  int all;
  int b;
  int *cellA;
  int *cellB;
  int code;
  int collision_id;
  int counter;
  int crem;
  int D;
  int i;
  int idx;
  int idy;
  int iframe = 0;
  int isx;
  int isy;
  int itA;
  int itB;
  int ix;
  int iy;
  int j;
  int k;
  int m;
  int nrem;
  int nsample = 0;
  int p;
  int *rem;
  int sizeA;
  int sizeB;
  int step;
  int steps_per_frame = 4 * 200;
  std::map<int, double[2]>::iterator it;
  unsigned long lImageSize;

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
    x[i] = -0.4 + (i % D + 0.5) * h;
    y[i] = -0.5 + (i / D + 0.5) * h;
  }
  for (i = 0; i < n; i++)
    if (x[i] > 0) {
      color[i][0] = 63. / 64;
      color[i][1] = 51. / 64;
      color[i][2] = 87. / 320;
    } else {
      color[i][0] = 363. / 640;
      color[i][1] = 183. / 640;
      color[i][2] = 0;
    }

  time = 0;
  for (step = 0; time < final_time; step++) {
    for (i = 0; i < n; i++)
      ay[i] -= 1;

    nrem = 0;
    crem = 0;
    rem = NULL;
    b = 0;
    for (it = collisions.begin(); it != collisions.end(); ++it) {
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
    for (it = nut_c2p.begin(); it != nut_c2p.end(); ++it) {
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
          fprintf(stderr, "main.cpp: error: realloc failed\n");
          exit(1);
        }
      }
      cells.data[k][cells.n[k]++] = i;
    }
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
                counter += add_collision(cellA[itA], cellB[itB]);
          else {
            m = sizeA / 2 + 1;
            for (i = 0; i < m; i++)
              for (j = i + 1; j < sizeA; j++, all++)
                counter += add_collision(cellA[i], cellA[j]);
          }
        }
      }

    for (i = 0; i < n; i++)
      if (pow(x[i] - nut.x, 2) + pow(y[i] - nut.y, 2) <= pow(radius + nut.r, 2))
        if (nut_c2p.find(i) == nut_c2p.end()) {
          nut_c2p[i][0] = 0;
          nut_c2p[i][1] = 0;
        }
    for (it = collisions.begin(); it != collisions.end(); ++it) {
      a = it->first % n;
      b = it->first / n;

      double x1[2] = {x[a], y[a]};
      double x2[2] = {x[b], y[b]};
      double r[2] = {x1[0] - x2[0], x1[1] - x2[1]};

      double v1[2] = {vx[a], vy[a]};
      double v2[2] = {vx[b], vy[b]};
      double f1[2];
      collision_compute(it->second[0], it->second[0], dt, radius, radius, r, v1,
                        v2, om[a], om[b], f1, &to[a], &to[b]);

      ax[a] += f1[0];
      ay[a] += f1[1];

      ax[b] -= f1[0];
      ay[b] -= f1[1];
    }

    double factor = pow(radius / nut.r, 2);
    for (it = nut_c2p.begin(); it != nut_c2p.end(); ++it) {
      a = it->first;

      double x1[2] = {x[a], y[a]};
      double r[2] = {x1[0] - nut.x, x1[1] - nut.y};
      double v1[2] = {vx[a], vy[a]};
      double v2[2] = {nut.u, nut.v};
      double f1[2];
      collision_compute(it->second[0], it->second[0], dt, nut.r, radius, r, v1,
                        v2, om[a], nut.omega, f1, &to[a], &nut.domegadt);
      ax[a] += f1[0];
      ay[a] += f1[1];
      nut.ax -= factor * f1[0];
      nut.ay -= factor * f1[1];
    }

    nrem = 0;
    crem = 0;
    rem = NULL;
    for (std::map<int, double[2]>::iterator it = boundary_collisions.begin();
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

    for (i = 0; i < n; i++) {
      for (p = 0; p < 4; p++) {
        d = box_distance_to_plane(x[i], y[i], p);
        if (fabs(d) <= radius) {
          collision_id = i + n * p;
          if (boundary_collisions.find(collision_id) ==
              boundary_collisions.end()) {
            boundary_collisions[collision_id][0] = 0;
            boundary_collisions[collision_id][0] = 0;
          }
        }
      }
    }
    for (p = 0; p < 4; p++) {
      d = box_distance_to_plane(nut.x, nut.y, p);
      if (fabs(d) <= nut.r) {
        collision_id = p;
        if (nut_c2b.find(collision_id) == nut_c2b.end()) {
          nut_c2b[collision_id][0] = 0;
          nut_c2b[collision_id][0] = 0;
        }
      }
    }
    for (it = boundary_collisions.begin(); it != boundary_collisions.end();
         ++it) {
      a = it->first % n;
      p = it->first / n;
      _update_bcollision(p, radius, x[a], y[a], vx[a], vy[a], om[a], it->second,
                         &ax[a], &ay[a], &to[a]);
    }
    for (it = nut_c2b.begin(); it != nut_c2b.end(); ++it) {
      p = it->first;
      _update_bcollision(p, nut.r, nut.x, nut.y, nut.u, nut.v, nut.omega,
                         it->second, &nut.ax, &nut.ay, &nut.domegadt);
    }

    box_update();
    nut.ay += -1.0;
    nut.u += dt * nut.ax;
    nut.v += dt * nut.ay;
    nut.x += dt * nut.u;
    nut.y += dt * nut.v;
    nut.omega += dt * nut.domegadt;
    nut.ax = 0;
    nut.ay = 0;
    nut.domegadt = 0;
    for (int i = 0; i < n; i++) {
      vx[i] += dt * ax[i];
      vy[i] += dt * ay[i];
    }
    memset(ax, 0, sizeof ax);
    memset(ay, 0, sizeof ay);
    for (int i = 0; i < n; i++) {
      om[i] += to[i] * inv_momOFinertia_dt;
      to[i] = 0;
    }
    for (int i = 0; i < n; i++) {
      x[i] += dt * vx[i];
      y[i] += dt * vy[i];
    }
    sum = 0;
    if (time > 0.95 * final_time) {
      sum += box_distance_to_plane(nut.x, nut.y, 1);
      nsample++;
    }

    if (step % steps_per_frame == 0) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glDisable(GL_LIGHTING);
      glLineWidth(2);
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
      glPopAttrib();
      glutSwapBuffers();
      if (bStoreImages) {
        sprintf(path, "%05d.ppm", iframe++);
        lImageSize = 3 * xsize * ysize;
        pBits = (GLbyte *)malloc(lImageSize);
        if (pBits == NULL)
          exit(1);
        glReadPixels(0, 0, xsize, ysize, GL_RGB, GL_UNSIGNED_BYTE, pBits);
        file = fopen(path, "wb");
        if (file == NULL) {
          free(pBits);
          exit(1);
        }
        fprintf(file, "P6\n%d %d\n255\n", xsize, ysize);
        fwrite(pBits, lImageSize, 1, file);
        free(pBits);
        fclose(file);
        printf("T=%2.2f\n", time);
      }
    }
    time += dt;
  }
  printf("%e\n", sum / nsample);
  exit(0);
}

int main(int argc, char **argv0) {
  argv = argv0;
  glutInit(&argc, argv);
  glutInitWindowSize(xsize, ysize);
  glutCreateWindow("Brazil Nut");
  glutIdleFunc(loop);
  glMatrixMode(GL_PROJECTION);
  glOrtho(-1.5, 1.5, -1.5, 1.5, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHT0);
  GLfloat lightPos[] = {0, 0, -4, 2};
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
  glutMainLoop();
}
