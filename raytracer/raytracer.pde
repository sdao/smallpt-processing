// Adapted from smallpt <http://www.kevinbeason.com/smallpt/>.

import java.util.Random;

static class Vec3d {
  static final Vec3d ZERO = new Vec3d(0.0);

  double x;
  double y;
  double z;

  Vec3d(double k) {
    this.x = k;
    this.y = k;
    this.z = k;
  }
  Vec3d(double x, double y, double z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  Vec3d gammaCorrected() {
    return new Vec3d(
      Math.pow(clamp(this.x), 1/2.2),
      Math.pow(clamp(this.y), 1/2.2),
      Math.pow(clamp(this.z), 1/2.2));
  }
  
  Vec3d plus(Vec3d rhs) {
    return new Vec3d(this.x + rhs.x, this.y + rhs.y, this.z + rhs.z);
  }
  Vec3d minus(Vec3d rhs) {
    return new Vec3d(this.x - rhs.x, this.y - rhs.y, this.z - rhs.z);
  }
  Vec3d mult(double k) {
    return new Vec3d(this.x * k, this.y * k, this.z * k);
  }
  Vec3d mult(Vec3d rhs) {
    return new Vec3d(this.x * rhs.x, this.y * rhs.y, this.z * rhs.z);
  }
  Vec3d div(double k) {
    return mult(1.0 / k);
  }
  double dot(Vec3d rhs) {
    return this.x * rhs.x + this.y * rhs.y + this.z * rhs.z;
  }
  double length() {
    return Math.sqrt(this.dot(this));
  }
  Vec3d normalized() {
    return this.div(length());
  }
  Vec3d cross(Vec3d rhs) {
    return new Vec3d(
      this.y * rhs.z - this.z * rhs.y,
      this.z * rhs.x - this.x * rhs.z,
      this.x * rhs.y - this.y * rhs.x);
  }
  void clear() {
    this.x = 0.0;
    this.y = 0.0;
    this.z = 0.0;
  }
  void accum(Vec3d other) {
    this.x += other.x;
    this.y += other.y;
    this.z += other.z;
  }
}

static class Ray {
  Vec3d origin;
  Vec3d direction;
  Ray(Vec3d origin, Vec3d direction) {
    this.origin = origin;
    this.direction = direction;
  }
  Vec3d atDistance(double distance) {
    return this.origin.plus(this.direction.mult(distance));
  }
}

static enum Material {
  DIFFUSE,
  SPECULAR,
  REFRACTIVE,
}

static final double NO_HIT = 0.0;

static class Sphere {
  double radius;
  Vec3d position;
  
  Material material;
  Vec3d surfaceColor;
  Vec3d emissionColor;
 
  Sphere(
      double radius,
      Vec3d position,
      Material material,
      Vec3d surfaceColor,
      Vec3d emissionColor) {
    this.radius = radius;
    this.position = position;
    this.material = material;
    this.surfaceColor = surfaceColor;
    this.emissionColor = emissionColor;
  }

  double intersect(Ray r) {
    // Ray-sphere intersection, see <https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection>.
    Vec3d op = this.position.minus(r.origin);
    final double epsilon = 1e-4;
    final double b = op.dot(r.direction);
    final double det = b * b - op.dot(op) + this.radius * this.radius;
    if (det < 0.0) {
      return NO_HIT; // behind
    }
    else {
      double sqrtDet = Math.sqrt(det);
      if (b - sqrtDet > epsilon) {
        return b - sqrtDet;
      }
      else if (b + sqrtDet > epsilon) {
        return b + sqrtDet;
      }
      else {
        return NO_HIT; // just glanced the sphere
      }
    }
  }
}

static double clamp(double x) {
  if (x < 0.0) {
    return 0.0;
  }
  else if (x > 1.0) {
    return 1.0;
  }
  else {
    return x;
  }
}

static class Intersection {
  Sphere sphere;
  double distance;
  Intersection(Sphere sphere, double distance) {
    this.sphere = sphere;
    this.distance = distance;
  }
}

static Intersection intersect(Sphere[] scene, Ray r) {
  Sphere bestSphere = null;
  double bestDist = 1e20;
  for (Sphere s : scene) {
    double dist = s.intersect(r);
    if (dist != NO_HIT && dist < bestDist) {
      bestSphere = s;
      bestDist = dist;
    }
  }

  if (bestSphere != null) {
    return new Intersection(bestSphere, bestDist);
  }
  return null;
}

static Vec3d radiance(Sphere[] scene, Random rand, Ray r, int recursionDepth) {
  Intersection isect = intersect(scene, r);
  if (isect == null || recursionDepth > 200) {
    return Vec3d.ZERO;
  }

  Sphere s = isect.sphere;
  Vec3d x = r.atDistance(isect.distance);
  Vec3d n = x.minus(s.position).normalized();
  Vec3d n1 = n.dot(r.direction) < 0.0 ? n : n.mult(-1.0);
  Vec3d f = s.surfaceColor;

  if (recursionDepth >= 5) {
    // Russian roulette.
    double p = Math.max(f.x, Math.max(f.y, f.z)); // max color component
    if (rand.nextFloat() < p) {
      // Ray survives.
      f = f.div(p);
    }
    else {
      // Ray dies.
      return s.emissionColor;
    }
  }

  // Compute ray transmission based on material type.
  if (s.material == Material.DIFFUSE) {
    // Ideal diffuser. Sample a random direction in the hemisphere.
    double r1 = TWO_PI * rand.nextFloat();
    double r2 = rand.nextFloat();
    double r2s = Math.sqrt(r2);
    Vec3d w = n1;
    Vec3d axis = Math.abs(w.x) > 0.1 ? new Vec3d(0.0, 1.0, 0.0) : new Vec3d(1.0, 0.0, 0.0);
    Vec3d u = axis.cross(w).normalized();
    Vec3d v = w.cross(u);
    Vec3d d =
      u.mult(Math.cos(r1) * r2s)
      .plus(v.mult(Math.sin(r1) * r2s))
      .plus(w.mult(Math.sqrt(1.0 - r2))).normalized();
    return s.emissionColor.plus(f.mult(radiance(scene, rand, new Ray(x, d), recursionDepth + 1)));
  }
  else if (s.material == Material.SPECULAR) {
    // Ideal specular (mirror) reflection.
    Vec3d reflDir = r.direction.minus(n.mult(2.0 * n.dot(r.direction)));
    return s.emissionColor.plus(f.mult(radiance(scene, rand, new Ray(x, reflDir), recursionDepth + 1))); 
  }
  else { // Material.REFRACTIVE
    // Ideal dielectric (non-metallic) refraction.
    Vec3d reflDir = r.direction.minus(n.mult(2.0 * n.dot(r.direction)));
    Ray reflRay = new Ray(x, reflDir);
    boolean into = n.dot(n1) > 0.0; // going outside in or vice versa?
    double nc = 1.0;
    double nt = 1.5;
    double nnt = into ? nc/nt : nt/nc; // Snell's law
    double ddn = r.direction.dot(n1);
    double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
    if (cos2t < 0.0) {
      // Total internal reflection.
      return s.emissionColor.plus(f.mult(radiance(scene, rand, reflRay, recursionDepth + 1)));
    }
    else {
      // Normal refraction.
      double intoFactor = into ? 1.0 : -1.0;
      Vec3d tdir = r.direction.mult(nnt).minus(n.mult(intoFactor * (ddn * nnt + Math.sqrt(cos2t)))).normalized();
      double a = nt - nc;
      double b = nt + nc;
      double R0 = (a * a) / (b * b);
      double c = 1.0 - (into ? -ddn : tdir.dot(n));
      double Re = R0 + (1.0 - R0) * c * c * c * c * c;
      double Tr = 1.0 - Re;
      double P = 0.25 + 0.5 * Re;
      double RP = Re / P;
      double TP = Tr / (1.0 - P);
      
      if (recursionDepth >= 2) {
        // Russian roulette for refraction.
        if (rand.nextFloat() < P) {
          return s.emissionColor.plus(f.mult(radiance(scene, rand, reflRay, recursionDepth + 1).mult(RP)));
        }
        else {
          return s.emissionColor.plus(f.mult(radiance(scene, rand, new Ray(x, tdir), recursionDepth + 1).mult(TP)));
        }
      }
      else {
        return s.emissionColor.plus(f.mult(
          radiance(scene, rand, reflRay, recursionDepth + 1).mult(Re)
          .plus(radiance(scene, rand, new Ray(x, tdir), recursionDepth + 1).mult(Tr))));
      }
    }
  }
}

Object globalLock = new Object();
Sphere[] scene;
PImage img;

Vec3d masterData[];

boolean threadResetFlags[] = new boolean[8];

void render1() {
  render(1, 8);
}
void render2() {
  render(2, 8);
}
void render3() {
  render(3, 8);
}
void render4() {
  render(4, 8);
}
void render5() {
  render(5, 8);
}
void render6() {
  render(6, 8);
}
void render7() {
  render(7, 8);
}
void render8() {
  render(8, 8);
}

void render(int threadId, int totalThreads) {
  assert(width * height == masterData.length);

  Random rand = new Random();
  Ray cam = new Ray(new Vec3d(50.0, 52.0, 295.6), new Vec3d(0.0, -0.042612, -1.0).normalized());
  Vec3d cx = new Vec3d(width * 0.5135 / height, 0.0, 0.0);
  Vec3d cy = cx.cross(cam.direction).normalized().mult(0.5135);
  Vec3d data[] = new Vec3d[masterData.length];
  for (int i = 0; i < masterData.length; ++i) {
    data[i] = new Vec3d(0.0);
  }

  int iteration = 0;
  while (true) {
    iteration++;
    for (int i = 0; i < masterData.length; ++i) {
      data[i].clear();
    }

    for (int y = threadId - 1; y < height; y += totalThreads) {
      if (threadResetFlags[threadId - 1]) {
        break;
      }

      for (int x = 0; x < width; ++x) {
        int i = (height - y - 1) * width + x;
  
        // Subpixel cols and rows (2x2).
        for (int sy = 0; sy < 2; sy++) {
          for (int sx = 0; sx < 2; sx++) {
            double r1 = 2.0 * rand.nextFloat();
            double dx = r1 < 1.0 ? Math.sqrt(r1) - 1.0 : 1.0 - Math.sqrt(2.0 - r1);
            double r2 = 2.0 * rand.nextFloat();
            double dy = r2 < 1.0 ? Math.sqrt(r2) - 1.0 : 1.0 - Math.sqrt(2.0 - r2);
            Vec3d d =
              cx.mult(((sx + 0.5 + dx) / 2.0 + x) / width - 0.5)
              .plus(cy.mult(((sy + 0.5 + dy) / 2.0 + y) / height - 0.5))
              .plus(cam.direction);
            Vec3d r = radiance(scene, rand, new Ray(cam.origin.plus(d.mult(140)), d.normalized()), 0);
            data[i].accum(r.mult(0.25)); // remember to divide by 4
          }
        }
      }
    }

    // Update once per iteration.
    synchronized (globalLock) {
      if (threadResetFlags[threadId - 1]) {
        threadResetFlags[threadId - 1] = false;
        iteration = 0;
        for (int y = threadId - 1; y < height; y += totalThreads) {
          for (int x = 0; x < width; ++x) {
            int i = (height - y - 1) * width + x;
            masterData[i].clear();
          }
        }
        continue;
      }

      print("Thread " + threadId + "; Iteration " + iteration + "\n");
      for (int y = threadId - 1; y < height; y += totalThreads) {
        for (int x = 0; x < width; ++x) {
          int i = (height - y - 1) * width + x;
          masterData[i].accum(data[i]);
          Vec3d rgb = masterData[i].div(iteration).gammaCorrected();
          img.pixels[i] = color((float) rgb.x, (float) rgb.y, (float) rgb.z);
        }
      }
      img.updatePixels();
    }
  }
}

void setup() {  
  size(800, 600);
  colorMode(RGB, 1.0);

  scene = new Sphere[] {
    //         radius position                        material             surfaceColor                 emissionColor
    new Sphere(1e5,   new Vec3d(1e5+1, 40.8, 81.6),   Material.DIFFUSE,    new Vec3d(0.75, 0.25, 0.25), new Vec3d(0.0)), // left
    new Sphere(1e5,   new Vec3d(-1e5+99, 40.8, 81.6), Material.DIFFUSE,    new Vec3d(0.25, 0.75, 0.25), new Vec3d(0.0)), // right
    new Sphere(1e5,   new Vec3d(50, 40.8, 1e5),       Material.DIFFUSE,    new Vec3d(0.75, 0.75, 0.75), new Vec3d(0.0)), // back
    new Sphere(1e5,   new Vec3d(50, 40.8, -1e5+170),  Material.DIFFUSE,    new Vec3d(0.0),              new Vec3d(0.0)), // front
    new Sphere(1e5,   new Vec3d(50, 1e5, 81.6),       Material.DIFFUSE,    new Vec3d(0.75, 0.75, 0.75), new Vec3d(0.0)), // bottom
    new Sphere(1e5,   new Vec3d(50, -1e5+81.6, 81.6), Material.DIFFUSE,    new Vec3d(0.75, 0.75, 0.75), new Vec3d(0.0)), // top
    new Sphere(16.5,  new Vec3d(27, 16.5, 47),        Material.SPECULAR,   new Vec3d(0.999),            new Vec3d(0.0)),
    new Sphere(16.5,  new Vec3d(73, 16.5, 78),        Material.REFRACTIVE, new Vec3d(0.999),            new Vec3d(0.0)),
    new Sphere(600,   new Vec3d(50, 681.6-.27, 81.6), Material.DIFFUSE,    new Vec3d(0.0),              new Vec3d(12.0)),
  };

  img = createImage(width, height, RGB);
  masterData = new Vec3d[width * height];
  for (int i = 0; i < masterData.length; ++i) {
    masterData[i] = new Vec3d(0.0);
  }

  thread("render1");
  thread("render2");
  thread("render3");
  thread("render4");
  thread("render5");
  thread("render6");
  thread("render7");
  thread("render8");
}

void draw() {
  image(img, 0, 0);
}

void mouseClicked() {
  Sphere sphereToIntersect;
  Sphere sphereToMove;
  if (mouseButton == LEFT) {
    sphereToIntersect = scene[4];
    sphereToMove = scene[6];
  }
  else if (mouseButton == RIGHT) {
    sphereToIntersect = scene[4];
    sphereToMove = scene[7];
  }
  else if (mouseButton == CENTER) {
    sphereToIntersect = scene[5];
    sphereToMove = scene[8];
  }
  else {
    return;
  }

  Ray cam = new Ray(new Vec3d(50.0, 52.0, 295.6), new Vec3d(0.0, -0.042612, -1.0).normalized());
  Vec3d cx = new Vec3d(width * 0.5135 / height, 0.0, 0.0);
  Vec3d cy = cx.cross(cam.direction).normalized().mult(0.5135);
  Vec3d d =
    cx.mult((float) mouseX / width - 0.5)
    .plus(cy.mult((float) (height - mouseY - 1) / height - 0.5))
    .plus(cam.direction);

  // Intersect with bottom surface.
  Ray r = new Ray(cam.origin.plus(d.mult(140)), d.normalized());
  double hitDist = sphereToIntersect.intersect(r);
  if (hitDist != NO_HIT) {
    Vec3d hitPoint = r.atDistance(hitDist);
    synchronized (globalLock) {
      print("Resetting render...\n");

      sphereToMove.position.x = hitPoint.x;
      sphereToMove.position.z = hitPoint.z;

      for (int i = 0; i < threadResetFlags.length; ++i) {
        threadResetFlags[i] = true;
      }
    }
  }
}
