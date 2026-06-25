
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Include ONLY your custom v8 Math Engine
#include "linear_algebra.h"
#include "fast_math.h"

// Use your custom Quake III Fast Inverse Sqrt to normalize vectors!
inline vec3 fast_normalize(vec3 v) {
    double mag2 = vec3_dot(v, v);
    if (mag2 < 1e-8) return (vec3){0, 0, 0};
    return vec3_scale(v, ml_fast_rsqrt(mag2));
}

typedef struct {
    vec3 center;
    double radius;
    vec3 color;
} Sphere;

// Ray-Sphere Intersection using your fast math
double intersect_sphere(vec3 ro, vec3 rd, Sphere s) {
    vec3 oc = vec3_sub(ro, s.center);
    double b = vec3_dot(oc, rd);
    double c = vec3_dot(oc, oc) - s.radius * s.radius;
    double h = b * b - c;
    if (h < 0.0) return -1.0;

    // Using ml_fast_rsqrt instead of standard sqrt!
    double sq = 1.0 / ml_fast_rsqrt(h);
    double t = -b - sq;
    if (t < 0.001) t = -b + sq;
    return t > 0.001 ? t : -1.0;
}

// Ray-Plane Intersection (for the floor)
double intersect_plane(vec3 ro, vec3 rd, double y) {
    if (fabs(rd.y) < 1e-6) return -1.0;
    double t = (y - ro.y) / rd.y;
    return t > 0.001 ? t : -1.0;
}

int main() {
    int W = 512, H = 384;
    FILE *f = fopen("render.ppm", "w");
    fprintf(f, "P3\n%d %d\n255\n", W, H);

    // Scene Setup
    Sphere spheres[3] = {
        {(vec3){-1.5, 0.0, -5.0}, 1.0, (vec3){0.9, 0.2, 0.2}}, // Red
        {(vec3){ 1.5, 0.0, -5.0}, 1.0, (vec3){0.2, 0.9, 0.2}}, // Green
        {(vec3){ 0.0, 1.5, -8.0}, 1.5, (vec3){0.2, 0.2, 0.9}}  // Blue
    };

    vec3 light_dir = fast_normalize((vec3){-1.0, 2.0, -3.0});
    vec3 cam_pos = {0.0, 1.0, 0.0};

    printf("Rendering %dx%d image using v8 math engine...\n", W, H);

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            // Generate Ray
            double u = (x - W / 2.0) / (W / 2.0);
            double v = (H / 2.0 - y) / (W / 2.0); // Aspect ratio correction
            vec3 rd = fast_normalize((vec3){u, v, -1.5});
            vec3 ro = cam_pos;

            double min_t = 1e9;
            vec3 hit_color = {0.15, 0.15, 0.25}; // Background sky
            vec3 hit_normal = {0, 0, 0};
            int hit_plane = 0;
            vec3 hit_pos;

            // 1. Check Spheres
            for (int i = 0; i < 3; i++) {
                double t = intersect_sphere(ro, rd, spheres[i]);
                if (t > 0 && t < min_t) {
                    min_t = t;
                    hit_color = spheres[i].color;
                    hit_pos = vec3_add(ro, vec3_scale(rd, t));
                    hit_normal = fast_normalize(vec3_sub(hit_pos, spheres[i].center));
                    hit_plane = 0;
                }
            }

            // 2. Check Floor
            double tp = intersect_plane(ro, rd, -1.0);
            if (tp > 0 && tp < min_t) {
                min_t = tp;
                hit_pos = vec3_add(ro, vec3_scale(rd, tp));
                hit_normal = (vec3){0, 1, 0};
                hit_plane = 1;

                // Checkerboard pattern
                int cx = (int)floor(hit_pos.x + 1000.0) % 2;
                int cz = (int)floor(hit_pos.z + 1000.0) % 2;
                hit_color = (cx ^ cz) ? (vec3){0.8, 0.8, 0.8} : (vec3){0.2, 0.2, 0.2};
            }

            // 3. Phong Lighting
            if (min_t < 1e8) {
                double diff = vec3_dot(hit_normal, light_dir);
                if (diff < 0) diff = 0;
                double amb = 0.2;
                double intensity = amb + diff * 0.8;
                hit_color = vec3_scale(hit_color, intensity);
            }

            // Output pixel
            int r = (int)(hit_color.x * 255); if(r>255)r=255; if(r<0)r=0;
            int g = (int)(hit_color.y * 255); if(g>255)g=255; if(g<0)g=0;
            int b = (int)(hit_color.z * 255); if(b>255)b=255; if(b<0)b=0;
            fprintf(f, "%d %d %d ", r, g, b);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    printf("✅ Render complete! Saved to render.ppm\n");
    return 0;
}
