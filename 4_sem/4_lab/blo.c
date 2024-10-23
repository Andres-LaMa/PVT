#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

const float G = 6.67e-11;
omp_lock_t *locks;

typedef struct Particle
{
    float x;
    float y;
    float z;
}particle;

double wtime(){
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1E-9;
}

void initial(particle *p, particle *f, particle *v, float *m, int n){
    for (int i = 0; i < n; i++){
        p[i].x = rand() / (float)RAND_MAX-0.5;
        p[i].y = rand() / (float)RAND_MAX-0.5;
        p[i].z = rand() / (float)RAND_MAX-0.5;

        v[i].x = rand() / (float)RAND_MAX-0.5;
        v[i].y = rand() / (float)RAND_MAX-0.5;
        v[i].z = rand() / (float)RAND_MAX-0.5;

        m[i] = rand() / (float)RAND_MAX*10+0.01;
        f[i].x = f[i].y = f[i].z = 0;
    }
}

void s_calc_forces(particle *p, particle *f, float *m, int n){
    for (int i = 0; i < n-1; i++){
        for (int j = i+1; j < n; j++){
            float dist = sqrtf(powf(p[i].x - p[j].x, 2)+powf(p[i].y - p[j].y, 2)+powf(p[i].z - p[j].z, 2));
            float mag = (G * m[i]*m[j])/powf(dist,2);
            particle dir = {
                .x = p[j].x - p[i].x,
                .y = p[j].y - p[i].y,
                .z = p[j].z - p[i].z,
            };

            f[i].x += mag*dir.x / dist;
            f[i].y += mag*dir.y / dist;
            f[i].z += mag*dir.z / dist;

            f[j].x += mag*dir.x / dist;
            f[j].y += mag*dir.y / dist;
            f[j].z += mag*dir.z / dist;
        }
    }
}

void n_block_calc(particle *p, particle *f, float *m, int n){
    #pragma omp for schedule(dynamic, 4) nowait
    for (int i = 0; i < n-1; i++){
        for (int j = i+1; j < n; j++){
            float dist = sqrtf(powf(p[i].x-p[j].x, 2)+
                                powf(p[i].y-p[j].y, 2)+
                                powf(p[i].z-p[j].z, 2));
            float mag = (G * m[i]*m[j])/powf(dist, 2);
            particle dir = {
                .x = p[j].x - p[i].x,
                .y = p[j].y - p[i].y,
                .z = p[j].z - p[i].z
            };
            omp_set_lock(&locks[i]);
            f[i].x += mag * dir.x / dist;
            f[i].y += mag * dir.y / dist;
            f[i].z += mag * dir.z / dist;
            omp_unset_lock(&locks[i]);
            omp_set_lock(&locks[j]);
            f[j].x -= mag * dir.x / dist;
            f[j].y -= mag * dir.y / dist;
            f[j].z -= mag * dir.z / dist;
            omp_unset_lock(&locks[j]);
        }
    }
    
}

void s_move_partic(particle *p, particle *f, particle *v, float *m, int n, double dt){
    for (int i = 0; i < n; i++){
        particle dv = {
            .x = f[i].x / m[i] * dt,
            .y = f[i].y / m[i] * dt,
            .z = f[i].z / m[i] * dt,
        };
        particle dp = {
            .x = (v[i].x + dv.x/2)*dt,
            .y = (v[i].y + dv.y/2)*dt,
            .z = (v[i].z + dv.z/2)*dt,
        };
        v[i].x += dv.x;
        v[i].y += dv.y;
        v[i].z += dv.z;

        p[i].x += dp.x;
        p[i].y += dp.y;
        p[i].z += dp.z;

        f[i].x = f[i].y = f[i].z = 0;
    }
}

void serial_v(particle *p, particle *f, particle *v, float *m, int n, double *time){
    double ttotal, tinit = 0, tforces = 0, tmove = 0;
    double dt = 1e-5;
    ttotal = omp_get_wtime();
    for (double t = 0; t <= 1; t+= dt){
        tforces -= omp_get_wtime();
        s_calc_forces(p, f, m, n);
        tforces += omp_get_wtime();
        
        tmove -= omp_get_wtime();
        s_move_partic(p, f, v, m, n, dt);
        tmove += omp_get_wtime();
    }
    ttotal = omp_get_wtime() - ttotal;
    *time = ttotal;
    printf("# NBody (n=%d)\n", n);
    printf("# Serial time (sec): ttotal %.6f, tinit %.6f, tforces %.6f, tmove %.6f\n", ttotal, tinit, tforces, tmove);
}

void n_block_v(particle *p, particle *f, particle *v, float *m, int n, double *time){
    double ttotal, tinit = 0, tforces = 0, tmove = 0;
    double dt = 1e-5;
    ttotal = omp_get_wtime();
    locks = malloc(sizeof(omp_lock_t)*n);
    for (int i = 0; i < n; i++){
        omp_init_lock(&locks[i]);
    }
    dt = 1e-5;
    for (double t = 0; t <= 1; t+= dt){
        tforces -= omp_get_wtime();
        n_block_calc(p, f, m, n);
        tforces += omp_get_wtime();
        #pragma omp barrier
        tmove -= omp_get_wtime();
        s_move_partic(p, f, v, m, n, dt);
        tmove += omp_get_wtime();
        #pragma omp barrier
    }
    ttotal = omp_get_wtime() - ttotal;
    *time = ttotal;
    printf("# NBody (n=%d)\n", n);
    printf("# N-block time (sec): ttotal %.6f, tinit %.6f, tforces %.6f, tmove %.6f\n", ttotal, tinit, tforces, tmove);
    free(locks);
}

int main (int argc, char *argv[]) {
    double tinit = 0;
    int n = (argc > 1) ? atoi(argv[1]) : 10;
    char *filename = "n_block.csv";

    tinit = -omp_get_wtime();
    particle *p = malloc(sizeof(*p)*n);
    particle *f = malloc(sizeof(*f)*n);
    particle *v = malloc(sizeof(*v)*n);
    float *m =  malloc(sizeof(*m)*n);

    initial(p, f, v, m, n);
    tinit += omp_get_wtime();
    

    double serial_time;
    // Serial version
    serial_v(p, f, v, m, n, &serial_time);

    FILE *fout = fopen(filename, "w");

    if(!fout){
        fprintf(stderr, "Can`t save file\n");
        exit(EXIT_FAILURE);
    }

    int threads[] = {2, 4, 6, 8};
    for (int i = 0; i < 4; i++) {
        omp_set_num_threads(threads[i]);

        double parallel_time;
        printf("%d: \n", threads[i]);
        //n_block version
        n_block_v(p, f, v, m, n, &parallel_time);

        fprintf(fout, "%d; %.6lf\n", threads[i], serial_time / parallel_time);
    }

    fclose(fout);

    free(m);
    free(v);
    free(f);
    free(p);

    return 0;
}