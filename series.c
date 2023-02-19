#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define WIDTH 800
#define HEIGHT 800
#define MAX_ITER 1000

int mandelbrot(double real, double imag) {
    double z_real = 0;
    double z_imag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        double z_real_squared = z_real * z_real;
        double z_imag_squared = z_imag * z_imag;

        if (z_real_squared + z_imag_squared > 4) {
            return i;
        }

        double new_z_real = z_real_squared - z_imag_squared + real;
        double new_z_imag = 2 * z_real * z_imag + imag;

        z_real = new_z_real;
        z_imag = new_z_imag;
    }
    return 0;
}

int main() {
    clock_t start, end;
    double time_used;

    start = clock();

    FILE* fp;
    char* filename = "mandelbrot.ppm";
    fp = fopen(filename, "wb");
    fprintf(fp, "P6\n%d %d\n255\n", WIDTH, HEIGHT);

    for (int row = 0; row < HEIGHT; row++) {
        for (int col = 0; col < WIDTH; col++) {
            double real = (col - WIDTH / 2.0) * 4.0 / WIDTH;
            double imag = (row - HEIGHT / 2.0) * 4.0 / WIDTH;

            int iter = mandelbrot(real, imag);

            unsigned char r, g, b;
            if (iter == 0) {
                r = g = b = 0;
            } else {
                r = (int) (sin(iter * 0.03) * 127 + 128);
                g = (int) (sin(iter * 0.05) * 127 + 128);
                b = (int) (sin(iter * 0.07) * 127 + 128);
            }

            fputc(r, fp);
            fputc(g, fp);
            fputc(b, fp);
        }
    }

    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Time taken: %f seconds\n", time_used);

    fclose(fp);
    return 0;
}


/**#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <time.h>
#define WIDTH 800
#define HEIGHT 800
#define MAX_ITER 1000
int mandelbrot(double real, double imag) {
    double z_real = 0;
    double z_imag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        double z_real_squared = z_real * z_real;
        double z_imag_squared = z_imag * z_imag;

        if (z_real_squared + z_imag_squared > 4) {
            return i;
        }

        double new_z_real = z_real_squared - z_imag_squared + real;
        double new_z_imag = 2 * z_real * z_imag + imag;

        z_real = new_z_real;
        z_imag = new_z_imag;
    }
    return 0;}
int main() {
   // time_t start, end;
    //double time_used;
   // start = time(NULL);

    FILE* fp;
    char* filename = "mandelbrot.ppm";
    fp = fopen(filename, "wb");
    fprintf(fp, "P6\n%d %d\n255\n", WIDTH, HEIGHT);


    for (int row = 0; row < HEIGHT; row++) {
        for (int col = 0; col < WIDTH; col++) {
            double real = (col - WIDTH / 2.0) * 4.0 / WIDTH;
            double imag = (row - HEIGHT / 2.0) * 4.0 / WIDTH;

            int iter = mandelbrot(real, imag);

            unsigned char r, g, b;
            if (iter == 0) {
                r = g = b = 0;
            } else {
                r = (int) (sin(iter * 0.03) * 127 + 128);
                g = (int) (sin(iter * 0.05) * 127 + 128);
                b = (int) (sin(iter * 0.07) * 127 + 128);
            }

            fputc(r, fp);
            fputc(g, fp);
            fputc(b, fp);
        }
    }


    //end = time(NULL);
   // time_used = difftime(end, start);
    printf("Time taken: %f seconds", time_used);

    fclose(fp);
    return 0;
}

**/

/**#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
int mandelbrot_point(double complex c, int max_iterations) {
    double complex z = 0;
    for (int i = 0; i < max_iterations; i++) {
        z = z*z + c;
        if (cabs(z) > 2) {
            return i;
        }
    }
    return max_iterations;
}
void mandelbrot_set_parallel(double xmin, double xmax, double ymin, double ymax, int width, int height, int max_iterations, int *mandelbrot) {
    double dx = (xmax - xmin) / width;
    double dy = (ymax - ymin) / height;
#pragma omp parallel for collapse(2)
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            double complex c = xmin + dx*j + (ymin + dy*i)*I;
            mandelbrot[i*width+j] = mandelbrot_point(c, max_iterations);
        }
    }
}
int main() {
    double xmin = -2;
    double xmax = 2;
    double ymin = -2;
    double ymax = 2;
    int width = 1000;
    int height = 1000;
    int max_iterations = 100;
    time_t start, end;
    double time_used;
    start = time(NULL);
   // double start_time, end_time, total_time;
    int *mandelbrot = malloc(width*height*sizeof(int));
   // start_time = clock();
  //  double start_time = omp_get_wtime();


    mandelbrot_set_parallel(xmin, xmax, ymin, ymax, width, height, max_iterations, mandelbrot);
  //  double end_time = omp_get_wtime();
  //  printf("Elapsed time: %.2f seconds\n", end_time - start_time);

    FILE *fp = fopen("seriesTest.pgm", "wb");
    fprintf(fp, "P5\n%d %d\n%d\n", width, height, 255);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            fputc(mandelbrot[i*width+j], fp);
        }
    }
    fclose(fp);
    free(mandelbrot);
    end = time(NULL);
    time_used = difftime(end, start);
    printf("Time taken: %f seconds", time_used);
    return 0;
}

**/

/** testing
 *
 * #include <stdio.h>
#include <complex.h>

#define WIDTH 800
#define HEIGHT 800
#define MAX_ITER 1000

int mandelbrot(complex double c)
{
    complex double z = 0;
    int i;
    for (i = 0; i < MAX_ITER; i++)
    {
        z = z * z + c;
        if (cabs(z) > 2)
        {
            return i;
        }
    }
    return MAX_ITER;
}

int main()
{
    int data[WIDTH][HEIGHT];
    int x, y;
    for (y = 0; y < HEIGHT; y++)
    {
        for (x = 0; x < WIDTH; x++)
        {
            double real = (x - WIDTH / 2) * 4.0 / WIDTH;
            double imag = (y - HEIGHT / 2) * 4.0 / HEIGHT;
            complex double c = real + imag * I;
            data[x][y] = mandelbrot(c);
        }
    }


        FILE *fp = fopen("series.ppm", "w");
        fprintf(fp, "P3\n%d %d\n255\n", WIDTH, HEIGHT);
        int i, j;
        for (i = 0; i < HEIGHT; i++)
        {
            for (j = 0; j < WIDTH; j++)
            {
                int index = i * WIDTH + j;
                int value = data[i][index];
                int r = value % 16;
                int g = (value / 16) % 16;
                int b = (value / 256) % 16;
                fprintf(fp, "%d %d %d ", r * 16, g * 16, b * 16);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);


    for (int i = 0; i <HEIGHT ; i++) {
        for (int j = 0; j < WIDTH; j++) {
            printf("%d ", data[i][j]);
        }
        printf("\n");
    }

    return 0;
}
**/