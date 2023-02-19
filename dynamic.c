#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define WIDTH 800
#define HEIGHT 800
#define MAX_ITER 10000
int mandelbrot(double x, double y)
{
    double zx = 0.0, zy = 0.0;
    int i;
    for (i = 0; i < MAX_ITER; i++)
    {
        double zx_new = zx * zx - zy * zy + x;
        double zy_new = 2 * zx * zy + y;
        zx = zx_new;
        zy = zy_new;
        if (zx * zx + zy * zy > 4.0)
        {
            return i;
        }
    }
    return MAX_ITER;
}
int main(int argc, char **argv)
{   int rank, size;

    clock_t start, end;
    double time_used;
    start = clock();


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    int rowStart = rank * HEIGHT / size;
    int rowEnd = (rank + 1) * HEIGHT / size;
    int *buffer = (int *)malloc(WIDTH * (rowEnd - rowStart) * sizeof(int));


    int i, j, k;
    for (i = rowStart, k = 0; i < rowEnd; i++, k += WIDTH)
    {
        for (j = 0; j < WIDTH; j++)
        {
            double x = ((double)j - WIDTH / 2) / (WIDTH / 4.0);
            double y = ((double)i - HEIGHT / 2) / (HEIGHT / 4.0);
            buffer[k + j] = mandelbrot(x, y);
        }
    }



    int *counts = (int *)malloc(size * sizeof(int));
    int *displs = (int *)malloc(size * sizeof(int));
    for (i = 0; i < size; i++)
    {
        counts[i] = WIDTH * (HEIGHT / size);
        displs[i] = i * counts[i];
    }

    int *recvbuf = (int *)malloc(WIDTH * HEIGHT * sizeof(int));
    MPI_Gatherv(buffer, counts[rank], MPI_INT, recvbuf, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);


    if (rank == 0)
    {
        FILE *fp = fopen("dynamic2.ppm", "w");
        fprintf(fp, "P3\n%d %d\n255\n", WIDTH, HEIGHT);
        int i, j;
        for (i = 0; i < HEIGHT; i++)
        {
            for (j = 0; j < WIDTH; j++)
            {
                int index = i * WIDTH + j;
                int value = recvbuf[index];
                int r = value % 16;
                int g = (value / 16) % 16;
                int b = (value / 256) % 16;
                fprintf(fp, "%d %d %d ", r * 16, g * 16, b * 16);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    free(buffer);
    free(counts);
    free(displs);
    free(recvbuf);

    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    if(rank==0) {
        printf("Time taken: %f seconds\n", time_used);
    }

    MPI_Finalize();
    return 0;
}

