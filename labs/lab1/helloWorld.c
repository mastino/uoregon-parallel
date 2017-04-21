  #include <omp.h>
#include <time.h>
#define NUM_ITERS 20

int main() {
  int i;
  #pragma omp parallel for
  for(i=0;i<NUM_ITERS;i++) {
    int thread_num=omp_get_thread_num();
    int num_threads=omp_get_num_threads();
    printf("thread %d of %d\n",thread_num,num_threads);
    sleep(1);
  }
}
