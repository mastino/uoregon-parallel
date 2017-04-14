#include <omp.h>
#include <time.h>
#define NUM_ITERS 10

int main() {
  int i;
  #pragma omp parallel for
  for(i=0;i<NUM_ITERS;i++) {
    int thread_num=omp_get_thread_num();
    printf("Alien #%d says hello, they're  #%d in line.\n",thread_num,i);
    sleep(1);
  }
}
