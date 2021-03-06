#define ARRAYSIZE 100

void init(int* a, int size){
  int i=0;
  for(i=0;i<size;i++){
    a[i]=i;
  }
}

int main() {
  int i=0;
  int* a=(int*)malloc(sizeof(int)*ARRAYSIZE);
  init(a,ARRAYSIZE);
  omp_set_num_threads(4);
  //intended behavior   : a[i] gets the sum of its value and its next (cyclic) neighbor's initial value.
  //observed behavior   : inconsistent quasi cumulative.
  //Luckily we have Computer Scientists here to fix this for us
  int a0 = a[0];
  #pragma omp parallel for
  for(i=0;i<ARRAYSIZE-1;i++){
    a[i]=a[i]+a[((i+1)%ARRAYSIZE)];
  }

  a[ARRAYSIZE-1] = a[ARRAYSIZE-1] + a0;

  for(i=0;i<ARRAYSIZE;i++){
    printf("Index: %d, Value: %d\n",i,a[i]);
  }
}
