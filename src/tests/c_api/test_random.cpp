#include <stdio.h>
#include <stdlib.h>
#include <sys/random.h>

int main(void) {

  printf("Size of random: %zu\n", sizeof(random()));
  printf("Size of long long: %zu\n", sizeof(long long));
  printf("Size of long: %zu\n", sizeof(long));
  printf("Size of int: %zu\n", sizeof(int));
  printf("Size of unsigned int (seed): %zu\n", sizeof(unsigned int));
  
  return 0;
}