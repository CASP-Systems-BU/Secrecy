#include "../../include/core/sharing.h"

// Sets random seed
void init_sharing() {
  if (sodium_init() == -1) {
    printf("Cannot initialize sodium library.\n");
    exit(-1);
  }
}

void arithmetic_to_boolean(AShare *x, BShare *s1, BShare *s2, BShare *s3) {
  generate_bool_share(*x, s1, s2, s3);
}

// Generate boolean shares
void generate_bool_share(Data x, BShare *x1, BShare *x2, BShare *x3) {
  // Generate two random numbers
  randombytes_buf(x1, sizeof(BShare));
  randombytes_buf(x2, sizeof(BShare));
  *x3 = x ^ (*x1 ^ *x2);
}

// Generate arithmetic shares
void generate_int_share(Data x, AShare *x1, AShare *x2, AShare *x3) {
  // Generate two random numbers
  randombytes_buf(x1, sizeof(AShare));
  randombytes_buf(x2, sizeof(AShare));
  *x3 = x - (*x1 + *x2);
}

// Generate a random seed
Seed generate_random_seed() {
  Seed s;
  randombytes_buf(&s, sizeof(s));
  return s;
}

// Generate boolean share tables
// Assumes preallocated BShareTables
void generate_bool_share_tables(Table *data, BShareTable *shares1,
                                BShareTable *shares2, BShareTable *shares3) {
  for (int i=0;i<data->numRows;i++) {
    int step = 0;
    for (int j=0;j<data->numCols;j++) {
       Data element = data->content[i][j];
       BShare s1, s2, s3;
       generate_bool_share(element, &s1, &s2, &s3);
       // Populate first table
       shares1->content[i][j + step] = s1;
       shares1->content[i][j + step + 1] = s2;
       // Populate second table
       shares2->content[i][j + step] = s2;
       shares2->content[i][j + step + 1] = s3;
       // Populate third table
       shares3->content[i][j + step] = s3;
       shares3->content[i][j + step + 1] = s1;
       step++;
    }
  }
}

// Generate arithmetic share tables
// Assumes preallocated AShareTables
void generate_int_share_tables(Table *data, AShareTable *shares1,
                               AShareTable *shares2, AShareTable *shares3) {
  for (int i=0;i<data->numRows;i++) {
    for (int j=0;j<data->numCols;j++) {
       Data element = data->content[i][j];
       AShare s1, s2, s3;
       generate_int_share(element, &s1, &s2, &s3);
       // Populate first table
       shares1->content[i][2 * j] = s1;
       shares1->content[i][2 * j + 1] = s2;
       // Populate second table
       shares2->content[i][2 * j] = s2;
       shares2->content[i][2 * j + 1] = s3;
       // Populate third table
       shares3->content[i][2 * j] = s3;
       shares3->content[i][2 * j + 1] = s1;
    }
  }
}

// Generate len random bits and populate arrays rb, ra with their
// arithmetic and binary shares.
void generate_rand_bit_shares(BShare *rb1, BShare *ra1, 
                              BShare *rb2, BShare *ra2,
                              BShare *rb3, BShare *ra3, int len) {

  // Generate len random bits
  unsigned int randombit[len];
  randombytes_buf(&randombit[0], sizeof(unsigned int) * len);

  // Generate corresponding AShares, BShares
  for (int i=0; i<len; i++) {
    generate_bool_share(randombit[i]&1, &rb1[i], &rb2[i], &rb3[i]);
    generate_int_share(randombit[i]&1, &ra1[i], &ra2[i], &ra3[i]);
  }
}
