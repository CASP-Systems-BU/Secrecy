#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 0

void test_bool_sharing(Data x);
void test_int_sharing(Data x);
void test_bool_share_tables(Table *t, BShareTable *s1, BShareTable *s2,
                            BShareTable *s3);
 void test_int_share_tables(Table *t, AShareTable *s1, AShareTable *s2,
                            AShareTable *s3);

int main(void) {

  init_sharing();

  #if DEBUG
  printf("Testing boolean sharing.\n");
  #endif

  // Test boolean sharing
  for (int i=0; i<100; i++) {
    Data x = random();
    test_bool_sharing(x);
  }

  #if DEBUG
  printf("Testing boolean share table generation.\n");
  #endif

  // Test BShareTable generation
  Table data;
  generate_random_table(&data, 10, 10);
  BShareTable t1 = {0,0,10,20};
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {0,1,10,20};
  allocate_bool_shares_table(&t2);
  BShareTable t3 = {0,2,10,20};
  allocate_bool_shares_table(&t3);
  // Generate boolean shares
  generate_bool_share_tables(&data, &t1, &t2, &t3);
  // Test sharing
  test_bool_share_tables(&data, &t1, &t2, &t3);

  #if DEBUG
  printf("Testing arithmetic sharing.\n");
  #endif

  // Test arithmetic sharing
  for (int i=0; i<100; i++) {
    Data x = random();
    test_int_sharing(x);
  }

  #if DEBUG
  printf("Testing arithmetic share table generation.\n");
  #endif

  // Test AShareTable generation
  AShareTable s1 = {0,0,10,20};
  allocate_int_shares_table(&s1);
  AShareTable s2 = {0,1,10,20};
  allocate_int_shares_table(&s2);
  AShareTable s3 = {0,2,10,20};
  allocate_int_shares_table(&s3);
  // Generate arithmetic shares
  generate_int_share_tables(&data, &s1, &s2, &s3);
  // Test sharing
  test_int_share_tables(&data, &s1, &s2, &s3);

  printf("TEST SHARING: OK.\n");
}

void test_bool_sharing(Data x) {
  BShare s1, s2, s3;
  generate_bool_share(x,&s1,&s2,&s3);
  assert((s1 ^ s2 ^ s3) == x);
}

void test_bool_share_tables(Table *data, BShareTable* t1, BShareTable* t2,
                                                          BShareTable *t3) {
  assert(data->numRows == t1->numRows);
  assert(data->numRows == t2->numRows);
  assert(data->numRows == t3->numRows);
  assert(2*data->numCols == t1->numCols);
  assert(2*data->numCols == t2->numCols);
  assert(2*data->numCols == t3->numCols);

  for (int i=0; i<data->numRows; i++) {
    for (int j=0; j<data->numCols; j++) {
      assert(t1->content[i][2 * j] == t3->content[i][2 * j + 1]);
      assert(t1->content[i][2 * j + 1] == t2->content[i][2 * j]);
      assert(t2->content[i][2 * j + 1] == t3->content[i][2 * j]);
      assert(data->content[i][j] == (t1->content[i][2 * j] ^
                                     t2->content[i][2 * j] ^
                                     t3->content[i][2 * j])
            );
    }
  }
}

void test_int_share_tables(Table *data, AShareTable* t1, AShareTable* t2,
                                                         AShareTable *t3) {
  assert(data->numRows == t1->numRows);
  assert(data->numRows == t2->numRows);
  assert(data->numRows == t3->numRows);
  assert(2*data->numCols == t1->numCols);
  assert(2*data->numCols == t2->numCols);
  assert(2*data->numCols == t3->numCols);

  for (int i=0; i<data->numRows; i++) {
    for (int j=0; j<data->numCols; j++) {
      assert(t1->content[i][2 * j] == t3->content[i][2 * j + 1]);
      assert(t1->content[i][2 * j + 1] == t2->content[i][2 * j]);
      assert(t2->content[i][2 * j + 1] == t3->content[i][2 * j]);
      assert(data->content[i][j] == (t1->content[i][2 * j] +
                                     t2->content[i][2 * j] +
                                     t3->content[i][2 * j])
            );
    }
  }
}

void test_int_sharing(Data x) {
  AShare s1, s2, s3;
  generate_int_share(x,&s1,&s2,&s3);
  assert((s1 + s2 + s3) == x);
}
