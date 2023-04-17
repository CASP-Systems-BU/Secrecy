#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "../../../include/core/primitives.h"

#define DEBUG 0

struct shares {
  AShare s1;
  AShare s2;
  AShare s3;
};

struct shares_b {
  BShare s1;
  BShare s2;
  BShare s3;
};

// arithmetic primitives
void test_add(struct shares, struct shares, long long);
void test_mul(struct shares, struct shares, struct shares, long long);
void test_eq(struct shares, struct shares, struct shares, struct shares,
             long long);

// binary primitives
void test_and_b(struct shares_b, struct shares_b, struct shares_b, long long);
void test_ltz(struct shares_b, long long);

int main(void) {

  long long x_in = 15;
  long long y_in = 15;
  long long z_in = 16;
  long long w_in = 2;

  struct shares x = {2, 10, 3}; // x = 15
  struct shares y = {-5, 5, 15}; // y = 15
  struct shares z = {-5, 5, 16}; // z = 16
  struct shares r = {4, 4, -8}; // r = 0
  struct shares w = {7, 9, -14}; // w = 2

  struct shares_b x_b = {2, 10, 7}; // x = 15
  struct shares_b y_b = {5, 5, 15}; // y = 15
  struct shares_b r_b = {4, 3, 7}; // r = 0
  struct shares_b xn_b = {785699, -784077, 42}; // x = -3014
  struct shares_b xn2_b = {-33355, -573328, -104827}; // x = -632000

  // test addition
  test_add(x, y, x_in + y_in);

  // test multiplication
  test_mul(x, y, r, x_in * y_in);

  // test equality
  test_eq(x, y, w, r, (x_in - y_in) * w_in);

  // test inequality
  test_eq(x, z, w, r, (x_in - z_in) * w_in);

  // test logical and
  test_and_b(x_b, y_b, r_b, x_in & y_in);

  // test ltz with positive number
  test_ltz(x_b, 0);

  // test ltz with 0
  test_ltz(r_b, 0);

  // test ltz with negative number
  test_ltz(xn_b, 1);

  // test ltz with negative number
  test_ltz(xn2_b, 1);

  printf("TEST PRIMITIVES: OK.\n");
  return 0;
}

/**
 * Tests the result of adding two shared integers.
 **/
void test_add(struct shares x, struct shares y, long long result) {
  assert((add(x.s1, y.s1) + add(x.s2, y.s2) + add(x.s3, y.s3)) == result);
}

/**
 * Tests the result of multiplying two shared integers.
 **/
void test_mul(struct shares x, struct shares y, struct shares r, long long result) {

  AShare z1 = mul(x.s1, x.s2, y.s1, y.s2, r.s1);
  AShare z2 = mul(x.s2, x.s3, y.s2, y.s3, r.s2);
  AShare z3 = mul(x.s3, x.s1, y.s3, y.s1, r.s3);

  #if DEBUG
  printf("%lld\n", z1);
  printf("%lld\n", z2);
  printf("%lld\n", z3);
  printf("%lld\n", z1+z2+z3);
  #endif

  assert((z1+z2+z3) == result);
}

/**
 * Tests the equality of two shared integers.
**/
void test_eq(struct shares x, struct shares y, struct shares w, struct shares r,
             long long result) {

  AShare z1 = eq(x.s1, x.s2, y.s1, y.s2, w.s1, w.s2, r.s1);
  AShare z2 = eq(x.s2, x.s3, y.s2, y.s3, w.s2, w.s3, r.s2);
  AShare z3 = eq(x.s3, x.s1, y.s3, y.s1, w.s3, w.s1, r.s3);

  #if DEBUG
  printf("%lld\n", z1);
  printf("%lld\n", z2);
  printf("%lld\n", z3);
  printf("%lld\n", z1+z2+z3);
  #endif

  assert((z1+z2+z3) == result);
}

/**
 * Tests the result of logical bitwise-and between two shared integers.
 **/
void test_and_b(struct shares_b x, struct shares_b y, struct shares_b r,
                long long result) {

  BShare z1 = and_b(x.s1, x.s2, y.s1, y.s2, r.s1);
  BShare z2 = and_b(x.s2, x.s3, y.s2, y.s3, r.s2);
  BShare z3 = and_b(x.s3, x.s1, y.s3, y.s1, r.s3);

  #if DEBUG
  printf("%lld\n", z1);
  printf("%lld\n", z2);
  printf("%lld\n", z3);
  printf("%lld\n", z1^z2^z3);
  #endif

  assert((z1^z2^z3) == result);
}

/**
 * Tests the result of < 0.
 */
void test_ltz(struct shares_b x, long long result) {

  BShare z1 = ltz_b(x.s1);
  BShare z2 = ltz_b(x.s2);
  BShare z3 = ltz_b(x.s3);

  #if DEBUG
  printf("%lld\n", z1);
  printf("%lld\n", z2);
  printf("%lld\n", z3);
  printf("%lld\n", z1^z2^z3);
  #endif

  assert((z1^z2^z3) == result);

}
