#include "../../include/core/party.h"

#include <sodium.h>
#include "mpi.h"

#include "../../include/core/comm.h"

#define PRIVATE static

#define RAND_STATE_SIZE 128
#define MSG_TAG 42

Seed seed_local = 0, seed_remote = 0;
char local_state[RAND_STATE_SIZE], remote_state[RAND_STATE_SIZE];

PRIVATE void check_init_seeds(const char*);

int exchange_rsz_seeds(int succ_rank, int pred_rank) {

  if (succ_rank == -1){
    succ_rank = get_succ();
  }

  if(pred_rank == -1){
    pred_rank = get_pred();
  }

  // initialize random number generator
  if (sodium_init() == -1) {
        return 1;
  }

  // generate local seed
  randombytes_buf(&seed_local, sizeof(seed_local));

  // send seed to successor
  MPI_Send(&seed_local, 1, MPI_LONG_LONG, succ_rank, MSG_TAG, MPI_COMM_WORLD);

  // receive remote seed
  MPI_Recv(&seed_remote, 1, MPI_LONG_LONG, pred_rank, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // init generator states
  initstate(seed_local, local_state, RAND_STATE_SIZE);
  initstate(seed_remote, remote_state, RAND_STATE_SIZE);

  return 0;
}

WSharePair get_next_w() {
    printf("NOT IMPLEMENTED: %s\n", __func__);
    exit(EXIT_FAILURE);
  /*  WSharePair wpair;
    wpair.first = 0;
    wpair.second = 0;
    return wpair; */
}

// Generate one random arithmetic share
AShare get_next_r() {
  check_init_seeds(__func__);
  setstate(local_state);
  long long r_local = random();
  setstate(remote_state);
  return r_local - random();
}

// Generate one random binary share
BShare get_next_rb() {
  check_init_seeds(__func__);

  setstate(local_state);
  long long r_local = random();
  setstate(remote_state);
  return r_local ^ random();
}

// Generate pairs of random shares
// one local and one remote
void get_next_rb_pair_array(BShare *r1, BShare *r2, int len) {
  check_init_seeds(__func__);

  // Generate len random shares using the local seed
  setstate(local_state);
  for (int i=0; i<len; i++) {
    r1[i] = random();
  }

  // Generate len random shares using the remote seed
  setstate(remote_state);
  for (int i=0; i<len; i++) {
    r2[i] = random();
  }
}

// Generate an array of random binary shares
void get_next_rb_array(BShare *rnum, int len) {
  check_init_seeds(__func__);

  // Generate len random shares using the local seed
  setstate(local_state);
  for (int i=0; i<len; i++) {
    rnum[i] = random();
  }

  // Generate len random shares using the remote seed
  // and xor them with the corresponding local ones
  setstate(remote_state);
  for (int i=0; i<len; i++) {
    rnum[i] ^= random();
  }
}

// Generate an array of random arithmetic shares
void get_next_array(BShare *rnum, int len) {
  check_init_seeds(__func__);
  setstate(local_state);
  for (int i=0; i<len; i++) {
    rnum[i] = random();
  }
  setstate(remote_state);
  for (int i=0; i<len; i++) {
    rnum[i] -= random();
  }
}

// check if seeds have been initialized
PRIVATE void check_init_seeds(const char* f) {
    if (seed_local == 0 || seed_remote == 0) {
        fprintf(stderr, "ERROR: exchange_rsz_seeds() must be called before %s\n", f);
        exit(EXIT_FAILURE);
    }
}