#ifndef SECRECY_CONSTANTS_H
#define SECRECY_CONSTANTS_H

struct Protocol {
    int NUM_PARTIES;
    int PARTY_REPLICATION;
};

static Protocol R_3PC{
        3,
        2
};

// Current protocol
static Protocol CR_P = R_3PC;

#endif //SECRECY_CONSTANTS_H
