// This file contains various functions to help manage the code.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<time.h>
#include<sys/time.h>
#include<inttypes.h>
#include "../include/utility.h"

// Function to check if string is empty:
int isEmpty(const char *s) {
  while (*s != '\0') {
    if (!isspace((unsigned char)*s))
      return 0;
    s++;
  }
  return 1;
}

// Get current time
uint64_t get_time()
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec * 1e9 + t.tv_nsec;
}

// Generate a uniformly distributed random number:
double randU(double low, double high)
{
	return (high - low) * drand48() + low;
}
