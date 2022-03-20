#pragma once

#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <assert.h>
#include <complex>
#include <ctime>
#include <experimental/filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi/mpi.h>
#include <random>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <vector>
// using namespace concurrency;
#define pi2 3.141592653589793 * 2
using namespace std;

extern bool energy_window;
extern int Rmax; // defined in compute_matrix_energy_window

extern bool detector_only;
extern bool Detector_Continue_Simulation;

extern double *coupling_strength_to_mode0;
extern double Coupling_between_electronic_state;
extern bool load_sampling_state;

// define function here
float ran2(long &idum);
void convert_dv(const vector<vector<int>> &vec_2d, vector<int> &vec_1d,
                vector<int> &displacement, vector<int> &element_size);
// used for cnostruct buffer for communication between process for matrix
// multiplication.
int construct_send_buffer_index(int *remoteVecCount, int *remoteVecPtr,
                                int *remoteVecIndex,
                                int *tosendVecCount_element,
                                int *tosendVecPtr_element,
                                int *&tosendVecIndex_ptr);

int compar(const void *a, const void *b);

int find_position_for_insert_binary(
    const vector<vector<int>> &vmode, const vector<int> &ndetector,
    bool &exist); // write in compute_matrix_energy_window.cpp

#endif // UTIL_H
