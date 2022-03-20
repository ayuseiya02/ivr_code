//
// Created by phyzch on 4/14/20.
// Used for putting non essential function in same file.
//
#include "System.h"
#include "util.h"
using namespace std;

void detector::shift_detector_Hamiltonian(ofstream &log) {
  int i;
  int irow_index, icol_index;
  // -----  compute state's energy and shift it before doing simulation
  // -------------
  vector<complex<double>> H_phi;
  H_phi.resize(dmatsize[0]);
  for (i = 0; i < dmatsize[0]; i++) {
    H_phi[i] = 0;
  }

  double de;
  double de_all;

  // use ground electronic state to compute energy and shift Hamiltonian
  // according to it. initial_state_in_sampling_state_index_list[0] is index for
  for (i = 0; i < dmatnum[0]; i++) {
    irow_index = local_dirow[0][i];
    icol_index = local_dicol[0][i]; // compute to point to colindex in
    H_phi[irow_index] =
        H_phi[irow_index] +
        dmat[0][i] *
            complex(
                xd[initial_state_in_sampling_state_index_list[0]][icol_index],
                yd[initial_state_in_sampling_state_index_list[0]][icol_index]);
  }
  de = 0;
  for (i = 0; i < dmatsize[0]; i++) {
    de = de +
         real(H_phi[i] *
              complex(xd[initial_state_in_sampling_state_index_list[0]][i],
                      -yd[initial_state_in_sampling_state_index_list[0]][i]));
  }
  MPI_Allreduce(&de, &de_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // shift Hamiltonian by detector energy de_all:
  if (my_id == 0) {
    cout << "Shift Hamiltonian by energy of system  " << de_all << endl;
    log << "Shift Hamiltonian by energy of system  " << de_all << endl;
  }
  for (i = 0; i < dmatsize[0]; i++) {
    dmat[0][i] = dmat[0][i] - de_all;
  }
}
