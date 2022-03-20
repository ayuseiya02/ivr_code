//
// Created by phyzch on 7/22/20.
//
#include "../System.h"
#include "../util.h"

using namespace std;
int Rmax; // maximum distance allowed in detector state space.

int state_distance(const vector<int> &ndetector, int *state1, int moddim) {
  // compute distance between two state : ndetector and state1. Dimension of mod
  // is moddim
  int i;
  int distance = 0;
  int initial_vibrational_state_index = 1;
  for (i = initial_vibrational_state_index; i < moddim; i++) {
    distance = distance + abs(state1[i] - ndetector[i]);
  }
  return distance;
}

int vmode_compare(const vector<int> &lhs, const vector<int> &rhs) {
  // compare function used for find_position_for_insert_binary.
  int i, j;
  int size = lhs.size();
  if (size != rhs.size()) {
    cout << "size of vector do not match." << endl;
    exit(-5);
  }
  for (i = 0; i < size; i++) {
    if (lhs[i] > rhs[i])
      return 1;
    else if (lhs[i] < rhs[i])
      return -1;
  }
  return 0;
}

// find position to insert the corresponding detector state, called in
// compute_matrix_size() function.  Binary search.
int find_position_for_insert_binary(const vector<vector<int>> &vmode,
                                    const vector<int> &ndetector, bool &exist) {
  // code used in compute_matrix_size() function
  // first compare first mode, 0 is in front of 1. So finally our mode is
  // |00000> , |00001>,|00010>, etc, |10000> This is just binary search
  // algorithm.
  if (vmode.size() == 0) {
    exist = false;
    return 0;
  }
  int left_flag = 0;
  int right_flag = vmode.size() - 1;
  int position = (left_flag + right_flag) / 2;
  exist = false;
  int mark;
  // binary search using compare function
  while (right_flag > left_flag) {
    mark = vmode_compare(ndetector, vmode[position]);
    if (mark == 1) {
      // ndetector should be after position
      left_flag = position + 1;
      position = (left_flag + right_flag) / 2;
    } else if (mark == -1) {
      right_flag = position - 1;
      position = (left_flag + right_flag) / 2;
    } else if (mark == 0) {
      exist = true;
      return position;
    }
  }
  // now left == right. Decide position to insert and return it.
  mark = vmode_compare(ndetector, vmode[position]);
  if (mark == 0) {
    exist = true;
    return position; // we will not insert.
  } else if (mark == -1) {
    exist = false;
    return position; // we will insert before this position.
  } else if (mark == 1) {
    exist = false;
    return position + 1; // we will insert after this position.
  }

  //
}

void full_system::compute_detector_matrix_size_MPI_sphere() {
  // use this function we construct detector state in a cube.
  if (my_id == 0) {
    int i, j, k;
    int i1;
    double energy;

    double energy_for_vibrational_state;
    // ndetector0 and ndetector1 indicate current detector mode quantum number
    // (0001001 for example)  we are considering.
    vector<int> ndetector0(d.nmodes[0]);

    int location;
    bool exist = false;

    int state_one_norm_distance;

    // ground_electronic_state_initial_energy_shift : energy shift for ground
    // state
    double ground_electronic_state_initial_energy_shift = 0;

    // excited_electronic_state_initial_energy_shift : energy shift for excited
    // state. this energy shift is due to different coupling strength in
    // different electronic state. See eq.(22 b ) in Logan's note.
    double excited_electronic_state_initial_energy_shift =
        pow(coupling_strength_to_mode0[0], 2) / d.mfreq[0][1] -
        pow(coupling_strength_to_mode0[1], 2) / d.mfreq[0][1];

    // energy for crossing region.  only for reference, not used in code.
    double crossing_region_energy =
        d.mfreq[0][1] / 2 *
        pow(d.mfreq[0][0] / (sqrt(2) * (coupling_strength_to_mode0[0] +
                                        coupling_strength_to_mode0[1])) +
                sqrt(2) * coupling_strength_to_mode0[0] / d.mfreq[0][1],
            2);

    cout << "Crossing region energy " << crossing_region_energy << endl;
    output << "Crossing region energy " << crossing_region_energy << endl;
    log << "Crossing region energy " << crossing_region_energy << endl;

    ndetector0[1] =
        -1; // this is for:  when we go into code: ndetector0[i]=
            // ndetector0[i]+1, our first vibrational  state is |000000>
    while (1) {
    label2:; // label2 is for detector0 to jump to beginning of while(1) loop
             // (this is inner layer of while(1))

      // ground electronic state
      ndetector0[0] = 0;
      // i1 begin from 1, this is because mode 0 is index to denote electronic
      // state. code below will let us go through all states (000000) - >
      // (100000) -> (nmax0, nmax1, .... , nmax_n) with constraint (energy
      // constrain + distance cutoff constraint) we add to sift states to use
      for (i1 = 1; i1 < d.nmodes[0]; i1++) { // loop through detector0
        // define the way we loop through vibrational quantum numbe   (000000) -
        // > (nmax0, nmax1, .... , nmax_n):
        ndetector0[i1] = ndetector0[i1] + 1;
        if (ndetector0[i1] <= d.nmax[0][i1])
          break;
        if (ndetector0[d.nmodes[0] - 1] > d.nmax[0][d.nmodes[0] - 1]) {
          ndetector0[d.nmodes[0] - 1] = 0;
          goto label1; // use goto to jump out of nested loop
        }
        ndetector0[i1] = 0;
      }

      energy = 0;
      energy_for_vibrational_state = 0;
      // calculate detector 0 energy.
      for (i = 0; i < d.nmodes[0]; i++) {
        energy = energy + ndetector0[i] * d.mfreq[0][i];
      }
      for (i = 1; i < d.nmodes[0]; i++) {
        energy_for_vibrational_state =
            energy_for_vibrational_state + ndetector0[i] * d.mfreq[0][i];
      }

      //--------------------------------------------------------------------------------------------
      // criteria below make sure detector 0 's energy is reasonable.
      // d.initial_Detector_energy[0] only include vibrational energy of
      // electronic state.
      if (energy_for_vibrational_state >
          d.initial_Detector_energy[0] + d.detector_energy_window_size) {
        // detector 0's energy can not be larger than its
        // initial_Detector_energy (This is energy of crossing region, defined
        // in function: construct_initial_state_MPI ) + energy window size we
        // impose. jump to next detector state whose energy may be lower. for
        // example (0031**) -> (0002**)  (we order state by their vibrational
        // mode, we first loop through mode 1 and then mode 2 , then mode 3 ...)

        // start with first vibrational mode
        k = 1;
        while (ndetector0[k] == 0) {
          ndetector0[k] = d.nmax[0][k];
          k++;
          if (k >= d.nmodes[0]) {
            break;
          }
        }
        if (k < d.nmodes[0]) {
          ndetector0[k] = d.nmax[0][k];
        }
        goto label2;
      }

      //------------------------------------------------------------------------------------------------
      // criteria below make sure detector 0 can not be too far away from
      // initial state this distance only apply to vibrational state space.
      // 1-norm distance applies.
      state_one_norm_distance =
          state_distance(ndetector0, d.initial_detector_state[0], d.nmodes[0]);
      if (state_one_norm_distance > Rmax) {
        goto label2;
      }

      // shift energy by ground state energy of each electronic state
      energy = ground_electronic_state_initial_energy_shift + energy;

      //--------------------------------------insert this state in detector's
      //state.-----------------------------------------------------------
      // --------------------------  data structure we use here is sorted linear
      // list. We sort state according to their vibrational quantum number (n0,
      // n1, n2, ...). Compare quanta in small mode index first
      location = find_position_for_insert_binary(
          vmode0, ndetector0,
          exist); // we check if this mode exist and the location we have to
                  // insert this state at the same time.
      if (!exist) {
        // when we push back we should consider arrange them in order. We
        // compute location to insert in find_position_for_insert_binary()
        // function:
        vmode0.insert(vmode0.begin() + location, ndetector0);
        dmat0.insert(dmat0.begin() + location, energy);
      }
    }
  label1:;

    for (i = 0; i < d.nmodes[0]; i++) {
      ndetector0[i] = 0;
    }

    ndetector0[1] =
        -1; // this is for:  when we go into code: ndetector0[i]=
            // ndetector0[i]+1, our first vibrational  state is |000000>
    while (1) {
    label4:; // label4 is for detector1 to jump to beginning of while(1) loop
             // (this is inner layer of while(1))

      // excited electronic state
      ndetector0[0] = 1;
      for (i1 = 1; i1 < d.nmodes[0];
           i1++) { // loop through detector1's vibrational state
        // define the way we loop through detector1:
        ndetector0[i1] = ndetector0[i1] + 1;
        if (ndetector0[i1] <= d.nmax[1][i1])
          break;
        if (ndetector0[d.nmodes[1] - 1] > d.nmax[1][d.nmodes[1] - 1]) {
          ndetector0[d.nmodes[1] - 1] = 0;
          goto label3; // use goto to jump out of nested loop
        }
        ndetector0[i1] = 0;
      }

      energy = 0;
      energy_for_vibrational_state = 0;
      // calculate detector 1 energy
      for (i = 0; i < d.nmodes[1]; i++) {
        energy = energy + ndetector0[i] * d.mfreq[1][i];
      }
      for (i = 1; i < d.nmodes[1]; i++) {
        energy_for_vibrational_state =
            energy_for_vibrational_state + ndetector0[i] * d.mfreq[1][i];
      }

      //--------------------------------------------------------------------------------------------
      // criteria below make sure detector 0 's energy is reasonable.
      // d.initial_Detector_energy[0] only include vibrational energy of
      // electronic state.
      if (energy_for_vibrational_state >
          d.initial_Detector_energy[1] + d.detector_energy_window_size) {
        // detector 0's energy can not be larger than its initial energy +
        // photon energy jump to next detector state.

        // start with first vibrational mode
        k = 1;
        while (ndetector0[k] == 0) {
          ndetector0[k] = d.nmax[0][k];
          k++;
          if (k >= d.nmodes[0]) {
            break;
          }
        }
        if (k < d.nmodes[0]) {
          ndetector0[k] = d.nmax[0][k];
        }
        goto label4;
      }

      //------------------------------------------------------------------------------------------------
      // criteria below make sure detector 1 can not be too far away from bright
      // state and lower bright state. this distance only apply to vibrational
      // state space.
      state_one_norm_distance =
          state_distance(ndetector0, d.initial_detector_state[1], d.nmodes[1]);
      if (state_one_norm_distance > Rmax) {
        goto label4;
      }

      // shift energy by ground state energy (energy of |00000> vibrational
      // state) of each electronic state
      energy = excited_electronic_state_initial_energy_shift + energy;

      // Notice, in code there are variables vmode1, dmat1. They are used in
      // TLS+ photon model (there Hamiltonian is product phase space of two
      // molecules. See our recent paper :
      // https://aip.scitation.org/doi/10.1063/5.0043665 .) In this program,
      // variable vmode1, dmat1 is not used and only vmode0 and ndetector0 is
      // used.
      //--------------------------------------insert this state in detector's
      //state.-----------------------------------------------------------
      location = find_position_for_insert_binary(
          vmode0, ndetector0,
          exist); // we check if this mode exist and the location we have to
                  // insert this state at the same time.
      if (!exist) {
        // when we push back we should consider arrange them in order. We
        // compute location to insert in find_position_for_insert_binary()
        // function:
        vmode0.insert(vmode0.begin() + location, ndetector0);
        dmat0.insert(dmat0.begin() + location, energy);
      }
    }
  label3:;
  }
}
