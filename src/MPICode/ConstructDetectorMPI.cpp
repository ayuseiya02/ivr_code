//
// Created by phyzch on 6/23/20.
// This file contain function that use our previously coded programm to
// construct matrix.
//
#include "../System.h"
#include "../util.h"
#include <random>

using namespace std;

int distance_cutoff_for_4_piont_corre;
double Energy_Range_4_point_corre_function_average =
    0; // average over bunch of states within energy window for computation of 4
       // point correlation function
int Distance_Range_4_point_corre_function_average =
    0; // average over bunch of states within distance cutoff for computation of
       // 4 point correlation function

void Broadcast_dmat_vmode(int stlnum, vector<double> &dmat0,
                          vector<double> &dmat1, vector<vector<int>> &vmode0,
                          vector<vector<int>> &vmode1);

// initialize parameters for detector (Molecules) and allocate space there.
void detector::allocate_space(int tlnum) {
  int i, j;
  dmatdim = detdim * detdim / fillfrac;

  nmodes = new int[tlnum]; //  number of modes in each detector

  proptime = new double[tlnum]; //  pre-coupling propagation time

  nmax = new int *[tlnum]; // maximum number of quanta in eqch mode.

  // if mode is bright or dark
  modtype = new int *[tlnum];

  mfreq = new double *[tlnum]; // mfreq: frequency of each mode here.

  aij = new double *[tlnum];

  modcoup = new double *[tlnum]; // coupling of matrix mod here.

  premodcoup = new double
      *[tlnum]; // coupling before system and detector contact with each other.

  dmat = new vector<double>[tlnum];

  for (i = 0; i < tlnum; i++) {
    vector<int> v1;
    dirow.push_back(v1);
    dicol.push_back(v1);
    vector<vector<int>> v2;
    dv.push_back(v2);
  }

  // matrix element number for detector matrix
  dmatnum = new int[tlnum];
  // off diagonal matrix element number for detector matrix
  doffnum = new int[tlnum];

  dmatsize = new int[tlnum]; // size of detector matrix

  // tell detector total matrix size..
  total_dmat_size.reserve(2);
  total_dmat_num.reserve(2);
  total_dmat_off_num.reserve(2);
  dmatsize_each_process = new int *[2];
  dmatsize_offset_each_process = new int *[2];
  doffnum_each_process = new int *[2];
  dmatnum_each_process = new int *[2];
  ; // record detector matrix element number in each process.
  dmat_offset_each_process = new int *[2];
  ; // record local first detector matrix's index in global matrix.
  for (i = 0; i < tlnum; i++) {
    dmatsize_each_process[i] = new int[num_proc];
    dmatsize_offset_each_process[i] = new int[num_proc];
    doffnum_each_process[i] = new int[num_proc];
    dmatnum_each_process[i] = new int[num_proc];
    dmat_offset_each_process[i] = new int[num_proc];
  }
};

void detector::allocate_space_single_detector(int detector_index) {
  int i = detector_index;
  nmax[i] = new int[nmodes[i]];
  modtype[i] = new int[nmodes[i]];
  mfreq[i] = new double[nmodes[i]];
  aij[i] = new double[nmodes[i]];
  modcoup[i] = new double[nmodes[i]];
  premodcoup[i] = new double[nmodes[i]];
}

// read parameters for  detector (Molecules) matrix construction.
void detector::read_MPI(ifstream &input, ofstream &output, ofstream &log,
                        int tlnum, int tldim, string path) {
  int i, j;

  stlnum = tlnum;
  stldim = tldim;
  allocate_space(tlnum);
  if (my_id == 0) {
    // maxdis: maximum allowed 1-norm distance in state space when constructing
    // coupling between states. for example distance between (1, 0, 1) and
    // (0,2,1) is 1 + 2 + 0 == 3 . cutoff: cutoff strength for intra-detector
    // coupling strength.  cutoff2: cutoff strength for inter-detector coupling.
    // kelvin: temperature used to initialize detector state. (currently no use
    // in program now. can be ignored)
    input >> maxdis >> cutoff >> kelvin;
    if (!Detector_Continue_Simulation) {
      output << "Detector  " << maxdis << " " << cutoff << " " << kelvin
             << endl;
    }
    for (i = 0; i < tlnum; i++) {
      input >> nmodes[i] >> proptime[i];
      allocate_space_single_detector(i);
      for (j = 0; j < nmodes[i]; j++) {
        // modtype , premodcoup and modcoup is currently not in use in our
        // program, can ignore them. mfreq: frequency of harmonic oscillator for
        // each dof. In this program. mfreq[i][0] which is first mode is chosen
        // to stands for ground state energy for electronic state. nmax: maximum
        // allowed quanta in each oscillator. here mfreq [ stlnum ] [ dof ] .
        // stlnum == 2 in our program, which represent two electronic state. dof
        // is vibrational degree of freedom in each electronic state
        input >> mfreq[i][j] >> nmax[i][j] >> modtype[i][j] >>
            premodcoup[i][j] >> modcoup[i][j];
      }
    }
    // -----------------------------------------------------------------------------------------------------
    if (!Detector_Continue_Simulation) {
      for (i = 0; i < tlnum; i++) {
        for (j = 0; j < nmodes[0]; j++) {
          // set precision of mfreq to 0.01. Convenient to restart simulation.
          mfreq[i][j] = floor(mfreq[i][j] * 100) / 100;
        }
      }
    } else {
      // this module is used to restart program. We can record wave function of
      // whole system into save_wavefunction.txt and save_Hamiltonian.txt then
      // read data there to restart program.  However, in current program we do
      // not use this module. Can ignore this part.
      ifstream load;
      load.open(path + "save_wavefunction.txt");
      if (load.is_open()) {
        for (i = 0; i < tlnum; i++) {
          for (j = 0; j < nmodes[0]; j++) {
            load >> mfreq[i][j];
          }
        }
        load.close();
      } else {
        log << "Can not load Detector frequency from save_wave function.txt"
            << endl;
        MPI_Abort(MPI_COMM_WORLD,
                  -6); // error code -6 : fail to load detector frequency
      }
    }
    //-----------------------------------------------------------------------------------------------------
    for (i = 0; i < tlnum; i++) {
      // output information to output.txt file.
      if (!Detector_Continue_Simulation) {
        output << nmodes[i] << " " << proptime[i] << endl;
      }

      for (j = 0; j < nmodes[i]; j++) {
        // a_ij is scaling factor. anharmonic coupling has the form:  V =
        // V^{(3)} a^{\Delta n} . here aij refer to scaling factor a. we
        // recently publish a paper with same scaling model, check
        // https://aip.scitation.org/doi/10.1063/5.0043665. eq(1)(2) is the
        // Hamiltonian we currently employ for anharmonic vibrational coupling
        // we choose mean frequency \Omega = 500
        aij[i][j] = a_intra * pow(double(mfreq[i][j]), 0.5) / pow(500, 0.5);
        // aij corresponding to scaling factor for f= f_{bright}/2 cm^{-1}.
        if (!Detector_Continue_Simulation) {
          output << mfreq[i][j] << " " << nmax[i][j] << " " << modtype[i][j]
                 << " " << premodcoup[i][j] << " " << modcoup[i][j] << endl;
        }
      }
    }
  }

  // Do same thing in other process use MPI function.
  MPI_Bcast(&nmodes[0], tlnum, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&proptime[0], tlnum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (my_id != 0) {
    for (i = 0; i < stlnum; i++) {
      allocate_space_single_detector(i);
    }
  }

  // variable used for constructing anharmonic coupling.
  deln = new int[max(nmodes[0], nmodes[1])];
  nbar = new double[max(nmodes[0], nmodes[1])];

  MPI_Bcast(&maxdis, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cutoff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&kelvin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (i = 0; i < stlnum; i++) {
    MPI_Bcast(&mfreq[i][0], nmodes[i], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nmax[i][0], nmodes[i], MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&modtype[i][0], nmodes[i], MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&premodcoup[i][0], nmodes[i], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&modcoup[i][0], nmodes[i], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&aij[i][0], nmodes[i], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
};

void detector::construct_dmatrix_MPI(ifstream &input, ofstream &output,
                                     ofstream &log, vector<double> &dmat0,
                                     vector<double> &dmat1,
                                     vector<vector<int>> &vmode0,
                                     vector<vector<int>> &vmode1) {
  int m;

  // previously information for state space is in process 0. Now we broadcast
  // this information to all procerss.
  construct_dv_dirow_dicol_dmatrix_MPI(log, dmat0, dmat1, vmode0, vmode1);

  // compute index for initial state (here state at crossing region) in
  // Hamiltonian.
  compute_important_state_index();

  // construct off-diagonal matrix element of Hamiltonian.
  compute_detector_offdiag_part_MPI(log, dmat0, dmat1, vmode0, vmode1);

  // record energy of vibrational state.
  dmat_energy = dmat0;

  // each process communicate with other process to tell them number of
  // anharmonic coupling (off diagonal element) they constructed
  broadcast_dmatnum_doffnum();

  // broadcast Hamiltonian: dmat  constructed to other process . (total_dmat).
  // same for dirow, dicol   ->  (total_dirow, total_dicol)
  broadcast_total_dmat();

  if (my_id == 0) {
    for (m = 0; m < 1; m++) {
      if (!Detector_Continue_Simulation) {
        output << "Matrix for detector " << m << " : " << total_dmat_size[m]
               << "  " << total_dmat_off_num[m] << endl;
      }
    }
  }

  if (not load_sampling_state) {
    // construct sampling_state_index which contain states we want to compute
    // IPR
    construct_sampling_state_index(dmat0);
  } else {
    // load sampling_state_index from file: sampling_state_info
    load_sampling_state_index(log);
  }
}

void detector::construct_sampling_state_index(vector<double> &dmat0) {
  // dmat0 is energy level for system
  int i, j;
  int state_number = dmat0.size();
  double energy;
  int index;

  int initial_state_index_in_all_process;
  bool initial_state_exist;
  int list_size;

  double min_energy = *min_element(dmat0.begin(), dmat0.end());
  double max_energy = *max_element(dmat0.begin(), dmat0.end());

  int block_number = 20;
  int element_num_in_each_block = 10;

  double energy_step = (max_energy - min_energy) / block_number;
  vector<vector<int>> state_in_energy_block_list;

  int sampling_state_num;
  int *sampling_state_array;

  // Code below divide vibrational states according to their energy into 20
  // blocks, each block we sample 10 states to compute IPR.
  if (my_id == 0) {

    for (i = 0; i < block_number; i++) {
      vector<int> v;
      state_in_energy_block_list.push_back(v);
    }

    for (i = 0; i < state_number; i++) {
      energy = dmat0[i];
      index = int((energy - min_energy) / energy_step);
      // in case we encounter max_energy element
      if (index == block_number) {
        index = block_number - 1;
      }
      state_in_energy_block_list[index].push_back(i);
    }

    for (i = 0; i < block_number; i++) {
      list_size = state_in_energy_block_list[i].size();
      if (element_num_in_each_block >= list_size) {
        for (j = 0; j < list_size; j++) {
          sampling_state_index.push_back(state_in_energy_block_list[i][j]);
        }
      } else {
        // Shuffle the list and choose first several element.
        std::shuffle(state_in_energy_block_list[i].begin(),
                     state_in_energy_block_list[i].end(),
                     std::mt19937(std::random_device()()));
        for (j = 0; j < element_num_in_each_block; j++) {
          sampling_state_index.push_back(state_in_energy_block_list[i][j]);
        }
      }
    }
  }

  // broadcasting element in sampling_state_index.
  if (my_id == 0) {
    sampling_state_num = sampling_state_index.size();
    sampling_state_array = new int[sampling_state_num];
    for (i = 0; i < sampling_state_num; i++) {
      sampling_state_array[i] = sampling_state_index[i];
    }
  }
  MPI_Bcast(&sampling_state_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (my_id != 0) {
    sampling_state_array = new int[sampling_state_num];
  }

  MPI_Bcast(&sampling_state_array[0], sampling_state_num, MPI_INT, 0,
            MPI_COMM_WORLD);
  if (my_id != 0) {
    for (i = 0; i < sampling_state_num; i++) {
      sampling_state_index.push_back(sampling_state_array[i]);
    }
  }

  delete[] sampling_state_array;

  // check if initial_state is in sampling state index. If not , add it. Also
  // find initial_state's index in sampling state index list
  for (i = 0; i < 2; i++) {
    initial_state_in_sampling_state_index_list.push_back(0);
  }
  sampling_state_num = sampling_state_index.size();
  for (i = 0; i < 2; i++) {
    // now initial_state_index record initial state in each electronic state
    initial_state_index_in_all_process =
        initial_state_index[i] +
        initial_state_pc_id[i] * int(state_number / num_proc);

    initial_state_exist = false;
    for (j = 0; j < sampling_state_num; j++) {
      if (initial_state_index_in_all_process == sampling_state_index[j]) {
        initial_state_exist = true;
        initial_state_in_sampling_state_index_list[i] = j;
        break;
      }
    }

    if (not initial_state_exist) {
      sampling_state_index.push_back(initial_state_index_in_all_process);
      initial_state_in_sampling_state_index_list[i] =
          sampling_state_index.size();
    }
  }

  sampling_state_num = sampling_state_index.size();

  // output vmode for state we selected as sampling state
  if (my_id == 0) {
    ofstream Sampling_state_mode_info(path + "sampling_state_info.txt");
    Sampling_state_mode_info << sampling_state_num << "  " << nmodes[0] << endl;
    for (i = 0; i < sampling_state_num; i++) {
      for (j = 0; j < nmodes[0]; j++) {
        Sampling_state_mode_info << dv_all[0][sampling_state_index[i]][j]
                                 << "  ";
      }
      Sampling_state_mode_info << endl;
    }

    Sampling_state_mode_info.close();
  }
}

void detector::load_sampling_state_index(ofstream &log) {
  int state_number = dv_all[0].size();
  if (my_id == 0) {
    cout << "Read sampling state from sampling_state_info.txt" << endl;
  }
  int i, j;
  int sampling_state_number;
  int number_of_mode;
  bool exist;
  int location;
  vector<int> mode_quanta;
  int quanta;
  string ss;
  if (my_id == 0) {

    ifstream load;
    load.open(path + "sampling_state_info.txt");
    if (load.is_open()) {
      load >> sampling_state_number >> number_of_mode;
      std::getline(load, ss);
      if (number_of_mode != nmodes[0]) {
        cout << "number of mode read from file not equal to nmodes[0]. Check "
                "sampling_state_info.txt. "
             << endl;
        log << "number of mode read from file not equal to nmodes[0]. Check "
               "sampling_state_info.txt. "
            << endl;
        MPI_Abort(MPI_COMM_WORLD, -23);
      }

      for (i = 0; i < sampling_state_number; i++) {
        mode_quanta.clear();
        for (j = 0; j < number_of_mode; j++) {
          load >> quanta;
          mode_quanta.push_back(quanta);
        }
        std::getline(load, ss);
        location =
            find_position_for_insert_binary(dv_all[0], mode_quanta, exist);
        if (not exist) {
          cout << "State with quanta: ";
          for (j = 0; j < number_of_mode; j++) {
            cout << mode_quanta[j] << " ";
          }
          cout << " Can not found in system. Check it. " << endl;
          MPI_Abort(MPI_COMM_WORLD, -23);
        } else {
          sampling_state_index.push_back(location);
        }
      }

    } else {
      cout << "Can not open sampling_state_info.txt. Check if this file exist "
              "in folder . "
           << endl;
      log << "Can not open sampling_state_info.txt. Check if this file exist "
             "in folder . "
          << endl;
      MPI_Abort(MPI_COMM_WORLD, -23);
    }
  }

  // Broadcast this to other process
  MPI_Bcast(&sampling_state_number, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int *sampling_state_array = new int[sampling_state_number];
  if (my_id == 0) {
    for (i = 0; i < sampling_state_number; i++) {
      sampling_state_array[i] = sampling_state_index[i];
    }
  }

  MPI_Bcast(&sampling_state_array[0], sampling_state_number, MPI_INT, 0,
            MPI_COMM_WORLD);
  if (my_id != 0) {
    for (i = 0; i < sampling_state_number; i++) {
      sampling_state_index.push_back(sampling_state_array[i]);
    }
  }

  // find initial state index
  int initial_state_index_in_all_process;
  bool initial_state_exist;
  for (i = 0; i < 2; i++) {
    initial_state_in_sampling_state_index_list.push_back(0);
  }
  sampling_state_number = sampling_state_index.size();
  for (i = 0; i < 2; i++) {
    // now initial_state_index record initial state in each electronic state
    initial_state_index_in_all_process =
        initial_state_index[i] +
        initial_state_pc_id[i] * int(state_number / num_proc);

    initial_state_exist = false;
    for (j = 0; j < sampling_state_number; j++) {
      if (initial_state_index_in_all_process == sampling_state_index[j]) {
        initial_state_exist = true;
        initial_state_in_sampling_state_index_list[i] = j;
        break;
      }
    }

    if (not initial_state_exist) {
      if (my_id == 0) {
        cout << "Initial state is not in sampling_state_info.txt. Check what's "
                "wrong.  "
             << endl;
        MPI_Abort(MPI_COMM_WORLD, -23);
      }
    }
  }

  delete[] sampling_state_array;
}

void detector::construct_dv_dirow_dicol_dmatrix_MPI(
    ofstream &log, vector<double> &dmat0, vector<double> &dmat1,
    vector<vector<int>> &vmode0, vector<vector<int>> &vmode1) {
  //  dv : vibrational quantum number for states in each process.
  // dv_all[0] = vmode0. vibrational quantum number for all states in this
  // program. total_dmat_size:  size of total Hamiltonian. (size of sparse
  // matrix or number of total vibrational states) dmatsize:  number of
  // vibrational states in each process. we only use dmat0 and vmode0, thus only
  // total_dmat_size[0] , dmat_size[0], dv_all[0], dv[0] is constructed.
  // total_dmat_size[1] , dmat_size[1] etc is not used.

  // fixme: Understanding detail of these MPI code is tedious, also writing them
  // is also tedious. (That's drawback of MPI). But their functions are pretty
  // clear.
  // fixme: I suggest not dig into details for MPICode part when you read the
  // code the first time.  Only understand what these function is about is
  // enough.

  int i, m;
  int vsize, vsize2;
  int displacement;
  if (my_id == 0) {
    total_dmat_size[0] = dmat0.size();
  }
  MPI_Bcast(&total_dmat_size[0], 2, MPI_INT, 0, MPI_COMM_WORLD);
  vector<int> *dirow_all, *dicol_all;
  vector<double> *dmat_all;
  int **vector_size, **displacement_list;

  for (i = 0; i < 1; i++) {
    vector<vector<int>> v1;
    dv_all.push_back(v1);
  }

  dmat_all = new vector<double>[2]; // do not confuse dmat_all with total_dmat.
                                    // dmat_all only contain diagonal term
  dirow_all = new vector<int>[2];
  dicol_all = new vector<int>[2];
  vector_size = new int *[2];
  displacement_list = new int *[2];

  Broadcast_dmat_vmode(
      stlnum, dmat0, dmat1, vmode0,
      vmode1); // broadcast dmat0, dmat1, vmode0, vmode1 to all process.
  dv_all[0] = vmode0;

  if (my_id == 0) {
    dmat_all[0] = dmat0;
    // prepare dirow, dicol for broadcasting.
    for (m = 0; m < 1; m++) {
      int size = dmat_all[m].size();
      for (i = 0; i < size; i++) {
        dirow_all[m].push_back(i);
        dicol_all[m].push_back(i);
      }
    }
    // prepare vector size and vector displacement:
    for (m = 0; m < 1; m++) {
      vsize = total_dmat_size[m] / num_proc;
      vsize2 = total_dmat_size[m] - (num_proc - 1) * vsize;
      vector_size[m] =
          new int[num_proc]; // size of vector to scatter to each process.
      displacement_list[m] = new int[num_proc]; // displacement of vector to
                                                // scatter to each process.
      displacement = 0;
      for (i = 0; i < num_proc - 1; i++) {
        vector_size[m][i] = vsize;
        displacement_list[m][i] = displacement;
        displacement = displacement + vsize;
      }
      vector_size[m][num_proc - 1] = vsize2;
      displacement_list[m][num_proc - 1] = displacement;
    }
  }

  // broadcast Hamiltonian diagonal element in dmat_all (dmat0) to all other
  // process. also for dirow_all (row index of element) and dicol_all (column
  // index for element)
  Scatter_dirow_dicol_dmatrix(dmat_all, dirow_all, dicol_all, vector_size,
                              displacement_list, log);

  // scatter dv_all to dv
  Scatter_dv(total_dmat_size);

  // construct detector matrix size for each process.
  if (my_id != num_proc - 1) {
    for (m = 0; m < 1; m++) {
      vsize = total_dmat_size[m] / num_proc;
      dmatsize[m] = vsize;
    }
  } else {
    for (m = 0; m < 1; m++) {
      vsize =
          total_dmat_size[m] - total_dmat_size[m] / num_proc * (num_proc - 1);
      dmatsize[m] = vsize;
    }
  }
  for (m = 0; m < 1; m++) {
    MPI_Allgather(&dmatsize[m], 1, MPI_INT, &dmatsize_each_process[m][0], 1,
                  MPI_INT, MPI_COMM_WORLD);
  }
  for (m = 0; m < 1; m++) {
    dmatsize_offset_each_process[m][0] = 0;
    for (i = 1; i < num_proc; i++) {
      dmatsize_offset_each_process[m][i] =
          dmatsize_offset_each_process[m][i - 1] +
          dmatsize_each_process[m][i - 1];
    }
  }

  if (my_id == 0) {
    for (m = 0; m < 1; m++) {
      delete[] vector_size[m];
      delete[] displacement_list[m];
    }
  }
  delete[] vector_size;
  delete[] displacement_list;
  delete[] dirow_all;
  delete[] dicol_all;
}

void Broadcast_dmat_vmode(int stlnum, vector<double> &dmat0,
                          vector<double> &dmat1, vector<vector<int>> &vmode0,
                          vector<vector<int>> &vmode1) {
  int i, j, m;
  int dmat0_size;
  // we broadcast dmat0, dmat1, .. to all other process. This is need for us to
  // compute off diagonal matrix first allocate space for dmat0 , dmat1.
  if (my_id == 0) {
    dmat0_size = dmat0.size();
  }
  MPI_Bcast(&dmat0_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (my_id != 0) {
    dmat0.resize(dmat0_size);
  }
  // Broadcast dmat0, dmat1
  MPI_Bcast((void *)dmat0.data(), dmat0_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // ----------------Broadcast vmode0, vmode1 -------------------------------
  vector<vector<int>> *v_ptr;
  // ------------- For sending vmode0,vmode1 --------------------------
  vector<int> vmode_1d;
  vector<int> displacement;
  vector<int> element_size;
  int vmode_1d_size;
  int element_number;
  //-------------------------------------------
  for (m = 0; m < 1; m++) {
    vmode_1d.clear();
    displacement.clear();
    element_size.clear();
    if (m == 0)
      v_ptr = &(vmode0);
    else
      v_ptr = &(vmode1);
    if (my_id == 0) {
      convert_dv(*v_ptr, vmode_1d, displacement, element_size);
      vmode_1d_size = vmode_1d.size();
      MPI_Bcast(&vmode_1d_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
      element_number = element_size.size();
      MPI_Bcast(&element_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast((void *)element_size.data(), element_number, MPI_INT, 0,
                MPI_COMM_WORLD); // send size of each 2d element
      MPI_Bcast((void *)vmode_1d.data(), vmode_1d_size, MPI_INT, 0,
                MPI_COMM_WORLD); // send converted 1d data to other process.
    } else {
      MPI_Bcast(&vmode_1d_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
      vmode_1d.reserve(vmode_1d_size);
      MPI_Bcast(&element_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
      element_size.reserve(element_number);
      MPI_Bcast((void *)element_size.data(), element_number, MPI_INT, 0,
                MPI_COMM_WORLD);
      MPI_Bcast((void *)vmode_1d.data(), vmode_1d_size, MPI_INT, 0,
                MPI_COMM_WORLD);

      int index = 0;
      // convert vmode_1d to vmode0 / vmode1
      for (i = 0; i < element_number; i++) {
        vector<int> dv_vmode;
        dv_vmode.reserve(element_size[i]);
        for (j = 0; j < element_size[i]; j++) {
          dv_vmode.push_back(vmode_1d[index]);
          index++;
        }
        (*v_ptr).push_back(dv_vmode);
      }
    }
  }
}

double factorial(int n) {
  double f = 1;
  int i;
  for (i = 2; i <= n; i++) {
    f = f * i;
  }
  return f;
}

void detector::compute_detector_offdiag_part_MPI(ofstream &log,
                                                 vector<double> &dmat0,
                                                 vector<double> &dmat1,
                                                 vector<vector<int>> &vmode0,
                                                 vector<vector<int>> &vmode1) {
  int i, j, m, k;
  int begin_index;
  int ntot;
  double value, lij;
  vector<vector<int>> *vmode_ptr;
  vector<double> *dmat_ptr;
  bool exist;
  int position;
  double random_number;

  bool Other_vibrational_state_same;
  int Minimum_Mode_1_index;
  double tunneling_matrix_element;
  double factorial_coefficient;
  int m_index, n_index;

  // See Logan's note for Alpha  (eq. 31(b) ) , Crossing_point_quanta_spin_down.
  // Crossing_point_quanta_spin_up. (eq.47) alpha = (\lambda_{up} +
  // \lambda_{down} ) / \hbar \omega_{0}.  \omega_{0} is frequency of
  // vibrational state couple to electronic state, here is mfreq[0][1]
  double Alpha =
      (coupling_strength_to_mode0[0] + coupling_strength_to_mode0[1]) /
      (mfreq[0][1]);

  double Crossing_point_quanta_spin_down =
      pow((mfreq[0][0]) / (coupling_strength_to_mode0[0] +
                           coupling_strength_to_mode0[1]) +
              2 * coupling_strength_to_mode0[0] / (mfreq[0][1]),
          2) /
      4;

  double Crossing_point_quanta_spin_up =
      pow((mfreq[0][0]) / (coupling_strength_to_mode0[0] +
                           coupling_strength_to_mode0[1]) -
              2 * coupling_strength_to_mode0[1] / (mfreq[0][1]),
          2) /
      4;

  if (my_id == 0) {
    cout << "ALpha :   " << Alpha << endl;
    log << "ALpha:  " << Alpha << endl;
  }

  double cutoff_for_coupling_between_electronic_state = 0.001;

  // different process do different amount of work.
  vmode_ptr = &(vmode0);
  dmat_ptr = &(dmat0);
  begin_index = total_dmat_size[0] / num_proc * my_id;
  // compute off diagonal matrix element, using dmat0
  for (i = begin_index; i < begin_index + dmatsize[0]; i++) {
    for (j = 0; j < total_dmat_size[0];
         j++) { // j is different from serial program. Now we record both
                // symmetric Hamiltonian element
      if (i == j)
        continue;

      // Coupling within same electronic state.
      if ((*vmode_ptr)[i][0] == (*vmode_ptr)[j][0]) {
        ntot = 0;
        // k begins with 1, because index 0 is to denote different electronic
        // state.
        for (k = 1; k < nmodes[0]; k++) {
          deln[k] = abs(
              (*vmode_ptr)[i][k] -
              (*vmode_ptr)[j][k]); //  deln[k] = abs(dv[m][i][k] - dv[m][j][k]);
          nbar[k] = sqrt(sqrt(
              double(max(1, (*vmode_ptr)[i][k])) *
              double(max(
                  1, (*vmode_ptr)[j][k])))); // sqrt( vi[k] * vj[k] ).  This
                                             // account for coefficient \sqrt{n
                                             // + 1} generated by a^{+}|n> etc.
          ntot = ntot + deln[k];
        }
        // this is because we don't have q1 or q1 * q2 in Harmonic oscillator's
        // Hamiltonian
        if (ntot == 2) {
          for (k = 1; k < nmodes[0]; k++) {
            if (deln[k] == 2)
              deln[k] = 4;
            if (deln[k] == 1)
              deln[k] = 2;
          }
        } else if (ntot == 1) {
          for (k = 1; k < nmodes[0]; k++) {
            if (deln[k] == 1)
              deln[k] = 3;
          }
        } else if (ntot == 0) {
          log << "Error! The off-diagonal element must differ in q.n." << endl;
          MPI_Abort(MPI_COMM_WORLD, -8);
        }
        // check ntot with maxdis in quantum number space:
        // maxdis: 1-norm distance cutoff in state space for coupling between
        // states.
        if (ntot < maxdis) {

          if (ntot % 2 == 0) {
            value = V_intra; // V=0.03 as requirement.
          } else {
            value = -V_intra;
          }

          for (k = 1; k < nmodes[0]; k++) {
            // aij is scaling factor for mode k.
            value = value * pow(aij[0][k] * nbar[k], deln[k]);
          }
          // lij == V / (\Delta E) . Here V is anharmonic coupling between
          // states, \Delta E is energy difference between states. Recall
          // time-independent perturbation in QM. only when lij > cutoff, this
          // anharmonic coupling will be important for dynamics.
          if ((*dmat_ptr)[i] != (*dmat_ptr)[j]) {
            lij = abs(value / ((*dmat_ptr)[i] - (*dmat_ptr)[j]));
            if (lij > cutoff) {
              dmat[0].push_back(value);
              dirow[0].push_back(i);
              dicol[0].push_back(j);
            }
          } else {
            dmat[0].push_back(value);
            dirow[0].push_back(i);
            dicol[0].push_back(j);
          }
        }
      }

      // Coupling between states with different electronic state.
      // Code below is for coupling between vibrational states with different
      // electronic state.
      // fixme: I think it's hard to understand code below before reading
      // Logan's note and my own derivation. You can first skip details here, we
      // can discuss them at meeting later.
      if ((*vmode_ptr)[i][0] != (*vmode_ptr)[j][0]) {
        // Check other vibrational state except vibrational mode with index: 1
        // is the same
        Other_vibrational_state_same = true;
        for (k = 2; k < nmodes[0]; k++) {
          if ((*vmode_ptr)[i][k] != (*vmode_ptr)[j][k]) {
            Other_vibrational_state_same = false;
            break;
          }
        }
        if (not Other_vibrational_state_same) {
          continue;
        }

        // include valid coupling
        // below we compute <spin_down, mode_1_index (m_index) | spin_up,
        // mode_1_index (n_index) > for <spin_up,  m_index , spin_down, n_index
        // > == < spin_down, n_index | spin_up, m_index>
        if ((*vmode_ptr)[i][0] == 0 and (*vmode_ptr)[j][0] == 1) {
          //  <spin_down | spin_up> case
          m_index = (*vmode_ptr)[i][1];
          n_index = (*vmode_ptr)[j][1];
        } else {
          // <spin_up | spin_down> case:
          m_index = (*vmode_ptr)[j][1];
          n_index = (*vmode_ptr)[i][1];
        }
        // below compute <spin_down, m_index | spin_up, n_index >
        Minimum_Mode_1_index = min(m_index, n_index);
        tunneling_matrix_element = 0;

        for (k = 0; k <= Minimum_Mode_1_index; k++) {
          // (-1)^{m_index - k}/ [(m_index - k)! * (n_index-k)! ] * sqrt(m! *
          // n!) / (k!)
          factorial_coefficient =
              1 / (factorial(m_index - k) * factorial(n_index - k)) *
              pow(-1, m_index - k) *
              sqrt(factorial(m_index) * factorial(n_index)) / (factorial(k));

          tunneling_matrix_element =
              tunneling_matrix_element +
              pow(Alpha, m_index + n_index - 2 * k) * factorial_coefficient;
        }
        // Coupling_between_electronic_state : t in Logan's note.
        tunneling_matrix_element = tunneling_matrix_element *
                                   exp(-pow(Alpha, 2) / 2) *
                                   Coupling_between_electronic_state;

        if (i == initial_state_index[0] and j == initial_state_index[1]) {
          cout << "Found " << endl;
        }

        // cutoff criteria. This equivalent to strong coupling is among states
        // in two electronic state which have similar energy. (near crossing
        // region)

        if ((*dmat_ptr)[i] != (*dmat_ptr)[j]) {
          lij =
              abs(tunneling_matrix_element / ((*dmat_ptr)[i] - (*dmat_ptr)[j]));
          if (lij > cutoff_for_coupling_between_electronic_state) {
            dmat[0].push_back(tunneling_matrix_element);
            dirow[0].push_back(i);
            dicol[0].push_back(j);
          }
        } else {
          dmat[0].push_back(tunneling_matrix_element);
          dirow[0].push_back(i);
          dicol[0].push_back(j);
        }
      }
    }
  }
}

void detector::broadcast_dmatnum_doffnum() {
  int m, i;
  // compute dmatnum and doffnum and dmatnum_each_process, total_dmatnum,
  // total_doffnum etc.
  for (m = 0; m < 1; m++) {
    dmatnum[m] = dmat[m].size();
    doffnum[m] = dmatnum[m] - dmatsize[m];
  }
  // compute toatl_dmatnum, total_dmatoffnum
  for (m = 0; m < 1; m++) {
    MPI_Allreduce(&dmatnum[m], &total_dmat_num[m], 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&doffnum[m], &total_dmat_off_num[m], 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
  }
  // broadcast detector matrix number and detector matrix off-diagonal number in
  // each process.
  for (m = 0; m < 1; m++) {
    MPI_Allgather(&dmatnum[m], 1, MPI_INT, &dmatnum_each_process[m][0], 1,
                  MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&doffnum[m], 1, MPI_INT, &doffnum_each_process[m][0], 1,
                  MPI_INT, MPI_COMM_WORLD);
  }
  // compute offset of detector matrix for each process.
  for (m = 0; m < 1; m++) {
    dmat_offset_each_process[m][0] = 0;
    for (i = 1; i < num_proc; i++) {
      dmat_offset_each_process[m][i] =
          dmat_offset_each_process[m][i - 1] + dmatnum_each_process[m][i - 1];
    }
  }
}

void detector::broadcast_total_dmat() {
  /*
   construct total_dmat, total_dirow, total_dicol
   use dmatnum_each_process,  dmat_offset_each_process\
  */
  total_dmat = new double *[stlnum];
  total_dirow = new int *[stlnum];
  total_dicol = new int *[stlnum];
  int m;
  for (m = 0; m < 1; m++) {
    total_dmat[m] = new double[total_dmat_num[m]];
    total_dirow[m] = new int[total_dmat_num[m]];
    total_dicol[m] = new int[total_dmat_num[m]];
    MPI_Allgatherv(&dmat[m][0], dmatnum[m], MPI_DOUBLE, &total_dmat[m][0],
                   dmatnum_each_process[m], dmat_offset_each_process[m],
                   MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(&dirow[m][0], dmatnum[m], MPI_INT, &total_dirow[m][0],
                   dmatnum_each_process[m], dmat_offset_each_process[m],
                   MPI_INT, MPI_COMM_WORLD);
    MPI_Allgatherv(&dicol[m][0], dmatnum[m], MPI_INT, &total_dicol[m][0],
                   dmatnum_each_process[m], dmat_offset_each_process[m],
                   MPI_INT, MPI_COMM_WORLD);
  }
}

//-------------------------------- Construct detector wave function
//------------------------------
void detector::initialize_detector_state_MPI(ofstream &log) {
  // this is the code for initializing detector state's wavefunction.
  int m, i, j;
  double norm;
  double total_norm;
  int local_index = 0;
  int dmat0_offset = my_id * int(total_dmat_size[0] / num_proc);
  int sampling_state_list_size = sampling_state_index.size();

  for (i = 0; i < sampling_state_list_size; i++) {
    vector<double> v1;
    v1.reserve(dmatsize[0]);
    xd.push_back(v1);
    yd.push_back(v1);
  }
  xd_all = new double *[sampling_state_list_size];
  yd_all = new double *[sampling_state_list_size];
  for (i = 0; i < sampling_state_list_size; i++) {
    xd_all[i] = new double[total_dmat_size[0]];
    yd_all[i] = new double[total_dmat_size[0]];
  }

  for (m = 0; m < sampling_state_list_size; m++) {
    norm = 0;
    for (i = 0; i < dmatsize[0]; i++) {
      xd[m].push_back(0);
      yd[m].push_back(0);
    }
    if (sampling_state_index[m] >= dmat0_offset and
        sampling_state_index[m] < dmat0_offset + dmatsize[0]) {
      // state m is in this process.
      local_index = sampling_state_index[m] - dmat0_offset;
      xd[m][local_index] = 1;
      norm = 1;
    }
    MPI_Allreduce(&norm, &total_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (total_norm == 0) {
      if (my_id == 0) {
        cout << "Norm for sampling state " << m << " is 0" << endl;
        log << "Norm for sampling state " << m << " is 0" << endl;
        MPI_Abort(MPI_COMM_WORLD, -10);
      }
    }
  }
}

void detector::construct_initial_state_MPI(ifstream &input, ofstream &output) {
  // MPI version of construct initial state
  /*
   *  Initial detector state: detector state populated at beginning of
   * simulation Bright  detector state: detector state with bright mode 1 quanta
   * higher in energy ,all other mode equal to initial detector state
   *  Initial_Detector_energy: Energy of detector at initial state.
   */
  int m, i, j;

  // The coupling is strong at crossing point for two electronic state.
  // (internal conversion) we choose mfreq[0][**] as frequency for lower
  // electronic state. mfreq[0][0] - mfreq[1][0] as  energy difference between
  // ground state for two electronic state. we choose mfreq[1][**] (** between 1
  // and dof) as frequency for higher electronic states. See Logan's note
  // eq.(47).
  double Crossing_point_quanta_spin_down =
      pow((mfreq[1][0]) / (coupling_strength_to_mode0[0] +
                           coupling_strength_to_mode0[1]) +
              2 * coupling_strength_to_mode0[0] / (mfreq[0][1]),
          2) /
      4;

  double Crossing_point_quanta_spin_up =
      pow((mfreq[1][0]) / (coupling_strength_to_mode0[0] +
                           coupling_strength_to_mode0[1]) -
              2 * coupling_strength_to_mode0[1] / (mfreq[0][1]),
          2) /
      4;

  // here initial_Detector_energy only count vibrational energy , not including
  // electronic ground state energy.
  initial_detector_state = new int *[stlnum];
  initial_Detector_energy = new double[stlnum];

  for (m = 0; m < stlnum; m++) {
    // nmodes[m] is number of modes (dof) for each electronic state.
    initial_detector_state[m] = new int[nmodes[m]];
    initial_Detector_energy[m] = 0;
  }

  if (my_id == 0) {

    // initialize our initial state from input.txt.  this is usually treated as
    // center point when constructing vibrational state in state space.
    for (m = 0; m < stlnum; m++) {
      for (i = 0; i < nmodes[m]; i++) {
        input >> initial_detector_state[m][i];
      }
    }

    // shift initial state to crossing point.
    // Initial_detector_state[0] is spin down state at crossing point.
    // Initial_detector_state[1] is spin up state and crossing point. which have
    // energy higher

    // first dof (index 0 ) is mode denote electronic state. second dof (index
    // 1) is mode denote vibrational harmonic oscillator couple to electronic
    // dof.
    initial_detector_state[0][0] = 0;
    initial_detector_state[0][1] = lround(Crossing_point_quanta_spin_down);

    initial_detector_state[1][0] = 1;
    initial_detector_state[1][1] = lround(Crossing_point_quanta_spin_up);
  }

  for (m = 0; m < stlnum; m++) {
    // Broad cast initial detector state to all process.
    MPI_Bcast(&initial_detector_state[m][0], nmodes[m], MPI_INT, 0,
              MPI_COMM_WORLD);
  }

  for (m = 0; m < stlnum; m++) {
    // initialize our initial detector(Molecule) state's vibrational energy.
    // This energy exclude mode 0 which represents electronic state
    for (i = 1; i < nmodes[m]; i++) {
      initial_Detector_energy[m] = initial_Detector_energy[m] +
                                   initial_detector_state[m][i] * mfreq[m][i];
    }
  }

  if (my_id == 0) { // output initial detector state to output.txt
    cout << "Initial detector state:" << endl;
    output << "Initial detector state:" << endl;
    for (m = 0; m < stlnum; m++) {
      for (i = 0; i < nmodes[m]; i++) {
        cout << initial_detector_state[m][i] << " ";
        output << initial_detector_state[m][i] << " ";
      }
      cout << endl;
      output << endl;
    }
  }
}

void detector::compute_important_state_index() {
  // compute initial state index for two detector and the process id it's in.
  // initial state index is index in local process.
  int m, i, j;
  int index;
  bool check_mark1, check_mark2;
  vector<int> special_state_vec;
  int **special_state;
  int *special_state_index;
  int *special_state_pc_id;
  int position;
  bool exist;
  bool *exist_bool_for_pc = new bool[num_proc];

  initial_state_index = new int[2];
  initial_state_pc_id = new int[2];
  // we initialize bright_state_index and initial_state_index.
  // We do this because sometimes we may not include bright state in our
  // simulation, then this index need to be initialized.
  for (m = 0; m < 2; m++) {
    initial_state_index[m] = 0;
    initial_state_pc_id[m] = 0;
  }

  special_state = initial_detector_state;
  special_state_index = initial_state_index;
  special_state_pc_id = initial_state_pc_id;

  // loop for state in two electronic state
  for (m = 0; m < 2; m++) {
    special_state_vec.clear();
    for (i = 0; i < nmodes[0]; i++) {
      special_state_vec.push_back(special_state[m][i]);
    }
    position = find_position_for_insert_binary(dv[0], special_state_vec, exist);

    MPI_Allgather(&exist, 1, MPI_C_BOOL, &exist_bool_for_pc[0], 1, MPI_C_BOOL,
                  MPI_COMM_WORLD);

    special_state_pc_id[m] = -1;

    for (i = 0; i < num_proc; i++) {
      if (exist_bool_for_pc[i]) {
        special_state_pc_id[m] = i;
      }
    }
    if (special_state_pc_id[m] == -1) {
      if (my_id == 0) {
        cout << " Can not find initial state or brigth state in all vmode. "
                "Must have bug here."
             << endl;
        MPI_Abort(MPI_COMM_WORLD, -7);
      }
    }
    MPI_Bcast(&position, 1, MPI_INT, special_state_pc_id[m], MPI_COMM_WORLD);
    special_state_index[m] = position;
  }

  delete[] exist_bool_for_pc;
}
