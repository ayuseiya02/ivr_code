//
// Created by phyzch on 6/23/20.
// Wrap up function to read data from input.txt and do output.
//
#include "../System.h"
#include "../util.h"

void full_system::read_input_with_MPI() {
  coupling_strength_to_mode0 = new double[2];

  if (my_id == 0) {
    double coupling_strength_to_mode_spin_up;
    double coupling_strength_to_mode_spin_down;

    input.open(path +
               "/input/InFile1.txt"); // information recorded in input.txt
    if (!input.is_open()) {
      cout << "THE INFILE FAILS TO OPEN!" << endl;
      log << "THE INFILE FAILS TO OPEN!" << endl;
      MPI_Abort(MPI_COMM_WORLD, -2); // Abort the process and return error code
                                     // -2. (input file can't open)
    }
    input >> energy_window >> detector_only >> Detector_Continue_Simulation;
    input >> Rmax >> d.V_intra >> d.a_intra >> d.detector_energy_window_size;

    // coupling strength of electronic state to mode 0 coordinate.
    input >> coupling_strength_to_mode_spin_up >>
        coupling_strength_to_mode_spin_down;

    coupling_strength_to_mode0[0] = coupling_strength_to_mode_spin_down;
    coupling_strength_to_mode0[1] = coupling_strength_to_mode_spin_up;

    input >> Coupling_between_electronic_state;

    // read time used for simulation.
    // delt: time step.
    // tstart: time start to turn on coupling between photon and detector. (No
    // use in current project) tmax: maximum time for simulation of photon +
    // detector. (No use in current project) tprint: time step to print
    // simulation result.
    input >> delt >> tstart >> tmax >> tprint;
    // check if input is valid
    if (!Detector_Continue_Simulation) {
      log.open(path +
               "log.txt"); // log is used to record information about program.
    } else {
      log.open(path + "log.txt", ios::app);
    }
    if (!log.is_open()) {
      cout << "THE LOG FILE FAILS TO OPEN" << endl;
      MPI_Abort(MPI_COMM_WORLD, -1); // Abort the process and return error code
                                     // -1.  (log file can't open)
    }
    log << "Number of Process running:  " << num_proc << endl;
    if (!Detector_Continue_Simulation) {
      // start from beginning.
      output.open(path + "/output/OutFile1.txt");
    } else {
      output.open(path + "/output/OutFile1.txt", ios::app);
    }
    if (!output.is_open()) {
      cout << "OUTPUT FILE FAILS TO OPEN." << endl;
      log << "OUTPUT FILE FAILS TO OPEN." << endl;
      MPI_Abort(MPI_COMM_WORLD, -3); // Abort the process and return error code
                                     // -3 (output file can't open).
    }

    if (delt <= 0 || tstart > tmax || delt > tstart) {
      cout << "Wrong time variable." << endl;
      log << "Wrong time variable."
          << "The simulation is cancelled." << endl;
      log.close();
      input.close();
      output.close();
      MPI_Abort(MPI_COMM_WORLD,
                -4); // Abort with error code -4: Wrong time variable.
    }
    if (!Detector_Continue_Simulation) {
      output << delt << " " << tstart << " " << tmax << " " << tprint << endl;
    }
  }
  // Broadcast  parameter from process 0  to all process.
  // function:  int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
  // int root, MPI_Comm comm)
  MPI_Bcast(&energy_window, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&detector_only, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Detector_Continue_Simulation, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

  MPI_Bcast(&Rmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&d.V_intra, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&d.detector_energy_window_size, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&coupling_strength_to_mode0[0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Coupling_between_electronic_state, 1, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);

  // Bcast delt tstart tmax tprint to other process.
  MPI_Bcast(&delt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tstart, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tprint, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // used to rescale the matrix element amplitude.
  cf = 0.0299792458 * delt * pi2;
  d.cf = cf;
}
