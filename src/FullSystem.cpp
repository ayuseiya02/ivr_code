#include "System.h"
#include "util.h"
using namespace std;
// advise: When you get lost in all these functions. Keep eye on
// fullsystem::fullsystem() and full_system::Quantum_evolution. Because
// they are functions do most of works and call other functions here.

// choice to only simulate detector
bool detector_only;

// choice to continue simulation of detector precoupling state: (Used to
// investigate decoherence of detector)
bool Detector_Continue_Simulation;

// energy window choice: We usually turn this on for large dof simulation. This
// will dramatically decrease computational cost for simulation
bool energy_window;

double *coupling_strength_to_mode0;
double Coupling_between_electronic_state;

// initialization of parameters and do some pre-coupling set up
full_system::full_system(string path1) {

  path = path1;
  d.path = path;

  // read parameter and time step from input.txt
  read_input_with_MPI();

  s.read_MPI(input, output, log);
  d.read_MPI(input, output, log, s.tlnum, s.tldim, path);
  d.construct_initial_state_MPI(input, output);

  // We construct vibrational state for each electronic state . vibrational
  // state is represented by their mode quanta , for example (1,0,0,0,2,0,1)
  // etc.
  compute_detector_matrix_size_MPI_sphere();

  // we construct anharmonic coupling between states, thus finish constructing
  // Hamiltonian.
  d.construct_dmatrix_MPI(input, output, log, dmat0, dmat1, vmode0, vmode1);
  if (my_id == 0) {
    cout << "Finish constructing Matrix" << endl;
  }
}

// Doing Quantum Simulation with SUR algorithm, parallelized version.
void full_system::Quantum_evolution() {
  // -----------------------------------------------------------------------------------
  // Now we construct our wavefunction /phi for our detector and full_system.
  // (For system it is already constructed in s.read())
  clock_t start_time, end_time, duration;

  start_time = clock();

  Evolve_full_system_MPI(); // Evolve the system.

  end_time = clock();
  duration = end_time - start_time;
  if (my_id == 0) {
    log << "The total run time for parallel computing is "
        << (double(duration) / CLOCKS_PER_SEC) / 60
        << " minutes  for simulation time  " << d.proptime[0] << endl;
    cout << "The total run time for parallel computing is "
         << (double(duration) / CLOCKS_PER_SEC) / 60
         << " minutes  for simulation time  " << d.proptime[0] << endl;
  }
}
