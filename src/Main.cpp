#include "System.h"
#include "util.h"

using namespace std;

namespace fs = std::experimental::filesystem;

void check_and_create_file(string parent_path, string path);

// Information for MPI program
int my_id;
int num_proc;

bool load_sampling_state = true;

int main(int argc, char *argv[]) {
  srand(time(0));
  // folder path to read input.txt and write output file.
  string path = "/Users/ayuseiya/research/ivr/code/src/Data";

  // Hamiltonian H = H_{0} + V
  // We have two options to construct anharmonic coupling part of Hamiltonian :
  // V , one is construct Hamiltonian using scaling law, based on paper: similar
  // to eq.4 ,5 in https://www.pnas.org/content/pnas/95/11/5960.full.pdf another
  // is read anharmonic coupling  V^{(3)} q_{1} q_{2} q_{3} from files.
  // cvpt_path is used for this.
  int Filenumber = 1;

  // MPI Command
  MPI_Init(&argc, &argv);
  // num_proc:  number of process run in parallel.
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  // my_id : id of current process.  my_id and num_proc will be global variable
  // in program.
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  for (int i = 0; i < Filenumber; i++) {

    { // the parenthese here let destructor called after we use this instance.
      // pay attention to destructor to avoid memory leak when we do 1000 or
      // more cases.
      full_system photon_entangled_system(
          path); // set parameter and construct Hamiltonian.
      photon_entangled_system
          .Quantum_evolution(); // creat initial state (or read from
                                // file(cvpt_path)  ). Then use SUR algorithm
                                // Evolve Schrodinger equation.
    }
  }
  MPI_Finalize();
}
