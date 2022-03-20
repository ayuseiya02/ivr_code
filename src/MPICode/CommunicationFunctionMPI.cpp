//
// Created by phyzch on 6/26/20.
//
#include "../System.h"
#include "../util.h"
using namespace std;

void convert_dv(const vector<vector<int>> &vec_2d, vector<int> &vec_1d,
                vector<int> &displacement, vector<int> &element_size) {
  // copy 2d vector data into 1d vector.
  // MPI can not work with 2 dimension vector<vector<int>>. we don't know the
  // size of vector<int> because it has container and many other part.
  // displacement will record the beginning of these buffer in vec_2d (in case
  // we have to send vec_2d with element vector with different size element size
  // will record size of element in 2d vector.
  int i, j;
  int loc = 0;
  int size_1 = vec_2d.size();
  int size_2;
  for (i = 0; i < size_1; i++) {
    size_2 = vec_2d[i].size();
    for (j = 0; j < size_2; j++) {
      vec_1d.push_back(vec_2d[i][j]);
    }
    displacement.push_back(loc);
    element_size.push_back(size_2);
    loc = loc + size_2;
  }
}

// this function is used in construct_dv_dirow_dicol_dmatrix_MPI.
void detector::Scatter_dirow_dicol_dmatrix(vector<double> *dmat_all,
                                           vector<int> *dirow_data,
                                           vector<int> *dicol_data,
                                           int **vector_size,
                                           int **vector_displs, ofstream &log) {
  // dmat_all , dirow_data, dicol_data contain data in all process. This
  // information is stored in master process (process id =0) vector_size
  // indicate size for each process. vector_displacement indicate displacement
  // for data in each process.
  int i;
  int v_size;
  // use MPI_Scatterv to scatter dmat to other process.
  int m;
  for (m = 0; m < 1; m++) {
    // tell each process the size of vector they will receive.
    MPI_Scatter(&vector_size[m][0], 1, MPI_INT, &v_size, 1, MPI_INT, 0,
                MPI_COMM_WORLD);
    // reserve space:
    dmat[m].reserve(v_size); // reserve space to receive buffer from Scatterv
    dirow[m].reserve(v_size);
    dicol[m].reserve(v_size);
    double *temp_dmat = new double[v_size];
    int *temp_dirow = new int[v_size];
    int *temp_dicol = new int[v_size];
    // use MPI_Scatterv to scatter dmat. for send buffer for std::vector, it
    // should be vector.data().
    MPI_Scatterv((void *)(dmat_all[m]).data(), vector_size[m], vector_displs[m],
                 MPI_DOUBLE, &temp_dmat[0], v_size, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    // for dirow, dicol.
    MPI_Scatterv((void *)(dirow_data[m]).data(), vector_size[m],
                 vector_displs[m], MPI_INT, &temp_dirow[0], v_size, MPI_INT, 0,
                 MPI_COMM_WORLD);
    MPI_Scatterv((void *)(dicol_data[m]).data(), vector_size[m],
                 vector_displs[m], MPI_INT, &temp_dicol[0], v_size, MPI_INT, 0,
                 MPI_COMM_WORLD);
    for (i = 0; i < v_size; i++) {
      dmat[m].push_back(temp_dmat[i]);
      dirow[m].push_back(temp_dirow[i]);
      dicol[m].push_back(temp_dicol[i]);
    };
    delete[] temp_dirow;
    delete[] temp_dicol;
    delete[] temp_dmat;
  }
}

// We will use this function to scatter dv to other process, used in
// construct_dirow_dicol_dv and load_detector_Hamiltonian()
void detector::Scatter_dv(vector<int> &total_mat_num) {
  // The number of state is evenly distributed for detector vmode in our
  // program, only know total_mat_num we could infer the size for each process
  // and corresponding displacement.
  int i, j;
  int *vector_2d_displs,
      *vector_2d_size; // vector_displs: displacement of vector, vector_size:
                       // size of vector to be sent. (For MPI_Scatterv)
  int vsize;
  int vsize_2d; // size of 2d vector to receive.
  //---------------------
  // for convert_dv() function
  vector<int> vmode_1d;      // 1 dimensional vmode vector
  vector<int> displacement;  // displacement of each 2d element (vmode)
  vector<int> element_size;  // size of each 2d element == moddim
  vector<int> receive_vmode; // vector to receive vmode_1d
  // ------------------
  vector_2d_displs =
      new int[num_proc]; // for sending 2d vector we need extra displacement.
  vector_2d_size =
      new int[num_proc]; // foor sending 2d vector we need new size of buffer.
  int m, p;
  for (m = 0; m < 1; m++) {
    if (my_id == 0) {
      //--------------------------------------
      vmode_1d.clear();
      displacement.clear();
      element_size.clear();
      receive_vmode.clear();
      // convert 2d vector into 1d vector to send.
      convert_dv(dv_all[m], vmode_1d, displacement, element_size);
      //----------------------------------------------------------
      // size to send for 1d array.
      vsize = total_mat_num[m] / num_proc;
      // vector_2d_displacement is displacement of 2d vector in corresponding 1d
      // vector vmode_1d. (1d vector is to be sent). vector_2d_size is size of
      // each process will receive for 2d vector
      for (p = 0; p < num_proc; p++) {
        vector_2d_displs[p] =
            displacement[vsize * p]; // displacement record the beginning of all
                                     // 1d vector element in 2d vector.
        // displacement[vector_displs[p]] contain beginning of buffer to send in
        // 1d vector : vmode_1d
        if (p == 0)
          ;
        else {
          vector_2d_size[p - 1] = vector_2d_displs[p] - vector_2d_displs[p - 1];
        }
      }
      vector_2d_size[num_proc - 1] =
          vmode_1d.size() - vector_2d_displs[num_proc - 1];
    } else {
      receive_vmode.clear();
      vsize = total_mat_num[m] / num_proc;
      if (my_id == num_proc - 1) {
        vsize =
            total_mat_num[m] - (total_mat_num[m] / num_proc) * (num_proc - 1);
      }
    }
    MPI_Scatter(&vector_2d_size[0], 1, MPI_INT, &vsize_2d, 1, MPI_INT, 0,
                MPI_COMM_WORLD);
    receive_vmode.reserve(vsize_2d);
    MPI_Scatterv((void *)vmode_1d.data(), vector_2d_size, vector_2d_displs,
                 MPI_INT, (void *)receive_vmode.data(), vsize_2d, MPI_INT, 0,
                 MPI_COMM_WORLD);
    // parse 1d array: vsize_2d to 2d dv.
    int index = 0;
    for (i = 0; i < vsize; i++) {
      vector<int> temp_vmode;
      for (j = 0; j < nmodes[m]; j++) {
        temp_vmode.push_back(receive_vmode[index]);
        index = index + 1;
      }
      dv[m].push_back(temp_vmode);
    }
  }
  delete[] vector_2d_size;
  delete[] vector_2d_displs;
}
