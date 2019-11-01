#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "mpi.h"

void adjtocsr(int** xadj,
              int** adjncy,
              int** vtxdist,
              int* nvtxs,
              int* nedges,
              std::string filepath,
              MPI_Comm comm){
    int id;
    int n_proc;
    int vtxs_per_proc;
    int total_vtxs;
    MPI_Comm_rank(comm, &id);
    MPI_Comm_size(comm, &n_proc);
    int* __vtxdist = new int[n_proc+1];
    if(id ==  n_proc-1) {
        std::ifstream adj_file;
        adj_file.open(filepath, std::ios::in);
        if(!adj_file) {
            std::cout<<"FILE CREATION FAILED"<<std::endl;
            MPI_Abort(comm, 1);
        }
        std::string s;
        getline(adj_file,s);
        total_vtxs = stoi(s);
        getline(adj_file,s);
        MPI_Bcast(&total_vtxs, 1, MPI_INT, n_proc-1, comm);
        vtxs_per_proc = total_vtxs/n_proc;
        int vtxs_read_sofar=0;
        std::vector<int> _xadj;
        std::vector<int> _adjncy;
        _xadj.push_back(0);
        int dest=0;
        while(getline(adj_file,s)){
            vtxs_read_sofar++;
            std::string delim = " ";
            auto start = 0U;
            auto end = s.find(delim);
            int count=0;
            while (end != std::string::npos)
            {
                start = end + delim.length();
                end = s.find(delim, start);
                _adjncy.push_back(stoi(s.substr(start,end-start)));
                count++;
            }
            _xadj.push_back(_xadj.back()+count);
            if(vtxs_read_sofar%vtxs_per_proc ==  0 && vtxs_read_sofar<total_vtxs - vtxs_per_proc){
                MPI_Send(_xadj.data(), _xadj.size(), MPI_INT, dest, 0, comm);
                int nlocal_edges = _adjncy.size();
                MPI_Send(&nlocal_edges, 1, MPI_INT, dest, 1, comm);
                MPI_Send(_adjncy.data(), _adjncy.size(), MPI_INT, dest++, 2, comm);
                _xadj.clear();
                _adjncy.clear();
                _xadj.push_back(0);
            } else if (vtxs_read_sofar == total_vtxs) {
                *nvtxs = _xadj.size()-1;
                *nedges = _adjncy.size();
                int* __xadj = new int[_xadj.size()];
                int* __adjncy = new int[_adjncy.size()];
                std::copy(_xadj.begin(), _xadj.end(), __xadj);
                std::copy(_adjncy.begin(), _adjncy.end(), __adjncy);
                _xadj.clear();
                _adjncy.clear();
                *adjncy = __adjncy;
                *xadj = __xadj;
            }
        }
        adj_file.close();
    } else {
        auto status = new MPI_Status();
        MPI_Bcast(&total_vtxs, 1, MPI_INT, n_proc-1, comm);
        *nvtxs = total_vtxs/n_proc;
        *xadj = new int[(*nvtxs)+1];
        MPI_Recv(*xadj, (*nvtxs)+1, MPI_INT, n_proc-1, 0, comm, status);
        MPI_Recv(nedges, 1, MPI_INT, n_proc-1, 1, comm, status);
        *adjncy = new int[*nedges];
        MPI_Recv(*adjncy, *nedges, MPI_INT, n_proc-1, 2, comm, status);
    }
    __vtxdist[0] = 0;
    vtxs_per_proc =  total_vtxs/n_proc;
    for(int i=1; i<n_proc; ++i) {
        __vtxdist[i] = i*vtxs_per_proc;
    }
    __vtxdist[n_proc] =  total_vtxs;
    *vtxdist = __vtxdist;
}
