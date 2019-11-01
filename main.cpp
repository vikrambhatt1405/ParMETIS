#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include <random>
#include <parmetis.h>
#include <mpi.h>

void adjtocsr(int** xadj,
              int** adjncy,
              int** vtxdist,
              int* nvtxs,
              int* nedges,
              std::string filepath,
              MPI_Comm comm);
int main(int argc, char* argv[]) {
    int id;
    int p;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    idx_t * adjncy;
    idx_t * xadj;
    idx_t * vtxdist;
    idx_t * adjwgt;
    idx_t nvtxs;
    idx_t nedges;
    idx_t wgtflag=1, numflag =0 , ncon = 1, nparts = 2;
    idx_t options[] = {0};
    idx_t edgecut=0;
    idx_t *part;
    auto tpwgts = new real_t[ncon*nparts];
    for(int i=0; i<ncon*nparts; ++i){
        tpwgts[i] = 1/((float_t )nparts);
    }
    auto ubvec = new real_t[ncon];
    for(int i=0; i<ncon; ++i){
        ubvec[i] = 2;
    }
    MPI_Comm comm = MPI_COMM_WORLD;
    adjtocsr(&xadj, &adjncy, &vtxdist, &nvtxs, &nedges, "/home/vikrambhatt/Academics/Research/Data/adjlist_small.txt" , MPI_COMM_WORLD);
    std::cout<<"Processor: "<<id<<std::endl;
    std::cout<<"Vertices: "<<nvtxs<<std::endl;
    std::cout<<"Edges: "<<nedges<<std::endl;
    std::cout<<std::endl;
    std::random_device r;
    std::seed_seq      seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937       eng(seed); // a source of random data
    std::uniform_int_distribution<int> dist(1,10);
    std::vector<int> v(nedges);
    generate(v.begin(), v.end(), bind(dist, eng));
    adjwgt = new int[nedges];
    std::copy(v.begin(), v.end(), adjwgt);
    part = new idx_t[nedges];
    int rval = ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, NULL, adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec,
                                options, &edgecut, part, &comm);
    if(rval==METIS_OK){
        std::cout<<"OK!"<<std::endl;
    }
    MPI_Finalize();
    return 0;
}