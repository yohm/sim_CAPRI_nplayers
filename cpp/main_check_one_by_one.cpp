#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <set>
#include <utility>
#include <deque>
#include <random>
#include "mpi.h"
#include "Strategy.hpp"
#include "MyLib.hpp"


struct Counts {
  uint64_t n_D_E;  // defensible & efficient
  uint64_t n_nD_E; // not defensible & efficient
  uint64_t n_D_nE; // defensible & not efficient
  uint64_t n_nD_nE; // not defensible & not efficient
  uint64_t n_Error;
  Counts() {
    n_D_E = 0;
    n_nD_E = 0;
    n_D_nE = 0;
    n_nD_nE = 0;
    n_Error = 0;
  }
  void Print(std::ostream& out) const {
    out << "D_E / D_nE / nD_E / nD_nE / n_Error : " << ToC(n_D_E) << " / " << ToC(n_D_nE) << " / " << ToC(n_nD_E) << " / " << ToC(n_nD_nE) << " / " << ToC(n_Error) << std::endl;
  }
  void Add(const Counts& x) {
    n_D_E += x.n_D_E;
    n_D_nE += x.n_D_nE;
    n_nD_E += x.n_nD_E;
    n_nD_nE += x.n_nD_nE;
    n_Error += x.n_Error;
  }
};

void CheckDFS(const Strategy& s, Counts& counter) {
  assert( s.NumU() == 0 );
  if( s.NumFixed() == 64 ) {
    bool d = s.IsDefensible2();
    bool e = s.IsEfficientTopo();
    bool l = s.IsEfficient();
    if( e && !l ) {
      l = s.IsEfficient(0.000001);  // check with a smaller threshold
      if(!l) {
        l = s.IsEfficient(0.0000001);  // check with an even smaller threshold
      }
    }
    if( e != l ) {
      std::cerr << s.ToString() << ' ' << d << ' ' << e << ' ' << l << std::endl;
      counter.n_Error++;
      return;
      // throw "must not happen";
    }
    if( d ) {
      if( e ) { counter.n_D_E++; }
      else { counter.n_D_nE++; }
    }
    else {
      if( e ) { counter.n_nD_E++; }
      else { counter.n_nD_nE++; }
    }
    return;
  }
  else {
    int idx = 0;
    for(; idx<64; idx++) {
      if(s.ActionAt(idx) == W || s.ActionAt(idx) == U) { break; }
    }
    for(int j=0; j<2; j++) {
      Strategy _s = s;
      _s.SetAction(idx, ((j==0)?C:D) );
      CheckDFS(_s, counter);
    }
  }
}

Counts CheckOneByOne(const Strategy& str) {
  assert( str.NumU() == 0 );
  Counts counter;
  CheckDFS(str, counter);
  return std::move(counter);
}

void testStrategy(const std::string& str, uint64_t D_E = 0, uint64_t D_nE = 0, uint64_t nD_E = 0, uint64_t nD_nE = 0) {
  Strategy s(str.c_str());
  auto c = CheckOneByOne(s);
  if(D_E+D_nE+nD_E+nD_nE > 0) {
    assert( c.n_D_E == D_E );
    assert( c.n_D_nE == D_nE );
    assert( c.n_nD_E == nD_E );
    assert( c.n_nD_nE == nD_nE );
  }
  c.Print(std::cout);
}

void test() {
  // testStrategy("ccdddcddccdcddcdccddddcdccccddddddcdccdddddddcdddcccddddcddcdddd", 1, 0);
  testStrategy("cdddcdcdccddddccccccdddddddccdddcdccdccdddcdddddcdcccdddccdccddd", 1);
  /*
  testStrategy("cdddcccdcdcdccdccccdddddcccdccddccddcd*cdccddcddddcd***c****ccdd", 0, 256);  // all inefficient
  testStrategy("ccddcccdccccddcdcdccddccdcddcccd*d*dccddddcdcc**c*ccddcc*c**cccd", 0, 256);  // all inefficient
  testStrategy("cdddddddddddccdcddcdddcddddddcddddddcdddcddd**cd**cd***cdc*cdcdd", 256, 0); // all efficient
  testStrategy("cdcdddcdddddccdcddcdddcdcddddcdddcdccdddcddd**cd**cd***cdc*cdcdd", 256, 0); // all efficient
  testStrategy("ccdd*cddc*ccdccdc*ddddccdcc**ccd*dc*ccdcdccddccdcccddd*cdccccccd", 0, 256); // all inefficient
  testStrategy("ccddddddcdccdccdccddddccdccccccddddcccdcdccddccd**cddd*cdc****cd", 0, 128); // inefficient
  testStrategy("ccddcd*dccccddcdccccddcddccc*dcdddc*ccc*ddcdddcdd*cddddddccdcddd", 0, 16, 0, 16); // inefficient
   */
}


int main(int argc, char** argv) {
#ifndef NDEBUG
  test();
  return 0;
#else
  MPI_Init(&argc, &argv);

  if( argc != 5 ) {
    std::cerr << "Error : invalid argument" << std::endl;
    std::cerr << "  Usage : " << argv[0] << " <in_format> <num_files> <out_format> <mode 0:all 1:random>" << std::endl;
    return 1;
  }

  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int num_procs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  const int N_FILES = std::strtol(argv[2],NULL,0);
  if(N_FILES <= 0) { throw std::runtime_error("invalid input"); }
  const int PROCS_PER_FILE = num_procs / N_FILES;
  char infile[256];
  sprintf(infile, argv[1], my_rank / PROCS_PER_FILE);
  std::cerr << "reading " << infile << " @ rank " << my_rank << std::endl;
  std::ifstream fin(infile);
  if( !fin.is_open() ) {
    std::cerr << "[Error] No input file " << infile << std::endl;
    throw std::runtime_error("no input file");
  }

  char outfile[256];
  sprintf(outfile, argv[3], my_rank);
  std::ofstream fout(outfile);

  const int mode = std::strtol(argv[4],NULL,0);
  std::mt19937 rnd( my_rank );

  Counts total;

  int count = 0;
  for( std::string line; fin >> line; count++) {
    if(count % 1000 == 0) {
      std::cerr << "step: " << count << " @ " << my_rank << "\t";
      total.Print(std::cerr);
    }

    if( count % PROCS_PER_FILE == my_rank%PROCS_PER_FILE ) {
      assert( line.size() == 64 );
      auto start = std::chrono::system_clock::now();
      if(mode > 0) {
        for(int i=0; i<64; i++) {
          if(line[i] == '_' || line[i] == '*') {
            line[i] = ((rnd() & 1) == 0) ? 'c' : 'd';
          }
        }
      }
      //std::cerr << "checking: " << line << std::endl;

      Strategy _str(line.c_str());

      Counts res = CheckOneByOne(_str);

      auto m1 = std::chrono::system_clock::now();
      double e1 = std::chrono::duration_cast<std::chrono::milliseconds>(m1-start).count();
      if(e1 > 3000.0) { std::cerr << "e1 > 3sec : " << line << std::endl; }

      fout << line << ' ' << res.n_D_E << ' ' << res.n_D_nE << ' ' << res.n_nD_E << ' ' << res.n_nD_nE << std::endl;
      total.Add(res);
    }
  }
  fout.close();

  Counts all_total;
  MPI_Reduce(&total.n_D_E, &all_total.n_D_E, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&total.n_D_nE, &all_total.n_D_nE, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&total.n_nD_E, &all_total.n_nD_E, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&total.n_nD_nE, &all_total.n_nD_nE, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&total.n_Error, &all_total.n_Error, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if(my_rank == 0) {
    all_total.Print(std::cerr);
  }

  MPI_Finalize();

  return 0;
#endif
}

