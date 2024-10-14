#ifndef CLIQUE_STACK_H
#define CLIQUE_STACK_H

#include <stack>
#include <vector>

// Class to manage the stack of cliques.
//
// The stack is used as follows:
// - use init to initialize the stack. If the parameters are correct, there will
//   be no reallocation of stack during its operation. The parameters indicate
//   the maximum size of the stack and the workspace.
// - use setup, with the size of the clique that is being computed, to
//   initialize enough elements in the workspace. This also returns a pointer to
//   write to the workspace.
// - use getChild to obtain information about the child clique that is at the
//   top of the stack. This also returns a pointer to read the data of the child
//   clique.
// - use popChild when the top child clique is no longer needed, to move on to
//   the next one.
// - use pushWork, with the ID of supernode that is being processed, when all
//   children clique have been assembled, the dense factorization is completed,
//   and the clique in the workspace is ready to be pushed onto the stack.

class CliqueStack {
  std::vector<double> stack_;
  std::vector<double> work_;
  int worksize_{};

  // pairs (sn, size) of supernodes that got pushed
  std::stack<std::pair<int, int>> sn_pushed_{};

  // entry corresponding to top of stack
  int top_{};

 public:
  void init(int stack_size, int max_clique_size);
  double* setup(int clique_size);
  const double* getChild(int& child_sn) const;
  void popChild();
  void pushWork(int sn);
};

// Example: sn 7 has children 5,6 and the stack looks like this:
//
// | - - - | - - - - - - - | - | - - - - - - - ...
// |   4   |       5       | 6 | t
// | - - - | - - - - - - - | - | - - - - - - - ...
//
// where t points to the top of the stack.
// The information about the supernodes pushed (in sn_pushed_) is:
// (sn 4, size 7), (sn 5, size 15), (sn 6, size 3)
//
// The clique of sn 7 is constructed in the workspace. After a child clique is
// assembled, it is popped from the stack, by moving t back. After all children
// are assembled and the dense factorization is done, the clique in workspace
// can be pushed onto the stack.
//
// At the end, the stack looks like this:
//
// | - - - | - - - - - - - - - - - | - - - - - ...
// |   4   |           7           | t
// | - - - | - - - - - - - - - - - | - - - - - ...
//
// The information in sn_pushed_ is:
// (sn 4, size 7), (sn 7, size 23)

#endif