#ifndef PRIORITIZEDQPSOLVER_H
#define PRIORITIZEDQPSOLVER_H

#include <OsqpEigen/OsqpEigen.h>
#include <Eigen/Dense>
#include <prioritized_qp_base/PrioritizedQPBaseSolver.h>

namespace prioritized_qp{
  class Task : public prioritized_qp_base::Task
  {
    /*
      Ax = b
      dl <= Cx <= du
     */
  public:
    // settingsの設定はユーザーが行うこと
    OsqpEigen::Solver& solver() { return solver_; }

    virtual bool isInitializeSolverRequired(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                      Eigen::SparseMatrix<double,Eigen::ColMajor>& A) override;
    virtual bool initializeSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                  Eigen::VectorXd& g,
                                  Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                                  Eigen::VectorXd& lowerBound,
                                  Eigen::VectorXd& upperBound) override;
    virtual bool updateSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                              Eigen::VectorXd& gradient,
                              Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                              Eigen::VectorXd& lowerBound,
                              Eigen::VectorXd& upperBound) override;
    virtual bool solve(bool forceColdStart=false)override;
    virtual Eigen::VectorXd getSolution()override;
  private:
    OsqpEigen::Solver solver_;
  };

};

#endif
