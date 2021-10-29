#ifndef PRIORITIZEDQPNASOQSOLVER_H
#define PRIORITIZEDQPNASOQSOLVER_H

#include <qpswifteigen/qpswifteigen.h>
#include <prioritized_qp_base/PrioritizedQPBaseSolver.h>

namespace prioritized_qp_qpswift{
  class Task : public prioritized_qp_base::Task
  {
  public:
    // settingsの設定はユーザーが行うこと
    qpswifteigen::solver& solver() { return solver_; }

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
    qpswifteigen::solver solver_;
    Eigen::SparseMatrix<double,Eigen::ColMajor> P_;
    Eigen::VectorXd c_;
    Eigen::SparseMatrix<double,Eigen::ColMajor> A_;
    Eigen::VectorXd b_;
    Eigen::SparseMatrix<double,Eigen::ColMajor> G_;
    Eigen::VectorXd h_;
  };

};

#endif
