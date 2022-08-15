#ifndef PRIORITIZEDQPSOLVERQPOASES_H
#define PRIORITIZEDQPSOLVERQPOASES_H

#include <Eigen/Dense>
#include <prioritized_qp_base/PrioritizedQPBaseSolver.h>
#include <qpOASES.hpp>

namespace prioritized_qp_qpoases{
  class Task : public prioritized_qp_base::Task
  {
    /*
      Ax = b
      dl <= Cx <= du
     */
  public:
    // optionsの設定はユーザーが行うこと
    qpOASES::Options& options() { return options_; }

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
    qpOASES::Options options_;
    qpOASES::SQProblem example_;
    Eigen::SparseMatrix<double,Eigen::ColMajor> H_;
    Eigen::VectorXd g_;
    Eigen::SparseMatrix<double,Eigen::ColMajor> A_;
    Eigen::VectorXd lowerBound_;
    Eigen::VectorXd upperBound_;
  };

};

#endif
