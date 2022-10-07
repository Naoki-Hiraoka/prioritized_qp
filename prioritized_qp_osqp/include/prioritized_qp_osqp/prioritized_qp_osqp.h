#ifndef PRIORITIZED_QP_OSQP_H
#define PRIORITIZED_QP_OSQP_H

#include <osqp++.h>
#include <prioritized_qp_base/PrioritizedQPBaseSolver.h>

namespace prioritized_qp_osqp{
  class Task : public prioritized_qp_base::Task
  {
    /*
      Ax = b
      dl <= Cx <= du
     */
  public:
    // settingsの設定はユーザーが行うこと (initializesolver前にやらないと反省されない)
    osqp::OsqpSettings& settings() { return settings_; }

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
    osqp::OsqpInstance instance_;
    osqp::OsqpSolver solver_;
    osqp::OsqpSettings settings_;
  };

};

#endif
