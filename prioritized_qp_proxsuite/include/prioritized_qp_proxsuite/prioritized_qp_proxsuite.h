#ifndef PRIORITIZED_QP_PROXSUITE_H
#define PRIORITIZED_QP_PROXSUITE_H

#include <iostream>
#include <proxsuite/proxqp/sparse/sparse.hpp>
#include <prioritized_qp_base/PrioritizedQPBaseSolver.h>

namespace prioritized_qp_proxsuite{
  class Task : public prioritized_qp_base::Task
  {
    /*
      Ax = b
      dl <= Cx <= du
     */
  public:
    proxsuite::proxqp::Settings<double>& settings() { return settings_; }

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
    std::shared_ptr<proxsuite::proxqp::sparse::QP<double, int>> solver_;
    proxsuite::proxqp::Settings<double> settings_;
  };

};

#endif
