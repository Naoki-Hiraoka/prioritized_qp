#ifndef PRIORITIZEDQPNASOQSOLVER_H
#define PRIORITIZEDQPNASOQSOLVER_H

#include <nasoq/nasoq.h>
#include <prioritized_qp_base/PrioritizedQPBaseSolver.h>

namespace prioritized_qp_nasoq{
  class Task : public prioritized_qp_base::Task
  {
    /*
      Ax = b
      dl <= Cx <= du
     */
  public:
    // settingsの設定はユーザーが行うこと
    // https://github.com/sympiler/nasoq/blob/master/include/nasoq/nasoq.h
    nasoq::QPSettings& settings() { return settings_; }

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
    nasoq::QPSettings settings_;
    Eigen::SparseMatrix<double,Eigen::ColMajor> H_;
    Eigen::VectorXd q_;
    Eigen::SparseMatrix<double,Eigen::ColMajor> A_;
    Eigen::VectorXd b_;
    Eigen::SparseMatrix<double,Eigen::ColMajor> C_;
    Eigen::VectorXd d_;
    Eigen::VectorXd x_;
    ;
  };

};

#endif
