#include <prioritized_qp_qpswift/PrioritizedQPQPSWIFTSolver.h>
#include <iostream>

namespace prioritized_qp_qpswift{

  bool Task::isInitializeSolverRequired(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                        Eigen::SparseMatrix<double,Eigen::ColMajor>& A){
    return true;
  }

  bool Task::initializeSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                              Eigen::VectorXd& gradient,
                              Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                              Eigen::VectorXd& lowerBound,
                              Eigen::VectorXd& upperBound){
    this->P_ = H;
    this->c_ = gradient;
    this->A_.resize(0,H.cols());
    this->b_.resize(0);
    Eigen::SparseMatrix<double,Eigen::RowMajor> G_rowmajor(A.rows()*2,A.cols());
    G_rowmajor.topRows(A.rows()) = -A;
    G_rowmajor.bottomRows(A.rows()) = A;
    this->G_ = G_rowmajor;
    this->h_.resize(A.rows()*2);
    this->h_.head(A.rows()) = -lowerBound;
    this->h_.tail(A.rows()) = upperBound;
    this->solver_.initialize(this->P_,this->c_,this->A_,this->b_,this->G_,this->h_);
    return true;
  }
  bool Task::updateSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                          Eigen::VectorXd& gradient,
                          Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                          Eigen::VectorXd& lowerBound,
                          Eigen::VectorXd& upperBound){
    this->P_ = H;
    this->c_ = gradient;
    this->A_.resize(0,H.cols());
    this->b_.resize(0);
    Eigen::SparseMatrix<double,Eigen::RowMajor> G_rowmajor(A.rows()*2,A.cols());
    G_rowmajor.topRows(A.rows()) = -A;
    G_rowmajor.bottomRows(A.rows()) = A;
    this->G_ = G_rowmajor;
    this->h_.resize(A.rows()*2);
    this->h_.head(A.rows()) = -lowerBound;
    this->h_.tail(A.rows()) = upperBound;
    this->solver_.initialize(this->P_,this->c_,this->A_,this->b_,this->G_,this->h_);
    return true;
  }

  bool Task::solve(bool forceColdStart){
    if(forceColdStart) this->solver_.initialize(this->P_,this->c_,this->A_,this->b_,this->G_,this->h_);
    return this->solver_.solve();
  }

  Eigen::VectorXd Task::getSolution(){
    Eigen::VectorXd ret;
    this->solver_.getSolution(ret);
    return ret;
  }
};
