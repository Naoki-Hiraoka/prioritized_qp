#include <prioritized_qp_nasoq/PrioritizedQPNASOQSolver.h>
#include <iostream>

namespace prioritized_qp_nasoq{

  bool Task::isInitializeSolverRequired(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                        Eigen::SparseMatrix<double,Eigen::ColMajor>& A){
    return false;
  }

  bool Task::initializeSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                              Eigen::VectorXd& gradient,
                              Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                              Eigen::VectorXd& lowerBound,
                              Eigen::VectorXd& upperBound){
    this->H_ = H;
    this->q_ = gradient;
    this->A_.resize(0,H.cols());
    this->b_.resize(0);
    Eigen::SparseMatrix<double,Eigen::RowMajor> C_rowmajor(A.rows()*2,A.cols());
    C_rowmajor.topRows(A.rows()) = -A;
    C_rowmajor.bottomRows(A.rows()) = A;
    this->C_ = C_rowmajor;
    this->d_.resize(A.rows()*2);
    this->d_.head(A.rows()) = -lowerBound;
    this->d_.tail(A.rows()) = upperBound;
    return true;
  }
  bool Task::updateSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                          Eigen::VectorXd& gradient,
                          Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                          Eigen::VectorXd& lowerBound,
                          Eigen::VectorXd& upperBound){
    this->H_ = H;
    this->q_ = gradient;
    this->A_.resize(0,H.cols());
    this->b_.resize(0);
    Eigen::SparseMatrix<double,Eigen::RowMajor> C_rowmajor(A.rows()*2,A.cols());
    C_rowmajor.topRows(A.rows()) = -A;
    C_rowmajor.bottomRows(A.rows()) = A;
    this->C_ = C_rowmajor;
    this->d_.resize(A.rows()*2);
    this->d_.head(A.rows()) = -lowerBound;
    this->d_.tail(A.rows()) = upperBound;
    return true;
  }

  bool Task::solve(bool forceColdStart){
    nasoq::Nasoq nasoq(H_.rows(),H_.outerIndexPtr(),H_.innerIndexPtr(),H_.valuePtr(),q_.data(),
                       A_.rows(),A_.cols(),A_.outerIndexPtr(),A_.innerIndexPtr(),
                       A_.valuePtr(),b_.data(),
                       C_.rows(),C_.cols(),C_.outerIndexPtr(),C_.innerIndexPtr(),
                       C_.valuePtr(),d_.data());
    /// Define solver settings if provided
    nasoq.max_iter_nas = settings_.max_iter_nas; // 4000
    nasoq.diag_perturb = settings_.diag_perturb; //1e-9
    nasoq.zero_thresh = settings_.diag_perturb; //1e-9
    nasoq.eps_abs = settings_.eps; // 1e-3
    nasoq.eps_rel = settings_.eps_rel; // no use 1e-3
    nasoq.warm_start = 0;
    nasoq.batch_size = settings_.batch_size;// 1
    nasoq.scaling = settings_.scaling; // 0
    nasoq.inner_iter_ref = settings_.inner_iter_ref; //0
    nasoq.outer_iter_ref = settings_.outer_iter_ref;//0
    nasoq.max_iter=settings_.max_iter;//0
    nasoq.stop_tol = settings_.stop_tol;//1e-15
    nasoq.auto_reg_en = 0;
    if (settings_.nasoq_variant == "fixed")//fixed
      nasoq.variant = nasoq::nasoq_mode::Fixed;
    else if (settings_.nasoq_variant == "tuned")
      nasoq.variant = nasoq::nasoq_mode::Tuned;
    else if(settings_.nasoq_variant == "auto")
      nasoq.variant = nasoq::nasoq_mode::AUTO;
    else
      nasoq.variant= nasoq::nasoq_mode::PREDET;
    //nasoq.sol_name = "NASOQ-Fixed";

    /// Solve Problem
    int exitflag = nasoq.solve();

    /// Copy output
    x_ = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> >(nasoq.primal_vars,H_.rows(),1);

    return exitflag;
  }

  Eigen::VectorXd Task::getSolution(){
    return this->x_;
  }
};
