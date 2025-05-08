#include <prioritized_qp_proxsuite/prioritized_qp_proxsuite.h>

namespace prioritized_qp_proxsuite{

  bool Task::isInitializeSolverRequired(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                        Eigen::SparseMatrix<double,Eigen::ColMajor>& A){
    return
      !this->solver_ ||
      !this->solver_->work.internal.is_initialized ||
      this->solver_->model.dim != H.rows() ||
      this->solver_->model.n_in != A.rows();
  }

  bool Task::initializeSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                              Eigen::VectorXd& gradient,
                              Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                              Eigen::VectorXd& lowerBound,
                              Eigen::VectorXd& upperBound){

    
    this->solver_ = std::make_shared<proxsuite::proxqp::sparse::QP<double, int>>(H.rows(), 0, A.rows());
    this->solver_->settings = this->settings_;
    this->solver_->init(H,
			gradient,
			proxsuite::nullopt,
			proxsuite::nullopt,
			A,
			lowerBound,
			upperBound);

    return this->solver_->work.internal.is_initialized;
  }
  bool Task::updateSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                          Eigen::VectorXd& gradient,
                          Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                          Eigen::VectorXd& lowerBound,
                          Eigen::VectorXd& upperBound){
    this->solver_->update(H,
			  gradient,
			  proxsuite::nullopt,
			  proxsuite::nullopt,
			  A,
			  lowerBound,
			  upperBound);

    return this->solver_->work.internal.is_initialized;
  }

  bool Task::solve(bool forceColdStart){
    if(forceColdStart) this->solver_->settings.initial_guess = proxsuite::proxqp::InitialGuessStatus::NO_INITIAL_GUESS;
    else this->solver_->settings.initial_guess = proxsuite::proxqp::InitialGuessStatus::EQUALITY_CONSTRAINED_INITIAL_GUESS;
    this->solver_->solve();
    return true;
  }

  Eigen::VectorXd Task::getSolution(){
    return this->solver_->results.x;
  }
};
