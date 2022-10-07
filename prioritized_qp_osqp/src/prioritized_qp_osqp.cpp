#include <prioritized_qp_osqp/prioritized_qp_osqp.h>

namespace prioritized_qp_osqp{

  bool Task::isInitializeSolverRequired(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                        Eigen::SparseMatrix<double,Eigen::ColMajor>& A){
    return
      !this->solver_.IsInitialized() ||
      this->instance_.num_variables() != H.rows() ||
      this->instance_.num_constraints() != A.rows();
  }

  bool Task::initializeSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                              Eigen::VectorXd& gradient,
                              Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                              Eigen::VectorXd& lowerBound,
                              Eigen::VectorXd& upperBound){
    this->instance_.objective_matrix = H;
    this->instance_.objective_vector = gradient;
    this->instance_.constraint_matrix = A;
    this->instance_.lower_bounds = lowerBound;
    this->instance_.upper_bounds = upperBound;

    if(!this->solver_.Init(this->instance_, this->settings_).ok()) return false;
    return this->solver_.IsInitialized();
  }
  bool Task::updateSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                          Eigen::VectorXd& gradient,
                          Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                          Eigen::VectorXd& lowerBound,
                          Eigen::VectorXd& upperBound){
    if(!this->solver_.UpdateObjectiveAndConstraintMatrices(H,A).ok()) return false;
    if(!this->solver_.SetObjectiveVector(gradient).ok()) return false;
    if(!this->solver_.SetBounds(lowerBound,upperBound).ok()) return false;

    return this->solver_.IsInitialized();
  }

  bool Task::solve(bool forceColdStart){
    if(forceColdStart) this->solver_.UpdateWarmStart(false);
    else this->solver_.UpdateWarmStart(true);
    return this->solver_.Solve() == osqp::OsqpExitCode::kOptimal;
  }

  Eigen::VectorXd Task::getSolution(){
    return this->solver_.primal_solution();
  }
};
