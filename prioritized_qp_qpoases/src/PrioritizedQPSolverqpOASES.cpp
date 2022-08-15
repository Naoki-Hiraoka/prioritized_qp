#include <prioritized_qp_qpoases/PrioritizedQPSolverqpOASES.h>

namespace prioritized_qp_qpoases{

  bool Task::isInitializeSolverRequired(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                        Eigen::SparseMatrix<double,Eigen::ColMajor>& A){
    return
      this->example_.getNV() != H.rows() ||
      this->example_.getNC() != A.rows();
  }

  bool Task::initializeSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                              Eigen::VectorXd& g,
                              Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                              Eigen::VectorXd& lowerBound,
                              Eigen::VectorXd& upperBound){
    this->example_ = qpOASES::SQProblem ( H.rows(),A.rows(), qpOASES::HST_UNKNOWN);
    this->H_ = H;
    this->g_ = g;
    this->A_ = A;
    this->lowerBound_ = lowerBound;
    this->upperBound_ = upperBound;
    return true;
  }
  bool Task::updateSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                          Eigen::VectorXd& g,
                          Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                          Eigen::VectorXd& lowerBound,
                          Eigen::VectorXd& upperBound){
    this->H_ = H;
    this->g_ = g;
    this->A_ = A;
    this->lowerBound_ = lowerBound;
    this->upperBound_ = upperBound;
    return true;
  }

  bool Task::solve(bool forceColdStart){
    int state_len = this->H_.rows();
    int inequality_len = this->A_.rows();

    qpOASES::real_t* H = new qpOASES::real_t[state_len*state_len];
    qpOASES::real_t* A = new qpOASES::real_t[inequality_len*state_len];
    qpOASES::real_t* g = new qpOASES::real_t[state_len];
    qpOASES::real_t* ub = new qpOASES::real_t[state_len];
    qpOASES::real_t* lb = new qpOASES::real_t[state_len];
    qpOASES::real_t* ubA = new qpOASES::real_t[inequality_len];
    qpOASES::real_t* lbA = new qpOASES::real_t[inequality_len];

    for (int k=0; k < this->H_.outerSize(); ++k){
      for (Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(this->H_,k); it; ++it){
        H[it.row()*state_len + it.col()] = it.value();
      }
    }
    for (size_t i = 0; i < state_len; i++) {
      g[i] = this->g_[i];
    }
    for (size_t i = 0; i < state_len; i++) {
      lb[i] = -1e20;
      ub[i] = 1e20;
    }
    for (int k=0; k < this->A_.outerSize(); ++k){
      for (Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(this->A_,k); it; ++it){
        A[it.row()*state_len + it.col()] = it.value();
      }
    }
    for (size_t i = 0; i < state_len; i++) {
      lbA[i] = this->lowerBound_[i];
      ubA[i] = this->upperBound_[i];
    }

    this->example_.setOptions(this->options_);
    qpOASES::returnValue _status;
    int nWSR = 1000;
    if(this->example_.isInitialised() && this->example_.isSolved() && !forceColdStart){
      _status = this->example_.hotstart( H,g,A,lb,ub,lbA,ubA, nWSR);
    }else{
      _status = this->example_.init( H,g,A,lb,ub,lbA,ubA, nWSR);
    }
    int status = qpOASES::getSimpleStatus(_status);

    delete[] H;
    delete[] A;
    delete[] g;
    delete[] ub;
    delete[] lb;
    delete[] ubA;
    delete[] lbA;

    return status==0;
  }

  Eigen::VectorXd Task::getSolution(){
    qpOASES::real_t* xOpt;
    this->example_.getPrimalSolution( xOpt );
    Eigen::VectorXd x(this->H_.rows());
    for(size_t i=0; i<this->H_.rows();i++){
      x[i]=xOpt[i];
    }
    return x;
  }
};
