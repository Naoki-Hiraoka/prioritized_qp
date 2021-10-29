#include <prioritized_qp/PrioritizedQPSolver.h>

namespace prioritized_qp{

  bool Task::isInitializeSolverRequired(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                        Eigen::SparseMatrix<double,Eigen::ColMajor>& A){
    return
      !this->solver_.isInitialized() ||
      this->solver_.workspace()->data->n != H.rows() ||
      this->solver_.workspace()->data->m != A.rows();
  }

  bool Task::initializeSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                              Eigen::VectorXd& gradient,
                              Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                              Eigen::VectorXd& lowerBound,
                              Eigen::VectorXd& upperBound){
    this->solver_.data()->clearHessianMatrix();
    this->solver_.data()->clearLinearConstraintsMatrix();
    this->solver_.clearSolver();

    this->solver_.data()->setNumberOfVariables(H.cols());
    this->solver_.data()->setNumberOfConstraints(A.rows());
    this->solver_.data()->setHessianMatrix(H);
    this->solver_.data()->setGradient(gradient);
    this->solver_.data()->setLinearConstraintsMatrix(A);
    this->solver_.data()->setLowerBound(lowerBound);
    this->solver_.data()->setUpperBound(upperBound);

    this->solver_.initSolver();
    return true;
  }
  bool Task::updateSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                          Eigen::VectorXd& gradient,
                          Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                          Eigen::VectorXd& lowerBound,
                          Eigen::VectorXd& upperBound){
    // OsqpEigenのsolver.data()はsetGradientとsetUpperBound, setLowerBound時にEigenの配列のポインタをそのまま保持する。solver.data()は、initSolver時にosqpにコピーされる.問題は、updateHessianMatrixやupdateLinearConstraintsMatrixの時に、配列の埋まっている部分が変化した場合にもinitSolverが呼ばれることである. このときにeigenの配列がデストラクトされていたらメモリ外にアクセスしてしまう)
    // solver.data()は、initSolver時にosqpにコピーされ, update関数はosqpにコピーされた方の値を直接操作する。結果,配列の埋まっている部分が変化した場合にinitSolverが呼ばれたときに、古いsolver.data()の値で初期化されてしまう。
    // これら2つのバグを避けるためには、solver.dataを毎回セットし直すしかない。これらの処理は、通常のupdate時には必要ないものである
    this->solver_.data()->clearHessianMatrix();
    this->solver_.data()->clearLinearConstraintsMatrix();
    this->solver_.data()->setHessianMatrix(H);
    this->solver_.data()->setGradient(gradient);
    this->solver_.data()->setLinearConstraintsMatrix(A);
    this->solver_.data()->setLowerBound(lowerBound);
    this->solver_.data()->setUpperBound(upperBound);

    this->solver_.updateHessianMatrix(Eigen::SparseMatrix<double,Eigen::ColMajor>(H)); // OsqpEigenの実装の都合上，ColMajorでないとサイレントにバグが起こる
    this->solver_.updateGradient(gradient);
    this->solver_.updateLinearConstraintsMatrix(Eigen::SparseMatrix<double,Eigen::ColMajor>(A)); // OsqpEigenの実装の都合上，ColMajorでないとサイレントにバグが起こる
    this->solver_.updateBounds(lowerBound,upperBound);//upperとlower同時にupdateしないと，一時的にupperがlowerを下回ってエラーになる
  }

  bool Task::solve(bool forceColdStart){
    if(forceColdStart) this->solver_.clearSolverVariables();
    return this->solver_.solve();
  }

  Eigen::VectorXd Task::getSolution(){
    return this->solver_.getSolution();
  }
};
