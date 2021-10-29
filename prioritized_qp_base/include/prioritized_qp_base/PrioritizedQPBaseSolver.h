#ifndef PRIORITIZEDQPBASESOLVER_H
#define PRIORITIZEDQPBASESOLVER_H

#include <Eigen/Eigen>
#include <memory>

namespace prioritized_qp_base{
  class Task
  {
    /*
      Ax = b
      dl <= Cx <= du
     */
  public:
    std::string name() const { return name_; }
    std::string& name() { return name_; }

    // A, C のcol数は全タスクで共通であるべき
    // w, wa, wcは正であること
    Eigen::SparseMatrix<double,Eigen::RowMajor> A() const { return A_; }
    Eigen::SparseMatrix<double,Eigen::RowMajor>& A() { return A_; }

    Eigen::VectorXd b() const { return b_; }
    Eigen::VectorXd& b() { return b_; }

    Eigen::VectorXd wa() const { return wa_; }
    Eigen::VectorXd& wa() { return wa_; }

    Eigen::SparseMatrix<double,Eigen::RowMajor> C() const { return C_; }
    Eigen::SparseMatrix<double,Eigen::RowMajor>& C() { return C_; }

    Eigen::VectorXd dl() const { return dl_; }
    Eigen::VectorXd& dl() { return dl_; }
    Eigen::VectorXd du() const { return du_; }
    Eigen::VectorXd& du() { return du_; }

    Eigen::VectorXd wc() const { return wc_; }
    Eigen::VectorXd& wc() { return wc_; }

    Eigen::VectorXd w() const { return w_; }
    Eigen::VectorXd& w() { return w_; }

    // 追加の変数に対応. 最大値最小化等で用いる. 設定しなくてもよい
    std::vector<std::string> id_ext() const { return id_ext_;} //各追加変数のタスク間の対応をあらわす
    std::vector<std::string>& id_ext() { return id_ext_;}

    Eigen::SparseMatrix<double,Eigen::RowMajor> A_ext() const { return A_ext_; }
    Eigen::SparseMatrix<double,Eigen::RowMajor>& A_ext() { return A_ext_; }

    Eigen::SparseMatrix<double,Eigen::RowMajor> C_ext() const { return C_ext_; }
    Eigen::SparseMatrix<double,Eigen::RowMajor>& C_ext() { return C_ext_; }

    Eigen::VectorXd w_ext() const { return w_ext_; }
    Eigen::VectorXd& w_ext() { return w_ext_; }

    // Trueにすると解く
    bool toSolve() const {return toSolve_; }
    bool& toSolve() {return toSolve_; }

    /*
      min 0.5 xHx + gx
      s.t lB <= Ax <= cB
    */
    virtual bool isInitializeSolverRequired(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                      Eigen::SparseMatrix<double,Eigen::ColMajor>& A) =0;
    virtual bool initializeSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                                  Eigen::VectorXd& g,
                                  Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                                  Eigen::VectorXd& lowerBound,
                                  Eigen::VectorXd& upperBound) =0;
    virtual bool updateSolver(Eigen::SparseMatrix<double,Eigen::ColMajor>& H,
                              Eigen::VectorXd& g,
                              Eigen::SparseMatrix<double,Eigen::ColMajor>& A,
                              Eigen::VectorXd& lowerBound,
                              Eigen::VectorXd& upperBound) =0;
    virtual bool solve(bool forceColdStart=false)=0;
    virtual Eigen::VectorXd getSolution()=0;
  private:
    std::string name_;

    Eigen::SparseMatrix<double,Eigen::RowMajor> A_;
    Eigen::VectorXd b_;
    Eigen::VectorXd wa_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> C_;
    Eigen::VectorXd dl_;
    Eigen::VectorXd du_;
    Eigen::VectorXd wc_;

    Eigen::VectorXd w_;

    std::vector<std::string> id_ext_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> A_ext_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> C_ext_;
    Eigen::VectorXd w_ext_;

    bool toSolve_;
  };

  bool solve(std::vector< std::shared_ptr<Task> >& tasks, Eigen::VectorXd& result, int debuglevel=0);

};

#endif
