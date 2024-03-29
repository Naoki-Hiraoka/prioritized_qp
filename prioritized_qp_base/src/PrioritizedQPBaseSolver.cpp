#include <prioritized_qp_base/PrioritizedQPBaseSolver.h>
#include <iostream>

namespace prioritized_qp_base{
  bool solve(std::vector< std::shared_ptr<Task> >& tasks, Eigen::VectorXd& result, int debuglevel){
    if(tasks.size()==0) {
      std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] no task is provided" << std::endl;
      return false;
    }

    // 次元をチェック
    int dim = -1;
    for(size_t i=0;i<tasks.size();i++){
      if(dim<0) dim = tasks[i]->A().cols();
      if(dim != tasks[i]->A().cols()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (dim != tasks[i]->A().cols())" << std::endl;
        return false;
      }
      if(dim != tasks[i]->C().cols()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (dim != tasks[i]->C().cols())" << std::endl;
        return false;
      }
      if(tasks[i]->toSolve() && (dim != tasks[i]->w().size())){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->toSolve() && (dim != tasks[i]->w().size()))" << std::endl;
        return false;
      }

      if(tasks[i]->A().rows() != tasks[i]->b().size()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->A().rows() != tasks[i]->b().size())" << std::endl;
        return false;
      }
      if(tasks[i]->toSolve() && tasks[i]->A().rows() != tasks[i]->wa().size()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->toSolve() && tasks[i]->A().rows() != tasks[i]->wa().size())" << std::endl;
        return false;
      }
      if(tasks[i]->C().rows() != tasks[i]->dl().size() || tasks[i]->C().rows() != tasks[i]->du().size()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->C().rows() != tasks[i]->dl().size() || tasks[i]->C().rows() != tasks[i]->du().size())" << std::endl;
        return false;
      }
      if(tasks[i]->toSolve() && tasks[i]->C().rows() != tasks[i]->wc().size()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->toSolve() && tasks[i]->C().rows() != tasks[i]->wc().size())" << std::endl;
        return false;
      }
      if(tasks[i]->A_ext().cols() != tasks[i]->C_ext().cols()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->A_ext().cols() != tasks[i]->C_ext().cols())" << std::endl;
        return false;
      }
      if(tasks[i]->id_ext().size() != 0 && tasks[i]->A_ext().cols() != tasks[i]->id_ext().size()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->id_ext().size() != 0 && tasks[i]->A_ext().cols() != tasks[i]->id_ext().size())" << std::endl;
        return false;
      }
      if(tasks[i]->A_ext().cols() != 0 && tasks[i]->A_ext().rows() != tasks[i]->A().rows()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->A_ext().cols() != 0 && tasks[i]->A_ext().rows() != tasks[i]->A().rows())" << std::endl;
        return false;
      }
      if(tasks[i]->C_ext().cols() != 0 && tasks[i]->C_ext().rows() != tasks[i]->C().rows()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->C_ext().cols() != 0 && tasks[i]->C_ext().rows() != tasks[i]->C().rows())" << std::endl;
        return false;
      }
      if(tasks[i]->toSolve() && tasks[i]->A_ext().cols() != tasks[i]->w_ext().size()){
        std::cerr << "[prioritized_qp_base::PriroritizedQPSolver::solve] dimension mismatch (tasks[i]->toSolve() && tasks[i]->A_ext().cols() != tasks[i]->w_ext().size())" << std::endl;
        return false;
      }
    }

    if (dim==0){
      result = Eigen::VectorXd(0);
      return true;
    }

    Eigen::VectorXd solution = Eigen::VectorXd::Zero(dim);
    Eigen::SparseMatrix<double,Eigen::RowMajor> As(0,dim);
    Eigen::VectorXd lBs(0);
    Eigen::VectorXd uBs(0);
    Eigen::VectorXd w_exts(0);
    std::map<std::string, size_t> id_exts;

    for(size_t i=0;i<tasks.size();i++){
      Eigen::SparseMatrix<double,Eigen::RowMajor> taskA = tasks[i]->A();
      Eigen::SparseMatrix<double,Eigen::RowMajor> taskC = tasks[i]->C();
      {
        size_t additionalCols = 0;
        if(tasks[i]->id_ext().size() == 0) additionalCols = tasks[i]->A_ext().cols();
        else{
          for(size_t j=0;j<tasks[i]->A_ext().cols();j++){
            if(id_exts.find(tasks[i]->id_ext()[j]) == id_exts.end()){
              id_exts[tasks[i]->id_ext()[j]] = As.cols() + additionalCols;
              additionalCols++;
            }
          }
        }
        Eigen::SparseMatrix<double,Eigen::ColMajor> taskA_ColMajor(taskA.rows(), As.cols() + additionalCols);
        taskA_ColMajor.leftCols(taskA.cols()) = taskA;
        if(tasks[i]->id_ext().size() == 0){
          taskA_ColMajor.rightCols(tasks[i]->A_ext().cols()) = tasks[i]->A_ext();
        }else{
          Eigen::SparseMatrix<double,Eigen::ColMajor> taskA_ext_ColMajor(tasks[i]->A_ext());
          for(size_t j=0;j<tasks[i]->A_ext().cols();j++){
            taskA_ColMajor.col(id_exts[tasks[i]->id_ext()[j]]) = taskA_ext_ColMajor.col(j);
          }
        }
        taskA = taskA_ColMajor;
        Eigen::SparseMatrix<double,Eigen::ColMajor> taskC_ColMajor(taskC.rows(), As.cols() + additionalCols);
        taskC_ColMajor.leftCols(taskC.cols()) = taskC;
        if(tasks[i]->id_ext().size() == 0){
          taskC_ColMajor.rightCols(tasks[i]->C_ext().cols()) = tasks[i]->C_ext();
        }else{
          Eigen::SparseMatrix<double,Eigen::ColMajor> taskC_ext_ColMajor(tasks[i]->C_ext());
          for(size_t j=0;j<tasks[i]->C_ext().cols();j++){
            taskC_ColMajor.col(id_exts[tasks[i]->id_ext()[j]]) = taskC_ext_ColMajor.col(j);
          }
        }
        taskC = taskC_ColMajor;
        w_exts.conservativeResize(w_exts.size() + additionalCols);
        if(tasks[i]->id_ext().size() == 0){
          w_exts.tail(tasks[i]->w_ext().size()) = tasks[i]->w_ext();
        }else{
          for(size_t j=0;j<tasks[i]->w_ext().size();j++){
            w_exts[id_exts[tasks[i]->id_ext()[j]] - dim] = tasks[i]->w_ext()[j];
          }
        }
      }

      As.conservativeResize(As.rows()+taskA.rows(),taskA.cols());
      As.bottomRows(taskA.rows()) = taskA;
      lBs.conservativeResize(lBs.size()+tasks[i]->b().size());
      uBs.conservativeResize(uBs.size()+tasks[i]->b().size());
      lBs.tail(tasks[i]->b().size()) = tasks[i]->b();
      uBs.tail(tasks[i]->b().size()) = tasks[i]->b();
      As.conservativeResize(As.rows()+taskC.rows(),taskC.cols());
      As.bottomRows(taskC.rows()) = taskC;
      lBs.conservativeResize(lBs.size()+tasks[i]->dl().size());
      uBs.conservativeResize(uBs.size()+tasks[i]->du().size());
      lBs.tail(tasks[i]->dl().size()) = tasks[i]->dl();
      uBs.tail(tasks[i]->du().size()) = tasks[i]->du();

      if(!tasks[i]->toSolve() ||
         (tasks[i]->A().rows()==0 && tasks[i]->C().rows()==0)){
        if(debuglevel > 1){
          std::cerr << tasks[i]->name() << std::endl;
          std::cerr << "A" << std::endl;
          std::cerr << taskA << std::endl;
          std::cerr << "b" << std::endl;
          std::cerr << tasks[i]->b() << std::endl;
          std::cerr << "C" << std::endl;
          std::cerr << taskC << std::endl;
          std::cerr << "dl du" << std::endl;
          for(size_t j=0;j<taskC.rows();j++){
            std::cerr << tasks[i]->dl()[j] << " " << tasks[i]->du()[j] << std::endl;
          }
          std::cerr << "w" << std::endl;
          std::cerr << tasks[i]->w() << std::endl;
          std::cerr << "w_exts" << std::endl;
          std::cerr << w_exts << std::endl;
        }
      }else{
        int num = As.cols() + tasks[i]->A().rows() + tasks[i]->C().rows();

        /*
          min 0.5 xHx + gx
          s.t lB <= Ax <= cB
         */
        Eigen::SparseMatrix<double,Eigen::ColMajor> H(num,num);
        Eigen::SparseMatrix<double,Eigen::ColMajor> A(As.rows(),num);
        Eigen::VectorXd gradient = Eigen::VectorXd::Zero(num);
        Eigen::VectorXd upperBound = uBs;
        Eigen::VectorXd lowerBound = lBs;

        for(size_t j=0;j<dim;j++) H.insert(j,j) = tasks[i]->w()(j);
        for(size_t j=0;j<w_exts.size();j++) H.insert(dim+j,dim+j) = w_exts(j);
        for(size_t j=0;j<tasks[i]->wa().size();j++) H.insert(As.cols()+j,As.cols()+j) = tasks[i]->wa()(j);
        for(size_t j=0;j<tasks[i]->wc().size();j++) H.insert(As.cols()+tasks[i]->wa().size()+j,As.cols()+tasks[i]->wa().size()+j) = tasks[i]->wc()(j);
        A.leftCols(As.cols()) = Eigen::SparseMatrix<double,Eigen::ColMajor>(As);
        for(size_t j=0;j<tasks[i]->A().rows() + tasks[i]->C().rows();j++){
          A.insert(As.rows() - tasks[i]->A().rows() - tasks[i]->C().rows() + j, As.cols() + j) = 1;
        }

        if(tasks[i]->isInitializeSolverRequired(H,A)){
          tasks[i]->initializeSolver(H,gradient,A,lowerBound,upperBound);
        }else{
          if(!tasks[i]->updateSolver(H,gradient,A,lowerBound,upperBound)){
            tasks[i]->initializeSolver(H,gradient,A,lowerBound,upperBound);
          }
        }

        bool solved = tasks[i]->solve();
        if(!solved) {
          solved = tasks[i]->solve(true);
        }
        if(!solved) {
          if(debuglevel > 1){
            std::cerr << tasks[i]->name() << std::endl;
            std::cerr << "A" << std::endl;
            std::cerr << taskA << std::endl;
            std::cerr << "b wa" << std::endl;
            for(size_t j=0;j<taskA.rows();j++){
              std::cerr << tasks[i]->b()[j] << " " << tasks[i]->wa()[j] << std::endl;
            }
            std::cerr << "C" << std::endl;
            std::cerr << taskC << std::endl;
            std::cerr << "dl du wc" << std::endl;
            for(size_t j=0;j<taskC.rows();j++){
              std::cerr << tasks[i]->dl()[j] << " " << tasks[i]->du()[j] << " " << tasks[i]->wc()[j] <<std::endl;
            }
            std::cerr << "w" << std::endl;
            std::cerr << tasks[i]->w() << std::endl;
            std::cerr << "w_exts" << std::endl;
            std::cerr << w_exts << std::endl;
            std::cerr << "!solved" << std::endl;
          }

          return false;
        }

        solution = tasks[i]->getSolution().head(As.cols());

        Eigen::VectorXd this_b = taskA * solution;
        Eigen::VectorXd this_d = taskC * solution;
        for(size_t j=0;j<tasks[i]->A().rows();j++){
          if(this_b(j)>tasks[i]->b()(j)) uBs(uBs.rows() - tasks[i]->A().rows() - tasks[i]->C().rows() + j) = this_b(j);
          if(this_b(j)<tasks[i]->b()(j)) lBs(lBs.rows() - tasks[i]->A().rows() - tasks[i]->C().rows() + j) = this_b(j);
        }
        for(size_t j=0;j<tasks[i]->C().rows();j++){
          if(this_d(j)>tasks[i]->du()(j)) uBs(uBs.rows() - tasks[i]->C().rows() + j) = this_d(j);
          if(this_d(j)<tasks[i]->dl()(j)) lBs(lBs.rows() - tasks[i]->C().rows() + j) = this_d(j);
        }

        if(debuglevel > 1){
          std::cerr << tasks[i]->name() << std::endl;
          std::cerr << "solution" << std::endl;
          std::cerr << solution << std::endl;
          std::cerr << "A" << std::endl;
          std::cerr << taskA << std::endl;
          std::cerr << "b this_b wa" << std::endl;
          for(size_t j=0;j<taskA.rows();j++){
            std::cerr << tasks[i]->b()[j] << " " << this_b[j] << " " << tasks[i]->wa()[j] << std::endl;
          }
          std::cerr << "C" << std::endl;
          std::cerr << taskC << std::endl;
          std::cerr << "dl this_d du wc" << std::endl;
          for(size_t j=0;j<taskC.rows();j++){
            std::cerr << tasks[i]->dl()[j] << " " << this_d[j] << " " << tasks[i]->du()[j] << " "  << tasks[i]->wc()[j] <<std::endl;
          }
          std::cerr << "w" << std::endl;
          std::cerr << tasks[i]->w() << std::endl;
          std::cerr << "w_exts" << std::endl;
          std::cerr << w_exts << std::endl;
        }
      }
    }

    result = solution.head(dim);

    if(debuglevel > 1){
      std::cerr << "result" << std::endl;
      std::cerr << solution << std::endl;
      std::cerr << "As" << std::endl;
      std::cerr << As << std::endl;
      std::cerr << "lBs result_b uBs" << std::endl;
      Eigen::VectorXd result_b = As * solution;
      for(size_t i=0;i<As.rows();i++){
        std::cerr << lBs[i] << " " << result_b[i] << " " << uBs[i] << std::endl;
      }
    }
    return true;
  }

};
