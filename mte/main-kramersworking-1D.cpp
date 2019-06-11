//matrix-coupledlog 23 June 2015
//this program should make a matrix of the coupled logistic equation transistion rates
//then "invert it" to solve the backwards Kolmogorov equation and the Kramers solution
//double = bad, double = better; but maybe this is the limit? - overflow errors!
//K from 2 to 50, a from 0.1 to 1.0 took 2600s
//K from 2^0 to 2^7, buffer 4, a 0.0 to 1.0 took 3934s (mostly from the last set); problem with 2^8, bad allocation?
//K from 2^0 to 2^7, buffer 5, a 0.0 to 1.0 took 9539s (~90% from the last set)
//K = 2^8 simply doesn't work, crashes after 10 minutes (without a useful error message)
//update 2016.04.14
//kay = double(8.+8.*jay) up to jay=16, so up to K=136 but it crashed but anyway gave negatives (ie didn't work)

#include <iostream>
#include <fstream>
//#include <typeinfo> //for: typeid - unnecessary?
#include <cmath> //for: pow
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include <complex>//unnecessary?
//#include <algorithm>//for: sort
#include <iomanip> //for setprecision

using namespace std;
//using namespace arma;
using namespace Eigen;


//define birth and death rates
double bn (int n,int m, double params[4]) {//should probably include "[6]" in all of them so no/less junk gets in
    return params[0]*double(n)*(1.0+params[1] - params[3]*double(n)/params[2]);} // ! should all the inputs be doubles? or at least the output specifically?
double dn (int n,int m, double params[4]) {
    return params[0]*n*(params[1] + (1.0-params[3])*n/params[2]);}
double gg (int n,int m, double params[]) {
    return bn(n,m,params) + dn(n,m,params);}
double gn (int n,int m, double params[]) {
    return dn(n,m,params);}

int vectorarg (int n,int m, int sumsize) {
    return int(m*sumsize-sumsize+n-1);}
double initconds (double params[6]) {//returns one double of K/(1+a), that is, where the fixed point is
    double tempini = params[2]/(1.+params[4]);
    if (tempini>1) {return tempini;}
    else {return 1.0;}
}

//really, should use Triplets
SparseMatrix<double> makesparsetrix1D(int interspec, double params[6]){    //makes the sparse matrix
//    could also use Triplets (http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling)
    int msize = interspec;
    SparseMatrix<double> mat(msize,msize);
    cout << "[square] matrix size = " << mat.rows() << endl;

    mat.reserve(VectorXi::Constant(msize,5));//each column will have at most five elements
    for (int j=0; j<interspec-1; j++){
                    mat.insert(j,j)=-gg(j+1,j+1,params);
                    mat.insert(j+1,j)=bn(j+1,j+1,params);
                    mat.insert(j,j+1)=dn(j+2,j+1,params);
    };
    int x=interspec-1;
    mat.insert(x,x)=-gn(interspec,interspec,params);

    mat.makeCompressed();

    return mat;
}

SparseMatrix<double> makesparsetrix1Dtranspose(int interspec, double params[6]){    //makes the sparse matrix
//    could also use Triplets (http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling)
    int msize = interspec;
    SparseMatrix<double> mat(msize,msize);
    cout << "[square] matrix size = " << mat.rows() << endl;

    mat.reserve(VectorXi::Constant(msize,5));//each column will have at most five elements
    for (int j=0; j<interspec-1; j++){
                    mat.insert(j,j)=-gg(j+1,j+1,params);
                    mat.insert(j+1,j)=bn(j+1,j+1,params);
                    mat.insert(j,j+1)=dn(j+2,j+1,params);
    };
    int x=interspec-1;
    mat.insert(x,x)=-gn(interspec,interspec,params);

    //I bet I have to compress before I transpose - "It is worth noting that most of our wrappers to external libraries requires compressed matrices as inputs."
    mat.makeCompressed();
//    mat.transpose();

//    mat.makeCompressed();

    return mat.transpose();
}

SparseMatrix<double> makesparsetrix1Dtranspose3(int interspec, double params[6]){    //makes the sparse matrix
//    could also use Triplets (http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling)
    int msize = interspec;
    SparseMatrix<double> mat(msize,msize);
    cout << "[square] matrix size = " << mat.rows() << endl;

    mat.reserve(VectorXi::Constant(msize,5));//each column will have at most five elements
    for (int j=0; j<interspec-1; j++){
                    mat.insert(j,j)=-gg(j+1,j+1,params);
                    mat.insert(j,j+1)=bn(j+1,j+1,params);
                    mat.insert(j+1,j)=dn(j+1,j+2,params);
    };
    int x=interspec-1;
    mat.insert(x,x)=-gn(interspec,interspec,params);

    mat.makeCompressed();

    return mat;
}


int main(){
    //define a file to write to
    ofstream filer; filer.open("kramers-1D-delta1E0-20170922.txt");

    double nK = 32.;  //number of different K values, up to 8*nK
//    double nK = 8.;  //number of different K values, up to 2^nK
    int buffer = 4; //how many times bigger than K is the artificial reflecting boundary (3 is okay, 4 is better, 5 is enough+)
    double rconst = 1.0; double dconst = 1.0;
    cout << "let us do " << nK << " of these..." << endl;

    for(double aaa=0.0; aaa<1.01; aaa+=0.1){
      for(double jay=0.;jay<nK+0.5;jay++){
        //define constants(x - c)
        double qconst = aaa;
        double Kconst = double(8.+8.*jay);//from 4 to 32*4 took 67430s (19h)
//        double kay = double(pow(2.,jay));
        double constants[4] = {rconst, dconst, Kconst, qconst};
        int sumsize = int(constants[2])*(buffer);
        int matsize = sumsize;//or could use pow

        //make and fill matrix
        SparseMatrix<double> mat = makesparsetrix1D(sumsize,constants);
        SparseMatrix<double> matT = makesparsetrix1Dtranspose(sumsize,constants);
//        SparseMatrix<double> matT3 = makesparsetrix1Dtranspose3(sumsize,constants);
//        SparseMatrix<double> matT2 = mat.transpose();//moretests - this is the good one
//        mat.makeCompressed(); matT.makeCompressed();//this line should be unnecessary as it's already done in initialization

        //solving Mx=b, really M_b T=-1, so define vectors - don't forget to transpose M!
        VectorXd waysout = VectorXd::Zero(matsize,1); waysout(0)=-dn(1,1,constants);
        VectorXd negones = VectorXd::Constant(matsize,1,double(-1.0));

        VectorXd probs, times, timesalt;
//        VectorXd probs2, times2, timesalt2;//moretests
//        VectorXd probs3, times3, timesalt3;//moretests3
//        VectorXd negonesvec = VectorXd::Constant(matsize,1,-1.0);

        //or possibly we need to sum all of one row of M^-l //but this may need transposes, depending on whether I want row or column //Kramer
        VectorXd sumthese;//Kramer
        VectorXd temp1 = VectorXd::Zero(int(constants[2]),1);
        VectorXd temp3 = VectorXd::Zero(matsize-int(constants[2])-1,1);
        VectorXd temp2(1); temp2 << double(1.0);
        VectorXd temp4(int(constants[2])+1); temp4 << temp1,temp2;
        VectorXd pickaline(matsize); pickaline << temp4, temp3;
//        VectorXd pickaline = VectorXd::Zero(matsize,1);//Kramer
//        pickaline(int(constants[2]))=1.0;

        //try solving it
        SparseLU<SparseMatrix<double>/**/, COLAMDOrdering<int>/**/ > solver;//not clear to me what the ordering is or why it's important
        solver.analyzePattern(mat);//this is also not clear to me why it's used, and/or if it's exclusive with the somment above
        solver.factorize(mat);//analyzePattern() and factorize() are done together in compute()
//        if(solver.info()!=Success){cout<<"could not factorize";}
        sumthese = solver.solve(pickaline);//Kramer
//        if(solver.info()!=Success){cout<<"could not solve";}

        SparseLU<SparseMatrix<double>/**/, COLAMDOrdering<int>/**/ > solver1;
        solver1.analyzePattern(matT);
        solver1.factorize(matT);
//        if(solver1.info()!=Success){cout<<"could not factorize";}
        probs = solver1.solve(waysout);
        times = -solver1.solve(probs);
        timesalt = solver1.solve(negones);
//        if(solver1.info()!=Success){cout<<"could not solve";}

///* moretests */
////        matT2.makeCompressed();
////        matT2.makeCompressed();
//        SparseLU<SparseMatrix<double>/**/, COLAMDOrdering<int>/**/ > solver2;
//        solver2.analyzePattern(matT2);
//        solver2.factorize(matT2);
////        if(solver2.info()!=Success){cout<<"could not factorize";}
//        probs2 = solver2.solve(waysout);    //this works
//        times2 = -solver2.solve(probs2);    //this works
//        timesalt2 = solver2.solve(negones); //this works
////        if(solver2.info()!=Success){cout<<"could not solve";}

///* moretests3 */
//        SparseLU<SparseMatrix<double>/**/, COLAMDOrdering<int>/**/ > solver3;
//        solver3.analyzePattern(matT3);
//        solver3.factorize(matT3);
//        probs3 = solver3.solve(waysout);
//        times3 = -solver3.solve(probs3);
//        timesalt3 = solver3.solve(negones); //THIS WORKS?! - I was using solver2

        //let's see
//        cout << "For K const " << kay << " and a " << aaa << " you get an extinction time of " << tims(int(kay*sumsize-sumsize+kay-1)) << endl;
//        cout << "For K const " << kay << " and a " << aaa << " you get an extinction time of " << -sumthese.sum() << " times " << fudgefactor << endl;
//        cout << endl;
        filer << Kconst << "\t" << dconst << "\t" << qconst << "\t" << probs(int(constants[2])) << "\t" << times(int(constants[2])) << "\n";
//        filer << std::setprecision(8) << Kconst << "\t" << dconst << "\t" << qconst << "\t\t" << -sumthese.sum() << "\t"
//         << times(int(constants[2])) << "\t" << timesalt(int(constants[2])) << "\t"
//         << times2(int(constants[2])) << "\t" << timesalt2(int(constants[2])) << "\t"
//         << times3(int(constants[2])) << "\t" << timesalt3(int(constants[2])) << "\t\t"
//         << probs(int(constants[2])) << "\t" << probs2(int(constants[2])) << "\t" << probs3(int(constants[2])) << "\n";

//        cout << matT.isApprox(matT) << " " << mat.isApprox(matT.transpose()) << " "  << mat.isApprox(matT2.transpose()) << " "  << matT2.isApprox(mat.transpose()) << "\n";
      }
    }

//    getchar();//to press enter - unnecessary
    filer.close();
//    return 0;
}

