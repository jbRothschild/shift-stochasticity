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

using namespace std;
//using namespace arma;
using namespace Eigen;


//define birth and death rates
double bn (int n,int m, double params[6]) {//should probably include "[6]" in all of them so no/less junk gets in
    return params[0]*double(n);} // ! should all the inputs be doubles? or at least the output specifically?
double bm (int n,int m, double params[]) {
    return params[1]*m;}
double dn (int n,int m, double params[]) {
    return params[0]*n*(n + params[4]*m)/params[2];}
double dm (int n,int m, double params[]) {
    return params[1]*m*(params[5]*n + m)/params[3];}
double gg (int n,int m, double params[]) {
    return bn(n,m,params) + bm(n,m,params) + dn(n,m,params) + dm(n,m,params);}
double gn (int n,int m, double params[]) {
    return bm(n,m,params) + dn(n,m,params) + dm(n,m,params);}
double gm (int n,int m, double params[]) {
    return bn(n,m,params) + dn(n,m,params) + dm(n,m,params);}
double gnm (int n,int m, double params[]) {
    return dn(n,m,params) + dm(n,m,params);}

int vectorarg (int n,int m, int sumsize) {
    return int(m*sumsize-sumsize+n-1);}
double initconds (double params[6]) {//returns one double of K/(1+a), that is, where the fixed point is
    double tempini = params[2]/(1.+params[4]);
    if (tempini>1) {return tempini;}
    else {return 1.0;}
}

double fudgefactor = 1.0e1; //included to prevent overflow error - grow the matrix to shrink its inverse

SparseMatrix<double> makesparsetrix(int interspec, double params[6]){    //makes the sparse matrix
//    could also use Triplets (http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling)
    int msize = interspec*interspec;
    SparseMatrix<double> mat(msize,msize);
    cout << "[square] matrix size = " << mat.rows() << endl;

    mat.reserve(VectorXi::Constant(msize,5));//each column will have at most five elements
    for (int j=0; j<interspec-1; j++){
          for (int l=0; l<interspec-1; l++){
                    int x=interspec*j+l; int y=interspec*j+l;
                    mat.insert(x,y)=-gg(l+1,j+1,params);
                    mat.insert(x+1,y)=bn(l+1,j+1,params);
                    mat.insert(x,y+1)=dn(l+2,j+1,params);
                    mat.insert(x+interspec,y)=bm(l+1,j+1,params);
                    mat.insert(x,y+interspec)=dm(l+1,j+2,params);
    }};
    for (int j=0; j<interspec-1; j++){
                int x=interspec*j+interspec-1;
                int y=interspec*(interspec-1)+j;
                mat.insert(x,x)=-gn(interspec,j+1,params);
                mat.insert(x+interspec,x)=bm(interspec,j+1,params);
                mat.insert(x,x+interspec)=dm(interspec,j+2,params);
                mat.insert(y,y)=-gm(j+1,interspec,params);
    };
    mat.insert(msize-1,msize-1)=-gnm(interspec,interspec,params);

    mat.makeCompressed();

    mat = mat*fudgefactor;

    return mat;
}


int main()
{
    //define a file to write to
    ofstream filer; filer.open("kramers-tauvKa-20180719.txt");

    double nK = 7.;  //number of different K values, up to 8*nK
//    double nK = 8.;  //number of different K values, up to 2^nK
    int buffer = 7; //how many times bigger than K is the artificial reflecting boundary (3 is okay, 4 is better, 5 is enough+)
    double rconst = 1.0;
    cout << "let us do " << nK << " of these..." << endl;

    for(double aaa=-0.1; aaa<1.21; aaa+=0.05){
      for(double jay=1.;jay<nK+0.5;jay++){
        //define constants(x - c)
        double kay = double(8.+8.*jay);//from 4 to 32*4 took 67430s (19h)
//        double kay = double(pow(2.,jay));
        double constants[6] = {rconst, rconst, kay, kay, aaa, aaa};
        int sumsize = int(constants[2])*(buffer);
        int matsize = sumsize*sumsize;//or could use pow

        //make and fill matrix
        SparseMatrix<double> mat = makesparsetrix(sumsize,constants);
        SparseMatrix<double> matT = mat.transpose();
        mat.makeCompressed(); matT.makeCompressed();

        //solving Mx=b, really M_b T=-1, so define vectors - don't forget to transpose M!
        VectorXd tims;
        VectorXd negonesvec = VectorXd::Constant(matsize,1,double(-1.0));
        VectorXd fluxout = VectorXd::Zero(matsize,1);
        //??THIS DOES NOT CURRENTLY WORK for(int eye=0;eye<matsize;eye++){fluxout(i)+=dn(i,1,constants);fluxout(i*sumsize)+=dm(1,i,constants);};
        //or possibly we need to sum all of one row of M^-l //but this may need transposes, depending on whether I want row or column //Kramer
        VectorXd sumthese;//Kramer
        VectorXd pickaline = VectorXd::Zero(matsize,1);//Kramer
        int startat = int(initconds(constants)); pickaline(vectorarg(startat,startat,sumsize))=1.0;//Kramer
//        pickaline(vectorarg(kay,kay,sumsize))=1.0;//Kramer

//        //check matrix elements for correctness - verified
//        MatrixXf densemat; densemat=MatrixXf(mat);
//        for (int i=0;i<mat.rows();i++){
//            for (int j=0;j<9;j++){
////            for (int j=9;j<mat.cols();j++){
//                cout << densemat(i,j) << "\t";
//            };
//            cout << endl;
//        };

        //try solving it
        SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;//not clear to me what the ordering is or why it's important
        solver.analyzePattern(mat);//this is also not clear to me why it's used, and/or if it's exclusive with the somment above
        solver.factorize(mat);
        sumthese = solver.solve(pickaline);//Kramer
//        solver.analyzePattern(mat.transpose());//this is also not clear to me why it's used, and/or if it's exclusive with the somment above
//        solver.factorize(mat.transpose());
//        tims = solver.solve(negonesvec);
        SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver2;
        solver2.analyzePattern(matT);
        solver2.factorize(matT);
        tims = -solver2.solve(negonesvec);

        //let's see
//        cout << "For K const " << kay << " and a " << aaa << " you get an extinction time of " << tims(int(kay*sumsize-sumsize+kay-1)) << endl;
//        cout << "For K const " << kay << " and a " << aaa << " you get an extinction time of " << -sumthese.sum() << " times " << fudgefactor << endl;
//        cout << endl;
        filer << kay << "\t" << aaa << "\t" << -sumthese.sum() << "\n";
      }
    }

//    getchar();//to press enter - unnecessary
    filer.close();
    return 0;
}

