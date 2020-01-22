#include <iostream>
#include <fstream>

#include <math.h>
#include <vector>
#include <unistd.h>

using namespace std;

double f(vector<double> x){
    // f(x1,x2)=x1+x2
    return x[0]+x[1];
}

vector<double> d_f(vector<double> x){
    // df(x1,x2)=[1 1]
    return {1.0,1.0};
}

double g_1(vector<double> x){
    // g1(x1,x2)=x1^2+x2^2-1
    return pow(x[0],2)+pow(x[1],2)-1;
}

vector<double> d_g_1(vector<double> x){
    // dg1(x1,x2)=[2*x1 2*x2]
    return {2*x[0],2*x[1]};
}

double g_2(vector<double> x){
    //g2(x1,x2)=-x1
    return -x[0];
}

vector<double> d_g_2(vector<double> x){
    //dg2(x1,x2)=[-1 0]
    return {-1.0,0.0};
}

double F(vector<double> x,double meu){
    return f(x)-meu*(log(-g_1(x))+log(-g_2(x)));
}

vector<double> d_F(vector<double> x,double meu){
    vector<double> df=d_f(x);
    double dF_1=df[0]-meu*(((1/g_1(x))*d_g_1(x)[0])+((1/g_2(x))*d_g_2(x)[0]));
    double dF_2=df[1]-meu*(((1/g_1(x))*d_g_1(x)[1])+((1/g_2(x))*d_g_2(x)[1]));
    return {dF_1,dF_2};
}



int main()
{
    // initial value for x
    vector<double> x{0.5,0.5};
    // log barrier method's hyperparamter
    double meu=1.0;
    bool stuck=false;
    // line search method's hyperparamters
    double tol=0.01;
    double wolfe_inc=1.2;
    double wolfe_dec=0.5;
    double wolf_par=0.01;
    double alpha=1.0;
    vector<double> delta;
    // Helper variables
    double wolf_term;
    double Fn;
    vector<double> d_Fn;
    unsigned int usecs=1000000;
    vector<double> x_opt;
    // lagrangian lamda parameter for inequalities
    vector<double> lamda{-(meu/g_1(x)),-(meu/g_2(x))};
    ofstream log_file_outer("path_outer.dat", ios_base::out | ios_base::app );
    ofstream log_file_inner("path_inner.dat", ios_base::out | ios_base::app );



        do{
            x_opt=x;
            cout<<"New Loop"<<endl;
            log_file_outer << x_opt[0] << " " << x_opt[1] << " " << F(x_opt,meu) <<endl;
        while (true){
            Fn=F(x,meu);
            log_file_inner << x[0] << " " << x[1] << " " << Fn <<endl;
            d_Fn=d_F(x,meu);
            delta={-d_Fn[0],-d_Fn[1]};
            wolf_term=wolf_par*(d_Fn[0]*alpha*delta[0]+d_Fn[1]*alpha*delta[1]);
            double f_next=F({x[0]+alpha*delta[0],x[1]+alpha*delta[1]},meu);
            while ((f_next>F(x,meu)+wolf_term) || isnan(f_next)){
                alpha=alpha*wolfe_dec;
                f_next=F({x[0]+alpha*delta[0],x[1]+alpha*delta[1]},meu);
                if (isnan(f_next)){
                    cout << "Function Value is Undefined "<<endl;
                }
//                usleep(usecs);
            }

            x[0]+=alpha*delta[0];
            x[1]+=alpha*delta[1];
            alpha=alpha*wolfe_inc;

            cout << "NEW SOLUTION : ("<< x[0] << ","<< x[1]<<")"<<endl;
            // terminaltion rule for the inner loop
            if(sqrt(pow(alpha*delta[0],2)+pow(alpha*delta[1],2))<tol)break;
            }

            lamda={-(meu/g_1(x)),-(meu/g_2(x))};
            meu*=0.5;
            cout<<"NEW MEU : "<<meu<<endl;
    // termination rule for the outer loop
    }while (sqrt(pow(x[0]-x_opt[0],2)+pow(x[1]-x_opt[1],2))>=0.001);

    x_opt=x;
    cout << "NEW SOLUTION : ("<< x_opt[0] << ","<< x_opt[1]<<")"<<endl;
    cout << "Lamda : ("<< lamda[0] << ","<< lamda[1]<<")"<<endl;
    log_file_inner << x[0] << " " << x[1] << " " << F(x,meu) <<endl;
    log_file_outer << x_opt[0] << " " << x_opt[1] << " " << F(x_opt,meu) <<endl;


    return 0;
}
