#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#define Fsun  1367.5           // sun constant
#define sigma   5.67E-8         //stefan-boltzmann constant
#define h   6.626007015E-34     //plank constant
//q   0.29+0.684*0.15    //phase integral  = 0.29+0.684*G( =0.15)
#define epsi   0.9              //radiance epsilon
#define kB   1.38064852E-23     //boltzmann constant  j/k
#define cl   299792458.0        //lightspeed m/s
#define au   149597870700.0
#define pi   3.1415926535

double flux(double eta,double D,double delta,double dd,double alpha,double lamda,double A,int Ndd){
    double summ = 0;
    double a = alpha-pi/2;
    double b = pi/2;
    double c = -pi/2 ;         
    double d = pi/2;
    double Tfit = pow((1-A)*Fsun/(eta*epsi*sigma*dd*dd),0.25);
    double dxdy = ((b-a)*(d-c)/(Ndd*Ndd)); 
    for (int j=0;j<Ndd;j++){ 
        for (int i=0;i<Ndd;i++){
            double phii = (c+(d-c)*j/Ndd);
            double thei = (a+(b-a)*i/Ndd);
            double T = Tfit*pow(cos(thei),0.25)*pow(cos(phii),0.25);
            double integral = cos(alpha)*cos(alpha)*cos(alpha-thei)/(exp(h*cl/(lamda*kB*T))-1);
            summ = summ + dxdy*integral;
        }
    }
    double F = (epsi*pow(D,2)*pi*h*cl*cl)*summ/(2*pow(delta,2)*pow(lamda,5));
    F = F*(lamda*1e-6)*(lamda*1e-6)/cl*1e29;
    return F;
}
int main(int argc,char*argv[]){
    /*for (int i=1;i<argc;i++){
        printf("%d ",jiechen(atoi(argv[i])));
    }
    //printf("\nargc=%d,argv=%c\n",argc,*argv[argc-1]);*/
    double eta,D,delta,d,alpha,lamda,A;
    int Ndd;
    eta = atof(argv[1]);
    D = atof(argv[2]);
    delta = atof(argv[3]);
    d = atof(argv[4]);
    alpha= atof(argv[5]);
    lamda = atof(argv[6]);
    A = atof(argv[7]);
    Ndd = atof(argv[8]);
    printf("%f",flux(eta,D,delta,d,alpha,lamda,A,Ndd));
    return 0;
}