//NOTE: BONUS INCLUDED
//#include<bits/stdc++.h>
//#include "stdc++.h"
#include <stdio.h>
#include "Vector2D.cpp"
using namespace std;
const double eps=1e-6;
//epsilon: if the relative speed of the two spheres at the point of contact is greater than epsilon, then there is friction. This prevents oscillatory behaviour from friction due to a non-zero time step
double m1=0.0254,m2=0.0254;//masses
double r=0.03021;//radius of both spheres
double I1=0.4*m1*r*r,I2=0.4*m2*r*r;//calculated moments of inertia
double u=0;//mu: coefficient of friction
Vector2D p1,p2,p;//momentums of spheres and total momentum
double L1,L2,L;//angular momentums
double E1,E2,E;//kinetic energies
double w1=0,w2=0;//angular velocities
Vector2D s1(0,0.01),s2(0.1,0);//displacement wrt the fiducial point 2D
Vector2D v1(0,0),v2(-2,0);//velocities 2D
//Vector2D s1(0,0),s2(1,0);//displacement wrt the fiducial point 1D
//Vector2D v1(1,0),v2(-2,0);//velocities 1D
double k=1000;//spring constant
const double dt=0.000001;//time step
void step(){
    Vector2D a1(0,0),a2(0,0);//accelerations
    double alpha1=0,alpha2=0;//angular accelerations
    double dis=(s1-s2).mag();//distance between two centers
    double PE1=0,PE2=0;//potential energies
    if(dis<2*r){//check for contact
        double x=(2*r-dis)*0.5;//compression amount for each sphere
        double Ff=u*k*x;//magnitude of frictional force
        //acceleration calculated from force which points perpendicular from the intersection of the two circles
        a1=(s1-s2).norm()*k*x/m1;
        a2=(s2-s1).norm()*k*x/m2;
        //calculates potential energy of the spring
        PE1=PE2=0.5*k*x*x;
        Vector2D perp=(s1-s2).norm().perp();//unit vector along the intersection of the two circles
        //speed along the direction of perp (adjusted for rotation at contact point
        double ev1=v1.dot(perp)-r*w1;
        double ev2=v2.dot(perp)+r*w2;

        double v12=ev1-ev2;//relative speed of object1 with respect to object 2
        if(v12>eps){//friction in one direction
            a1=a1-perp*Ff/m1;//adds frictional forces
            a2=a2+perp*Ff/m2;
            alpha1=((s2-s1)*0.5).cross(perp*-Ff)/I1;//calculates the angular acceleration
            alpha2=((s1-s2)*0.5).cross(perp*Ff)/I2;
        }
        if(v12<-eps){//friction in the other direction
            //same as above
            a1=a1+perp*Ff/m1;
            a2=a2-perp*Ff/m2;
            alpha1=((s2-s1)*0.5).cross(perp*Ff)/I1;
            alpha2=((s1-s2)*0.5).cross(perp*-Ff)/I2;
        }
    }
    p1=v1*m1,p2=v2*m2,p=p1+p2;//calculates momentums
    //calculates angular momentums
    L1=s1.cross(p1)+I1*w1;
    L2=s2.cross(p2)+I2*w2;
    L=L1+L2;
    //calculates energies
    E1=0.5*m1*v1.dot(v1)+0.5*I1*w1*w1+PE1;
    E2=0.5*m2*v2.dot(v2)+0.5*I2*w2*w2+PE2;
    E=E1+E2;
    //adjusts angular velocity, velocity and displacement for next frame
    w1+=alpha1*dt;
    w2+=alpha2*dt;
    s1=s1+v1*dt+a1*dt*dt*0.5;
    s2=s2+v2*dt+a2*dt*dt*0.5;
    v1=v1+a1*dt;
    v2=v2+a2*dt;
}
int main(){
    //freopen("//Users//kevinwan//Desktop//sim1.csv","w",stdout);
    freopen("2D-NONMOVINGTARGETSimulationData-x-t.csv","w",stdout);//prepares the csv file
    printf("time,s1.x,s1.y,s2.x,s2.y,v1.x,v1.y,v2.x,v2.y,p1.x,p1.y,p2.x,p2.y,p.x,p.y,w1,w2,L1,L2,L,E1,E2,E\n");//writes the header
    double time = 0;
    for(int i=0;i<=500000;i++){
        step();//a single timestep
        //prints the relavant information
        if(i%1000==0)
            printf("%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",time,s1.x,s1.y,s2.x,s2.y,v1.x,v1.y,v2.x,v2.y,p1.x,p1.y,p2.x,p2.y,p.x,p.y,w1,w2,L1,L2,L,E1,E2,E);
        time += dt;
    }
}