#include<math.h>
class Vector3D{
public:
    double x,y,z;
    Vector3D(double x,double y,double z){this->x=x,this->y=y,this->z=z;}
    Vector3D(){this->x=this->y=this->z=0;}
    //Vector3D polar(double r,double t,double p){return {r*cos(t),r*sin(t)};}
    Vector3D operator+(Vector3D ot){return {this->x+ot.x,this->y+ot.y,this->z+ot.z};}
    Vector3D operator-(Vector3D ot){return {this->x-ot.x,this->y-ot.y,this->z-ot.z};}
    double dot(Vector3D ot){return this->x*ot.x+this->y*ot.y+this->z*ot.z;}
    Vector3D cross(Vector3D ot){return {this->y*ot.z-this->z*ot.y,this->z*ot.x-this->x*ot.z,this->x*ot.y-this->y*ot.x};}
    double mag(){return sqrt(this->x*this->x+this->y*this->y);}
    //double dir(){return atan2(this->y,this->x);}
    Vector3D norm(){return *this/mag()}
    Vector3D operator*(double c){return {this->x*c,this->y*c,this->z*c};}
    Vector3D operator/(double c){return {this->x/c,this->y/c,this->z*c};}
    double proj(Vector3D ot){return this->dot(ot)/ot.mag();}
};