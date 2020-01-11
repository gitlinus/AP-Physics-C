#include<math.h>
class Vector2D{
public:
    double x,y;
    Vector2D(double x,double y){this->x=x,this->y=y;}
    Vector2D(){this->x=this->y=0;}
    Vector2D polar(double r,double t){return {r*cos(t),r*sin(t)};}
    Vector2D operator+(Vector2D ot){return {this->x+ot.x,this->y+ot.y};}
    Vector2D operator-(Vector2D ot){return {this->x-ot.x,this->y-ot.y};}
    double dot(Vector2D ot){return this->x*ot.x+this->y*ot.y;}
    double cross(Vector2D ot){return this->x*ot.y-this->y*ot.x;}
    double mag(){return sqrt(this->x*this->x+this->y*this->y);}
    double dir(){return atan2(this->y,this->x);}
    Vector2D norm(){double t=dir();return {cos(t),sin(t)};}
    Vector2D operator*(double c){return {this->x*c,this->y*c};}
    Vector2D operator/(double c){return {this->x/c,this->y/c};}
    double proj(Vector2D ot){return this->dot(ot)/ot.mag();}
    Vector2D perp(){return {y,-x};}
};
