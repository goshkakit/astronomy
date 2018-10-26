#include "mytypes.h"
#include "math.h"
#include <fstream>
#include "gaussJordan_3d.h"
//////////////////////////////////////////////////////////////////////////
// SCoordinate

void SCoordinate::move(double _x, double _y)
{
	x += _x;
	y += _y;
}

void SCoordinate::print()
{
	cout << "\n x: " << x << " y : " << y << endl;
}

bool SCoordinate::operator==(SCoordinate& aRight) const
{
	return x == aRight.x && y == aRight.y;
}
SCoordinate SCoordinate::operator +(SCoordinate &r){
	return SCoordinate(x + r.x, y + r.y);
}

SCoordinate SCoordinate::operator -(SCoordinate &r){
	return SCoordinate(x - r.x, y - r.y);
}
SCoordinate SCoordinate::operator *(double r){
	return SCoordinate(x*r, y*r);
}
SCoordinate SCoordinate::operator /(double r){
	return SCoordinate(x/r, y/r);
}


bool SCoordinate::isZero()
{
	return x == 0 && y == 0 ? true : false;
}

//////////////////////////////////////////////////////////////////////////
// S3DCoordinate

void S3DCoordinate::print()
{
	cout << "\n x: " << x << " y : " << y << " z : " << z << endl;
}

double S3DCoordinate::norm()
{
	return sqrt(x*x + y*y + z*z);
}

bool S3DCoordinate::operator==(S3DCoordinate& aRight) const
{
	return x == aRight.x && y == aRight.y && z == aRight.z;
}

/*S3DCoordinate S3DCoordinate::operator=(S3DCoordinate &r)
{
x = r.x;
y = r.y;
z = r.z;
return *this;
}*/

S3DCoordinate S3DCoordinate::operator=(const S3DCoordinate &r)
{
	x = r.x;
	y = r.y;
	z = r.z;
	return *this;
}

S3DCoordinate S3DCoordinate::operator+(S3DCoordinate &r)
{
	double x_,y_,z_;
	x_ = x + r.x;
	y_ = y + r.y;
	z_ = z + r.z;
	return S3DCoordinate(x_, y_, z_);
}
S3DCoordinate S3DCoordinate::operator-(S3DCoordinate &r)
{
	double x_,y_,z_;
	x_ = x - r.x;
	y_ = y - r.y;
	z_ = z - r.z;
	return S3DCoordinate(x_, y_, z_);
}
S3DCoordinate S3DCoordinate::operator*(double r)
{
	double x_,y_,z_;
	x_ = x*r;
	y_ = y*r;
	z_ = z*r;
	return S3DCoordinate(x_, y_, z_);
}
double S3DCoordinate::operator*(S3DCoordinate &r)
{
	double s;
	s = x*r.x + y*r.y + z*r.z;
	return s;
}

S3DCoordinate S3DCoordinate::operator/(double r)
{
	double x_,y_,z_;
	x_ = x/r;
	y_ = y/r;
	z_ = z/r;
	return S3DCoordinate(x_, y_, z_);
}

S3DCoordinate S3DCoordinate::operator^(S3DCoordinate &r){
	return S3DCoordinate(y*r.z-z*r.y, z*r.x-x*r.z, x*r.y-y*r.x);
}

S3DCoordinate S3DCoordinate::direction(){
	return this->operator *(1/this->norm());
}
SCoordinate SCoordinate::direction(){
	double s = sqrt(x*x+y*y);
	return SCoordinate(x/s, y/s);
}
double SCoordinate::length(){
	return sqrt(x*x+y*y);
}
//////////////////////////////////////////////////////////////////////////
// SRect

bool SRect::isVoid()
{
	return xLeft == xRight && yBottom == yTop;
}
int SRect::width()
{
	return abs(xRight - xLeft);
}
int SRect::height()
{
	return abs(yTop - yBottom);
}
SRect::SRect(const char *FileName){
	int centerX, centerY, w, h;
	ifstream fst;
	fst.open(FileName);
	if (fst.is_open()){
//		fst>>centerX>>centerY>>w>>h; 
		fst>>xLeft>>yTop>>xRight>>yBottom;
	}
	fst.close(); 
	fst.clear();
/*	xLeft = centerX - w/2;
	yTop = centerY-h/2;
	xRight = centerX + w/2;
	yBottom = centerY + h/2;
	*/
}

//////////////////////////////////////////////////////////////////////////
// SSize

int SSize::Area()
{
	return width*height;
}

//////////////////////////////////////////////////////////////////////////
// STriangle

STriangle::STriangle(int _n1, int _n2, int _n3)
{
	n1 = _n1;
	n2 = _n2;
	n3 = _n3;
}

SCTriangle::SCTriangle(S3DCoordinate _n1, S3DCoordinate _n2, S3DCoordinate _n3)
{
	n1 = _n1;
	n2 = _n2;
	n3 = _n3;
}

S3DCoordinate SCTriangle::centre()
{
	S3DCoordinate centre;
	centre.x = (n1.x + n2.x + n3.x)/3;
	centre.y = (n1.y + n2.y + n3.y)/3;
	centre.z = (n1.z + n2.z + n3.z)/3;
	return centre;
}

double SCTriangle::angle(S3DCoordinate aVertex)
{
	S3DCoordinate A, B, C;
	if(n1 == aVertex) {A = n1; B = n2; C = n3;}
	if(n2 == aVertex) {A = n2; B = n1; C = n3;}
	if(n3 == aVertex) {A = n3; B = n1; C = n2;}
	double a,b,c;
	a = distance3d(B,C); 
	b = distance3d(A,C); 
	c = distance3d(B,A);
	double cosa = (b*b+c*c-a*a)/(2*b*c);
	return acos(cosa);


}

//////////////////////////////////////////////////////////////////////////

double distance3d(S3DCoordinate aPoint1, S3DCoordinate aPoint2)
{
	double distance = pow((float)(aPoint1.x - aPoint2.x), 2) + pow((float)(aPoint1.y - aPoint2.y), 2) + pow((float)(aPoint1.z - aPoint2.z), 2);	
	distance = sqrt(distance);
	return distance;
}

double distance2d(SCoordinate aPoint1, SCoordinate aPoint2)
{
	double distance = sqrt(pow((float)(aPoint1.x - aPoint2.x), 2) + pow((float)(aPoint1.y - aPoint2.y), 2));	
	return distance;
}

S3DCoordinate operator-(const S3DCoordinate &r){
	return S3DCoordinate(-r.x, -r.y,-r.z);
}


S3DMatrix::S3DMatrix()
{ 
	memset(a, 0, sizeof(a));
}


S3DMatrix S3DMatrix::operator+(S3DMatrix r)
{
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			r.a[i][j]+=a[i][j];	
	return r;
}

S3DMatrix S3DMatrix::operator-(S3DMatrix r)
{
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			r.a[i][j]+=a[i][j];	
	return r;
}

S3DCoordinate S3DMatrix::operator*(S3DCoordinate r){
	return S3DCoordinate(	a[0][0]*r.x+a[0][1]*r.y+a[0][2]*r.z,
		a[1][0]*r.x+a[1][1]*r.y+a[1][2]*r.z,
		a[2][0]*r.x+a[2][1]*r.y+a[2][2]*r.z);
}

S3DMatrix S3DMatrix::operator*(S3DMatrix r){
	S3DMatrix res;
	for (int i=0; i<3; i++){
		for (int k=0; k<3; k++){
			for (int j=0; j<3; j++){
				res.a[i][k] += a[i][j]*r.a[j][k];
			}
		}
	}
	return res;
}

S3DMatrix S3DMatrix::operator*(double alfa){
	S3DMatrix Z;
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			Z.a[i][j] = a[i][j]*alfa;
		}
	}
	return Z;
}
S3DMatrix S3DMatrix::inverse(){
	S3DMatrix in = E();
/*	long double det=a[0][0]*(a[1][1]*a[2][2]-a[2][1]*a[1][2])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);

	in.a[0][0]=(a[1][1]*a[2][2]-a[2][1]*a[1][2])/det;
	in.a[0][1]=-(a[1][0]*a[2][2]-a[1][2]*a[2][0])/det;
	in.a[0][2]=(a[1][0]*a[2][1]-a[2][0]*a[1][1])/det;
	in.a[1][0]=-(a[0][1]*a[2][2]-a[0][2]*a[2][1])/det;
	in.a[1][1]=(a[0][0]*a[2][2]-a[0][2]*a[2][0])/det;
	in.a[1][2]=-(a[0][0]*a[2][1]-a[2][0]*a[0][1])/det;
	in.a[2][0]=(a[0][1]*a[1][2]-a[0][2]*a[1][1])/det;
	in.a[2][1]=-(a[0][0]*a[1][2]-a[1][0]*a[0][2])/det;
	in.a[2][2]=(a[0][0]*a[1][1]-a[1][0]*a[0][1])/det;
	*/

	gjelim3d(this, &in, 3, 3);
	return in;
}

CameraProjector::CameraProjector(const S3DMatrix _R, const S3DCoordinate _t) : centerCalculated(false){
	R = _R;
	t = _t;
}
CameraProjector::CameraProjector(const char* FileName) : centerCalculated(false){
	ifstream fst;
	fst.open(FileName);
	if (fst.is_open()){
		fst>>R.a[0][0]>>R.a[0][1]>>R.a[0][2]>>t.x;
		fst>>R.a[1][0]>>R.a[1][1]>>R.a[1][2]>>t.y;
		fst>>R.a[2][0]>>R.a[2][1]>>R.a[2][2]>>t.z;
	}
	fst.close(); 
	fst.clear();
}
CameraProjector::CameraProjector(const double a[3][4]) : centerCalculated(false){
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			R.a[i][j] = a[i][j];
		}
	}
	t.x=a[0][3];
	t.y=a[1][3];
	t.z=a[2][3];
}
SCoordinate CameraProjector::operator *(S3DCoordinate r){
	S3DCoordinate h = R*r+t;
	return  SCoordinate(h.x/h.z, h.y/h.z);
}


CameraProjector CameraProjector::operator*(double alfa){
	CameraProjector p(R, t);
	p.R.a[2][0]/=alfa;
	p.R.a[2][1]/=alfa;
	p.R.a[2][2]/=alfa;
	p.t.z/=alfa;
	return p;
}
CameraProjector CameraProjector::inverse(){
	S3DMatrix Rin=R.inverse();
	CameraProjector p(Rin, - (Rin*t));
	return p;
}
S3DCoordinate CameraProjector::getCenter(){
	if (centerCalculated)
		return center;
	else {
		CameraProjector inP = this->inverse();
		center = inP.R*S3DCoordinate(0,0,0)+inP.t;
		centerCalculated = true;
	}
	return center;
}

void CameraProjector::MoveArea(double dx, double dy){
	R.a[0][0]-=dx*R.a[2][0];
	R.a[0][1]-=dx*R.a[2][1];
	R.a[0][2]-=dx*R.a[2][2];
	t.x -=dx*t.z;

	R.a[1][0]-=dy*R.a[2][0];
	R.a[1][1]-=dy*R.a[2][1];
	R.a[1][2]-=dy*R.a[2][2];
	t.y -=dy*t.z;
}