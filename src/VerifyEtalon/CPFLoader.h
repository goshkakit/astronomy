
class CPFLoader
{
private:

public:
	CPFLoader();
	~CPFLoader();

	void LoadCpffile( char *Fname );
	void GetPos( double jd, double secinday, double *outvec );
	int GetNoradId();
};