//-------------------------------------------------------------------------------------------------//
// cpu pointer
//-------------------------------------------------------------------------------------------------//
#ifndef _CPU_PTR_H_
#define _CPU_PTR_H_
#include <stddef.h>

template < typename T > class CPU_Ptr
{
public:

	CPU_Ptr() : ptr(0) 
	{
		size = 0;
    }

	CPU_Ptr(size_t s) : ptr(0)
	{
		size = s;
		ptr = new T[s];
    }
    
    ~CPU_Ptr() 
	{
		if ( ptr ) 
		{
            delete[] ptr;
			size = 0;
			ptr = 0;
		}
    }
    
	void setValueZero()
	{
		for (size_t i = 0; i < size; i++)
		{
			ptr[i] = 0;
		}
	}

	void resize(size_t s) 
	{
		if (s == size)
			return;

		if ( ptr ) {
            delete[] ptr;
			size = 0;
			ptr = 0;
		}
		size = s;
		ptr = new T[s];
    }

	operator T * () const
	{
		return ptr;
    }

	int count()
	{
		return size;
    }

	int memorySize()
	{
		return size * sizeof(T);
	}
private:
	T * ptr;
	size_t size;
};

template < typename T > class CPU_SimdPtr
{
public:

	CPU_SimdPtr() : ptr(0) 
	{
		size = 0;
    }

	CPU_SimdPtr(size_t s) : ptr(0)
	{
		size = s;
		ptr = (T*)simd_malloc(sizeof(T) * s);
    }
    
    ~CPU_SimdPtr() 
	{
		if ( ptr ) 
		{
            free(ptr);
			size = 0;
			ptr = 0;
		}
    }
    
	void setValueZero()
	{
		for (size_t i = 0; i < size; i++)
		{
			ptr[i] = 0;
		}
	}

	void resize(size_t s) 
	{
		if (s == size)
			return;

		if ( ptr ) {
            free( ptr);
			size = 0;
			ptr = 0;
		}
		size = s;
		ptr = (T*)simd_malloc(sizeof(T) * s);
    }

	operator T * () const
	{
		return ptr;
    }

	int count()
	{
		return size;
    }

	int memorySize()
	{
		return size * sizeof(T);
	}
private:
	T * ptr;
	size_t size;
};

#endif 

