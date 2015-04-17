#pragma once
#include <cstdint>

template<typename T, std::uint_fast8_t N>
class LookUpTable
{
	public:
		typedef std::size_t uint_t;
		
		LookUpTable(uint_t dim)
		{}
};

/*
template<typename T>
class LookUpTable<T, 2>
{
public:
	typedef std::size_t uint_t;
	
	LookUpTable()
		: dimX(0), dimY(0)
	{
		table = 0;
	}

	LookUpTable(uint_t dimX, uint_t dimY)
	{
		AllocateTable(dimX, dimY);
	}
	virtual ~LookUpTable()
	{
		DeallocateTable();
	}

	LookUpTable(const LookUpTable& rhs)
	{
		AllocateTable(rhs.dimX, rhs.dimY);
		for (uint_t i = 0; i < dimX; ++i)
		{
			for (uint_t j = 0; j < dimY; ++j)
			{
				table[i][j] = rhs.table[i][j];
			}
		}
	}
	LookUpTable& operator=(const LookUpTable& rhs)
	{
		if (*this == rhs)
			return *this;
		if (dimX != rhs.dimX || dimY != rhs.dimY)
		{
			DeallocateTable();
			AllocateTable(rhs.dimX, rhs.dimY);
		}
		for (uint_t i = 0; i < dimX; ++i)
		{
			for (uint_t j = 0; j < dimY; ++j)
			{
				table[i][j] = rhs.table[i][j];
			}
		}
		return *this;
	}

	T*& operator[](uint_t index)
	{
		return table[index];
	}

	const T*& operator[](uint_t index) const
	{
		return table[index];
	}

	void AllocateTable(uint_t dimX, uint_t dimY)
	{
		this->dimX = dimX;
		this->dimY = dimY;
		table = new T*[dimX];
		for (uint_t i = 0; i < dimX; ++i)
			table[i] = new T[dimY];
	}

	void DeallocateTable()
	{
		for (uint_t i = 0; i < dimX; ++i)
			delete[] table[i];
		delete[] table;
	}
private:
	T** table;
	uint_t dimX;
	uint_t dimY;
};
*/

template<typename T>
class LookUpTable<T, 2>
{
public:
	typedef std::size_t uint_t;
	
	LookUpTable()
		: dimX(0), dimY(0)
	{
		table = 0;
	}

	LookUpTable(uint_t dimX, uint_t dimY)
	{
		AllocateTable(dimX, dimY);
	}
	virtual ~LookUpTable()
	{
		DeallocateTable();
	}

	LookUpTable(const LookUpTable& rhs)
	{
		AllocateTable(rhs.dimX, rhs.dimY);
		for (uint_t i = 0; i < dimX; ++i)
		{
			for (uint_t j = 0; j < dimY; ++j)
			{
				table[i*dimY + j] = rhs.table[i];
			}
		}
	}
	LookUpTable& operator=(const LookUpTable& rhs)
	{
		if (*this == rhs)
			return *this;
		if (dimX != rhs.dimX || dimY != rhs.dimY)
		{
			DeallocateTable();
			AllocateTable(rhs.dimX, rhs.dimY);
		}
		for (uint_t i = 0; i < dimX; ++i)
		{
			for (uint_t j = 0; j < dimY; ++j)
			{
				table[i*dimY + j] = rhs.table[i];
			}
		}
		return *this;
	}

	inline T& operator()(uint_t i1, uint_t i2)
	{
		return table[i1*dimY + i2];
	}

	inline const T& operator()(uint_t i1, uint_t i2) const
	{
		return table[i1*dimY + i2];
	}

	void AllocateTable(uint_t dimX, uint_t dimY)
	{
		this->dimX = dimX;
		this->dimY = dimY;
		table = new T[dimX*dimY];
	}

	void DeallocateTable()
	{
		delete[] table;
	}
private:
	T* table;
	uint_t dimX;
	uint_t dimY;
};

template<typename T>
class LookUpTable<T, 3>
{
	public:
		typedef std::size_t uint_t;
		
		LookUpTable()
			: dimX(0), dimY(0), dimZ(0)
		{
			table = 0;
		}

		LookUpTable(uint_t dimX, uint_t dimY, uint_t dimZ)
		{
			AllocateTable(dimX, dimY, dimZ);
		}
		virtual ~LookUpTable()
		{
			DeallocateTable();
		}

		LookUpTable(const LookUpTable& rhs)
		{
			AllocateTable(rhs.dimX, rhs.dimY, rhs.dimZ);
			for (uint_t i = 0; i < dimX; ++i)
			{
				for (uint_t j = 0; j < dimY; ++j)
				{
					for (uint_t k = 0; k < dimZ; ++k)
					{
						table[i][j][k] = rhs.table[i][j][k];
					}
				}
			}
		}
		LookUpTable& operator=(const LookUpTable& rhs)
		{
			if (*this == rhs)
				return *this;
			if (dimX != rhs.dimX || dimY != rhs.dimY || dimZ != rhs.dimZ)
			{
				DeallocateTable();
				AllocateTable(rhs.dimX, rhs.dimY, rhs.dimZ);
			}
			for (uint_t i = 0; i < dimX; ++i)
			{
				for (uint_t j = 0; j < dimY; ++j)
				{
					for (uint_t k = 0; k < dimZ; ++k)
					{
						table[i][j][k] = rhs.table[i][j][k];
					}
				}
			}
			return *this;
		}

		inline T**& operator[](uint_t index)
		{
			return table[index];
		}
		
		inline const T**& operator[](uint_t index) const
		{
			return table[index];
		}

		void AllocateTable(uint_t dimX, uint_t dimY, uint_t dimZ)
		{
			this->dimX = dimX;
			this->dimY = dimY;
			this->dimZ = dimZ;
			table = new T**[dimX];
			for (uint_t i = 0; i < dimX; ++i)
				table[i] = new T*[dimY];
			for (uint_t i = 0; i < dimX; ++i)
				for (uint_t j = 0; j < dimY; ++j)
					table[i][j] = new T[dimZ];
		}

		void DeallocateTable()
		{
			for (uint_t i = 0; i < dimX; ++i)
				for (uint_t j = 0; j < dimY; ++j)
					delete[] table[i][j];
			for (uint_t i = 0; i < dimX; ++i)
				delete[] table[i];
			delete[] table;
		}
	private:
		T*** table;
		uint_t dimX;
		uint_t dimY;
		uint_t dimZ;
};