#pragma once

#include "mpi.h"
#include "bounding_box.h"

/// TPointType should provide access operator [] to its coordinate and deep copy operator=
template <typename TPointType>
class DistributedBoundingBox : public BoundingBox<TPointType> {
	static constexpr int Dimension = 3;

	int mpi_size;
  int mpi_rank;

public:
	DistributedBoundingBox(TPointType const& MinPoint, TPointType const& MaxPoint) :
			BoundingBox<TPointType>(MinPoint, MaxPoint) {

		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

		CalculateGlobalBoundingBox();
	}

	template<typename TIteratorType>
	DistributedBoundingBox(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd) : 
			BoundingBox<TPointType>(PointsBegin, PointsEnd) {

		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

		CalculateGlobalBoundingBox();
	}

	void CalculateGlobalBoundingBox() {
		TPointType recvMinPoint;
		TPointType recvMaxPoint;

		MPI_Allreduce(&this->mMinPoint[0], &recvMinPoint[0], 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&this->mMaxPoint[0], &recvMaxPoint[0], 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		for(std::size_t i = 0; i < Dimension; i++) {
			this->mMinPoint[i] = recvMinPoint[i];
			this->mMaxPoint[i] = recvMaxPoint[i];
		}
	}
};