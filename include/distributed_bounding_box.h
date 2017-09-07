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

		ReduceBoundingBox();
	}

	template<typename TIteratorType>
	DistributedBoundingBox(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd) : 
			BoundingBox<TPointType>(PointsBegin, PointsEnd) {

		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

		ReduceBoundingBox();
	}

	void ReduceBoundingBox() {
		double * sendMinPoint = new double[Dimension];
		double * sendMaxPoint = new double[Dimension];
		double * recvMinPoint = new double[Dimension];
		double * recvMaxPoint = new double[Dimension];

		for(std::size_t i = 0; i < Dimension; i++) {
			sendMinPoint[i] = this->mMinPoint[i];
			sendMaxPoint[i] = this->mMaxPoint[i];
		}

		MPI_Allreduce(sendMinPoint, recvMinPoint, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(sendMaxPoint, recvMaxPoint, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		for(std::size_t i = 0; i < Dimension; i++) {
			this->mMinPoint[i] = recvMinPoint[i];
			this->mMaxPoint[i] = recvMaxPoint[i];
		}

		delete[] sendMinPoint;
		delete[] sendMaxPoint;
		delete[] recvMinPoint;
		delete[] recvMaxPoint;
	}
};