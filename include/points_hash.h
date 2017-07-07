#pragma once

#include <limits>
#include <cmath>

#include "bounding_box.h"
#include "spatial_search_result.h"
#include "ParallelCoherentHash.h"

template <typename TObjectType>
class PointsHash {
	static constexpr int Dimension = 3;
public:
	using InternalPointType = std::array<double, Dimension>;
	using PointerType = TObjectType*;
	using ResultType = SpatialSearchResult<TObjectType>;

	template<typename TIteratorType>
	PointsHash(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd) : mBoundingBox(PointsBegin, PointsEnd) {

		mNumberOfPoints = std::distance(PointsBegin, PointsEnd);
		if (mNumberOfPoints == 0) {
			mNumberOfCells = { 1,1,1 };
			mpPoints = nullptr;
			return;
		}
		CalculateCellSize();
		for (auto number_of_cells : mNumberOfCells)
			std::cout << number_of_cells << std::endl;

	}

	void SearchInRadius(TObjectType const& ThePoint, double Radius, std::vector<ResultType>& rResults) {
		//InternalPointType min_point;
		//std::array<std::size_t, Dimension> length;
		//const double radius2 = Radius * Radius;

		//for (int i = 0; i < Dimension; i++) {
		//	min_point[i] = ThePoint[i] - Radius;
		//	length[i] = CalculatePosition(ThePoint[i] + Radius, i) - CalculatePosition(ThePoint[i] - Radius, i) + 1;
		//}
		//auto min_cell = CalculateCellIndex(min_point);

		//for (std::size_t i_z = 0; i_z < length[2]; i_z++) {
		//	auto y_position = min_cell + i_z * mNumberOfCells[0] * mNumberOfCells[1];
		//	for (std::size_t i_y = 0; i_y < length[1]; i_y++) {
		//		for (std::size_t offset = mCellsOffsets[y_position]; offset < mCellsOffsets[y_position + length[0]]; offset++) {
		//			TObjectType* p_point = mpPoints[offset];
		//			if (Distance2(*p_point, ThePoint) <= radius2) {
		//				rResults.push_back(ResultType(p_point));
		//			}
		//		}
		//		y_position += mNumberOfCells[0];
		//	}
		//}
	}

	ResultType SearchNearest(TObjectType const& ThePoint) {
		//auto cell_index = CalculateCellIndex(ThePoint);
		ResultType current_result;

		//if (mNumberOfPoints == 0)
		//	return current_result;

		//current_result.SetDistance2(std::numeric_limits<double>::max());
		//double radius = std::max(mCellSize[0], mCellSize[1]);
		//radius = std::max(radius, mCellSize[2]) * .5;

		//while (!current_result.IsObjectFound()) {
		//	InternalPointType min_point;
		//	std::array<std::size_t, Dimension> length;
		//	const double radius2 = radius * radius;

		//	for (int i = 0; i < Dimension; i++) {
		//		min_point[i] = ThePoint[i] - radius;
		//		length[i] = CalculatePosition(ThePoint[i] + radius, i) - CalculatePosition(ThePoint[i] - radius, i) + 1;
		//	}
		//	auto min_cell = CalculateCellIndex(min_point);

		//	for (std::size_t i_z = 0; i_z < length[2]; i_z++) {
		//		auto y_position = min_cell + i_z * mNumberOfCells[0] * mNumberOfCells[1];
		//		for (std::size_t i_y = 0; i_y < length[1]; i_y++) {
		//			for (std::size_t offset = mCellsOffsets[y_position]; offset < mCellsOffsets[y_position + length[0]]; offset++) {
		//				TObjectType* p_point = mpPoints[offset];
		//				double distance_2 = Distance2(*p_point, ThePoint);
		//				if (distance_2 < current_result.GetDistance2()) {
		//					current_result.Set(p_point);
		//					current_result.SetDistance2(distance_2);
		//				}
		//			}
		//			y_position += mNumberOfCells[0];
		//		}
		//	}
		//	radius *= 2.00;
		//}
		return current_result;
	}


private:
	std::size_t mNumberOfPoints;
	std::array<std::size_t, Dimension> mNumberOfCells;
	BoundingBox<InternalPointType> mBoundingBox;
	InternalPointType  mCellSize;
	InternalPointType  mInverseOfCellSize;
	std::vector<std::size_t> mCellsOffsets;
	TObjectType** mpPoints;
	ParallelCoherentHash<std::size_t> mCells;


	template <typename TPointType>
	std::size_t CalculateCellIndex(TPointType const& ThePoint) {
		std::size_t result = 0;
		for (std::size_t i_dim = Dimension - 1; i_dim > 0; i_dim--)
		{
			result += CalculatePosition(ThePoint[i_dim], i_dim) ;
			result *= mNumberOfCells[i_dim - 1];
		}
		result += CalculatePosition(ThePoint[0], 0);
		return result;
	}

	std::size_t CalculatePosition(double Coordinate, int ThisDimension) {
		auto distance = Coordinate - mBoundingBox.GetMinPoint()[ThisDimension];
		distance = (distance < 0.00) ? 0.00 : distance;
	    std:size_t position = static_cast<std::size_t>(distance * mInverseOfCellSize[ThisDimension]);
		return (position > mNumberOfCells[ThisDimension] - 1) ? mNumberOfCells[ThisDimension] - 1 : position;
	}


	void CalculateCellSize() {
		std::size_t avarage_number_of_cells = static_cast<std::size_t>(std::pow(static_cast<double>(mNumberOfPoints), 1.00 / Dimension));
		std::array<double, 3> lengths;
		double avarage_length = 0.00;
		for (int i = 0; i < Dimension; i++) {
			lengths[i] = mBoundingBox.GetMaxPoint()[i] - mBoundingBox.GetMinPoint()[i];
			avarage_length += lengths[i];
		}
		avarage_length *= 1.00 / 3.00;

		if (avarage_length < std::numeric_limits<double>::epsilon()) {
			mNumberOfCells = { 1,1,1 };
			return;
		}

		for (int i = 0; i < Dimension; i++) {
			std::size_t number_of_cells_in_direction = static_cast<std::size_t>(lengths[i] / avarage_length * avarage_number_of_cells) + 1;
			std::size_t power_of_two = 1;
			while (power_of_two < number_of_cells_in_direction)
				power_of_two <<= 1;
			mNumberOfCells[i] = power_of_two;
			if (mNumberOfCells[i] > 1)
				mCellSize[i] = lengths[i] / mNumberOfCells[i];
			else
				mCellSize[i] = avarage_length;

			mInverseOfCellSize[i] = 1.00 / mCellSize[i];
		}

	}


	template<typename TIteratorType>
	void InitializeCellsOffsets(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd) {
		mCellsOffsets.resize(mNumberOfCells[0] * mNumberOfCells[1] * mNumberOfCells[2]+1, 0);
		for (auto i_point = PointsBegin; i_point != PointsEnd; i_point++) {
			mCellsOffsets[CalculateCellIndex(*i_point)+1]++;
		}

		for (std::size_t i_cell_offset = 1; i_cell_offset < mCellsOffsets.size(); i_cell_offset++) {
			mCellsOffsets[i_cell_offset] += mCellsOffsets[i_cell_offset - 1];
		}
	}

	template<typename TIteratorType>
	void AssignPointsToCells(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd) {
		for (auto i_point = PointsBegin; i_point != PointsEnd; i_point++) {
			auto index = CalculateCellIndex(*i_point);
			for (std::size_t offset = mCellsOffsets[index]; offset < mCellsOffsets[index + 1]; offset++)
				if (mpPoints[offset] == nullptr) {
					mpPoints[offset] = &(*i_point);
					break;
				}
		}

	}

	void SearchNearestInCell(std::size_t CellIndex, TObjectType const& ThePoint, ResultType& rCurrentResult) {
		for (std::size_t offset = mCellsOffsets[CellIndex]; offset < mCellsOffsets[CellIndex + 1]; offset++) {
			TObjectType* p_point = mpPoints[offset];
			auto distance_2 = Distance2(*p_point, ThePoint);
			if ( distance_2 <= rCurrentResult.GetDistance2()) {
				rCurrentResult.Set(p_point);
				rCurrentResult.SetDistance2(distance_2);
			}
		}
	}

	double Distance2(TObjectType const& FirstPoint, TObjectType const& SecondPoint) {
		double result = double();
		for (int i = 0; i < Dimension; i++) {
			auto distance_i = FirstPoint[i] - SecondPoint[i];
			result += distance_i * distance_i;
		}
		return result;
	}
};