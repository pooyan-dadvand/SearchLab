#pragma once

#include <array>	
#include "bounding_box.h"

template <typename TObjectType>
class SpatialSearchResult {
	TObjectType* mpObject;
	double mDistance2;
	bool mIsObjectFound;
	bool mIsDistanceCalculated;
public:
	SpatialSearchResult() : mpObject(nullptr), mDistance2(0.00), mIsObjectFound(false), mIsDistanceCalculated(false) {}
	SpatialSearchResult(TObjectType* pObject) : mpObject(pObject), mDistance2(0.00), mIsObjectFound(false), mIsDistanceCalculated(false) {
		if (mpObject != nullptr)
			mIsObjectFound = true;
	}

	SpatialSearchResult(SpatialSearchResult const& Other) = default;

	SpatialSearchResult(SpatialSearchResult&& Other) = default;

	TObjectType* const Get() const { return mpObject; }
	void Set(TObjectType* pObject) {
		mpObject = pObject;
		mIsObjectFound = true;
	}
	bool IsObjectFound() const { return mIsObjectFound; }

	double GetDistance2() const { return mDistance2; }
	void SetDistance2(double TheDistance2) { 
		mDistance2 = TheDistance2; 
		mIsDistanceCalculated = true;
	}
	bool IsDistanceCalculated() const { return mIsDistanceCalculated; }

	void Reset() {
		mpObject = nullptr;
		mDistance = 0.00;
		mIsObjectFound = false;
		mIsDistanceCalculated = false;
	}

	SpatialSearchResult& operator=(SpatialSearchResult const& Other) = default;

};

template <typename TObjectType>
class PointsBins {
	static constexpr int Dimension = 3;
public:
	using InternalPointType = std::array<double, Dimension>;
	using PointerType = TObjectType*;
	using ResultType = SpatialSearchResult<TObjectType>;

	template<typename TIteratorType>
	PointsBins(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd) : mBoundingBox(PointsBegin, PointsEnd) {
		mNumberOfPoints = std::distance(PointsBegin, PointsEnd);
		if (mNumberOfPoints == 0) {
			mNumberOfCells = { 1,1,1 };
			mpPoints = nullptr;
			return;
		}
		CalculateCellSize();
		mpPoints = new PointerType[mNumberOfPoints];
		for (std::size_t i = 0; i < mNumberOfPoints; i++)
			mpPoints[i] = nullptr;

		InitializeCellsOffsets(PointsBegin, PointsEnd);
		AssignPointsToCells(PointsBegin, PointsEnd);
	}

	void SearchInRadius(TObjectType const& ThePoint, double Radius, std::vector<ResultType>& rResults) {
		InternalPointType min_point;
		std::array<std::size_t, Dimension> length;
		const double radius2 = Radius * Radius;

		for (int i = 0; i < Dimension; i++) {
			min_point[i] = ThePoint[i] - Radius;
			length[i] = CalculatePosition(ThePoint[i] + Radius, i) - CalculatePosition(ThePoint[i] - Radius, i) + 1;
		}
		auto min_cell = CalculateCellIndex(min_point);

		for (std::size_t i_z = 0; i_z < length[2]; i_z++) {
			auto y_position = min_cell + i_z * mNumberOfCells[0] * mNumberOfCells[1];
			for (std::size_t i_y = 0; i_y < length[1]; i_y++) {
				for (std::size_t offset = mCellsOffsets[y_position]; offset < mCellsOffsets[y_position + length[0]]; offset++) {
					TObjectType* p_point = mpPoints[offset];
					if (Distance2(*p_point, ThePoint) <= radius2) {
						rResults.push_back(ResultType(p_point));
					}
				}
				y_position += mNumberOfCells[0];
			}
		}
	}

	ResultType SearchNearest(TObjectType const& ThePoint) {
		auto cell_index = CalculateCellIndex(ThePoint);
		ResultType current_result;

		if (mNumberOfPoints == 0)
			return current_result;

		current_result.SetDistance2(std::numeric_limits<double>::max());
		double radius = std::max(mCellSize[0], mCellSize[1]);
		radius = std::max(radius, mCellSize[2]) * .5;

		while (!current_result.IsObjectFound()) {
			InternalPointType min_point;
			std::array<std::size_t, Dimension> length;
			const double radius2 = radius * radius;

			for (int i = 0; i < Dimension; i++) {
				min_point[i] = ThePoint[i] - radius;
				length[i] = CalculatePosition(ThePoint[i] + radius, i) - CalculatePosition(ThePoint[i] - radius, i) + 1;
			}
			auto min_cell = CalculateCellIndex(min_point);

			for (std::size_t i_z = 0; i_z < length[2]; i_z++) {
				auto y_position = min_cell + i_z * mNumberOfCells[0] * mNumberOfCells[1];
				for (std::size_t i_y = 0; i_y < length[1]; i_y++) {
					for (std::size_t offset = mCellsOffsets[y_position]; offset < mCellsOffsets[y_position + length[0]]; offset++) {
						TObjectType* p_point = mpPoints[offset];
						double distance_2 = Distance2(*p_point, ThePoint);
						if (distance_2 < current_result.GetDistance2()) {
							current_result.Set(p_point);
							current_result.SetDistance2(distance_2);
						}
					}
					y_position += mNumberOfCells[0];
				}
			}
			radius *= 2.00;
		}
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
			mNumberOfCells[i] = static_cast<std::size_t>(lengths[i] / avarage_length * avarage_number_of_cells) + 1;
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
		//double number_of_empty_cells = 0;
		//double number_of_single_point_cells = 0;
		//double number_of_multi_point_cells = 0;
		//std::size_t max_cell_occupation = 0;
		//for (std::size_t i_cell_offset = 1; i_cell_offset < mCellsOffsets.size(); i_cell_offset++) {
		//	if (mCellsOffsets[i_cell_offset] == 0)
		//		number_of_empty_cells++;
		//	else if (mCellsOffsets[i_cell_offset] == 1)
		//		number_of_single_point_cells++;
		//	else
		//		number_of_multi_point_cells++;
		//	
		//	if (mCellsOffsets[i_cell_offset] > max_cell_occupation)
		//		max_cell_occupation = mCellsOffsets[i_cell_offset];
		//}

		//std::cout << mCellsOffsets.size() << " cells with " << 100.00* number_of_empty_cells / mCellsOffsets.size() << " empty cells ";
		//std::cout << 100.00* number_of_single_point_cells / mCellsOffsets.size() << " single object cells ";
		//std::cout << 100.00* number_of_multi_point_cells / mCellsOffsets.size() << " multi object cells ";
		//std::cout << "(max occupation = " << max_cell_occupation << ")";

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