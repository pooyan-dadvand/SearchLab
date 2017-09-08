#pragma once

#include "bounding_box.h"

class BinsCellsContainer {
	static constexpr int Dimension = 3;

protected:

	std::array<std::size_t, Dimension> mNumberOfCells;
	using InternalPointType = std::array<double, Dimension>;
	BoundingBox<InternalPointType> mBoundingBox;
	InternalPointType  mCellSize;
	InternalPointType  mInverseOfCellSize;
	std::vector<std::size_t> mCellsBeginIndices;


public:

	template<typename TIteratorType>
	BinsCellsContainer(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd)
		: mBoundingBox(PointsBegin, PointsEnd), mNumberOfCells({ 1,1,1 }) {

		std::size_t approximated_number_of_cells = std::distance(PointsBegin, PointsEnd);

		if (approximated_number_of_cells == 0)
			return;

		CalculateCellSize(approximated_number_of_cells);

		InitializeCellsBeginIndices(PointsBegin, PointsEnd);
	}
	
	void SetNumberOfCells(std::array<std::size_t, Dimension> const& TheNumberOfCells) { 
		mNumberOfCells = TheNumberOfCells; 
	}

	void SetNumberOfCells(std::size_t Axis, std::size_t TheNumberOfCells) { 
		mNumberOfCells[Axis] = TheNumberOfCells; 
	}

	void SetCellBeginIndex(std::size_t Index, std::size_t Value) {
		mCellsBeginIndices[Index] = Value;
	}

	std::size_t GetNumberOfCells(std::size_t Axis) const { 
		return mNumberOfCells[Axis];
	}

	std::size_t GetTotalNumberOfCells() const { 
		return mCellsBeginIndices.size();
	}

	double GetCellSize(std::size_t Axis) {
		return mCellSize[Axis];
	}

	std::size_t GetCellBeginIndex(std::size_t Index) {
		return mCellsBeginIndices[Index];
	}

	template <typename TPointType>
	std::size_t CalculateCellIndex(TPointType const& ThePoint) {
		std::size_t result = 0;
		for (std::size_t i_dim = Dimension - 1; i_dim > 0; i_dim--)
		{
			result += CalculatePosition(ThePoint[i_dim], i_dim);
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

private:

	void CalculateCellSize(std::size_t ApproximatedSize) {
		std::size_t avarage_number_of_cells = static_cast<std::size_t>(std::pow(static_cast<double>(ApproximatedSize), 1.00 / Dimension));
		std::array<double, 3> lengths;
		double avarage_length = 0.00;
		for (int i = 0; i < Dimension; i++) {
			lengths[i] = mBoundingBox.GetMaxPoint()[i] - mBoundingBox.GetMinPoint()[i];
			avarage_length += lengths[i];
		}
		avarage_length *= 1.00 / 3.00;

		if (avarage_length < std::numeric_limits<double>::epsilon()) {
			SetNumberOfCells({ 1,1,1 });
			return;
		}

		for (int i = 0; i < Dimension; i++) {
			SetNumberOfCells(i, static_cast<std::size_t>(lengths[i] / avarage_length * avarage_number_of_cells) + 1);
			if (mNumberOfCells[i] > 1)
				mCellSize[i] = lengths[i] / mNumberOfCells[i];
			else
				mCellSize[i] = avarage_length;

			mInverseOfCellSize[i] = 1.00 / mCellSize[i];
		}

	}

protected:

	template<typename TIteratorType>
	void InitializeCellsBeginIndices(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd) {
		mCellsBeginIndices.resize(mNumberOfCells[0] *mNumberOfCells[1] * mNumberOfCells[2] + 1, 0);
		for (auto i_point = PointsBegin; i_point != PointsEnd; i_point++) {
			mCellsBeginIndices[CalculateCellIndex(*i_point) + 1]++;
		}

		for (std::size_t i_cell_begin = 1; i_cell_begin < mCellsBeginIndices.size(); i_cell_begin++) {
			mCellsBeginIndices[i_cell_begin] += mCellsBeginIndices[i_cell_begin - 1];
		}
	}



};