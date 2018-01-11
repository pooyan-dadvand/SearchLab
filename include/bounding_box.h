#pragma once

/// TPointType should provide access operator [] to its coordinate and deep copy operator=
template <typename TPointType>
class BoundingBox {
	static constexpr int Dimension = 3;
protected:
	TPointType mMinPoint;
	TPointType mMaxPoint;
public:
	BoundingBox(TPointType const& MinPoint, TPointType const& MaxPoint) :
		mMinPoint(MinPoint), mMaxPoint(MaxPoint) {}
	BoundingBox( const BoundingBox &BBox) :
		mMinPoint( BBox.mMinPoint), mMaxPoint( BBox.mMaxPoint) {}

	template<typename TIteratorType>
	BoundingBox(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd) {
		if (PointsBegin == PointsEnd) {
			for (int i = 0; i < Dimension; i++)
			{
				mMinPoint[i] = 0.00;
				mMaxPoint[i] = 0.00;
			}
			return;
		}

		for (int i = 0; i < Dimension; i++)
		{
			mMinPoint[i] = (*PointsBegin)[i];
			mMaxPoint[i] = (*PointsBegin)[i];
		}

		for (TIteratorType Point = PointsBegin; Point != PointsEnd; Point++)
			for (int i = 0; i < Dimension; i++)
			{
				if ((*Point)[i] < mMinPoint[i]) mMinPoint[i] = (*Point)[i];
				if ((*Point)[i] > mMaxPoint[i]) mMaxPoint[i] = (*Point)[i];
			}
	}

	TPointType& GetMinPoint() { return mMinPoint; }
	TPointType const& GetMinPoint() const { return mMinPoint; }

	TPointType& GetMaxPoint() { return mMaxPoint; }
	TPointType const& GetMaxPoint() const { return mMaxPoint; }
};
