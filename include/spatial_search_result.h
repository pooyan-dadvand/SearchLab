#pragma once

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

	//SpatialSearchResult(SpatialSearchResult&& Other) = default;

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
		mDistance2 = 0.00;
		mIsObjectFound = false;
		mIsDistanceCalculated = false;
	}

	SpatialSearchResult& operator=(SpatialSearchResult const& Other) = default;

};
