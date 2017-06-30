#pragma once

#include <vector>
#include <array>	
#include <random>

template <typename TDataType>
class ParallelCoherentHash {
	static constexpr int MaximumAge = 63;
	static constexpr int Dimension = 3;
	static constexpr double LoadFactor = 0.9;

public:
	ParallelCoherentHash() {
		FillRandomOffsets();
	}

	template<typename TIteratorType>
	ParallelCoherentHash(TIteratorType const& DataBegin, TIteratorType const& DataEnd) {
		FillRandomOffsets();
		Insert(DataBegin, DataEnd);
	}

	template<typename TIteratorType>
	void Insert(TIteratorType	DataBegin, TIteratorType const& DataEnd) {
		for (auto i_data = DataBegin; i_data DataEnd; i_data++)
			Insert(*i_data);
	}

	void Insert(TDataType& rData) {
		HashData new_hash_data(GetHash(rData, 0), &rData);

	}

	std::size_t Size() {
		return mSize;
	}

	std::size_t Capacity() {
		return static_cast<std::size_t>(mHashTable.size()*LoadFactor);
	}

private:
	class HashData {
		char mMaximumAge;
		std::uint64_t mKey;
		TDataType* mpObject;
	public:
		HashData(std::uint64_t Key, TDataType* pObject) : mMaximumAge(0), mKey(Key), mpObject(pObject) {}
		HashData(HashData const& Other) = default;

		void SetMaximumAge(char MaximumAge) { mMaximumAge = MaximumAge; }
		char GetMaximumAge() { return mMaximumAge; }

		std::uint64_t GetKey() { return mKey; }
		TDataType* GetObject() { return mpObject; }
	};

	std::size_t mSize;
	std::vector<HashData> mHashTable;
	std::array<std::size_t, MaximumAge> mHashOffsets;

	void FillRandomOffsets() {
		std::random_device rd;
		for (auto& offset : mHashOffsets) {
			offset = rd();
		}
	}

	std::uint64_t GetHash(TDataType& rObject, char Age) {
		return CalculateCellIndex(rObject) + mHashOffsets[Age];
	}


};
