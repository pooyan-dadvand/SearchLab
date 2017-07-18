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
	using KeyType = std::uint64_t;

	ParallelCoherentHash() {
		FillRandomOffsets();
	}

	void Insert(KeyType Key, TDataType& rData) {
		for (int age = 0; age < MaximumAge; age++) {
			auto hash_position = GetHash(Key, age);
			auto& hash_data = mHashTable[hash_position];
			if (hash_data.IsEmpty()) {
				hash_data.SetKey(Key);
				hash_data.SetData(rData);
				return;
			}
		}
		throw std::runtime_error("Inserting failed for all hash functions");
	}

	std::size_t Size() {
		return mSize;
	}

	std::size_t Capacity() {
		return static_cast<std::size_t>(mHashTable.size()*LoadFactor);
	}

	void Reserve(std::size_t NewCapacity) {
		if (Size() == 0)
			mHashTable.resize(static_cast<std::size_t>((NewCapacity / LoadFactor) + 1));
		else
			throw std::runtime_error("Resizing a non empty hash table is not supported yet");
	}

private:
	class HashData {
		char mMaximumAge;
		std::uint64_t mKey;
		TDataType mData;
		bool mIsEmpty;
	public:
		HashData() : mMaximumAge(0), mKey(0), mData(), mIsEmpty(true) {}
		HashData(std::uint64_t Key, TDataType& rData) : mMaximumAge(0), mKey(Key), mData(rData), mIsEmpty(false) {}
		HashData(HashData const& Other) = default;

		void SetMaximumAge(char MaximumAge) { mMaximumAge = MaximumAge; }
		char GetMaximumAge() { return mMaximumAge; }

		void SetKey(KeyType Key) { mKey = Key; }
		std::uint64_t GetKey() { return mKey; }

		void SetData(TDataType& rData) { 
			mIsEmpty = false;
			mData = rData; 
		}

		TDataType& GetData() { return mData; }
		bool IsEmpty() { return mIsEmpty; }
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

	std::uint64_t GetHash(KeyType TheKey, char Age, std::size_t HashTableSize) {
		return (TheKey + mHashOffsets[Age])%HashTableSize;
	}

};
