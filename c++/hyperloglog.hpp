//#############################
//# Copyright 2016 Otmar Ertl #
//#############################

#ifndef _HYPERLOGLOG_HPP_
#define _HYPERLOGLOG_HPP_

#include <vector>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <random>
#include <iostream>

class HyperLogLog {
    
    unsigned char p;
    unsigned char q;
    
    std::vector<unsigned char> registers;
    
    std::vector<int> counts;
    
    unsigned char minCount;
    uint_fast64_t registerValueFilter;
    
    HyperLogLog(unsigned char p_, unsigned char q_, std::vector<unsigned char> registers_) : p(p_), q(q_), registers(registers_), counts(q_+2), minCount(0) {
        assert((size_t(1) << p_ ) == registers.size());
        for (unsigned char r : registers) {
            assert(r < q+2);
            counts[r] += 1;
        }
        while(counts[minCount]>0) {
            minCount += 1;
        }
        registerValueFilter = (~uint_fast64_t(0)) << minCount;
    } 

public:
    
    HyperLogLog(unsigned char p_, unsigned char q_) : p(p_), q(q_), registers(std::size_t(1) << p_), counts(q_+2), minCount(0), registerValueFilter(~uint_fast64_t(0)) {
        counts[0] = registers.size();
    }

    static HyperLogLog createFromCounts(const std::vector<int>& counts, std::mt19937_64& rng) {
        
        assert(counts.size() >= 2);
        unsigned char q = counts.size() - 2;
        size_t sum = 0;
        for (int c : counts) {
            sum += c;
        }
        int p;
        std::frexp(sum, &p);
        p -= 1;
        assert(sum == size_t(1) << p);

        std::vector<unsigned char> registers(sum);
        
        size_t idx = 0;
        for (size_t j = 0; j < counts.size(); ++j) {
            for (int k = 0; k < counts[j]; ++k) {
                registers[idx] = j;
                idx += 1;
            }
        }
        
        std::shuffle(registers.begin(), registers.end(), rng);
        
        return HyperLogLog(p, q, registers);
    }
    
    static HyperLogLog merge(const HyperLogLog& hll1, const HyperLogLog& hll2) {
        assert(hll1.p == hll2.p);
        assert(hll1.q == hll2.q);
        
        unsigned char p = hll1.p;
        unsigned char q = hll1.q;
        
        std::vector<unsigned char> registers(hll1.registers.size());
        for (size_t i = 0; i < registers.size(); ++i) {
            registers[i] = std::max(hll1.registers[i], hll2.registers[i]);
        }
        
        return HyperLogLog(p, q, registers);
    }

    static std::vector<int> getJointStatistic(const HyperLogLog& hll1, const HyperLogLog& hll2) {
        assert(hll1.p == hll2.p);
        assert(hll1.q == hll2.q);
        
        size_t stride = hll1.q+2;
        std::vector<int> data(stride*5);
        for (size_t i = 0; i < hll1.registers.size(); ++i) {
            
            unsigned char val1 = hll1.registers[i];
            unsigned char val2 = hll2.registers[i];
            
            if (val1 < val2) {
                data[stride*1 + val2] += 1; // up
                data[stride*2 + val1] += 1; // right
            }
            else if (val1 > val2) {
                data[stride*3 + val2] += 1; // down
                data[stride*4 + val1] += 1; // left
            }
            else {
                data[stride*0 + val1] += 1; // center
            }
        }
        return data;
    }
    
    // example with q = 1, p = 2
    // register values of hll1 : {0, 1, 2, 0}
    // register values of hll2 : {2, 1, 0, 2}
    //
    // the corresponding  count matrix is
    //
    //                / 0 0 2 \
    // count matrix = | 0 1 0 |
    //                \ 1 0 0 /
    //
    // or as array representation
    //
    // count matrix = [0, 0, 2, 0, 1, 0, 1, 0, 0]
    static std::vector<int> getCountMatrix(const HyperLogLog& hll1, const HyperLogLog& hll2) {
        assert(hll1.p == hll2.p);
        assert(hll1.q == hll2.q);
        
        unsigned char q = hll1.q;
        size_t stride = q+2;
        
        std::vector<int> countMatrix(stride*stride);
        
        for (size_t i = 0; i < hll1.registers.size(); ++i) {
            countMatrix[hll1.registers[i] + stride*hll2.registers[i]] += 1;
        }
        
        return countMatrix;
    }

        
    void add(std::uint_fast64_t hashValue) {
        
        if ((registerValueFilter | hashValue) == (~uint_fast64_t(0))) {
            
            std::size_t registerIdx = static_cast<std::size_t>(hashValue >> (64-p));
            unsigned char runLength = 1;
            while(runLength <= q && (hashValue & 1)) {
                runLength += 1;
                hashValue >>= 1;
            }
            
            unsigned char oldRunLength = registers[registerIdx];
            if (oldRunLength < runLength) {
                counts[oldRunLength] -= 1;
                registers[registerIdx] = runLength;
                counts[runLength] += 1;
                if (counts[oldRunLength] == 0 && oldRunLength==minCount) {
                    while(counts[minCount] == 0) {
                        minCount += 1;
                    }
                    registerValueFilter = (~uint_fast64_t(0)) << minCount;
                }
            }
        }
    }
    
    const std::vector<int>& getCounts() const {
        return counts;
    }
    
    HyperLogLog reduce(unsigned char newP, unsigned char newQ) const {
        
        assert(newP <= p);
        assert(newP + newQ <= p + q);
        
        std::size_t oldIdx = 0;
        
        std::vector<unsigned char> newRegisters(std::size_t(1) << newP);
        
        for(std::size_t newIdx = 0; newIdx < (std::size_t(1) << newP); ++newIdx) {
            for(std::size_t subIdx = 0; subIdx < (std::size_t(1) << (p-newP)); ++subIdx) {
                unsigned char runLength = registers[oldIdx];
                std::uint_fast64_t hashValue = oldIdx;
                if (runLength == q+1) {
                    while(runLength <= newQ && (hashValue & 1)) {
                        runLength += 1;
                        hashValue >>= 1;
                    }
                }
                newRegisters[newIdx] = std::max(newRegisters[newIdx], runLength);
                ++oldIdx;
            }
            
            newRegisters[newIdx] = std::min(newRegisters[newIdx], static_cast<unsigned char>(newQ+1));
        }
        
        return HyperLogLog(newP, newQ, newRegisters);
    }
    
};

#endif // _HYPERLOGLOG_HPP_




