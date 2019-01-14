/* The MIT License

   Copyright (c) 2011 Erik Garrison <erik.garrison@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#ifndef __GTF_INTERVAL_TREE_H
#define __GTF_INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>
#include <cassert>

template <class Scalar, typename Value>
class gtfInterval {
public:
    Scalar start;
    Scalar stop;
    Value value;
    gtfInterval(const Scalar& s, const Scalar& e, const Value& v)
    : start(std::min(s, e))
    , stop(std::max(s, e))
    , value(v) 
    {}
};

template <class Scalar, typename Value>
Value intervalStart(const gtfInterval<Scalar,Value>& i) {
    return i.start;
}

template <class Scalar, typename Value>
Value intervalStop(const gtfInterval<Scalar, Value>& i) {
    return i.stop;
}

template <class Scalar, typename Value>
std::ostream& operator<<(std::ostream& out, const gtfInterval<Scalar, Value>& i) {
    out << "gtfInterval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class Scalar, class Value>
class gtfIntervalTree {
public:
    typedef gtfInterval<Scalar, Value> gtf_interval;
    typedef std::vector<gtf_interval> gtf_interval_vector;


    struct gtfIntervalStartCmp {
        bool operator()(const gtf_interval& a, const gtf_interval& b) {
            return a.start < b.start;
        }
    };

    struct gtfIntervalStopCmp {
        bool operator()(const gtf_interval& a, const gtf_interval& b) {
            return a.stop < b.stop;
        }
    };

    gtfIntervalTree()
        : left(nullptr)
        , right(nullptr)
        , center(0)
    {}

    ~gtfIntervalTree() = default;

    std::unique_ptr<gtfIntervalTree> clone() const {
        return std::unique_ptr<gtfIntervalTree>(new gtfIntervalTree(*this));
    }

    gtfIntervalTree(const gtfIntervalTree& other)
    :   gtf_intervals(other.gtf_intervals),
        left(other.left ? other.left->clone() : nullptr),
        right(other.right ? other.right->clone() : nullptr),
        center(other.center)
    {}

    gtfIntervalTree& operator=(gtfIntervalTree&&) = default;
    gtfIntervalTree(gtfIntervalTree&&) = default;

    gtfIntervalTree& operator=(const gtfIntervalTree& other) {
        center = other.center;
        gtf_intervals = other.gtf_intervals;
        left = other.left ? other.left->clone() : nullptr;
        right = other.right ? other.right->clone() : nullptr;
        return *this;
    }

    gtfIntervalTree(
            gtf_interval_vector&& ivals,
            std::size_t depth = 16,
            std::size_t minbucket = 64,
            std::size_t maxbucket = 512, 
            Scalar leftextent = 0,
            Scalar rightextent = 0)
      : left(nullptr)
      , right(nullptr)
    {
        --depth;
        const auto minmaxStop = std::minmax_element(ivals.begin(), ivals.end(), 
                                                    gtfIntervalStopCmp());
        const auto minmaxStart = std::minmax_element(ivals.begin(), ivals.end(), 
                                                     gtfIntervalStartCmp());
        if (!ivals.empty()) {
            center = (minmaxStart.first->start + minmaxStop.second->stop) / 2;
        }
        if (leftextent == 0 && rightextent == 0) {
            // sort intervals by start
            std::sort(ivals.begin(), ivals.end(), gtfIntervalStartCmp());
        } else {
            assert(std::is_sorted(ivals.begin(), ivals.end(), gtfIntervalStartCmp()));
        }
        if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
            std::sort(ivals.begin(), ivals.end(), gtfIntervalStartCmp());
            gtf_intervals = std::move(ivals);
            assert(is_valid().first);
            return;
        } else {
            Scalar leftp = 0;
            Scalar rightp = 0;

            if (leftextent || rightextent) {
                leftp = leftextent;
                rightp = rightextent;
            } else {
                leftp = ivals.front().start;
                rightp = std::max_element(ivals.begin(), ivals.end(),
                                          gtfIntervalStopCmp())->stop;
            }

            gtf_interval_vector lefts;
            gtf_interval_vector rights;

            for (typename gtf_interval_vector::const_iterator i = ivals.begin(); 
                 i != ivals.end(); ++i) {
                const gtf_interval& iv = *i;
                if (iv.stop < center) {
                    lefts.push_back(iv);
                } else if (iv.start > center) {
                    rights.push_back(iv);
                } else {
                    assert(iv.start <= center);
                    assert(center <= iv.stop);
                    gtf_intervals.push_back(iv);
                }
            }

            if (!lefts.empty()) {
                left.reset(new gtfIntervalTree(std::move(lefts), 
                                            depth, minbucket, maxbucket,
                                            leftp, center));
            }
            if (!rights.empty()) {
                right.reset(new gtfIntervalTree(std::move(rights), 
                                             depth, minbucket, maxbucket, 
                                             center, rightp));
            }
        }
        assert(is_valid().first);
    }

    // Call f on all intervals near the range [start, stop]:
    template <class UnaryFunction>
    void visit_near(const Scalar& start, const Scalar& stop, UnaryFunction f) const {
        if (!gtf_intervals.empty() && ! (stop < gtf_intervals.front().start)) {
            for (auto & i : gtf_intervals) {
              f(i);
            }
        }
        if (left && start <= center) {
            left->visit_near(start, stop, f);
        }
        if (right && stop >= center) {
            right->visit_near(start, stop, f);
        }
    }

    // Call f on all gtf_intervals crossing pos
    template <class UnaryFunction>
    void visit_overlapping(const Scalar& pos, UnaryFunction f) const {
        visit_overlapping(pos, pos, f);
    }

    // Call f on all intervals overlapping [start, stop]
    template <class UnaryFunction>
    void visit_overlapping(const Scalar& start, const Scalar& stop, UnaryFunction f) const {
        auto filterF = [&](const gtf_interval& interval) {
            if (interval.stop >= start && interval.start <= stop) {
                // Only apply f if overlapping
                f(interval);
            }
        };
        visit_near(start, stop, filterF);
    }

    // Call f on all intervals contained within [start, stop]
    template <class UnaryFunction>
    void visit_contained(const Scalar& start, const Scalar& stop, UnaryFunction f) const {
        auto filterF = [&](const gtf_interval& iv) {
            if (start <= iv.start && iv.stop <= stop) {
                f(iv);
            }
        };
        visit_near(start, stop, filterF);
    }

    gtf_interval_vector findOverlapping(const Scalar& start, const Scalar& stop) const {
        gtf_interval_vector result;
        visit_overlapping(start, stop,
                          [&](const gtf_interval& interval) { 
                            result.emplace_back(interval); 
                          });
        return result;
    }

    gtf_interval_vector findContained(const Scalar& start, const Scalar& stop) const {
        gtf_interval_vector result;
        visit_contained(start, stop,
                        [&](const gtf_interval& interval) { 
                          result.push_back(interval); 
                        });
        return result;
    }
    
    bool empty() const {
        if (left && !left->empty()) {
            return false;
        }
        if (!gtf_intervals.empty()) { 
            return false;
        }
        if (right && !right->empty()) {
            return false;
        }
        return true;
    }

    template <class UnaryFunction>
    void visit_all(UnaryFunction f) const {
        if (left) {
            left->visit_all(f);
        }
        std::for_each(gtf_intervals.begin(), gtf_intervals.end(), f);
        if (right) {
            right->visit_all(f);
        }
    }

    std::pair<Scalar, Scalar> extentBruitForce() const {
        struct Extent {
            std::pair<Scalar, Scalar> x = {std::numeric_limits<Scalar>::max(),
                                                       std::numeric_limits<Scalar>::min() };
            void operator()(const gtf_interval & iv) {
                x.first  = std::min(x.first,  iv.start);
                x.second = std::max(x.second, iv.stop);
            }
                                                                };
	Extent extent;
					    
        visit_all([&](const gtf_interval & interval) { extent(interval); });
        return extent.x;
                                            }

    // Check all constraints.
    // If first is false, second is invalid.
    std::pair<bool, std::pair<Scalar, Scalar>> is_valid() const {
        const auto minmaxStop = std::minmax_element(gtf_intervals.begin(), gtf_intervals.end(), 
                                                    gtfIntervalStopCmp());
        const auto minmaxStart = std::minmax_element(gtf_intervals.begin(), gtf_intervals.end(), 
                                                     gtfIntervalStartCmp());
        
        std::pair<bool, std::pair<Scalar, Scalar>> result = {true, { std::numeric_limits<Scalar>::max(),
                                                                     std::numeric_limits<Scalar>::min() }};
        if (!gtf_intervals.empty()) {
            result.second.first   = std::min(result.second.first,  minmaxStart.first->start);
            result.second.second  = std::min(result.second.second, minmaxStop.second->stop);
        }
        if (left) {
            auto valid = left->is_valid();
            result.first &= valid.first;
            result.second.first   = std::min(result.second.first,  valid.second.first);
            result.second.second  = std::min(result.second.second, valid.second.second);
            if (!result.first) { return result; }
            if (valid.second.second >= center) {
                result.first = false;
                return result;
            }
        }
        if (right) {
            auto valid = right->is_valid();
            result.first &= valid.first;
            result.second.first   = std::min(result.second.first,  valid.second.first);
            result.second.second  = std::min(result.second.second, valid.second.second);
            if (!result.first) { return result; }
            if (valid.second.first <= center) { 
                result.first = false;
                return result;
            }
        }
        if (!std::is_sorted(gtf_intervals.begin(), gtf_intervals.end(), gtfIntervalStartCmp())) {
            result.first = false;
        }
        return result;        
    }

    friend std::ostream& operator<<(std::ostream& os, const gtfIntervalTree& itree) {
        return writeOut(os, itree);
    }

    friend std::ostream& writeOut(std::ostream& os, const gtfIntervalTree& itree, 
                                  std::size_t depth = 0) {
        auto pad = [&]() { for (std::size_t i = 0; i != depth; ++i) { os << ' '; } };
        pad(); os << "center: " << itree.center << '\n';
        for (const gtf_interval & inter : itree.gtf_intervals) {
            pad(); os << inter << '\n';
        }
        if (itree.left) {
            pad(); os << "left:\n";
            writeOut(os, *itree.left, depth + 1);
        } else {
            pad(); os << "left: nullptr\n";
        }
        if (itree.right) {
            pad(); os << "right:\n";
            writeOut(os, *itree.right, depth + 1);
        } else {
            pad(); os << "right: nullptr\n";
        }
        return os;
    }

private:
    gtf_interval_vector gtf_intervals;
    std::unique_ptr<gtfIntervalTree> left;
    std::unique_ptr<gtfIntervalTree> right;
    Scalar center;
};

#endif
