/*
 * stacked_iterator.h
 * Copyright 2015 John Lawson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once
#ifndef PTOPE_STACKED_ITERATOR_H_
#define PTOPE_STACKED_ITERATOR_H_

#include <utility>

namespace ptope {
template <class FeedIt, class EatIt, class Output>
class StackedIterator {
	public:
		StackedIterator(FeedIt && it) : _it(std::move(it)), _it2(_it.next()) {}
		bool has_next() {
			return _it2.has_next() || _it.has_next();
		}
		const Output & next() {
			/* TODO Here we never check that the newly created iterator actually contains
			 * anything. It is possible that _it2.has_next would return false
			 * immediately after construction. */
			if(!_it2.has_next()) {
				_it2 = EatIt(_it.next());
			}
			return _it2.next();
		}
	private:
		FeedIt _it;
		EatIt _it2;
};
}
#endif


