/*
 * filtered_iterator.h
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
#ifndef PTOPE_FILTERED_ITERATOR_H_
#define PTOPE_FILTERED_ITERATOR_H_
namespace ptope {
template <class It, class Output, class Filter, bool positive>
class FilteredIterator {
	public:
		FilteredIterator(It && it) : _it(it) {}
		bool has_next() {
			return _it.has_next();
		}
		Output next() {
			get_next();
			return _next;
		}
	private:
		It _it;
		Output _next;
		Filter _filter;

		void get_next() {
			if(!_it.has_next()) {
				return;
			}
			_next = _it.next();
			while(_filter(_next) != positive && _it.has_next()) {
				_next = _it.next();
			}
		}
};
}
#endif

