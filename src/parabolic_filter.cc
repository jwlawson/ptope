/*
 * parabolic_filter.h
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
#include "parabolic_filter.h"

namespace ptope {
bool ParabolicFilter::operator()(const arma::mat m) {
	_unvisited.clear();
	_components.clear();
	while(!_queue.empty()) {
		_queue.pop();
	}
	for(arma::uword i = 0; i < m.n_cols; ++i) {
		_unvisited.push_back(i);
	}
	arma::uword current_component = 0;
	while(!_unvisited.empty()) {
		_components.push_back(std::vector<arma::uword>());
		_queue.push(_unvisited.back());
		_unvisited.pop_back();
		while(!_queue.empty()) {
			arma::uword col = _queue.front();
			_queue.pop();
			for(auto unv_it = _unvisited.end(); unv_it > _unvisited.begin();) {
				/* Here we need a reverse iterator, wbut one which can be used to erase
				 * elements. Hence this slightly hacky version. */
				--unv_it;
				arma::uword row = *unv_it;
				if(m(row,col) != 0) {
					/* Vertex row is connected to vertex col */
					_unvisited.erase(unv_it);
					_components[current_component].push_back(row);
					_queue.push(row);
				}
			}
		}
	}
	return _components.size() == 1;
}
}

