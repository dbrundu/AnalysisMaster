/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2016 Antonio Augusto Alves Junior e Andrea Contu
 *
 *   This file is part of MassMeasurement software.
 *
 *   MassMeasurement is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   MassMeasurement is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Hydra.  If not, see <http://www.gnu.org/licenses/>.
 *
 *---------------------------------------------------------------------------*/

/*
 * utils.h
 *
 *  Created on: 30/04/2018
 *      Author: augalves
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>

#define VERBOSE_LINE(flag, message)\
if(flag){\
	std::cout<< "\033[1;35mVerbose: \033[0m"\
			 << message << std::endl;\
}


#define INFO_LINE(message)\
	std::cout<< "\033[1;34mInfo: \033[0m"\
			 << message << std::endl;\

#define WARNING_LINE(message)\
	std::cout<< "\033[1;34mWarning: \033[0m"\
			 << message << std::endl;\




#endif /* UTILS_H_ */
